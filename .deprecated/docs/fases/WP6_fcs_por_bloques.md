# WP6. Fase 4: FCS por bloques de columnas y etapa `GENE_P`

Estado: INTEGRACIÓN IMPLEMENTADA Y VALIDADA CON PRUEBAS DIRIGIDAS (2026-10-23); aceptación integrada y medidas reales a escala pendientes. Depende de: WP1 (clasificación nueva, RER paramétrico) y WP5 (matrices por bloque).
Conversación independiente. Antes de empezar, leer `docs/fases/README.md`.

## Objetivo

Que el FCS (CAAS y RER) procese bloques de columnas del nulo y acumule contadores, con equivalencia al contrato AUC vigente sobre los mismos rankings. A tamaño de bloque fijo, el recuento de rutas no retiene todas las columnas; GENE_P sí retiene todos los ciclos de un lote de genes y su memoria crece con N.

## Decisiones aplicables (del usuario)

- Ranking del FCS de CAAS = −log10 del p de permulación del gen, construido sobre el parámetro `gene_stat` (el estadístico base lo decide WP9; sin interruptor de ranking crudo en el código: la comparación con el crudo es de validación).
- GENE_P por transposición previa: un paso por bloque escribe trozos por lote de genes; cada GENE_P lee solo sus trozos.
- RER: `RER_CONT` escribe `corStat` por lote, con semilla por lote como parámetro (desviación de RERconverge a justificar); sustituye a `RER_PERMS_SPLIT`.
- RER: estadístico paramétrico (WP1); permulaciones a 10.000.
- El coste del FCS es asumible si también lo es para RER; aun así se implementa por bloques.

## Estado de partida

- Medido (FCS de neoplasia completo con el ranking calibrado; por encima de 1000 ciclos, columnas recicladas: los recursos valen y los p-valores no): 250, 500, 1000, 2000, 4000 y 10.000 ciclos → 2,8, 3,3, 5,2, 7,9, 14,4 y 32,1 GB; 7,5, 7,9, 12,0, 16,5, 27,1 y 58,5 min. Crecimiento lineal de ~3 MB y ~0,3 s por ciclo.
- `fcs_run_all` (`fcs_enrich.R:~620-808`) no separa el estadístico observado, el recuento nulo y el cierre: el observado vive en `enrich_df` y `p.perm` se escribe en su sitio; el nulo se carga entero (`readRDS`, todas las matrices) y se hacen copias densas. `FCS_COMPUTE_BATCHED` ya reparte por GMT, pero cada tarea lee el RDS completo. `fcs_align_null` es independiente por columna. El p de permulación es un recuento por conjunto, aditivo entre bloques: `p.perm = (Σ recuento + 1)/(Σ columnas válidas + 1)`.
- RER: `RER_CONT` es una tarea serie de R (8 CPU, 32 GB × intento) que guarda tres data.frames genes × N (`corP`, `corRho`, `corStat`) y los lee todos cada tarea de FCS; a 10.000 ciclos son ~3,6 GB residentes más copias densas (~6 GB antes del trabajo de conjuntos × N).

## Hallazgo que fija el diseño

Con el ranking por p de gen el nulo **no es independiente por columna**: el −log10 p de un gen en el ciclo j se puntúa contra los otros N−1 ciclos (`ge/N`, como `loo_log10p` de los controles y `fcs_empirical_corstat` en RER antiguo), de modo que necesita las puntuaciones del gen en todos los ciclos. Hace falta una etapa `GENE_P` antes de los bloques del FCS: dispersión por lote de genes, que lee las filas de esos genes de todos los bloques y escribe los p por bloque. Pseudo-cuenta decidida una vez y aplicada igual al observado y al nulo (observado: `(1 + c)/(N + 1)`, con c = nº de ciclos nulos con score ≥ obs; nulo: `(1 + ge)/N`, con ge = nº de otros ciclos con score ≥ el del ciclo; la regla es `(1 + recuento)/(M + 1)` con M = nº de comparaciones, así que `ge/N` daría p = 0 cuando el ciclo tiene el máximo estricto del gen), con comparaciones por `TIE_TOL`. En RER no hace falta (paramétrico).

## Alcance (ficheros)

`subworkflows/ENRICHMENT/local/src/fcs_enrich.R`, `fcs_compute.R`, `subworkflows/ENRICHMENT/fcs.nf` (`FCS_COMPUTE_BATCHED` ~194, `FCS_CONCAT` ~238), `workflows/enrichment.nf` (~212 a 230), `workflows/rerconverge.nf` (~196 a 207), `subworkflows/RERCONVERGE/continuous_rer.R` y `rer_cont.nf`, `conf/resources.config` (labels `process_fcs_batched`, `RER_CONT`), GUI/plantillas si se añaden parámetros.

## Cambios

1. Partir `fcs_run_all` en `fcs_observed()` (tabla observada y `realenrich` por ranking), `fcs_null_count(bloque)` (por conjunto: recuento de estadísticos nulos ≥ observado y columnas válidas) y `fcs_close()` (`p.perm`, BH, clasificación de WP1). `fcs_run_all` pasa a ser la composición de un solo bloque.
2. Procesos: `FCS_OBSERVED` por lote de GMT; `GENE_P` (CAAS); `FCS_NULL_BLOCK` por lote de GMT × bloque; `FCS_CLOSE`; `FCS_CONCAT` sin cambio.
3. RER: un paso único `RER_PERMS_SPLIT` que escribe columnas de `corStat` por bloques (RDS reducido de solo `corStat`), para que cada tarea no lea el objeto completo. Si se reparten los lotes de permulaciones de RER en tareas paralelas con semillas por lote, no es reproducible bit a bit y exige su propia justificación como desviación de RERconverge.

## Interfaz hacia WP7 y WP8

Ficheros de bloque de `corStat` y de p de gen; `fcs_observed()`, `fcs_null_count()` y `fcs_close()` como API.

## Validación (criterios de paso)

- La composición nueva frente a `legacy_fcs/` (WP0) sobre marine toy, PEPC y un conjunto sintético de 400 genes y 25 columnas: recuentos exactamente iguales y `stat` a 1e-12.
- Resultado invariante al tamaño de bloque, probado con tamaños {1, 7, N}.
- Ampliar `PhyloPhere_validation/unification/check_fcs_null_equivalence.R`.
- Medir RSS y tiempo a 1000, 4000 y 10.000 ciclos con bloques: memoria por tarea independiente de N (criterio de paso de memoria: pocos GB).

## Fuera de alcance

Construcción de las matrices por bloque (WP5), formato del RDS único (WP7).

## Riesgos y puntos a verificar

- Ranking firmado (RER): dirección por `sign(stat)`; el recuento bilateral usa `|nulo| ≥ |obs|`.
- Si `GENE_P` lee las filas de un lote de genes de todos los bloques, el patrón de E/S puede ser caro; elegir el diseño de almacenamiento por bloque con eso en mente (por ejemplo, filas agrupadas por lote de genes).

## Registro de la fase

```
- Conversación y fechas:
- Commit de partida / commit de cierre:
- Qué se hizo (ficheros y cambios):
- Comprobaciones ejecutadas (comando, conjunto de datos, tolerancia, resultado):
- Mediciones (recurso antes y después):
- Desviaciones del plan y decisiones tomadas (con su justificación):
- Pendiente que pasa a otra fase:
```

## Avance temático guardado (2026-10-22)

Se versionan `gene_p.R`, `caas_blocks.R`, `fcs_blocks.R` y `fcs_blocks.nf`: transposición, calibración por gen, observado, recuentos por bloque y cierre. La prueba de contrato AUC compartida con WP1 se incluye en `tests/test_fcs_auc_contract.R` y utiliza el helper `gene_p.R` de esta fase. El camino productivo sigue leyendo el RDS único; este commit no activa el subflujo de bloques.

Pendientes: cableado productivo y parámetros/GUI, RER por lotes, calibración amplia, pruebas de memoria/E/S y aceptación integrada. Las referencias al comparador anterior de p analíticos de CAAS deben actualizarse para usar la decisión AUC de WP1. Por decisión del usuario no se hace una corrida completa por WP; el avance se guarda para continuar la integración y ejecutar después la aceptación conjunta.

Pruebas dirigidas antes del commit parcial: `test_gene_p.py` (oráculo de conteo exhaustivo) y `test_caas_blocks.py::test_batches_over_blocks_equal_the_whole_null`: **2 passed**. Las comprobaciones que comparan con el comparador anterior de p analíticos no son todavía un criterio vigente tras la decisión AUC de WP1; deben actualizarse al continuar esta fase.


## Integración productiva y revisión del 2026-10-23 (vigente)

Este registro sustituye los criterios históricos de comparar p analíticos permulados y de identidad con el ranking crudo. No cierra la aceptación integrada ni la escala real. El estadístico base sigue siendo `gene_stat`; WP9 decide su elección.

### Cableado e interfaces

- `main.nf` transmite a ENRICHMENT las carpetas de CAAS_CORE_MERGE (sin b0); los precomputados pueden usar `caas_null_dir/blocks`. `FCS_CAAS` resuelve estas matrices, transpone una sola vez cada bloque para los tres rankings y calibra con GENE_P. Las columnas observadas reales son `score_global`, `score_top`, `score_bottom` de `fcs_stats.tsv`.
- GENE_P recibe únicamente los trozos de su ranking/lote de genes; cada recuento recibe únicamente los trozos calibrados de su ranking/bloque. Se rechazan ciclos solapados y IDs duplicados. Un bloque de otro `gene_stat` se rechaza antes de calibrar: no puede producirse un ranking calibrado válido ni sustituirse silenciosamente por el crudo. El lector antiguo de validación conserva su comportamiento de p.perm NA ante estampas incompatibles.
- `FCS_RER` conecta tanto el FCS central de ENRICHMENT como el FCS independiente de RER_MAIN. Conserva NA/NaN e infinitos y la evaluabilidad por ciclo: el nulo RER no se limita a los genes con score observado disponible. Los procesos RER interrumpen el flujo si fallan, evitando perder silenciosamente un nulo solicitado.
- RER_CONT/RER_BIN generan y escriben un lote cada vez, sin combinePermData global. La ruta `<trait>.{continuous,binary}.perms.rds` es ahora **un directorio** de `blk_*.rds`, cada uno con `corStat`, identificadores únicos de ciclo, semilla y ajustes de correlación. Los lectores productivos y la GUI aceptan ese directorio o el archivo RDS histórico. No se cambia la ruta de los precomputados.
- El p de gen RER sigue siendo el nativo de `permpvalcor` (mediana de corRho); no se sustituye por el p bilateral del AUC de rutas. Se escriben trozos temporales de corRho por lote de genes, se concatenan todos sus ciclos al calibrar y se eliminan al terminar. El directorio público solo lleva corStat y metadatos. La generación de lotes es serie dentro de la tarea RER, evitando replicar el bosque de árboles entre tareas.
- Se mantienen `caas_perms.rds`, su ensamblado y sus consumidores. Para nulos CAAS precomputados sin bloques, un adaptador lee ese RDS una vez y escribe bloques del estadístico; después aplica la misma calibración. Esta importación sigue cargando el objeto completo. Sin nulo CAAS utilizable, o con la permulación de enriquecimiento desactivada, se emite una tabla vacía y un aviso: no se publica un ranking crudo como sustituto. RER sin nulo conserva el observado analítico con p.perm NA.
- El observado interno entre procesos es RDS; el resultado público sigue siendo TSV. Se reprodujo que el roundtrip TSV cambiaba AUC hasta 2,78e-17 y nueve recuentos marinos (global), con diferencias de p.perm hasta 0,0297. RDS preserva los empates inclusivos sin añadir tolerancias al comparador de rutas.
- Parámetros nuevos: `fcs_block_batch_genes=1000` y `rer_batch_seed=1998`, con configuración, modelos/campos GUI, generación de scripts, ayuda y plantillas sincronizados. Los cambios previos de validation/ no se incluyen en esta entrega; en sus plantillas solo se versionan los dos campos nuevos. Las plantillas locales ignoradas también reciben esos campos, conservando el resto.

### Semillas y reproducibilidad RER

El lote i ejecuta `set.seed(rer_batch_seed + i - 1)`; el contexto continuo (incluida resolución del árbol de simulación) se prepara con la semilla base. Con entradas, versiones de R/paquetes, ajustes y tamaño de lote fijos, el resultado es reproducible y añadir lotes conserva los anteriores. **Cambiar `rer_perms_per_batch` cambia las asignaciones simuladas**: no se exige identidad con la secuencia antigua sin semillas por lote ni con otra partición. Es una adaptación explícita del simulador de RERconverge; la generación nativa de fenotipos y los ajustes armonizados no cambian. No se paraleliza este paso.

Se conservan los presupuestos explícitos de las plantillas y el default de 10 x 100. Para la aceptación a 10.000 se configura 100 x 100 y se guarda la semilla junto con las versiones. No se ha ejecutado ese trabajo real aquí.

### Validación dirigida y evidencia

Base de logs, snapshots estables, trabajo Nextflow, oráculos y medidas: `/home/miguel/IBE-UPF/PhD/PhyloPhere_validation/.tmp/wp6_integration/`. PATH de R/Python: `~/micromamba/envs/phylophere/bin`. No se ha hecho una corrida completa, ni se ha enviado ningún trabajo al cluster en esta revisión.

- Suite versionada final: `python -m pytest -q tests/test_wp6.py tests/test_wp6_cli.py`, **12 passed en 39,69 s** (`pytest_release.log`), incluido el contrato AUC de WP1.
- Oráculo versionado `test_wp6_blocks.R`: 400 genes, 25 columnas, 30 conjuntos; conteo exhaustivo de p de gen y leave-one-out, AUC nativo `fastwilcoxGMTall` y recuentos con bloques {1,7,N}, CAAS unilateral y RER bilateral, NA e infinitos. Calibración física por lotes de genes {1,7,G}, IDs y solapamientos.
- Regresión CLI CAAS: observado RDS, GENE_P, bloques y cierre frente al contrato en memoria; **recuentos enteros idénticos**, p/q a 1e-12. Regresión CLI RER: genes con observado NA que son evaluables en nulos, lateralidad y cierre. Input resolver: ausencia de nulo, RDS histórico, estampa incompatible y ciclos solapados.
- Subflujo real Nextflow sobre marine WP5 b25: **399 filas CAAS**, 100 ciclos, 3 rankings, 17 GMT, 2 lotes de GMT; p.perm frente a memoria completa sobre exactamente los mismos rankings calibrados con delta máximo 1,11e-16, NA y significación idénticos (`compare_nf_caas_rds.log`). Importación del RDS histórico en un bloque frente a los cuatro bloques vivos: todas las columnas stat/p/q a 1e-12 (`compare_compat.log`). Rama CAAS sin nulo: produce tabla vacía sin lanzar transposición/FCS (`run_no_null.log`).
- Subflujo Nextflow RER sobre los **mismos 1.000 ciclos guardados**, 2 bloques de 500, 17 GMT: **133 filas**, stat, pval, p.adj, p.perm, q.perm y sig_perm coinciden con memoria completa a 1e-12 (`compare_nf_rer_final.log`). Un primer ensayo detectó y permitió corregir la restricción errónea del universo nulo a genes evaluables observados; no se clasifica como fallo histórico.
- RER: oráculo nativo con semillas fijas (`test_rer_permulations.R`) aprobado para continuo y binario, clades all/terminal/ancestral, min.sp/min.pos y winsorizaciones no predeterminadas. CLI real sobre árboles/matriz marinos para binario y continuo: observado idéntico a la correlación directa con esos ajustes, dos lotes de dos ciclos, bloques reducidos y p de gen presentes (`rer_cli.log`). Además, los procesos Nextflow reales RER_CONT y RER_BIN generan y publican correctamente los directorios de bloques sobre esas entradas guardadas (`run_producers.log`, dos procesos aprobados); cada tarea limpia únicamente sus salidas nulas propias al reintentar un fallo. El helper de lotes reproduce permpvalcor sobre el nulo completo para lotes de genes {1,2,G}, semillas 1998..2000, NA y P=0.
- GUI: **19 passed** en el control de parámetros/modelo y JSON de plantillas (`gui.log`). WP5 versionado: **19 passed** (en la primera suite combinada hubo 22 PASS y un error de nombre en el test nuevo, después corregido; no se reejecuta WP5 sin cambios). Suite externa GENE_P/CAAS/FCS: cinco aprobadas; la CLI tuvo un fallo al leer código que se estaba editando, y se repitió sobre snapshot estable: **1 passed** (`external_cli_final.log`). Las referencias externas a la copia anterior se actualizan al contrato AUC actual; los ficheros legacy congelados no se modifican.
- `check_fcs_null_equivalence.R` externo ampliado para el ranking calibrado leave-one-out. PEPC guardado tiene **un gen x 1.000 ciclos**: no hay rutas evaluables con los GMT/min_g=10, y el test versionado ahora reconoce ese caso sin intentar agrupar una tabla sin columnas; no se cuenta como equivalencia de rutas positiva.
- Parseo del pipeline real con `nextflow ... main.nf --help` aprobado (`main_help.log`); no es aceptación extremo a extremo. Preparación dirigida de cluster en `tests/cluster/prepare_wp6_snapshot.sh` y copia/README en `PhyloPhere_validation/wp6_check/cluster/`, `bash -n` aprobado. El usuario lo ejecuta dentro de una asignación tras sincronizar el commit; exige temporales bajo `$HOME/scratch`.

### Memoria y E/S: medido frente a leído en el código

Benchmark **sintético local**, una ejecución por punto, 1.000 genes, 30 conjuntos, un ranking, bloque FCS=250. GNU time MaxRSS (KiB) y tiempo de pared (s), logs `bench_*`:

| Etapa | N | MaxRSS KiB | Pared s |
|---|---:|---:|---:|
| Recuento FCS por bloques | 1.000 | 412.856 | 3,78 |
| Recuento FCS por bloques | 4.000 | 419.736 | 4,09 |
| Recuento FCS por bloques | 10.000 | 423.192 | 7,43 |
| GENE_P, 1.000 genes | 1.000 | 129.688 | 0,67 |
| GENE_P, 1.000 genes | 4.000 | 190.100 | 1,22 |
| GENE_P, 1.000 genes | 10.000 | 383.372 | 4,64 |

Estos puntos no miden la transposición completa ni el bosque RER real, ni demuestran memoria constante del pipeline. Por código, el FCS retiene genes x B y matrices de rutas del lote de GMT; GENE_P retiene G x N con copias y resultados. La calibración del p de gen RER tiene una retención análoga para hasta 1.000 genes, además del bosque y la matriz de RER. La reducción CAAS que produce el RDS único sigue creciendo con N hasta WP7.

Por código (sin medir ancho de banda): transposición CAAS lee una vez el RDS de cada bloque y escribe trozos por ranking/gen; GENE_P lee esos trozos una vez por lote, escribe p por bloque, y cada lote GMT vuelve a leer los p de su bloque. RER escribe corStat público y corRho por lote de genes, y lee cada trozo de corRho al calibrar el p de gen. La validación de directorios RER hace una pasada adicional de los bloques. La adaptación de un RDS histórico añade una lectura completa y escritura de sus bloques. Se conserva el almacén WP5 y el RDS único, de modo que el espacio total no se reduce todavía como en WP7.

### Registro de entrega y pendientes

- Fechas: 2026-10-22/23. HEAD inicial `50eb0eb`; durante el trabajo el checkout recibió `6ae3a4f` (corrección WP4 de otra revisión). Sus cambios de documentación y los cambios previos de validation/ y prompts se conservan fuera del commit WP6. Cierre temático: commit de integración WP6 en el historial; no acredita aceptación integrada.
- Pendientes: corrida integrada neoplasia/maligno y strict, calibración amplia/WP9 y 10.000 reales, perfil de transposición/GENE_P/FCS/RER y E/S en cluster, aceptación conjunta de los WP. No se repiten aquí los siete fallos de cableado general ya registrados por WP5 ni se consideran verdes.
- Entrada a WP7: FCS productivo ya consume matrices por bloque y el importador mantiene los RDS históricos; se puede retirar el ensamblado CAAS único solo en WP7, auditando los consumidores restantes (SCORING/informes/GUI).
