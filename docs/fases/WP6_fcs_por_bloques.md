# WP6. Fase 4: FCS por bloques de columnas y etapa `GENE_P`

Estado: IMPLEMENTACIÓN PARCIAL; avance temático sin cierre. Depende de: WP1 (clasificación nueva, RER paramétrico) y WP5 (matrices por bloque).
Conversación independiente. Antes de empezar, leer `docs/fases/README.md`.

## Objetivo

Que el FCS (CAAS y RER) use memoria constante en N procesando bloques de columnas del nulo y acumulando contadores, con el mismo resultado que la versión actual.

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
