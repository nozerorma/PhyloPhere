# WP5. Fase 3: tablas por (gen, ciclo), merge por bloque, SCORING y acumulación

Estado: implementación local revisada (2026-10-22); validación a escala real pendiente. Estado vigente de esta revisión al final de la página. Los resultados del registro original se conservan como historial. Depende de: WP4. Serial con WP4 (formato de tesela).
Conversación independiente. Antes de empezar, leer `docs/fases/README.md`.

## Objetivo

Sustituir las ~80 millones de filas de texto del nulo por tablas compactas por (gen, ciclo) y un merge por bloque de ciclos, sin cambiar ningún número, y adaptar los consumidores (SCORING, acumulación, informe 11).

## Decisiones aplicables (del usuario)

- Eje ASR fuera (3 de las 6 matrices, columnas del cuantil 90): sin lectores en el repositorio (solo copias en `~/.cc/r11.R` y un `.history`); corregir el comentario obsoleto de `fcs_enrich.R:661-663`.
- Todas las variantes del estadístico de gen (máximo crudo, `F(max)`, `F(max)^n`) las calcula el merge; la elección es un parámetro (`gene_stat`) y la estampa de `fcs_enrich.R` y `scoring_compute.R` debe seguirla (hoy `size_adj_max`). La decisión de cuál usar la toma WP9.
- POSENRICH desactivado (WP2): ya no se necesita el fichero por posición para él. La acumulación pasa a un indicador por posición.
- Recalibrar los recursos del merge y de las teselas desde la traza de 1000 ciclos.

## Estado de partida

El detalle es gz de texto por gen; el paso B (`_finalize_perm_scores`, `gene_wrapper.py:753`) lo lee dos veces (pools por ciclo y genes) y reescribe `perm_pos_cycle_caas.tsv.gz` con todas las filas (80,9 millones a 1000 ciclos en cáncer; el 28 % con score de US). El merge actual: 32,5 min y 4,1 GB (neoplasia y maligno completo). `scoring_caas_perms.R` lee entera la tabla ancha: ~3,4 MB por ciclo (~34 GB a 10.000, extrapolado). Consumidores del fichero por posición: `scoring_compute.R:431` (`p.emp`, `p.emp_fact`), `accumulation_test.py`, informe 11 (lee `perm_pos_sample.tsv` y `perm_pos_quantiles.tsv`), POSENRICH y su informe 14.

## Hechos del código que fijan el diseño (verificados)

- `core.postproc.gene_removal` agrupa por ciclo (`pools` por `u.labeling`): los umbrales (cuartiles del recuento por gen, percentil de la densidad por longitud) se calculan dentro de cada ciclo. Por eso el merge es por bloque de ciclos (todos los genes del ciclo juntos) y no por lote de genes. Como la decisión de eliminar una unidad (ciclo, gen) necesita todos los genes del ciclo, la tarea **no puede emitir un pool agregado**: emite recuentos y scores por (gen, ciclo) y el pool se forma en el merge, sin las unidades eliminadas.
- El pool descarta además filas en cluster una a una (`_is_clustered`), mientras que el filtro de genes usa una bandera de cluster por (gen, ciclo): la tarea emite ambas cosas (los recuentos incluyen filas sin score de US).
- `pool_labelings(strict=True)` descarta hipótesis antes de agrupar por ciclo base: la tesela lleva una columna `conserved` por hipótesis o un segundo juego de tablas, y se mantienen `CAAS_CORE_MERGE_STRICT`, `CAAS_CORE_OBSERVED_STRICT` y `b0_observed_strict`.
- `p.emp_fact` (`scoring_compute.R:129`) se reduce de forma exacta en dos etapas: `p = min(1, (nd + 1)/(N + 1) × S)`, `S = (1 + nº de detecciones de la clase con score > obs − TIE_TOL)/(1 + tamaño de la clase)`; las clases dependen solo del histograma de `nd` por posición (suma de contadores); basta, por posición observada, un recuento por clase aditivo entre bloques. La curva de saturación sale del mínimo entre bloques del primer ciclo en que aparece cada posición. El filtro de universo es estático por fila.
- `perm_pos_sample.tsv` y `perm_pos_quantiles.tsv` usan un RNG único en orden de genes (reservoir): no pueden ser idénticos bit a bit; sembrar por ciclo y comparar estadísticamente.

## Alcance (ficheros)

`subworkflows/CT/caas_permulation.nf` (`CAAS_CORE_BATCHED`, `CAAS_CORE_MERGE` ~428 a 508, `CAAS_CORE_MERGE_STRICT`, `CAAS_CORE_OBSERVED(_STRICT)`), `subworkflows/CT_DISAMBIGUATION/local/reaggregate_perm_scores.py` y `src/utils/gene_wrapper.py` (`_cycle_gene_removal_from_detail`, `_build_cycle_score_pools`, `_finalize_perm_scores`), `src/core/scores.py`, `src/core/postproc.py`, `subworkflows/SCORING/local/src/scoring_caas_perms.R`, `scoring_compute.R`, `subworkflows/SCORING/local/11.Scoring_report.Rmd`, `subworkflows/CT_ACCUMULATION/local/accumulation_test.py`, `workflows/ct_accumulation.nf`, `main.nf` (canales `scoring_caas_*_ch` ~559 a 563), `conf/resources.config` (~227 a 245), `fcs_enrich.R:661-663`.

## Cambios

1. **Esquema de tesela** (f64, sin rebajar precisión): `units` (gen, ciclo, posiciones detectadas distintas, bandera de cluster, nº puntuadas por dirección) y `scores` (gen, ciclo, posición, lado, score, `conserved`; ya filtrado por `remove_clusters`). Formato binario columnar (decidir en la conversación: Arrow/Parquet u otro; debe poder leerse desde R y Python).
2. **`CAAS_CORE_MERGE` por bloque** (dispersión): une los lotes de genes del bloque, aplica `gene_removal` por ciclo, forma los pools exactos por (ciclo, dirección), calcula las tres variantes por (gen, ciclo) y escribe la matriz genes × B de los tres rankings, `removed_units`, filas de roster y cuantiles por ciclo (RNG por ciclo). `CAAS_CORE_MERGE_STRICT` ejecuta el mismo código sobre el juego estricto. Reducción final pequeña: roster y universo.
3. **SCORING:** `p.emp` con contadores por bloque (`nd`, recuento con score ≥ obs − `TIE_TOL`, primer ciclo), sumados en una reducción; `p.emp_fact` en dos etapas; los scores observados salen de la tesela `b_0`; `p.adj_bh` sobre la familia U (posiciones del universo).
4. **Acumulación:** indicador por posición con score positivo en algún ciclo.
5. Dejar de construir el eje ASR y de escribir `perm_pos_cycle_caas.tsv.gz` (no hay consumidores tras WP2).
6. Recalibrar `conf/resources.config` (merge del orden de 2 GB + k·B, con k medido, en escalera de intentos; tesela reajustada; ~12 GB para `GENE_P` en WP6).

## Interfaz hacia WP6 y WP7

Por bloque: `blk_k.matrices.rds` (o equivalente) con las matrices genes × B de los tres rankings del estadístico elegido, y los datos por (gen, ciclo) necesarios para derivar las tres variantes. WP6 (`GENE_P`) lee de ahí las filas de un lote de genes de todos los bloques.

## Validación (criterios de paso)

- Sobre marine toy, PEPC y cancer_b0_toy, por ciclo: `global_caas`, `top_caas`, `bottom_caas` iguales a las nuevas matrices (máx. |dif| ≤ 1e-12, mismos conjuntos de claves); `removed_units`, `position_scores.tsv`, `p.emp` y `p.emp_fact` (1e-12), roster y `b0/` iguales a la referencia; contadores enteros exactamente iguales; `strict/` con los mismos controles.
- A escala real, a 1000 ciclos en neoplasia y maligno completo frente a los directorios de referencia del cluster (por `sbatch`), incluido `strict/` y las tablas de SCORING.
- Portar las pruebas que asumen `gene_cycle_scores.tsv` (`test_null_consumers`, `test_scores`, `test_compare_null`, `test_null_roster`, `test_detail_only`, `test_union_shards`).
- Memoria por tarea independiente de N (medir a 1000, 4000 y 10.000 con las trazas de Nextflow).

## Fuera de alcance

FCS por bloques (WP6), formato del RDS único (WP7), POSENRICH.

## Riesgos y puntos a verificar

- Roster, denominador N y familia de BH reconstruidos desde las sumas con los mismos números que hoy.
- Orden de sumas y empates (`TIE_TOL`): mantener `math.fsum` y las comparaciones por `bisect_right(pool, max + TIE_TOL)`.
- El esquema equivocado descubierto tarde: probar en toy, luego a 1000 ciclos en el cluster y solo después a 10.000.

## Diseño (lo que produce esta fase)

Formato: Parquet (pyarrow, zstd, f64). `pyarrow` está en `environment/phylophere.yml` (25.0.0 en el entorno local, 23.0.1 en el lock); R no tiene `arrow` ni `nanoparquet`, ni en el entorno local ni en la especificación del contenedor, así que solo procesos Python leen tablas de tesela y de bloque, y lo que lee R es pequeño (matrices en RDS, contadores en TSV). Decidido con el usuario: no se añade `r-arrow`.

**Tesela** (`null_tables/`, y `null_tables_strict/` con el pooling estricto; un par de ficheros por lote de genes y bloque; `src/utils/null_tables.py`):

- `units.parquet`: una fila por (gen, ciclo) con detección: `gene, cycle, n_pos` (posiciones detectadas distintas, todas las filas: el recuento de `gene_removal`), `clust` (alguna posición en tren), `n_rows` (registros que quedan tras quitar trenes), `n_top, n_bottom, n_all` (posiciones con score > 0 por dirección o en cualquiera) y `max_top, max_bottom, max_all` (NaN si n = 0).
- `scores.parquet`: una fila por (gen, ciclo, posición) con algún score > 0, tras quitar trenes: `gene, cycle, position, top, bottom` (NaN si ese lado no puntúa; `all = max(top, bottom)`).
- Las filas sin score (47 % en marine, 58 % en el toy de cáncer, medido) no se guardan: solo su recuento (`n_rows - n_top - n_bottom`).
- `meta.json`: `by_labeling` (el grano de los ciclos: ciclo base con FOP, etiqueta sin él).

**Merge por bloque** (`CAAS_CORE_MERGE`, una tarea por bloque y una para `b0`; `merge_null_block.py`, `src/utils/null_merge.py`): `gene_removal` por ciclo con `core.postproc`, pools exactos por (ciclo, dirección) con las unidades no eliminadas, y por (gen, ciclo, dirección) `n`, `max`, `fmax = bisect_right(pool, max + TIE_TOL) / len(pool)` y `size_adj_max = fmax ** n` (potencia en Python puro, la misma expresión de `core.scores`). Salida en `caas_permulation/blocks/<bloque>/` (y `caas_permulation/b0/`):

- `matrices.rds`: `list(genes, cycles, gene_stat, n_units, stat = list(global, top, bottom), parts = list(max = list(...), fmax = list(...), n = list(...)))`. Matrices genes x ciclos del bloque, 0 donde la unidad no tiene valor. `stat` es el estadístico elegido con `gene_stat`; `parts` conserva lo que hace falta para cualquier variante (también `F(max)^(n^β)` de WP9). **Esta es la interfaz hacia WP6 y WP7.**
- `positions.parquet` (`gene, cycle, position, all`: las posiciones puntuadas de las unidades no eliminadas, sin trenes si se quitan), `detections.parquet` (por posición: `nd` y el primer ciclo en orden de reproducción), `removed_units.tsv`, `cycles.txt` (el roster del bloque), `perm_pos_sample.tsv`, `perm_pos_quantiles.tsv` (RNG sembrado con la semilla y el ciclo), `cycle_summary.tsv` y `cycle_hist.tsv` (por ciclo y dirección: entradas con score, genes con una, histograma en bins de 0,01; los lee el informe 11), `counts.json` (unidades, eliminadas, posiciones, ciclos con fila, `score_rule`).

**Reducción** (`CAAS_CORE_MERGE_REDUCE`; `reduce_null_blocks.py`, `scoring_caas_perms.R`): `cycle_roster.txt` (N de todos los consumidores), `removed_units.tsv`, `perm_pos_sample.tsv`, `perm_pos_quantiles.tsv`, `cycle_summary.tsv`, `cycle_hist.tsv`, `counts.json`, `null_positions.parquet` (por posición que puntúa el nulo: `nd`, `first_rank`) y `caas_perms.rds` (`list(caas_corStat_byrank = list(global, top, bottom), gene_stat)`, las matrices de todos los bloques con las columnas `sort(unique(roster))` de siempre; lo sustituye WP7). Con `--caas_full_perms 0` solo se fusiona `b0/` y no hay reducción.

**Contadores de SCORING** (`subworkflows/SCORING/scoring_null_positions.nf`, `null_position_stats.py`; cada contador es entero y suma entre bloques):

- `CAAS_NULL_UNIVERSE` (una vez): `nd` (suma de bloques) -> universo U (conjunto testado, genes del fondo, regla de tasa >= 0,5; más toda posición observada) -> clases de propensión de `p.emp_fact` (percentiles de detecciones: peso `v x nº de posiciones con nd = v`) -> `null_classes.parquet`, `null_fit.tsv`, `null_universe_base.tsv`.
- `CAAS_NULL_COUNT` (por bloque): `k_emp` (ciclos con `all >= obs - TIE_TOL`) y `c_k` (detecciones de las posiciones de la clase k con score >= obs - TIE_TOL) por posición observada.
- `CAAS_NULL_COUNTERS` (una vez): suma -> `null_counters.tsv`. `null_position_p.R` calcula `p.emp`, `p.emp_fact` y las familias BH (las posiciones solo del nulo entran con p = 1, solo por el tamaño de la familia) con las mismas fórmulas y guardas de antes; `scoring_compute.R` ya no lee ninguna tabla por posición.

**Parámetros nuevos:** `gene_stat` (`size_adj_max` por defecto, `fmax`, `max`; campo en la pestaña CAAS de la GUI y en las 14 plantillas de `gui/templates`; sigue el merge, `observed_core_scores.py` y la estampa de `fcs_enrich.R`) y `caas_null_dir` (directorio `caas_permulation/` de una corrida anterior, sustituye a `caas_pos_cycle_caas_file`, `caas_pos_sample_file`, `caas_pos_quantiles_file` y `caas_cycle_roster_file`). **Desaparecen:** el eje ASR, `perm_pos_cycle_caas.tsv.gz`, `gene_cycle_scores.tsv`, `union_shards.sh`, `merge_block_shards.py`, `reaggregate_perm_scores.py` y `n_detected`; POSENRICH deja de ejecutarse de hecho (su código no se ha tocado, ver `docs/POSENRICH_PORTE_EVALUACION.md`).

## Registro de la fase

```
- Conversación y fechas: WP5, 2026-10-22, conversación independiente, local, rama remove_caap.
- Commit de partida / commit de cierre: 90c06a4 (cierre de WP4, hecho por Claude con la autorización del usuario: código y docs/, sin validation/) / 0834081 (cierre de la implementación inicial; la revisión actual añade cambios locales sin commit).
- Qué se hizo (ficheros y cambios): ver "Diseño" arriba. Nuevos: src/utils/null_tables.py, src/utils/null_merge.py, merge_null_block.py, reduce_null_blocks.py (CT_DISAMBIGUATION/local), SCORING/local/src/{null_position_stats.py, null_position_p.R, assemble_block_matrices.R}, SCORING/scoring_null_positions.nf. Reescritos: scoring_caas_perms.R, gene_wrapper.py (réplica por gen -> tablas), disambiguation_perms_main.py, caas_permulation.nf (CAAS_CORE_MERGE por bloque + CAAS_CORE_MERGE_REDUCE), main.nf, workflows/{scoring,enrichment,ct_accumulation}.nf, scoring_compute.R/.nf, scoring_report.nf, 11.Scoring_report.Rmd (lee cycle_summary/cycle_hist), accumulation_test.py (lee null_positions.parquet), core/scores.py (gene_stat, cycle_order_key), observed_core_scores.py (--gene-stat), fcs_enrich.R/fcs_compute.R/fcs.nf/12.FCS_general_report.Rmd (estampa sigue gene_stat; comentario del eje ASR corregido), conf/{scoring,resources}.config, help.nf, README, GUI (campo gene_stat, caas_null_dir en Precomputed Run, report_registry, run_single/sbatch_array, 14 plantillas de gui/templates). Eliminados: reaggregate_perm_scores.py, union_shards.sh, merge_block_shards.py, el eje ASR, perm_pos_cycle_caas.tsv.gz, gene_cycle_scores.tsv, n_detected y los parámetros caas_pos_cycle_caas_file, caas_pos_sample_file, caas_pos_quantiles_file, caas_cycle_roster_file. POSENRICH: bloque retirado de enrichment.nf con aviso (su código no se tocó).
- Comprobaciones ejecutadas (comando, conjunto de datos, tolerancia, resultado):
  * Paso 0 (marine toy de WP3 y WP4): ver el registro de WP4. Referencia a92fd2d completada (17 min), el commit 90c06a4 con bloques 25 y 40 igual a ella en el lado CAAS, con delta 0.0; golden_marine_wp4/.
  * Pruebas nuevas contra una definición directa (oráculo escrito en la prueba) o contra el código original de 581b999 (unification/legacy_scoring/, sin cambios): test_null_tables (matrices bit a bit con 1 y 3 teselas, sin y con trenes, posiciones, muestra, resumen), test_null_position_stats (p.emp, p.emp_fact, BH, universo y curva de saturación idénticos al código viejo sobre una tabla en formato viejo, con 1, 3 y 6 bloques, con y sin conjunto testado, y las guardas), test_null_blocks_cli, test_gene_stat, test_accumulation_universe, test_compare_null (modos nuevos), test_tile_driver y las reescritas test_chunk_failures, test_dotted_gene_keys, test_postproc_wiring, test_train_columns, test_factorized_p_calibration. Mutación (quitar TIE_TOL de la comparación del pool): la prueba de matrices falla.
  * compare_null.py --tol 1e-12 contra la referencia vieja (formato tablas) de cada conjunto: PEPC c4 a 1000 ciclos con bloque 250: 8 passed, 0 failed (caas_perms.rds eje CAAS all.equal TRUE, posiciones por (gen, posición, ciclo) delta 0.0, b0, roster, scoring/position_scores.tsv con p.emp y p.emp_fact delta 0.0, fcs_stats); strict 6 passed, 0 failed. Toy de cáncer (bloques 25 y 40): matrices, posiciones, removed_units, roster, b0, position_scores, fcs_stats y leading edge de CAAS idénticos; muestra y cuantiles estadísticamente iguales (filas por ciclo idénticas, |z| máx 0,52); los 2 fallos son las tablas del FCS con esquema anterior a WP1; strict 7 passed. Marine toy (bloque 25) contra la referencia vieja y contra golden_marine_wp4: lado CAAS idéntico; difieren solo el lado RER (WP1 y permulación de RER no reproducible, ver abajo) y el esquema del FCS; strict 6 passed.
  * Bloque 25 contra bloque 40 (toy de cáncer): 14 passed, 0 failed, incluidas muestra y cuantiles (el RNG es por ciclo: no dependen del bloque); strict 8 passed.
  * Suites: unification y resources con RUN_SLOW=1 --continue-on-collection-errors: 235 failed, 354 passed, 17 errors (WP4: 250, 314, 17); ningún fallo fuera de la lista de la línea base. wiring sin Nextflow (GUI, auditoría de parámetros, nombres retirados): 55 failed, 110 passed, los mismos 55 de la lista. NO ejecutadas tras el cierre: las pruebas de wiring con Nextflow (test_core_batched, test_stream_vs_export, test_wiring, test_subset_b0), reescritas pero sin correr; la corrida real de la tubería en los tres conjuntos sí las ejercita.
- Mediciones (recurso antes y después):
  * Toy (medido en trazas de Nextflow): CAAS_CORE_MERGE de 4 a 192 MB y 1 a 2 s por bloque; CAAS_NULL_* de 4 a 45 MB; teselas de 0,26 a 2,4 GB (como en WP4, no cambian).
  * Escala de cáncer, SINTÉTICA (teselas generadas con las densidades del nulo de 1000 ciclos de cáncer: 11.800 genes, 9000 unidades y 22.400 posiciones puntuadas por ciclo, bloque de 500 ciclos; wp5_check/scale_merge.py y scale_bench.sh; una ejecución): merge_null_block.py 43 s y 2,9 GB de RSS; ensamblado del RDS del bloque 13 s y 0,6 GB; matrices en bruto 473 MB, matrices.rds 198 MB, positions.parquet 41 MB y las 590 teselas 161 MB. Reducción (reduce_null_blocks.py + scoring_caas_perms.R, bloques copia del primero con ciclos propios): 2 bloques (1000 ciclos) 11 + 12 s, 2,4 y 1,3 GB; 8 bloques (4000) 23 + 44 s, 2,7 y 4,9 GB; 20 bloques (10.000) 49 + 160 s, 2,7 y 7,6 GB. El merge de un bloque no depende de N (el tamaño del bloque lo fija); la reducción final en R sí crece con N porque caas_perms.rds sigue siendo un objeto único (lo sustituye WP7). Antes: merge de 32,5 min y 4,1 GB (neoplasia y maligno completo, 1000 ciclos), y scoring_caas_perms.R leía la tabla ancha entera.
  * No medido: los contadores de SCORING a escala (solo toy, 4 a 45 MB); RP1 (RSS de la tesela del gen mayor con el flujo de WP4); nada a 4000 o 10.000 ciclos con datos reales.
- Desviaciones del plan y decisiones tomadas (con su justificación):
  * Lado en columnas (top, bottom) en vez de filas por lado, y strict como un segundo par de tablas en vez de una columna conserved (el pooling estricto cambia el score, no solo filtra filas). Confirmadas con el usuario.
  * Formato Parquet solo para Python; R lee RDS y TSV pequeños (R no tiene arrow; no se añade r-arrow). Confirmado con el usuario.
  * perm_pos_sample.tsv: las filas sin score muestreadas llevan Position NA (no se guardan). Sin lector en el repositorio. Confirmado.
  * caas_perms.rds sigue existiendo como un objeto único (3 matrices del estadístico elegido), porque FCS (WP6) y el RDS único (WP7) no están hechos.
  * Sin permuladas (--caas_full_perms 0) no hay reducción ni caas_perms.rds (antes salía un RDS vacío): los consumidores toman el centinela y p.perm queda NA; CT_ACCUMULATION con permulación sigue parando con su error.
  * Se corrige un límite latente del código viejo: en R, K * cumsum(n) de las clases de p.emp_fact es entero y se desborda (NA) con más de 107 millones de detecciones (K = 20); el código nuevo usa enteros de 64 bits.
  * La estampa de la regla de score (SCORE_RULE) pasa de la columna score_aggregation de la tabla por posición a counts.json (score_rule); null_position_stats.py rechaza un nulo con otra regla.
  * remove_intermediates ya no actúa sobre nada (las teselas no se duplican); su destino lo decide WP8.
  * POSENRICH: el nulo por posición ya no existe como tabla, así que el módulo deja de ejecutarse con aviso; su código queda para el porte (docs/POSENRICH_PORTE_EVALUACION.md).
  * Las corridas de validación usan una copia del árbol (.tmp/code_wp5_*), porque las tareas copian el código al arrancar: editar el repositorio durante una corrida la rompe (ocurrió una vez).
  * validation/tier*/templates no se han tocado (las gestiona el usuario): sin gene_stat toman size_adj_max.
- Pendiente que pasa a otra fase:
  * Validación a 1000 ciclos en neoplasia y maligno completo en el cluster (y RP1 de WP4): scripts en PhyloPhere_validation/wp5_check/cluster/ (README.md, prepare_wp5_cluster.py); el usuario sincroniza y lanza.
  * Hallazgo (no de WP5): la permulación de RER no es reproducible entre dos corridas del mismo código (p.perm de RER difiere hasta 0,07 y q.perm hasta 0,15 en marine toy, bloques 25 y 40); stat y pval son idénticos.
  * test_scoring_compute y otras pruebas de la línea base siguen rotas por el código con CAAP (no se arreglan aquí).
```

## Qué debe comprobar primero WP6 (y WP7)

1. **Interfaz de entrada:** `caas_permulation/blocks/<bloque>/matrices.rds` (estructura arriba). Las filas de cada bloque son el universo (cleaned_background_main.txt) más los genes con alguna unidad en el bloque: WP6 debe alinear por nombre de gen y tomar 0 donde falta. Cada bloque lleva su estampa `gene_stat`; el control de `fcs_enrich.R` (nulo construido con otro estadístico: p.perm NA) debe conservarse al leer por bloques. Orden de columnas: `sort(unique(roster))` de R en `caas_perms.rds`; los bloques traen sus ciclos y el orden numérico es el de réplica.
2. **Lo que WP5 no hizo:** `caas_perms.rds` se sigue construyendo entero en `CAAS_CORE_MERGE_REDUCE` (7,6 GB en R con 20 bloques, sintético): WP7 lo sustituye, y mientras tanto la petición de memoria de ese proceso crece con N (`conf/resources.config`). Las variantes `max`, `fmax` y `n` están en `parts` de cada bloque; `F(max)^(n^β)` de WP9 se construye de `fmax` y `n` (el exponente en R no es bit a bit el de Python: solo la variante elegida lo es).
3. **Validación pendiente:** no hay medidas con datos reales a 1000, 4000 ni 10.000 ciclos de los procesos nuevos; la corrida del cluster (scripts preparados) da las trazas, y con ellas hay que revisar las peticiones de `conf/resources.config`, que hoy salen de datos sintéticos.
4. **Trampas:** la permulación de RER no es reproducible entre corridas (comparar solo stat y pval); corridas de validación desde una copia del árbol; `compare_null.py` acepta ahora ambos formatos (tablas y bloques).

## Revisión de WP5 (2026-10-22, posterior a la revisión de WP1)

### Fallo reproducido y corrección

El reductor aceptaba ciclos repetidos dentro de un bloque o compartidos por dos bloques. Al pasar dos veces un bloque de dos ciclos, `cycle_roster.txt` contenía dos ciclos pero `counts.json` sumaba cuatro; la tabla de detecciones también se sumaba dos veces. El ensamblador R aceptaba esas mismas columnas y las sobrescribía sin error. Esto se reprodujo ejecutando ambos scripts, no se detectó únicamente por inspección. No se ha demostrado que las corridas guardadas contengan ese solapamiento.

- `reduce_null_blocks.py` rechaza ciclos repetidos, solapamiento entre bloques y discrepancias entre `counts.json` y el roster de cada bloque, antes de publicar el roster final.
- `scoring_caas_perms.R` rechaza ciclos repetidos dentro o entre matrices, también si el ciclo no tiene detecciones. Los bloques válidos conservan sus resultados.
- Se añadieron seis casos en `PhyloPhere_validation/unification/test_null_blocks_cli.py`: seis fallos antes de la corrección; todos incluidos en la regresión posterior.

### Evidencia nueva

- Entorno: `micromamba/envs/phylophere`; `TMPDIR` y `--basetemp` bajo `PhyloPhere_validation/.tmp/wp5_review_20261022/`. Logs y JSON de comparación en ese directorio.
- Suite dirigida (pytest, siete archivos: `test_null_tables`, `test_null_position_stats`, `test_null_blocks_cli`, `test_gene_stat`, `test_accumulation_universe`, `test_compare_null`, `test_tile_driver`): **62 passed antes** y **68 passed después** de añadir las guardas y los seis controles. No sustituye la suite completa de unification/resources.
- Comparaciones repetidas de las corridas guardadas, tolerancia 1e-12: PEPC completo 8 PASS, 0 FAIL; estricto 6/0. Cáncer completo 11 PASS y 2 FAIL en tablas FCS; estricto 7/0. Marino completo 10 PASS y 6 FAIL en tablas FCS/RER; estricto 6/0. Cáncer bloques 25 frente a 40: completo 14/0, estricto 8/0. Los fallos FCS/RER se registran como tales; no se cuentan como equivalencia global. Los controles de matrices, posiciones, roster y SCORING de WP5 pasan. Son salidas guardadas de la implementación anterior, no una nueva ejecución completa del pipeline corregido.
- Reconstrucción con el código corregido a partir de bloques guardados: PEPC (1000 ciclos), cáncer y marino (100 cada uno), completo y strict. Las seis reducciones y ensamblados reproducen `caas_perms.rds` (`all.equal TRUE`, tolerancia 1e-12), roster, `null_positions.parquet` y `counts.json`: 24 controles aprobados. No se regeneró el replay de las teselas.
- Cableado real con Nextflow (`test_core_batched`, `test_stream_vs_export`, `test_subset_b0`, `test_wiring`): **38 passed, 7 failed**. Los 14 controles de core/replay pasan. Los siete casos fallidos ya aparecen entre los nueve del log previo `.tmp/wp4_tree_wiring.log`; no se declara verde la suite general. Fallan aserciones sobre la ruta de acumulación, un informe ASR no activado en el caso, una ruta de discovery con la configuración actual, exports retirados y el parser estático de `publishDir`. Las comprobaciones funcionales de publicación sí pasan. Estos tests generales siguen pendientes de adaptar; no se cambió código para satisfacer contratos retirados.
- Preparación de cluster: el generador actual renderiza correctamente neoplasia y maligno; `bash -n` pasa en los cuatro scripts. No se ha lanzado ningún trabajo ni sincronizado código al cluster.

### Pendiente vigente

- Corridas reales de 1000 ciclos en neoplasia y maligno completo, con strict y la misma caché ASR que las referencias; comparación de los resultados y trazas de recursos.
- Medición real de las teselas (incluido RP1), merge y contadores; los números de 4000/10.000 del registro siguen siendo sintéticos. Las medidas reales a esas escalas siguen pendientes en el plan de escalado.
- Recalibrar recursos a partir de esas trazas y conservar explícitamente la reducción R que crece con N hasta WP7.
- Actualizar la memoria de aceptación cuando se satisfagan los criterios integrados. La entrega temática de esta revisión queda autorizada por el usuario y no acredita la validación de escala.

Casos pendientes de `test_wiring.py` en esta revisión:

- `test_the_live_null_replays_and_disambiguates_in_one_process_family`
- `test_the_live_observed_labeling_is_the_b0_slice_of_the_core_and_feeds_everything_downstream`
- `test_a_discovery_tab_that_already_exists_is_scored_by_caas_observed_and_the_core_writes_no_observed_files`
- `test_perm_discovery_exports_that_already_exist_go_to_the_same_process`
- `test_perm_discovery_exports_made_without_the_real_labeling_are_rejected[b_1~H1\ta\tb\n]`
- `test_perm_discovery_exports_made_without_the_real_labeling_are_rejected[b_1\ta\tb\n]`
- `test_the_observed_processes_publish_their_output_directories_by_name`

### Entrega temática y aceptación integrada

Por decisión del usuario (2026-10-22), se entregan las guardas de ciclos y las pruebas locales en un commit de WP5. La corrida completa de neoplasia/maligno, la comparación a escala y los recursos se integran en una validación conjunta tras los WP; no se repite una corrida completa para cada fase. Los siete tests generales de cableado registrados siguen pendientes de adaptar. El commit no acredita validación de escala. Las seis regresiones de ciclos se incluyen también en `tests/test_wp5_null_blocks_cli.py`, dentro del repositorio, para acompañar el cambio versionado.

Prueba versionada antes del commit: `python -m pytest -q tests/test_wp5_null_blocks_cli.py` con TMPDIR en IBE: **19 passed**. Esta prueba se ejecuta sin el checkout externo de validación.
