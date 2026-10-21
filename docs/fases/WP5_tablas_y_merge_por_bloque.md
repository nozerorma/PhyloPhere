# WP5. Fase 3: tablas por (gen, ciclo), merge por bloque, SCORING y acumulación

Estado: pendiente. Depende de: WP4. Serial con WP4 (formato de tesela).
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
