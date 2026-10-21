# Plan: significancia del FCS, nulo de CAAS a 10.000 ciclos y limpieza

## Contexto

Tres líneas de trabajo ya discutidas, recogidas en `docs/FCS_REVISION_RERCONVERGE.md` (puntos P1 a P15) y `docs/PLAN_NULO_10K_PERMULACIONES.md` (diagnóstico, diseño y mediciones):

1. **Capa de significación del FCS**: regla basada en el p de permulación, sin Lachenbruch, con RER solo global y firmado.
2. **Nulo de CAAS a 10.000 ciclos con RAM pequeña por tarea**: bloques de ciclos, replay en flujo, estadísticos por (gen, ciclo) en lugar de filas, FCS por bloques. Medido sin fragmentar a 10.000 ciclos: FCS 32 GB y 58 min; `scoring_caas_perms.R` ~34 GB (extrapolado); RESAMPLE 4,4 GB y 11 min; etiquetados ~3,3 GB por proceso de replay.
3. **Limpieza**: separar `remove_intermediates` y `remove_work`.

Regla de trabajo: **ningún número cambia sin comparación con lo actual**. El árbol de trabajo está sucio y el usuario hace los commits: los valores de referencia ("golden") se atan a un commit suyo. Cómputo en cluster solo por `sbatch`; nada se edita en el cluster; sin `/tmp` (usar `~/.cc`); comentarios de código solo descriptivos del comportamiento actual; todo parámetro nuevo con campo en la GUI y plantillas sincronizadas (`gui/templates` y `validation/tier*/templates`).

## Decisiones del usuario incorporadas

- **FCS (sin cambios respecto a lo decidido):** `fcs_min_genes` 10; `fcs_max_genes` sin límite (0); voto bilateral en rankings firmados; criterio principal = BH sobre el p de permulación por ranking y base de datos, q < 0,15; Lachenbruch fuera del FCS; el FCS se presenta como resultado.
- **Segundo umbral:** `fcs_pperm_thr` 0,025 se conserva como segundo umbral. Es lo que hace Redlich et al. 2024 (suplemento, Methods): p ajustado < 0,15 y p de permulación < 0,025, justificado porque el p de permulación corrige la no independencia de los rangos de genes y de la filogenia; ellos reconocen que parte del refinamiento viene de ese filtro. Fig. S5 usa 0,10 y 0,05; Fig. S7 usa 0,15 y 0,025; 10.000 permulaciones.
- **RER:** solo ranking global y firmado; estadístico paramétrico `sign(Rho) · -log10(P)` en observado y nulo (como `getStat`); permulaciones a 10.000.
- **Gen de CAAS:** el ranking del FCS es −log10 del p de permulación del gen (decidido). Está **por decidir** el estadístico base que construye ese p: máximo crudo, `F(max)` o `F(max)^n`. No se confirma `F(max)^n`: hay que verificarlo y medir el peso de `^n`.
- **Eje ASR:** dejar de construirlo y corregir el comentario obsoleto; confirmado que no hay lectores en el repositorio (solo copias en `~/.cc/r11.R` y un `.history`).
- **Escalado:** eliminar la reutilización de exports; implementar fases 1, 2, 3 y 4; bloque de ciclos de tamaño fijo; recalibrar el merge; todo validado en profundidad contra lo actual.
- **`caas_perms.rds`:** bloques canónicos + exportación opcional del RDS único.
- **POSENRICH:** se desactiva con aviso, y se redacta un documento sobre cómo y si portarlo.
- **Informe 15:** el peso por evidencia pasa a ser el número de clústeres significativos distintos (tope 2); `evidence_score` pasa a −log10 del menor q de permulación.
- **Limpieza:** separar las dos opciones; decidir si `remove_intermediates` conserva sentido una vez cambiado el nulo (puerta de decisión tras el paquete de estadísticos por (gen, ciclo)).
- **Mediciones:** se hacen ligadas a cada cambio, no por separado.

## Hallazgos de la exploración que cambian el diseño

1. **El FCS por bloques no es independiente por columna con el ranking por p de gen.** −log10 p de un gen en el ciclo j necesita sus puntuaciones en los N ciclos (el nulo se puntúa contra los otros N−1, como hace hoy `fcs_empirical_corstat` en RER). Hace falta una etapa `GENE_P` antes de los bloques del FCS (dispersión por lote de genes, que lee las filas de esos genes de todos los bloques). Pseudo-cuenta decidida una vez y aplicada igual al observado y al nulo; comparaciones con `TIE_TOL`. En RER no hay problema (paramétrico).
2. **Filas por posición del nulo.** `perm_pos_cycle_caas.tsv.gz` lo leen `scoring_compute.R`, `posenrich_prep_caas_null.py`, el informe 14 y `accumulation_test.py`. Sin ese fichero: SCORING pasa a contadores por bloque, la acumulación a un indicador por posición, POSENRICH se desactiva con aviso; el informe 11 lee `perm_pos_sample.tsv` y `perm_pos_quantiles.tsv`.
3. **Muestra y cuantiles no pueden ser idénticos bit a bit** (el reservoir usa un RNG único en orden de genes): se siembra por ciclo y se compara estadísticamente.
4. **Filtro de genes y de clústeres.** El filtro de genes (`core/postproc.gene_removal`) se calibra por ciclo, con el recuento de posiciones detectadas por (gen, ciclo) y la bandera de cluster por unidad; el pool descarta además filas en cluster una a una. La tarea debe emitir ambos: bandera por (gen, ciclo) y scores ya filtrados por `remove_clusters`; los recuentos incluyen filas sin score de US.
5. **Salidas estrictas** (`perm_pos_detail_strict`, `b0_observed_strict`, `CAAS_CORE_MERGE_STRICT`, `CAAS_CORE_OBSERVED_STRICT`): `pool_labelings(strict=True)` descarta hipótesis antes de agrupar por ciclo base, no es un filtro posterior. La tarea emite una columna `conserved` por hipótesis o un segundo juego de tablas.
6. **Estampa `gene_stat` = `size_adj_max`** comprobada en `fcs_enrich.R` y `scoring_compute.R`: debe seguir al estadístico elegido.
7. **RER binario** (`rer_bin`) comparte `rer_perm_batches`: subir a 10.000 le afecta. `continuous_rer.R:165` usa el método por defecto de `correlateWithContinuousPhenotype` y `getPermsContinuous` usa `method = "p"`: comprobar que coinciden.
8. Los recuentos de validación de los documentos (1, 0, 6, 6 conjuntos) se calcularon sin la puerta 0,025: comparar primero con la puerta apagada.

## Paquetes de trabajo

| WP | Contenido | Depende de |
|---|---|---|
| WP0 | Congelar línea base | commit del usuario |
| WP1 | Estadísticos del FCS y RER, informes, reducción de redundancia (a posteriori) | WP0 |
| WP2 | Eliminar reutilización de exports; desactivar POSENRICH; separar opciones de limpieza | WP0 |
| WP3 | Fase 1: bloques de ciclos | WP2 |
| WP4 | Fase 2: replay en flujo y fusión | WP3 |
| WP5 | Fase 3: tablas por (gen, ciclo), merge por bloque, SCORING y acumulación | WP4 |
| WP6 | Fase 4: FCS por bloques y `GENE_P` | WP1, WP5 |
| WP7 | `caas_perms` por bloques y exportación opcional | WP5 |
| WP8 | Decisión final de limpieza, documentación, sincronía GUI y plantillas | WP5 |
| WP9 | Experimento del estadístico base del gen | A1 ya; A2 tras WP5 |

Paralelismo: WP1 y WP2/WP3 tocan ficheros disjuntos; WP9-A1 corre aparte; WP4 y WP5 son seriales (comparten el formato de tabla).

### WP0. Línea base

- El usuario hace commit y etiqueta. Salidas de referencia en `~/.cc/golden_<etiqueta>/` con el código sin modificar:
  - `validation/tier2` marine toy (121 genes, 100 ciclos, script `run_tier2_marine_toy_local_complete.sh`), `validation/tier1` pepc (1 gen, 1000 ciclos, 100 hipótesis) y `validation/unification/cancer_b0_toy` (45 genes).
  - Salidas completas del cluster de neoplasia y maligno completo (incluido `strict/`), con sumas de comprobación (por `sbatch`).
  - `fcs_stats.tsv`, `caas_perms.rds`, un `.continuous.perms.rds` de RER y las tablas FCS antiguas.
- Copia congelada de `fcs_run_all` y `fcs_empirical_corstat` como `PhyloPhere_validation/unification/legacy_fcs.R`.
- Ejecutar `RUN_SLOW=1 python3 -m pytest -q unification resources` (desde `PhyloPhere_validation/`) y anotar qué pasa.
- Ampliar `compare_null.py` (tablas y `--extra`) a `scoring/position_scores.tsv` y a las tablas FCS.
- Riesgo oculto: `test_null_consumers`, `test_scores` y `test_compare_null` asumen `gene_cycle_scores.tsv`; se portan en WP5.

### WP1. Estadísticos del FCS, RER e informes

Ficheros: `subworkflows/ENRICHMENT/local/src/fcs_enrich.R` (funciones de Lachenbruch 233 a 281, `fcs_classify_evidence` 564, esquema de columnas 538 y 635 a 644, `fcs_empirical_corstat` 593), `fcs_compute.R`, `subworkflows/ENRICHMENT/fcs.nf`, `conf/enrichment.config`, `conf/rerconverge.config`, `subworkflows/RERCONVERGE/rer_gene_lists.nf`, `continuous_rer.R`, `scoring_compute.R:1176-1186` (escribe `fcs_stats_rer.tsv` desde el p de permulación; hace falta `rer_P` paramétrico aguas arriba), `12.FCS_general_report.Rmd`, `15.Comparison_report.Rmd`, pestaña de enriquecimiento de la GUI, `gui/i18n.py`, `gui/templates/*.json` y `validation/tier*/templates`, `README.md`, `workflows/help.nf`.

1. **Parámetros:** `fcs_min_genes` 10 y `fcs_max_genes` 0 en config, GUI, plantillas, `help.nf`, README y valores por defecto de R (hoy 500 en varias funciones y 5 y 500 en `run_single.sh.j2`).
2. **Lachenbruch fuera:** funciones, columnas `lach_*`, `sig_*`, `evidence_*` y parámetro `fcs_fdr_lachenbruch` (GUI, i18n, plantillas). POSENRICH (copia propia del test) queda intacto y cae a `fcs_fdr`.
3. **Clasificación nueva** (`fcs_classify_significance`): BH sobre `p.perm` por ranking × base → `q.perm`; `sig_perm = q.perm < fcs_fdr`; `sig_top = sig_perm & p.perm < fcs_pperm_thr`; dirección = `sign(stat)` en rankings firmados; se conserva el q analítico como columna.
4. **RER:** ranking global únicamente (se eliminan acelerado y decelerado: `rer_gene_lists.nf:81-82` y las pestañas RER de top y bottom del informe 15); `score_global = sign(Rho) · -log10(P)`; el nulo es `corStat` tal cual (se elimina `fcs_empirical_corstat` y el camino de `corRho`); regla para el leading edge de un ranking firmado (hoy `12.Rmd:318-322` descarta `drv_score <= 0`). Permulaciones a 10.000 con medición previa de `RER_CONT` (`combinePermData` duplica memoria de forma transitoria).
5. **Informe 15:** se sustituye el peso Supported = 2 y Exploratory = 1 (`.evidence_weight_of`, líneas ~139-146 y 1046-1105) por el número de clústeres significativos distintos que contienen el gen (tope 2), y `evidence_score` por −log10 del menor `q.perm`.
6. **Reducción de redundancia (a posteriori, solo para informar):** nuevo `fcs_cluster.R` llamado en el informe 12 tras `sig`: asignación voraz por significancia (el término con menor `q.perm` es representante; cada término siguiente entra en el primer representante con Jaccard de leading edge ≥ umbral, o crea uno nuevo; determinista y sin encadenamiento), salida `cluster_id` y `is_representative`, parámetros de similitud y umbral expuestos en la GUI. No modifica ningún q.
7. **Validación:** `p.perm` de CAAS igual al de referencia hasta 1e-14 relativo; la clasificación igual a una recomputación independiente con `p.adjust` desde la tabla antigua (recuentos 1, 0, 6, 6 con la puerta apagada); estadístico observado de RER frente a `RERconverge::getStat` del commit fijado 2bd328f hasta 1e-12; identidad numérica `corStat == sign(corRho) · -log10(corP)` en el objeto de permulaciones; agrupación probada con conjuntos sintéticos de solapamiento conocido, estabilidad al orden de filas y partición determinista.

### WP2. Reutilización de exports, POSENRICH y limpieza

- **Eliminar la reutilización:** `permDiscFiles`, `perm_disc_in`, `prune_exports`, la publicación de `perm_disc`, el canal `reuse_disc` y la segunda vía de `CAAS_CORE_BATCHED` (`caas_permulation.nf:163-192, 225-292, 345-374`); `main.nf:321-391`, `scoring.nf:183-213`, `enrichment.nf:184`; parámetros `caas_perms_file` y `caas_pos_detail_file` (`conf/scoring.config:46,52`, `help.nf:343-344`, README); campos de `gui/models/precomputed.py`, `run_single.sh.j2:127-150, 243-259, 646, 652`, `sbatch_array.sh.j2` y `report_registry.py:323`. Claves retiradas por el mecanismo `_RETIRED_*` de `serialization.migrate()` sin subir `SCHEMA_VERSION`.
- **POSENRICH:** el módulo se omite con un aviso explícito cuando falta el nulo por posición. Entregable: `docs/POSENRICH_PORTE_EVALUACION.md` (qué necesita del nulo, alternativas sin Lachenbruch, coste, si compensa portarlo).
- **Limpieza:** `cleanup_policy` se sustituye por dos booleanos, `remove_intermediates` y `remove_work` (`conf/common.config:63-65`, `help.nf:165`, `caas_permulation.nf:48-53`, `gui/models/runtime.py:51`, `runtime_tab.py:43-100, 317`, `validate.py:66`, `run_single.sh.j2:366, 513, 797-812`, `sbatch_array.sh.j2:72`, todas las plantillas JSON). `migrate()`: `none` → (F, F), `prune_intermediates` → (T, F), `remove_work` → (T, T). Se añaden pruebas de `migrate()` y de los scripts generados (hoy no existen).
- **Validación:** es borrado puro. Las tablas de marine toy y pepc idénticas a la referencia (`compare_null.py`, tolerancia 1e-12) y DAG de Nextflow sin cambios.

### WP3. Fase 1: bloques de ciclos

- Parámetro fijo `caas_cycle_block_size` (por defecto 500; entrada de config en `conf/ct.config` junto a `ct_core_batch_size`, campo de GUI en `caas_tab.py`, `modules.py`, `i18n.py`, ambos juegos de plantillas).
- `SUBSET_RESAMPLE_PERMS` y `permulations.R` escriben `blocks/blk_NNNN/{fop_labelings.tab, fop_pairs.tsv}` y `b0/` aparte. `CAAS_CORE_BATCHED` pasa a teselas lote de genes × bloque (`combine`); `b_0` en su propia tesela; etiquetados cargados una vez por tarea (se acaba el `ct perm-replay` que reparsea el fichero por gen).
- Interfaz hacia WP4: directorio de tesela con `perm_pos_detail/<Gen>.tsv.gz` (filas de ese bloque), `b0_observed/` y los trozos estrictos; `CAAS_CORE_MERGE` aún los une.
- **Validación:** `compare_null.py` sobre marine toy (bloque de 25) y pepc (1000 ciclos, bloque de 250): todas las tablas idénticas a 1e-12 y el mismo orden de filas por gen (el orden entra en las sumas en coma flotante). Además `test_union_shards`, `test_gene_chunks` y `test_null_roster`.

### WP4. Fase 2: replay en flujo y fusión

- Módulo nuevo `perm_stream.py` en `subworkflows/CT_DISAMBIGUATION/local/src/utils/`: llama a `VectorizedPermReplay` (`subworkflows/CT/local/modules/perm_replay.py`) en proceso, un ciclo base cada vez, y alimenta `driver.score_labelings` y `pool_labelings`. En `gene_wrapper.process_all_genes_perms` se elimina el troceo por memoria (`_REPLAY_MB_PER_MROW`, `_chunk_gene_cycles`, `pending_pooled`). La vía de export de `ct perm-replay` queda solo para el observado.
- **Validación:** prueba nueva `test_stream_vs_export`: ambas vías sobre los genes de la referencia, con los trozos de detalle iguales byte a byte tras descomprimir; `test_kernel_fuzz` y `test_discovery_kernel` repetidos; `test_chunk_failures` reescrito; RSS del gen mayor (RP1) con bloque de 500: ya no debe depender del número de filas del gen; tiempo de carga del contexto de ASR por (gen, bloque), que decide si 500 es el tamaño adecuado.

### WP5. Fase 3: tablas por (gen, ciclo) y merge por bloque

- **Esquema de tesela** (f64, sin rebajar precisión): `units` (gen, ciclo, posiciones detectadas distintas, bandera de cluster, nº puntuadas por dirección) y `scores` (gen, ciclo, posición, lado, score, `conserved`; ya filtrado por `remove_clusters`).
- **`CAAS_CORE_MERGE` por bloque** (dispersión): une los lotes de genes del bloque, aplica `gene_removal` por ciclo, forma los pools exactos por (ciclo, dirección) y calcula las tres variantes (máximo crudo, `F(max)`, `F(max)^n`), con la salida elegida por parámetro `gene_stat`; escribe la matriz genes × B de los tres rankings, `removed_units`, filas de roster, cuantiles por ciclo (RNG por ciclo). `CAAS_CORE_MERGE_STRICT` ejecuta el mismo código sobre el juego estricto; una reducción final pequeña fija roster y universo.
- **SCORING:** `p.emp` con contadores por bloque (`nd`, recuento con score ≥ obs − `TIE_TOL`, primer ciclo), sumados en una reducción; `p.emp_fact` en dos etapas (suma del histograma de `nd` y definición de clases; recuento por clase de detecciones por encima de cada score observado); los scores observados salen de la tesela `b_0`. Acumulación: indicador por posición. Se deja de construir el eje ASR (3 matrices, columnas q90) y se corrige el comentario de `fcs_enrich.R:661-663`.
- **Recursos** recalibrados desde la traza de 1000 ciclos (`conf/resources.config:227-245`): merge del orden de 2 GB + k · B con k medido y escalera de intentos; tesela reajustada; `GENE_P` ~12 GB.
- **Validación:** sobre marine toy, pepc y cancer_b0_toy, por ciclo: `global_caas`, `top_caas` y `bottom_caas` iguales a las nuevas matrices (máx. |dif| ≤ 1e-12 y mismos conjuntos de claves); `removed_units`, `perm_pos_cycle_caas`, `position_scores.tsv`, `p.emp`, `p.emp_fact` (1e-12), roster y `b0/` iguales a la referencia; contadores enteros exactamente iguales; `strict/` repetido con los mismos controles. A escala real, a 1000 ciclos en neoplasia y maligno completo frente a los directorios de referencia, por `sbatch`.

### WP6. Fase 4: FCS por bloques y `GENE_P`

- `fcs_run_all` se parte en `fcs_observed()` (tabla observada y `realenrich`), `fcs_null_count(bloque)` (por conjunto, recuento de estadísticos nulos ≥ observado y columnas válidas) y `fcs_close()` (`p.perm = (Σ recuento + 1) / (Σ válidas + 1)`, BH y clasificación); `fcs_run_all` pasa a ser la composición de un solo bloque.
- Procesos: `FCS_OBSERVED` por lote de GMT, `GENE_P` (dispersión por lote de genes; lee las filas de esos genes de todos los bloques) para CAAS, `FCS_NULL_BLOCK` por lote de GMT × bloque, `FCS_CLOSE`, y `FCS_CONCAT` sin cambio. Para RER, un `RER_PERMS_SPLIT` único que escribe columnas de `corStat` por bloques (RDS reducido de solo `corStat`), de modo que cada tarea no lea el objeto completo.
- **Validación:** la composición nueva frente a `legacy_fcs.R` sobre marine toy, pepc y un conjunto sintético de 400 genes (recuentos exactamente iguales y `stat` a 1e-12); resultado invariante al tamaño de bloque probado en {1, 7, N}; ampliación de `check_fcs_null_equivalence.R`.

### WP7. `caas_perms` por bloques

- Canónico: `blocks/blk_k.matrices.rds` + roster + sello; `scoring_caas_perms.R` lee la tabla por trozos y rellena solo las matrices de CAAS; `--export-single` reconstruye el RDS único bajo demanda. Los consumidores (`scoring_compute.R:520` solo lee el roster, que sale de `cycle_roster.txt`; el FCS lee bloques).

### WP8. Cierre

- Decidir si `remove_intermediates` se retira (si el nulo no deja nada que podar) o pasa a significar "no publicar el almacén masivo por posición". Los `rm` residuales de las plantillas son higiene sin condición.
- README, `help.nf`, docs (`FCS_REVISION_RERCONVERGE.md` con los controles del ranking calibrado, el máximo crudo, las tres pruebas y el BH por gen; `PLAN_NULO_10K_PERMULACIONES.md`), y sincronía GUI y plantillas.

### WP9. Estadístico base del p de gen

- **A1 (ahora, sin fase nueva):** el `perm_pos_cycle_caas.tsv.gz` de referencia ya trae los scores por (gen, ciclo) tras el filtro; un `sbatch` calcula n, máximo y `F(max)` por (gen, ciclo, dirección) y las tres variantes sin repetir el replay. Se añade el exponente `β` con `score = F(max)^(n^β)`, β ∈ {0; 0,25; 0,5; 0,75; 1} (β = 0 es `F(max)`, β = 1 es `F(max)^n`) y el máximo crudo como tercer brazo.
- **Métricas:** uniformidad dentro del ciclo (KS), Spearman con n y con la propensión del gen, ICC del efecto de gen; p de gen contra los otros ciclos (genes con p ≤ 0,01 y 0,05 por ciclo nulo, percentil del observado entre nulos de nivel parecido); sensibilidad al nivel (Spearman entre genes puntuados y genes llamados en ciclos nulos); control "ciclo nulo como observado" del FCS con la regla nueva para el mejor candidato y el máximo crudo; fracción de genes cuyo p cambia más de un decimal al mover β de 0 a 1 y efecto con n > 8.
- **A2 (tras WP5):** repetir a 10.000 ciclos, porque el suelo del p (1/1001) limita a 1000.
- **Criterio:** `F(max)^n` solo si su calibración no es peor que `F(max)` y la dependencia de n que quita es relevante; en otro caso `F(max)`; máximo crudo solo si está mejor calibrado (los controles previos no lo sugieren). Las tres variantes salen del merge, así que la elección es un parámetro.

## Redundancia de términos: análisis

Hay dos usos distintos y la diferencia importa para tu duda.

- **A posteriori (informar, WP1.6):** agrupar los términos ya significativos no cambia ningún p ni q. Solo organiza el resultado y evita contar dos veces términos redundantes en el informe 15. Es seguro.
- **A priori (podar la colección antes de contrastar):** reduce la familia m de BH y por tanto puede subir el número de significativos. Es válido solo si la poda no depende de los datos (solo de la composición de los conjuntos); un filtro que mire significancia o leading edge antes de BH es selección y no controla el FDR. BH controla el FDR bajo dependencia positiva (PRDS; [Benjamini y Yekutieli 2001](https://projecteuclid.org/euclid.aos/1013699998)), de modo que la redundancia no infla el FDR: desperdicia potencia. El marco general de filtrar de forma independiente del estadístico es [Bourgon, Gentleman y Huber 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2906865), que exige independencia bajo la nula y avisa de que ciertos pares filtro-prueba pierden el control del error de tipo I.

**Medido con los GMT del repositorio** (conjuntos con al menos 10 genes en el universo de 12.473 genes de CAAS, solapamiento sobre la composición completa, sin datos): la redundancia está concentrada.

| Base de datos | Conjuntos | Copias exactas | Con Jaccard ≥ 0,9 | Con Jaccard ≥ 0,5 |
|---|---|---|---|---|
| Phenotype Human Phenotype Ontology | 3113 | 162 | 752 | 2079 |
| DisGeNET | 1874 | 321 | 611 | 783 |
| Reactome | 1236 | 88 | 317 | 892 |
| KEGG, Hallmark, c6, c8, GO CC, Wikipathways, cinasas | 50 a 808 | 0 a 1 | 0 a 8 | < 10 % |

- Quitar solo las copias exactas reduce m en un 5 a 17 % en tres bases: efecto despreciable sobre q. Con N = 10.000 el suelo del p baja a 1e-4 y la resolución, no la redundancia, es lo que más limita la significación en las bases grandes (en HPO con m = 3113 hacen falta ~21 conjuntos en el suelo a q < 0,15 con 1000 ciclos y ~2 con 10.000).
- El coeficiente de solapamiento (el que recomienda EnrichmentMap para GO) con asignación voraz por tamaño colapsa Reactome a 26 y HPO a 1 por anidamiento jerárquico: no es adecuado como valor por defecto en bases jerárquicas. Con Jaccard ≥ 0,5 salen 657 y 1858.
- **Conclusión del plan:** la inferencia principal se mantiene sobre las colecciones completas por base de datos, como RERconverge y Redlich; la reducción de redundancia es a posteriori y de presentación (similitud Jaccard por defecto, umbral calibrado con nuestros datos). Como comprobación barata y sin recomputar, WP1 recalcula BH sobre las tablas existentes quitando copias exactas y casi copias (Jaccard ≥ 0,9, ambas independientes de los datos) para cuantificar cuánto cambiaría el número de significativos; no entra como valor por defecto. Un umbral y una similitud distintos requieren calibración, porque las recomendaciones de EnrichmentMap (solapamiento 0,5, Jaccard 0,25) vienen de su documentación, no de un benchmark independiente.

## Riesgos y cómo los mitiga el orden

- **Orden de sumas y empates en coma flotante:** WP3 conserva el orden de filas, WP4 compara los trozos byte a byte, WP5 compara tablas a 1e-12.
- **Confusión de causas:** WP2 es borrado puro con tablas idénticas; WP1 se valida contra `p.perm` congelado antes de tocar el nulo, de modo que un número que se mueva se asigna a un solo cambio.
- **Consumidores ocultos de las filas antiguas** (POSENRICH, informes 11 y 14, acumulación): la auditoría arriba, desactivación explícita de POSENRICH y documento de evaluación.
- **Esquema de tesela equivocado descubierto tarde:** WP4 y WP5 comparten un formato; WP5 se prueba primero en el toy, luego a 1000 ciclos en el cluster y solo después a 10.000.
- **Dependencia entre ciclos del p de gen:** etapa `GENE_P` explícita.

## Verificación extremo a extremo

1. Suite existente: `RUN_SLOW=1 python3 -m pytest -q unification resources` en `PhyloPhere_validation/` tras cada WP, más las pruebas nuevas (`migrate()`, scripts generados, `test_stream_vs_export`, invarianza al bloque).
2. `compare_null.py` y `compare_b0.py` contra `~/.cc/golden_<etiqueta>/` en marine toy, pepc y cancer_b0_toy, tolerancia 1e-12.
3. Escala real por `sbatch`: neoplasia y maligno completo a 1000 ciclos frente a los directorios de referencia (matrices, `removed_units`, `position_scores.tsv`, FCS), incluido `strict/`.
4. Recursos: tras WP5 y WP6, medir RSS y tiempo a 1000, 4000 y 10.000 ciclos con las trazas de Nextflow y `sacct`; criterio: memoria por tarea independiente de N.
5. Objetivo final: una corrida a 10.000 ciclos de neoplasia completo con memoria por tarea de pocos GB y los mismos resultados que la versión actual donde ambas son comparables.

## Decisiones aún abiertas (puertas)

- Tamaño de bloque (propuesta 500, fijo): se confirma con el tiempo de carga del contexto de ASR por (gen, bloque) medido en WP4.
- Estadístico base del p de gen: lo decide WP9.
- Destino de `remove_intermediates`: tras WP5.
- Umbral y similitud por defecto de la reducción de redundancia: calibración con tablas existentes.
- Si POSENRICH se porta: tras leer `docs/POSENRICH_PORTE_EVALUACION.md`.
