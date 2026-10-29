# WP7. Bloques CAAS canónicos y exportación opcional

Estado: IMPLEMENTADO Y VALIDADO CON PRUEBAS DIRIGIDAS (2026-10-23). Aceptación integrada y recursos reales a escala pendientes. Depende de WP5; conserva el FCS calibrado de WP6 y la integración RER.

## Contrato vigente

El nulo canónico es `caas_permulation/blocks/blk_*/matrices.rds`, `cycle_roster.txt`, `null_metadata.json` y las tablas/contadores de WP5. Se conserva exactamente el nombre y contenido de los bloques: `genes`, `cycles`, `gene_stat`, `n_units`, `stat` por ranking y `parts` (`max`, `fmax`, `n`). El índice JSON `caas_blocks_v1` recoge universo, ciclos, estadístico, unidades e identidad MD5 del contenido de cada bloque. `counts.json` conserva además la regla del score por posición.

`--caas_export_single` es booleano y **false por defecto**. Solo activándolo se reconstruye y publica `caas_perms.rds`, tanto en completo como en strict, con el esquema histórico `list(caas_corStat_byrank = list(global, top, bottom), gene_stat)`. Se conserva también la lista vacía histórica cuando no hay unidades; en el formato canónico el roster sigue contando todos los ciclos sin detecciones.

La reducción normal lee una matriz de bloque cada vez para validarla y no asigna matrices genes × N. Retiene IDs y metadatos; la reducción Python retiene contadores del universo de posiciones. La exportación explícita sí asigna tres matrices globales y su memoria crece con el nulo completo. Las reservas de recursos distinguen ambos casos; no son nuevas medidas de escala. Tampoco se afirma memoria constante del pipeline: GENE_P retiene todos los ciclos de un lote de genes, y el bosque RER tiene costes separados.

## Auditoría de productores y consumidores

| Componente | Uso vigente y cambio |
|---|---|
| `CAAS_CORE_MERGE`, alias strict | Siguen produciendo y publicando los bloques WP5, sin cambios numéricos ni de formato. |
| `CAAS_CORE_MERGE_REDUCE`, alias strict | Reducen tablas, posiciones y contadores; validan bloques y escriben índice. El RDS único es output opcional. Reciben el fichero completo de etiquetados de `CAAS_CORE` para detectar un bloque entero ausente. |
| `main.nf` | Pasa `core.labelings` a ambas reducciones; entrega bloques a consumidores tras validarse el índice. Reutilización: entrega el directorio canónico completo al resolver FCS, con roster e índice. |
| `SCORING_COMPUTE` / `scoring_compute.R` | No tienen argumento RDS CAAS ni lo leen. N procede del roster y los p de contadores WP5. Se rechaza roster duplicado. |
| `SCORING_NULL_POSITIONS` | Valida matrices frente al roster también al reutilizar un directorio guardado; los contadores dependen de esa validación. Comprueba que `counts.json` conserve el mismo N. |
| `SCORING` / `ENRICHMENT` | Conservan el paso compatible del RDS opcional con `ifEmpty(NO_FILE)`; su ausencia no vacía las entradas del cálculo. |
| `FCS_CAAS` / `fcs_prepare_blocks.R` | Prefieren bloques vivos o directorio canónico guardado. Validan IDs, ciclos, estadístico, universo e índice cuando existe. Conservan importación de RDS histórico, incluido el vacío. |
| `SCORING_REPORT` / informe 11 | Ya no leían el RDS. Se elimina el slot y parámetro obsoletos del generador de informes de GUI. |
| `SCORING_FCS_REPORT` / informe 12 | El proceso recibe el enriquecimiento precomputado de WP6; se elimina su entrada RDS y su paso al renderer. El Rmd conserva su interfaz histórica para uso independiente. |
| GUI/configuración/ayuda/plantillas | Nuevo booleano en modelo, campo y traducciones, contexto, JSON generado, plantillas GUI y de validación, README y ayuda. Sigue el mecanismo de renderizado de los booleanos CAAS existentes. |
| Cuatro controles `adhoc_scripts/fcs*.sbatch` | Reconstruyen explícitamente desde bloques en un temporal, o leen la referencia histórica si no hay almacén por bloques. Los calibrados usan `gene_p.R`; el control crudo obtiene `parts$max` y `b0/matrices.rds`, con fallback al detalle histórico solo en referencias antiguas. Son controles de memoria completa, no el camino productivo por bloques. |
| `scoring_perms_scaling.sbatch` | Sustituida la entrada retirada `gene_cycle_scores.tsv` por un almacén canónico; modo explícito `validate`/`export`. |
| `~/.cc` | Contiene copias extraídas de informes y lanzadores históricos (`r11.R`, `file*.R`, `carn_render_test.sh`) con interfaces retiradas. No son fuentes productivas ni se han reescrito los artefactos de scratch. Las referencias antiguas siguen legibles por exportación/importación explícita. |

`gene_stat` sigue siendo un parámetro; WP9 decide cuál elegir. No se alteran ranking CAAS calibrado, AUC centrado, lateralidad, empates, pseudo-cuenta, BH, Redlich, observado interno RDS ni disponibilidad por ciclo de RER. No se adelantan decisiones de limpieza de WP8.

## Integridad

El contrato compartido `caas_block_contract.R` valida secuencialmente existencia de bloques; unicidad y validez de IDs; disjunción de ciclos; dimensiones, nombres y orden de matrices; rankings; `gene_stat` común y consistencia con `parts`; universo; unidades; roster local y correspondencia exacta con el roster global. Un ciclo sin detecciones requiere su columna cero, no se elimina. La reducción contrasta además los ciclos contra el fichero de replay completo, excluyendo b_0 y respetando el grano de hipótesis: así detecta la pérdida de un bloque entero aunque los restantes parezcan coherentes.

Los almacenes WP5 guardados sin índice nuevo se aceptan con su roster obligatorio y las validaciones de matrices. En un almacén nuevo se contrasta también la identidad del contenido con el índice; la ausencia total de bloques no se convierte en ausencia opcional de nulo. La compatibilidad RDS histórica es una vía explícita distinta. No puede detectarse una pérdida histórica si tanto el propio roster como todos sus contadores ya fueron truncados de forma coherente y no se dispone del replay original.

## Registro de implementación y validación (2026-10-23)

- Rama de partida: `remove_caap`, HEAD `3b19e24`, confirmado sin cambios ajenos de HEAD durante el trabajo. Implementación y pruebas guardadas en `955e0b5`; registro y tabla de fases en el commit documental siguiente. Se preservan los cambios ajenos de validation, WP4, su fila y prompts; los cambios WP7 en archivos compartidos se guardan por hunks.
- Evidencia, scripts, snapshots estables, comparaciones y logs: `/home/miguel/IBE-UPF/PhD/PhyloPhere_validation/.tmp/wp7/`. Entorno `~/micromamba/envs/phylophere/bin`, TMPDIR y basetemp en ese directorio. No se ha ejecutado código en el cluster.
- `pytest tests/test_wp5_null_blocks_cli.py tests/test_wp6_cli.py tests/test_wp6.py`: **31 PASS**. Incluye los contratos de FCS/RER de WP6 afectados por el resolver compartido. Una primera ejecución detectó el tipo del vector vacío en la comprobación de `parts`; se corrigió y la regresión pasó.
- `pytest tests/test_wp7.py tests/test_wp5_null_blocks_cli.py tests/test_wp6_cli.py`: **40 PASS**; después se añadieron dos esquinas de almacén guardado y se repitió su selección: **4 PASS**. La ejecución final de WP7/resolver: **23 PASS** (`pytest_release.log`), incluidos los 18 casos WP7. En total, 49 casos dirigidos distintos de WP5/WP6/WP7 pasan.
- Exportación de los bloques guardados frente al RDS guardado, `all.equal(tolerance=1e-12)` del objeto completo: **6 PASS**, PEPC completo/strict (1 gen × 1000 ciclos), marine completo/strict (121 × 100), cáncer completo (45 × 100) y strict (41 × 100). Mismos rankings, filas, columnas, orden, ciclos y `gene_stat`. Script `compare.R`, log `compare.log`.
- Nextflow **dos ejecuciones secuenciales**, snapshot `code/`, `run_nf.sh`: reducción completo y strict, `SCORING`, `FCS_CAAS` y ambos informes. Exportación false y true completadas. En false no se publica RDS; en true se publican ambos RDS; en ambos se publica índice y terminan los informes 11/12. La prueba parte de bloques guardados: no regenera teselas ni vuelve a probar el replay completo. Se corrigió la ubicación inicial del harness, que hacía resolver `baseDir` dentro de `tests`; se ejecuta desde la raíz del snapshot.
- Exportación on/off: igualdad exacta de `position_scores.tsv` (3825 filas), `gene_scores.tsv` y `fcs_stats.tsv` (121), resumen del universo (18) y FCS (399 rutas/rankings). Script/log `compare_outputs`.
- Publicación Nextflow: los dos RDS marinos publicados también coinciden con la referencia a 1e-12 (`compare_published.log`). Las cuatro reducciones (exportación on/off × completo/strict) conservan exactamente roster, JSON de contadores, posiciones Parquet y valores de tablas; las tablas TSV se comparan ordenadas por todas sus columnas (`compare_reductions.log`). El sample guardado tenía otro orden de filas, sin diferencias de contenido.
- Frente a controles guardados: SCORING WP5 conserva claves y `p.emp`, `p.emp_fact`, BH, `nd` y scores génicos (1e-12); FCS conserva los 399 resultados calibrados de WP6, estadísticas, p/q y significancia (1e-12). Script/log `compare_prior`. PEPC no tiene rutas evaluables con min_g=10 y no se cuenta como equivalencia positiva de rutas.
- GUI ampliada: `test_gui_params_json.py` + `test_gui_core_batch.py`: **49 PASS / 10 FAIL**. Se reprodujo exactamente el mismo resultado contra una copia de HEAD `3b19e24` (`gui_baseline_final.log`): dos plantillas mantienen selectores antiguos y ocho casos requieren `cancer_multi_final.json`, que ya no existe. No se declara verde esa suite ni se atribuyen sus fallos sin comparación. No se hace limpieza WP8 para satisfacerlos.
- Los cinco controles de `adhoc_scripts/` estaban bajo una regla ignore del directorio; se versionan explícitamente solo estos cinco archivos para que las correcciones viajen al cluster.
- La ayuda de `main.nf --help` desde snapshot estable compila y muestra el nuevo parámetro (log `main_help_stdout.log`).
- `bash -n` de los cinco controles y del launcher WP7: PASS. Los cuerpos R de los cuatro controles FCS pasan `parse()`. `git diff --check`: PASS.
- Launcher dirigido versionado `tests/cluster/prepare_wp7_snapshot.sh`: crea snapshot con `git archive HEAD` y ejecuta `tests/wp7_validation.nf` desde su raíz. El harness está preparado para la referencia marina con completo/strict; el usuario puede sincronizarlo y lanzarlo dentro de una asignación. No se envían trabajos ni se modifica el checkout remoto.

## Pendientes de aceptación conjunta

Corridas integradas, regeneración de nulos y medidas reales de RSS/tiempo/E/S a 1000/4000/10000 ciclos, completo y strict, siguen pendientes de la aceptación WP6/WP7 y del plan general. Las trazas locales de este WP corresponden al toy marino y no acreditan escalabilidad ni memoria constante. La exportación tiene un coste global deliberado y debe medirse aparte del camino normal. Se conservan abiertos los fallos generales documentados, WP1 en revisión, las decisiones WP8 y la elección WP9.
