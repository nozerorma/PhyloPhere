# WP3. Fase 1: bloques de ciclos

Estado: HECHO (2026-10-21), validado con PEPC c4 (1000 ciclos, bloque de 250) y el toy de cáncer (100 ciclos, bloques de 25 y de 40); marine toy pendiente (disco externo sin montar). Depende de: WP2 (una sola vía en `CAAS_CORE_BATCHED`).
Conversación independiente. Antes de empezar, leer `docs/fases/README.md`.

## Objetivo

Fragmentar el nulo por ciclos con un tamaño de bloque fijo, de modo que cada tarea reciba solo las etiquetas de su bloque y las cargue una vez. Es la base de la RAM independiente de N.

## Decisiones aplicables (del usuario)

- Tamaño de bloque fijo (propuesta 500). Se confirma con el tiempo de carga del contexto de ASR por (gen, bloque) que se mida en WP4.
- RESAMPLE escribe etiquetas por bloques.

## Estado de partida

Medido: cada tarea de `CAAS_CORE_BATCHED` recibe el `resample_perms.tab` completo (~100 hipótesis por ciclo; 11,7 MB a 1000 ciclos, ~117 MB a 10.000) y el `fop_pairs.tsv` (19 MB a 1000, ~197 MB a 10.000). `run_ct_perm_replay_batch.sh` lanza un `ct perm-replay` por gen y `ct` (línea ~248) llama a `simtrait_revive` sobre el fichero entero cada vez: 2,1 s y 328 MB por proceso a 1000 ciclos (~3,3 GB y ~24 s a 10.000, extrapolado), con 4 workers por tarea. `process_all_genes_perms` reparte `cycle_labelings` y `fop_pairs` completos en cada tarea del pool (`gene_wrapper.py:~1151`). RESAMPLE a 10.000 ciclos: 10,8 min y 4,4 GB (no limita).

## Alcance (ficheros)

- Parámetro nuevo `caas_cycle_block_size` (fijo; por defecto 500): `conf/ct.config` junto a `ct_core_batch_size` (~línea 46), `gui/models/modules.py`, `gui/widgets/tabs/caas_tab.py`, `gui/i18n.py`, `gui/generation/templates/sbatch_array.sh.j2` y `run_single.sh.j2`, `gui/models/serialization.py`, todas las plantillas JSON (`gui/templates`, `validation/tier1/templates`, `validation/tier2/templates`; la plantilla tier 1 y `gui/templates/tier1_pepc_c4.json` deben quedar idénticas), `workflows/help.nf`, `README.md`.
- `subworkflows/CT/caas_permulation.nf`: `SUBSET_RESAMPLE_PERMS` (~72 a 160) y `CAAS_CORE_BATCHED` (~171 a 317); flujo de `CAAS_CORE` (~323 a 380) con `combine` de lotes de genes × bloques.
- `subworkflows/CT/local/scripts/permulations.R` (escritura de `fop_labelings.tab` y `fop_pairs.tsv`), `subworkflows/CT/ct_resample.nf`, `subworkflows/CT/local/scripts/run_ct_perm_replay_batch.sh`, `subworkflows/CT/local/ct` (invocación `perm-replay`), `build_b0_labelings.py`.

## Cambios

1. `SUBSET_RESAMPLE_PERMS` y `permulations.R` escriben `blocks/blk_NNNN/{fop_labelings.tab, fop_pairs.tsv}` y `b0/` (diseño y etiquetas de `b_0`) aparte.
2. `CAAS_CORE_BATCHED` pasa a teselas lote de genes × bloque; `b_0` en su propia tesela (una sola vez).
3. Las etiquetas se cargan una vez por tarea (se acaba el reparseo por gen).
4. El replay en dos pasos sigue dentro de cada tarea (exports locales a la tarea).

## Interfaz hacia WP4

Directorio de tesela con `perm_pos_detail/<Gen>.tsv.gz` (filas de ese bloque), `b0_observed/` y los trozos estrictos (`perm_pos_detail_strict`, `b0_observed_strict`); `CAAS_CORE_MERGE` aún los une como hoy. Documentar el diseño de directorios de bloque.

## Validación (criterios de paso)

- `compare_null.py` sobre marine toy (bloque de 25) y PEPC (1000 ciclos, bloque de 250): todas las tablas idénticas a la referencia, tolerancia 1e-12 y mismo orden de filas por gen (el orden entra en las sumas en coma flotante).
- `test_union_shards`, `test_gene_chunks` y `test_null_roster` pasan; la suite no empeora respecto al registro de WP0.
- Resultado independiente de B: probar al menos dos tamaños de bloque en marine toy.

## Fuera de alcance

Replay en flujo (WP4), formato de tabla nuevo (WP5).

## Riesgos y puntos a verificar

- El orden de las hipótesis y de las filas dentro del gen entra en las sumas: no reordenar al partir en bloques.
- Roster y N: los ciclos sin detección deben contarse igual (el roster lista todos los ciclos repetidos).

## Diseño de los directorios de bloque (lo que produce esta fase)

`SUBSET_RESAMPLE_PERMS` toma los N primeros ciclos permutados por orden numérico del identificador (`b_1`, `b_2`, ..., no el orden del texto) y los reparte en bloques de `caas_cycle_block_size` ciclos consecutivos de esa lista. Un ciclo base con todas sus hipótesis (`b_12~H1`, `b_12~H2`, ...) cae siempre en un solo bloque. Salidas:

```
resample_perms.tab            b_0 y los N ciclos, completo (lista de ciclos del nulo: N y roster; se publica como antes)
blocks/b0/                    el etiquetado real b_0, nunca parte de un bloque
    resample_perms.tab          sus etiquetas (`b_0` o `b_0~H<m>`)
    fop_pairs.tsv               PSS de b_0 (solo con hipótesis FOP)
blocks/blk_0001/              primer bloque (blk_0002, ... los siguientes)
    resample_perms.tab          etiquetas del bloque, en el orden numérico de los ciclos
    fop_pairs.tsv               filas de PSS de los ciclos del bloque, con cabecera (solo FOP)
```

Con N = 0 solo existe `blocks/b0/`. Los bloques crecen al aumentar N sin cambiar los anteriores (la lista de ciclos es un prefijo estable), lo que permitirá reutilizarlos en WP5 a WP7.

`CAAS_CORE_BATCHED` ejecuta una tesela por (lote de genes, directorio de bloque): la tesela `b0` reproduce el etiquetado real y es la única que escribe `b0_observed/` (y `b0_observed_strict/`); las teselas de bloque solo escriben `perm_pos_detail/` (y `perm_pos_detail_strict/`). Cada tesela deja `perm_pos_detail/<Gen>.tsv.gz` con las filas de su bloque (la tesela `b0` deja `perm_pos_detail/b0/<Gen>.tsv.gz`).

`CAAS_CORE_MERGE` une las teselas con `union_shards.sh`: el trozo de un gen que está en una sola tesela se enlaza como antes; el de un gen presente en varias teselas se une con `merge_block_shards.py` en un trozo único con el orden de filas de una réplica sin bloques (ciclos en orden ascendente del texto del identificador, `b_100` antes de `b_11` antes de `b_3`; las filas de un ciclo seguidas y en el orden en que se produjeron) y con `n_detected` sumado entre bloques (cada bloque cuenta solo sus ciclos). Las lecturas posteriores (`reaggregate_perm_scores.py`, `CT_ACCUMULATION`) ven el mismo directorio de trozos por gen que antes.

## Registro de la fase

```
- Conversación y fechas: WP3, 2026-10-21 (noche), conversación independiente, local, rama remove_caap.
- Commit de partida / commit de cierre: ebec513 / pendiente (lo anota el usuario al hacer el commit de cierre).
- Qué se hizo (ficheros y cambios):
  * Parámetro caas_cycle_block_size (por defecto 500): conf/ct.config, workflows/help.nf, README.md (texto y fila de la tabla), gui/models/modules.py, gui/widgets/tabs/caas_tab.py (campo), gui/i18n.py (etiqueta en los seis idiomas), run_single.sh.j2 (params.json), sbatch_array.sh.j2 (export), y las 14 plantillas de gui/templates (valor "500"). gui/models/serialization.py no cambia: un proyecto sin la clave toma el valor por defecto del modelo. validation/tier*/templates no se han tocado (las gestiona el usuario): sin la clave toman 500; tier1_pepc_c4.json de validation/ y de gui/templates ya difieren en otras claves por cambios del usuario.
  * subworkflows/CT/caas_permulation.nf. SUBSET_RESAMPLE_PERMS recibe el tamaño de bloque y escribe blocks/b0/ y blocks/blk_NNNN/ (diseño en la sección anterior); conserva resample_perms.tab completo (lista de ciclos del nulo y salida publicada) y deja de sacar el fop_pairs.tsv global. CAAS_CORE_BATCHED es una tesela (lote de genes x directorio de bloque): replay con ct perm-replay sobre las etiquetas del bloque, ASR del mismo bloque con --resample-dir <bloque> y --fop-pairs <bloque>/fop_pairs.tsv; la tesela b0 es la única que exporta y puntúa b_0 (b0_observed pasa a salida opcional). CAAS_CORE hace live.combine(blocks.flatten()); CAAS_CORE_MERGE une las teselas con union_shards.sh. main.nf y workflows/ct.nf: el canal caas_fop_pairs se sustituye por caas_blocks (CAAS_CORE recibe subset, blocks, tree, lengths, taxid).
  * subworkflows/CT/local/ct: perm-replay admite --manifest (id y alignment por línea), --alignment_dir y --claim_dir; con manifest las rutas de exportación llevan {id}. El proceso parsea las etiquetas una vez y toma genes del manifest reclamando <claim_dir>/<id>. run_ct_perm_replay_batch.sh lanza min(workers, genes) procesos ct, cada uno con su watchdog y con el estado de salida propagado como antes.
  * subworkflows/CT/local/scripts/union_shards.sh (detecta los trozos presentes en varias teselas) y merge_block_shards.py (nuevo): unen los trozos de un gen en un trozo único, con las filas en el orden de una réplica sin bloques y n_detected sumado entre bloques.
  * permulations.R y ct_resample.nf no cambian (ver desviaciones). build_b0_labelings.py no cambia.
  * Pruebas (PhyloPhere_validation): unification/test_cycle_blocks.py (nuevo, 8 pruebas: unión de trozos de varios bloques, orden y n_detected, trozo fuera de orden rechazado, enlaces y fuentes intactas, ct con manifest y varios workers frente a la llamada directa); unification/test_perm_replay_runner.py adaptado al lanzador de un proceso por worker; wiring/test_subset_b0.py adaptado al nuevo SUBSET (+4 pruebas de bloques: tamaños 1, 2, 3 y 500, FOP con pares por bloque, orden numérico y no textual, tamaño 0 rechazado); wiring/test_core_batched.py y mini_core*.nf adaptados a las teselas (prueba nueva: trozos con bloques de 2+1 ciclos, lotes de 1 y 2 genes y réplica directa sin bloques, idénticos texto a texto) y con constantes de la época CAAP corregidas donde bloqueaban la cadena; wiring/test_wiring.py comprueba el cuarto argumento (blocks) de CAAS_CORE.
- Comprobaciones ejecutadas (comando, conjunto de datos, tolerancia, resultado):
  * PEPC c4, 1000 ciclos, bloque 250 (4 bloques + b0: 5 teselas de un gen). Plantilla validation/tier1/templates/tier1_pepc_c4.json con caas_cycle_block_size=250 (copia en wp3_check/pepc_b250/template.json), render_tier1_scripts.py, scripts copiados a la raíz del repo y retirados después. compare_null.py --a golden_581b999/pepc_c4_tolerance/results/c4_complete --b <corrida> --tol 1e-12: 11 passed, 0 failed, 10 skipped (ausentes en ambos); delta máximo 0.0 y 0 celdas distintas bit a bit en gene_cycle_scores, perm_pos_sample, perm_pos_quantiles, perm_pos_cycle_caas, caas_perms.rds (all.equal TRUE), los mismos de b0/, position_scores y fcs_stats. strict/: 9 passed, 0 failed. Salvedad: esta corrida usó una caché de ASR vacía (el renderizador coloca asr_cache bajo <outdir>/scripts, no donde yo copié la de la referencia), así que codeml se recalculó; con un solo gen dio el mismo resultado que la referencia (delta 0.0), pero la comparación no se apoya en la caché de la referencia.
  * Toy de cáncer (project.json de la referencia: 45 genes muestreados, 31 tras filtros, 100 ciclos, remove_intermediates, fcs_min_genes 5, fcs_max_genes 1000) con la caché de ASR de la referencia copiada a <outdir>/scripts/asr_cache (cp -a; 0 mensajes "computing it": la caché solo se lee, se reescribe tree_paml.nwk.node_map.json). Bloque 25 (4 bloques + b0, 15 teselas) y bloque 40 (3 bloques + b0, 12 teselas). compare_null.py contra golden_581b999/cancer_b0_toy/results/neoplasia_prevalence_toy_complete, 1e-12: 15 passed, 2 failed, 4 skipped en cada tamaño; delta 0.0 y 0 celdas distintas en las tablas del nulo (gene_cycle_scores, perm_pos_*, removed_units, caas_perms.rds), en b0/, en position_scores, fcs_stats, fcs_stats_accum y fcs_leading_edge. Los dos fallos son tablas del FCS cuyo esquema cambió en WP1 (fcs_all_results.tsv: columnas lach_* sustituidas por q.perm, sig_perm, sig_top, redlich_level, cluster_id; fcs_leading_edge_composition.tsv no existía en la referencia); en las 183 filas comunes stat, pval, p.adj, p.perm y num.genes coinciden con delta 0.0. strict/: 11 passed, 0 failed (ambos tamaños). Bloque 25 contra bloque 40, ambos entre sí: 17 passed, 0 failed (resultado independiente del tamaño de bloque).
  * Primer intento del toy de cáncer (caché vacía por el mismo motivo): 14 tablas distintas, delta máximo 0,31, con tres genes (MAS1L, OR7G3, SHROOM4) de ASR distinto; es la no reproducibilidad de codeml de WP0 y desapareció al usar la caché, no es un efecto del código.
  * Primer intento de PEPC con bloques: CAAS_CORE_MERGE falló (filas "no en orden de ciclo"): yo había supuesto el orden de las etiquetas ordenadas; el orden real es el del texto del ciclo (lo fija _merge_gene_chunks). Corregido y repetido con -resume (las cinco teselas estaban en caché).
  * Marine toy: NO ejecutado. La corrida de referencia de golden_581b999/marine_toy se interrumpió con el reinicio del equipo (21:07), results/ solo tiene las primeras etapas, y la reanudación con RESUME=1 falló porque /media/miguel/external (alineamientos toy) no está montado (dispositivo sda1, etiqueta external). Anotado en su PROVENANCE.txt. En su lugar se probaron dos tamaños de bloque (25 y 40) en el toy de cáncer.
  * Pruebas: unification y resources con RUN_SLOW=1: 252 failed, 311 passed, 17 errors (referencia de WP0: 252 failed, 302 passed, 17 errors; los mismos fallos, +9 pasadas = pruebas nuevas). test_union_shards 4/4, test_gene_chunks 2/2, test_null_roster 7/8 (el fallo, test_p_emp_counts_the_cycles_of_the_roster..., ya fallaba en la referencia). Los cuatro fallos de test_perm_replay_cli/runner que ejercen el ct real siguen siendo los de la referencia (usan --caap_mode, retirado); las pruebas nuevas cubren ese camino sin la opción. wiring sin Nextflow (11 ficheros, con test_cleanup_options): 55 failed y 110 passed, los mismos 55 fallos de la lista de referencia (98 passed de la referencia más las 12 pruebas de WP2); WP3 no añade pruebas en esos ficheros. wiring con Nextflow: test_subset_b0 14/14; test_core_batched 12/14 (fallan test_the_b0_master_..._frozen_one, que compara con un master congelado con caap_group, y test_a_discovery_tab_that_exists..., cuyo mini_caas_observed.nf llama a CAAS_OBSERVED con 4 de 5 entradas: los dos son anteriores a WP3 y no pasan por las teselas); test_wiring.py 15/26: los mismos 11 fallos que en un árbol limpio de ebec513 (git archive en .tmp, retirado), la prueba del argumento blocks pasa ahora y en ebec513 no. nextflow lint: caas_permulation.nf y ct.nf sin errores; main.nf 4 errores previos (WorkflowMap no resuelta por el linter, procesos dentro de cierre).
- Mediciones (recurso antes y después; una sola ejecución cada una, con distinta concurrencia, no son referencias de rendimiento):
  * PEPC, un gen, 1000 ciclos: referencia, una tarea de 1 m 3 s y 4,4 GB de RSS pico (8 GB virtual). WP3: cuatro teselas de bloque de 250 ciclos de 50 a 54 s y 1,5 a 1,6 GB de RSS pico cada una (las cuatro a la vez en 12 CPU), y la tesela b0 de 9,5 s y 0,6 GB.
  * Toy de cáncer, 100 ciclos, lotes de 20, 20 y 5 genes: referencia (una tarea por lote): lotes de 20, 20 y 5 genes con 1,6 GB, 1,6 GB y 1,9 GB de RSS pico y 2 m 39 s, 3 m 12 s y 2 m 2 s. Bloque 25: teselas de 12,6 a 28 s y de 0,73 a 1,2 GB; bloque 40: 8,3 a 28,4 s y de 0,92 a 1,3 GB.
  * Lectura de las etiquetas (simtrait_revive, proceso nuevo, 1000 ciclos x 100 hipótesis): fichero completo 20,8 MB, 100.100 etiquetas, 3,5 s y 443 MB de RSS pico; un bloque de 250 ciclos, 5,2 MB, 25.000 etiquetas, 0,64 s y 118 MB; el de b_0, 21 kB y 10 MB. Antes se leía el fichero completo en cada gen; ahora una vez por proceso ct (hasta task.cpus procesos por tesela, 12 en process_resample).
  * No medido: el tiempo de carga del contexto de ASR por (gen, bloque), que decide si 500 es un buen tamaño (WP4). Con bloques cada gen carga su contexto una vez por bloque, de modo que el coste fijo por gen se multiplica por el número de bloques.
- Desviaciones del plan y decisiones tomadas (con su justificación):
  * permulations.R y ct_resample.nf no escriben los bloques: SUBSET_RESAMPLE_PERMS ya elige los N primeros ciclos por orden numérico y ve el fichero completo una vez (a 10.000 ciclos RESAMPLE y su fichero no limitan, medido en el plan), y partir ahí evita tocar el cosechador de RESAMPLE y mantiene la elección de los N primeros en un solo sitio. Coste: el subset lee los ficheros enteros; es una tarea process_low de una sola vez.
  * Nombre del fichero de etiquetas del bloque: resample_perms.tab (no fop_labelings.tab), el mismo en modo FOP y en modo simple: es el nombre que ya leen read_labelings (resample_*.tab) y el lanzador, y evita un camino por modo.
  * "Etiquetas cargadas una vez por tarea" se cumple por proceso ct, no por tarea: cada tesela lanza min(workers, genes) procesos ct, cada uno con una lectura y varios genes (reclamo por directorio). Una sola lectura por tarea exigiría compartir el objeto entre procesos (fork tras la carga o hilos), a costa del aislamiento por proceso que da la salida 137 del OOM y el watchdog de CPU. Con 12 workers y bloques de 500 ciclos son 12 lecturas por tesela (de unos 0,2 GB cada una, extrapolado desde los 118 MB medidos a 250 ciclos), en lugar de una por gen.
  * El orden de las filas de un gen sin bloques es el del texto del ciclo, no el numérico; los bloques son por orden numérico (prefijo estable al crecer N, base de WP5 a WP7). La unión reordena por el texto del ciclo (mezcla de flujos ya ordenados, sin cargar el trozo entero) para que el trozo unido coincida con el de una réplica sin bloques; esto es un artefacto transitorio de WP3 que WP5 elimina al sustituir el trozo por tablas por (gen, ciclo).
  * n_detected del trozo se suma entre bloques en la unión. Nadie lo lee hoy (solo se escribe en el detalle), pero se conserva el valor de una réplica sin bloques.
  * El fop_pairs.tsv global deja de existir como entrada de CAAS_CORE_BATCHED y como salida de SUBSET (cada tesela lee las filas de PSS de su bloque; a 10.000 ciclos el global era de unos 197 MB por tarea, cifra extrapolada del plan).
  * serialization.py no cambia (el valor por defecto del modelo basta para proyectos sin la clave).
  * Pruebas de wiring: se corrigieron constantes obsoletas que impedían ejercer la cadena de teselas (caap_mode en la réplica directa, tax_id repetido en el fixture, global_meta_caas, 217 filas del master). No se arreglaron las dos pruebas que dependen del master congelado con caap_group ni la de CAAS_OBSERVED.
- Pendiente que pasa a otra fase:
  * WP4 (réplica en flujo). Qué comprobar primero: (1) medir en la tesela el tiempo de carga del contexto de ASR por (gen, bloque) y el reparto entre ct perm-replay y la réplica de ASR (en PEPC, 38 s de 50 s de una tesela son el ASR); decide el tamaño de bloque (hoy 500 por propuesta). (2) El intercambio ct -> ASR sigue por disco (perm_disc/<gen>.perm_replay.discovery.output): sigue siendo el objeto que WP4 sustituye; el consumidor es disambiguation_perms_main.py (_genes_from_export, _parse_discovery_entries). (3) Cada tesela sigue sacando perm_pos_detail/<Gen>.tsv.gz (+ strict y b0/) y CAAS_CORE_MERGE aún los une; WP5 cambia ese formato. (4) La caché de ASR la leen a la vez varias teselas del mismo gen: el cierre de exclusión (_gene_cache_lock) existe, pero con la caché fría varias teselas esperan a la primera que calcula; WP4 debe contar con ello. (5) Con 12 CPU por tesela el replay de ct y la réplica de ASR corren seguidos, no a la vez.
  * Marine toy: ejecutar la comparación (bloque 25 y otro tamaño) cuando el usuario monte /media/miguel/external y se regenere la referencia (golden_581b999/marine_toy/launch_resume.sh con RESUME=1).
  * Sincronizar validation/tier*/templates con caas_cycle_block_size (usuario).
```
