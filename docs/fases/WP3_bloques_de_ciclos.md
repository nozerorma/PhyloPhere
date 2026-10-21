# WP3. Fase 1: bloques de ciclos

Estado: pendiente. Depende de: WP2 (una sola vía en `CAAS_CORE_BATCHED`).
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
