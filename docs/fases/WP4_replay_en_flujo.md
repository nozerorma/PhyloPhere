# WP4. Fase 2: replay en flujo y fusión de descubrimiento y replay

Estado: pendiente. Depende de: WP3. Serial con WP5 (comparten el formato de la tesela).
Conversación independiente. Antes de empezar, leer `docs/fases/README.md`.

## Qué hay que comprobar primero (anotado al cerrar WP3)

WP3 queda en el registro de `WP3_bloques_de_ciclos.md` (diseño de bloques en su sección propia). Antes de cambiar nada:

1. **Tiempo del contexto de ASR por (gen, bloque).** Cada tesela (lote de genes x bloque) vuelve a cargar el contexto de ASR de cada gen, de modo que el coste fijo por gen se multiplica por el número de bloques. No está medido: en una tesela de PEPC (bloque de 250 ciclos, un gen) la réplica de ASR ocupó 38 s de los 50 s de la tarea. Medirlo separando carga del contexto y réplica; decide si 500 es un buen tamaño de bloque.
2. **Qué sustituye WP4.** El intercambio entre `ct perm-replay` y la réplica de ASR sigue siendo el fichero `perm_disc/<gen>.perm_replay.discovery.output` (lo lee `disambiguation_perms_main.py`, `_genes_from_export`). Los procesos `ct` ya reciben un manifest con varios genes y parsean las etiquetas del bloque una vez (opciones `--manifest`, `--alignment_dir`, `--claim_dir`).
3. **Interfaz de la tesela que no debe cambiar (WP5 la cambia después):** `perm_pos_detail/<Gen>.tsv.gz` con las filas del bloque (más `perm_pos_detail_strict/` y `b0/`), `b0_observed/` solo en la tesela `b0`. La unión de los trozos de un gen entre bloques (`union_shards.sh` + `merge_block_shards.py`) asume que cada trozo trae sus ciclos en orden ascendente del texto del ciclo y suma `n_detected` entre bloques: si WP4 cambia el orden de filas dentro del trozo, la unión lo rechaza con un error.
4. **Caché de ASR compartida.** Varias teselas del mismo gen leen y, con la caché fría, calculan la misma entrada; el cierre de exclusión por gen (`_gene_cache_lock`) hace que las demás esperen a la primera. Con la caché fría y muchos bloques esa espera crece; compararlo con la caché de la referencia.
5. **Referencias y trampas de validación.** Comparar con `compare_null.py` contra `golden_581b999/{pepc_c4_tolerance,cancer_b0_toy}`; el renderizador de scripts coloca `asr_cache`, `work` y `results` bajo `<outdir>/scripts`, así que la caché de la referencia debe copiarse a `<outdir>/scripts/asr_cache` (con `cp -a`); si no, codeml se recalcula y tres genes del toy de cáncer (MAS1L, OR7G3, SHROOM4) difieren. Las tablas `fcs_all_results.tsv` y `fcs_leading_edge_composition.tsv` de la referencia tienen el esquema anterior a WP1 y no son comparables columna a columna (las columnas comunes sí lo son). El marine toy sigue sin referencia completa (disco externo).

## Lo que WP3 no hizo o hizo distinto (tenerlo en cuenta)

- **Marine toy sin validar.** La referencia de `golden_581b999/marine_toy` está incompleta (murió con el reinicio de las 21:07) y no se puede regenerar sin montar `/media/miguel/external`. WP3 se validó con PEPC c4 (bloque 250) y el toy de cáncer (bloques 25 y 40). Si el disco está montado, ejecutar `bash golden_581b999/marine_toy/launch_resume.sh` y comparar también WP4 en marine.
- **`permulations.R` y `ct_resample.nf` sin cambios:** los bloques se cortan en `SUBSET_RESAMPLE_PERMS`, no en el cosechador de RESAMPLE.
- **Etiquetas parseadas una vez por proceso `ct`, no por tarea:** una tesela lanza hasta `task.cpus` procesos `ct` (12 en `process_resample`), cada uno con su propia lectura. No hay objeto compartido entre procesos.
- **Fichero de etiquetas del bloque:** `resample_perms.tab` (no `fop_labelings.tab`) en modo FOP y simple.
- **`fop_pairs.tsv` global eliminado:** cada tesela lee solo el de su bloque.
- **Orden de filas transitorio:** la unión de trozos entre bloques reordena por el texto del ciclo para igualar una réplica sin bloques, y suma `n_detected`. Es un artefacto que WP5 elimina.
- **No medido:** tiempo de carga del contexto de ASR por (gen, bloque); sin él el tamaño 500 sigue siendo solo una propuesta.
- **Pruebas ya rotas antes de WP3 (no son de WP4):** `test_core_batched` (master congelado con `caap_group`; `CAAS_OBSERVED` con 4 de 5 entradas), 11 de `test_wiring.py`, y 4 de `test_perm_replay_cli/runner` que usan `--caap_mode`. Una fase no debe añadir fallos a las listas de `golden_581b999/pytest_baseline/` ni a este conjunto.
- **`validation/tier*/templates` sin la clave `caas_cycle_block_size`** (las gestiona el usuario): toman 500.
- **Hallazgo abierto de otra fase:** DOMINO debe saltarse el paso en silencio sin aristas (README de fases).

## Objetivo

Que la memoria del worker no dependa del número de filas del gen ni de N: procesar un ciclo base cada vez dentro del worker y fusionar el descubrimiento con el replay de ASR, sin export en disco.

## Estado de partida

- Memoria del replay ∝ filas del gen: `_parse_discovery_entries` (`gene_wrapper.py:181`) carga como objetos todas las filas de los ciclos del trozo; `score_labelings` (`driver.py:185`) devuelve registros de todas las etiquetas; `pool_labelings` (`driver.py:282`) agrupa por ciclo base después. El código estima ~4400 MB por millón de filas (`_REPLAY_MB_PER_MROW`, `gene_wrapper.py:306`) y RP1 llega a 14 GB por worker. El agrupado ya es por ciclo base.
- El descubrimiento (`subworkflows/CT/local/modules/perm_replay.py`, Python) recorre trozos de etiquetados × posiciones (`count`, bucle `for start in range(0, self.B, b_chunk)`), mantiene en listas todos los aciertos del gen y escribe el export en orden de posición, como texto. El replay lo relee y lo reagrupa por ciclo. El export solo lo lee el replay de la misma tarea (exports locales tras WP2).
- `n_detected` solo rellena una columna del detalle; no tiene consumidor.
- Trazas de las pasadas finales de 1000 ciclos (malignant completo): `CAAS_CORE_BATCHED` 83 tareas, mediana 16,6 GB y máximo 30,4 GB (petición 20 GB), CPU ~255 % de 400 %; neoplasia llega a 46,8 GB con 6 fallos.

## Alcance (ficheros)

- Módulo nuevo `subworkflows/CT_DISAMBIGUATION/local/src/utils/perm_stream.py`: llama a `VectorizedPermReplay` en proceso, un ciclo base cada vez (los trozos de etiquetados deben alinearse con ciclos base completos: las hipótesis de un ciclo base se agrupan juntas), y alimenta `driver.score_labelings` y `pool_labelings`.
- `gene_wrapper.process_all_genes_perms`: eliminar el troceo por memoria (`_REPLAY_MB_PER_MROW`, `replay_chunk_count`, `_chunk_gene_cycles`, `pending_pooled`) y las variables de entorno asociadas (`CAAS_PERMS_WORKER_MEM_MB`, `CT_PERM_REPLAY_CHUNK_MEM_MB` en `caas_permulation.nf`).
- `subworkflows/CT/caas_permulation.nf`: `CAAS_CORE_BATCHED` sin los pasos intermedios de export.
- La vía de export de `ct perm-replay` permanece solo para el observado (`b_0`).

## Interfaz hacia WP5

El worker produce, por (gen, ciclo base), los registros agrupados (`pool_labelings`) listos para que WP5 los reduzca a la tabla de tesela. Mantener, de momento, el mismo detalle de salida (`perm_pos_detail/<Gen>.tsv.gz`), para comparar byte a byte con la vía anterior.

## Validación (criterios de paso)

- Prueba nueva `test_stream_vs_export`: ambas vías sobre los genes de la referencia con los trozos de detalle iguales byte a byte tras descomprimir.
- Repetir `test_kernel_fuzz` y `test_discovery_kernel`; reescribir `test_chunk_failures` y `test_gene_chunks` para el nuevo reparto.
- Medir el RSS del gen mayor (RP1) con bloque de 500: ya no debe depender del número de filas del gen. Anotar RSS y tiempo antes y después.
- Medir el tiempo de carga del contexto de ASR por (gen, bloque): decide si B = 500 es el tamaño adecuado (confirmar o ajustar `caas_cycle_block_size`).
- `compare_null.py` sobre marine toy y PEPC: tablas idénticas a la referencia (tolerancia 1e-12).

## Fuera de alcance

Cambiar el formato de la salida de la tesela (WP5), FCS (WP6).

## Riesgos y puntos a verificar

- Gene filas contiguas y orden de hipótesis (`labelings.py`, lectura de ficheros por nombre en orden alfabético: H1, H10, H100, H11…): el orden entra en sumas en coma flotante (`test_pool_domains_order.py`).
- `strict`: `pool_labelings(strict=True)` descarta hipótesis antes de agrupar; el flujo debe seguir produciendo ambos juegos de registros.
- Genes sin detecciones y genes cuyo contexto de ASR no carga: conservar el comportamiento actual (el gen queda fuera del nulo como del observado).

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
