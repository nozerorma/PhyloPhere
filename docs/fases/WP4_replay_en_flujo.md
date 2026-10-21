# WP4. Fase 2: replay en flujo y fusión de descubrimiento y replay

Estado: pendiente. Depende de: WP3. Serial con WP5 (comparten el formato de la tesela).
Conversación independiente. Antes de empezar, leer `docs/fases/README.md`.

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
