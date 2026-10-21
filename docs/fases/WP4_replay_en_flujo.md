# WP4. Fase 2: replay en flujo y fusión de descubrimiento y replay

Estado: HECHO (2026-10-21/22), validado con PEPC c4 (1000 ciclos, bloque 250) y el toy de cáncer (bloques 25 y 40); marine toy pendiente (disco externo sin montar). Depende de: WP3. Serial con WP5 (comparten el formato de la tesela).
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
- Conversación y fechas: WP4, 2026-10-21 (noche) y 2026-10-22, conversación independiente, local, rama remove_caap.
- Commit de partida / commit de cierre: 8169f5c / (lo anota el usuario al hacer el commit; ver el final de este registro).
- Qué se hizo (ficheros y cambios):
  * subworkflows/CT/local/modules/perm_replay.py (kernel de ct): `VectorizedPermReplay.count` pasa a ser el total de un generador nuevo, `count_chunks`, que codifica las posiciones una sola vez y emite (inicio, fin, recuentos, aciertos) por rango de etiquetados; `label_chunks` parte los etiquetados en rangos de ciclos base completos (rechaza hipótesis no contiguas). Las filas del export salen de un generador compartido, `perm_discovery_rows`, que usan el export de ct y el flujo (misma reconstrucción de cada acierto); `KernelDisagreement` (subclase de RuntimeError) es la señal de desacuerdo kernel/escalar. El export de ct sigue existiendo, y en la tubería solo lo usa la tesela b0 (exports de b_0).
  * subworkflows/CT_DISAMBIGUATION/local/src/utils/perm_stream.py (nuevo): `ReplaySettings` (umbrales de `ct perm-replay` leídos del mismo fichero de argumentos que arma la tarea), `ReplayContext` (etiquetados del bloque y kernel, por proceso), `stream_gene` (por rango de ciclos base: aciertos -> entradas -> `score_labelings` -> `pool_labelings` -> registros reducidos `SlimRecord`; el contexto de ASR se carga solo al primer ciclo con acierto).
  * src/utils/gene_wrapper.py: `process_all_genes_perms` sin export en disco ni troceo por memoria; cada gen es una tarea del pool (`_replay_gene_task`) que replica, reduce y escribe su trozo; el pool recibe el estado del bloque por inicializador (`_init_replay_worker`). Eliminados `_parse_discovery_entries`, `_scan_perm_discovery_dir`, `_REPLAY_MB_PER_MROW` y constantes, `replay_chunk_count`, `_chunk_gene_cycles`, `_perms_worker_replay`, `_merge_gene_chunks`, `_require_consistent_chunks`, `pending_pooled`, `build_cycle_inputs`, `_read_fop_pairs` y las variables CAAS_PERMS_WORKER_MEM_MB, CAAS_PERMS_MB_PER_MROW, CAAS_PERMS_BASE_MB.
  * disambiguation_perms_main.py: entradas nuevas --manifest, --manifest-alignment-dir, --design, --ct-dir, --replay-args-file, --ali-format, --chunk-labelings; ya no existe --perm-discovery ni `_genes_from_export` (los genes son los del manifest, el mayor alineamiento primero).
  * subworkflows/CT/caas_permulation.nf: CAAS_CORE_BATCHED sin pasos intermedios: la tesela de bloque no lanza ct ni escribe perm_disc; la tesela b0 lanza ct solo con los exports de b_0 (discovery y background); desaparecen CAAS_PERMS_WORKER_MEM_MB y CT_PERM_REPLAY_CHUNK_MEM_MB del .nf (el kernel conserva su presupuesto propio de 512 MB por llamada, con esa variable de entorno como ajuste opcional de ct). nextflow lint sin errores.
  * Sin parámetros nuevos de la tubería (no hay campo de GUI ni plantilla que sincronizar); `--chunk-labelings` es solo una entrada de pruebas del worker.
  * Pruebas (PhyloPhere_validation): wiring/test_stream_vs_export.py (nueva: vía de export reconstruida con ct + lector de texto frente al flujo, FOP con 8 ciclos x 4 hipótesis y etiquetados simples, con rangos por defecto, de un ciclo y de tres; trozos normal y strict iguales byte a byte tras descomprimir; y un gen sin detecciones no carga su ASR); unification/test_count_chunks.py (nueva, 7 pruebas del kernel por rangos); unification/test_chunk_failures.py reescrita (gen no replicable, gen sin detecciones, error, desacuerdo del kernel detiene la corrida); test_detail_only.py y test_dotted_gene_keys.py adaptadas a la tarea por gen (antes fallaban por caap_group); wiring/test_core_batched.py: la réplica directa pasa a ser el flujo sobre el fichero entero. Eliminadas test_gene_chunks.py (la función de fusión desaparece; el orden se cubre en test_stream_vs_export con rangos distintos) y ct_replay_chunking/test_chunk_count.py (la función desaparece).
- Comprobaciones ejecutadas (comando, conjunto de datos, tolerancia, resultado):
  * test_stream_vs_export (3 pruebas): pasan; mutación (orden de hipótesis invertido dentro del ciclo) la hace fallar.
  * PEPC c4, 1000 ciclos, bloque 250 (wp4_check/run_pepc.sh 250; plantilla de validation/ copiada y editada, caché de ASR de la referencia copiada a <outdir>/scripts/asr_cache con cp -a, 0 mensajes "computing it"). compare_null.py --a golden_581b999/pepc_c4_tolerance/results/c4_complete --b <corrida> --tol 1e-12: 11 passed, 0 failed, 10 skipped (ausentes en ambos); delta 0.0 en gene_cycle_scores, perm_pos_sample, perm_pos_quantiles, perm_pos_cycle_caas, b0/, position_scores y fcs_stats; caas_perms.rds all.equal TRUE. strict/: 9 passed, 0 failed.
  * Toy de cáncer (project.json de la referencia; wp4_check/run_cancer.sh), bloques 25 y 40, caché de la referencia, 0 "computing it": 15 passed, 2 failed, 4 skipped en cada tamaño; delta 0.0 y 0 celdas distintas en el nulo, b0/, removed_units, position_scores, fcs_stats, fcs_stats_accum y fcs_leading_edge. Los 2 fallos son los de esquema de WP1 (fcs_all_results.tsv y fcs_leading_edge_composition.tsv, la referencia es anterior). strict/: 11 passed. Bloque 25 contra bloque 40: 17 passed, 0 failed.
  * test_kernel_fuzz y test_discovery_kernel: se repitieron y siguen fallando igual que en la línea base (NameError SCHEMES/US en el escalar de referencia, reference_discovery.py, que sigue escrito para el código con CAAP/GS); no validan nada. En su lugar: test_count_chunks (kernel por rangos frente al recuento de una pasada) y test_stream_vs_export (flujo frente a la vía de export).
  * Suite: unification y resources con RUN_SLOW=1 --continue-on-collection-errors: 250 failed, 314 passed, 17 errors (WP3: 252, 311, 17). Frente a la lista de WP3 no hay fallos nuevos: 8 tests dejan de fallar (chunk_failures, detail_only, dotted_gene_keys, un test de ASR en caché fría); aparecen 6 en test_fcs_no_null solo porque mi lanzador usó el Rscript del sistema (sin RERconverge); con el PATH del entorno pasan 6/6. wiring con Nextflow (test_wiring, test_subset_b0, test_cleanup_options, test_retired_names, test_gui_core_batch, test_stream_vs_export): 23 failed, 78 passed; los 11 de test_wiring.py y los 12 de test_gui_core_batch/test_retired_names ya estaban en las listas de WP0/WP3, ninguno nuevo. test_core_batched: 12 de 14 (los dos de siempre: master congelado con caap_group y CAAS_OBSERVED con 4 de 5 entradas); la prueba principal (trozos independientes del lote y de los bloques e iguales a la réplica directa por flujo) pasa.
- Mediciones (recurso antes y después; una ejecución cada una):
  * Carga del contexto de ASR por (gen, bloque) y reparto de la tesela (registro del propio flujo, PEPC, caché de la referencia): carga del contexto 0,07 a 0,5 s por tesela de 250 ciclos y 0,2 s en la de 1000; núcleo (kernel) 3 a 5,5 s por 250 ciclos; construcción de entradas 1,3 a 2,9 s; puntuación (ASR) 10,7 a 20 s. En el toy de cáncer la carga es de 1,4 a 1,5 s por gen frente a 0,03 s de puntuación (ciclos de 25 a 40 por bloque, genes pequeños). Conclusión: con PEPC la carga del contexto es del 1 al 4 % de una tesela de 250 ciclos, así que 500 ciclos por bloque no es un mal tamaño y no se cambia; en genes pequeños (toy) la carga domina, pero la tesela dura segundos. El 7,3 s de ctx load de WP3 en PEPC no se reproduce (esta vez con la caché ya escrita y sin esperar el cierre de exclusión); no se midió con caché fría.
  * Tesela de PEPC de 250 ciclos, 1 gen (traza de Nextflow, suma del árbol de procesos; 4 workers): WP3 50 a 54 s y 1,5 a 1,6 GB; WP4 18 a 34 s y 1,5 GB (tesela b0: 7,7 s y 0,59 GB). Toy de cáncer: teselas de 10 a 29 s y 0,74 a 2,3 GB (bloque 25 y 40), WP3 de 8 a 28 s y 0,73 a 1,3 GB: sin mejora de RSS de la suma de procesos (RSS de la tesela dominado por el pool y las etiquetas por proceso), el tiempo baja en PEPC.
  * RSS máximo de un solo proceso (/usr/bin/time -v, 1 worker, PEPC, mismos bloques de 250/500/1000 ciclos = 96k/202k/413k filas de descubrimiento): flujo 153/184/246 MB y 11,8/24,7/49,8 s; vía anterior (HEAD 8169f5c reconstruido en .tmp) ct 564 MB y 1,97 GB, réplica 257 MB y 666 MB, tiempo ct+réplica 13,8 y 51,9 s (250 y 1000 ciclos). Pendiente en el flujo con las filas: ~0,29 GB por millón de filas frente a ~4,4 GB por millón de la estimación del código (medida en el cluster, no repetida aquí). El RSS del gen mayor (RP1, 3,55 M filas a 1000 ciclos) NO se midió: no hay RP1 local; hay que medirlo en el cluster (ver pendiente).
- Desviaciones del plan y decisiones tomadas (con su justificación):
  * Rangos de "ciclos base completos" y no de un ciclo base: el kernel hace un bucle Python por posición dentro de cada rango; con un ciclo base (~100 etiquetados) por rango ese bucle se repetiría 500 veces por tesela. El rango es el que da el presupuesto de memoria del kernel (512 MB por llamada), extendido hasta cerrar el ciclo base; el escalado y la reducción sí son por ciclo base.
  * El trozo del gen sigue ocupando memoria ∝ a sus filas reducidas (SlimRecord de 5 campos por registro): el criterio "RSS independiente de las filas del gen" se cumple para las entradas y los registros completos (lo que dominaba), no para el trozo de salida, que WP5 sustituye por tablas por (gen, ciclo). Medido arriba: 0,29 GB/Mrow.
  * El gen entero es una tarea (sin repartir un gen grande entre workers): una tesela con un solo gen usa un worker (en PEPC las teselas bajan de tiempo igualmente). El paralelismo es entre genes de la tesela y entre teselas.
  * Los genes sin detección ya no cargan su contexto de ASR (antes el export vacío de un gen sin acierto lo cargaba); un gen cuyo contexto no carga queda fuera del nulo con aviso, como antes.
  * Un desacuerdo kernel/escalar detiene la corrida (excepción propagada), no deja el gen fuera del nulo.
  * Las etiquetas del bloque se leen una vez por proceso del pool (el pool arranca por forkserver en Python 3.14, así que no se comparten entre procesos): hasta task.cpus lecturas por tesela, como WP3 con ct; el proceso padre las lee otra vez solo para la lista de ciclos.
  * No se tocó caas_cycle_block_size (500) ni validation/tier*/templates.
- Pendiente que pasa a otra fase:
  * Medir en el cluster (sbatch) el RSS y el tiempo de la tesela del gen mayor (RP1) con bloque de 500 y con el código nuevo; hasta entonces no hay cifra de RP1.
  * Marine toy sin validar (disco externo sin montar).
  * Pendiente sin cambiar de WP3: el orden de filas por texto del ciclo y la suma de n_detected siguen en la unión; WP5 los elimina.
  * test_kernel_fuzz/test_discovery_kernel/reference_discovery.py siguen rotos desde la línea base (referencia escalar escrita para CAAP/GS).
```

## Qué debe comprobar primero WP5 (anotado al cerrar WP4)

1. **El trozo de la tesela no ha cambiado de formato:** `perm_pos_detail/<Gen>.tsv.gz` (+ `perm_pos_detail_strict/` y `b0/`), con las filas del bloque en orden de texto del ciclo y `n_detected` por bloque. Lo escribe ahora `gene_wrapper._replay_gene_task` en el worker (ya no el padre) tras `_perms_worker_finalize`; WP5 puede sustituir esa reducción sin tocar el flujo.
2. **Lo que produce el flujo hacia WP5:** `perm_stream.stream_gene` devuelve, por gen, `[(ciclo o ciclo base, [SlimRecord(position, asr_path_score, side, n_hypotheses, n_species)])]` en orden de texto del ciclo, y lo mismo con el pooling strict. Es el punto natural para reducir a tablas por (gen, ciclo) sin pasar por el detalle por fila.
3. **Lo que WP4 no hizo o hizo distinto:** (a) RP1 sin medir (cluster); (b) el trozo de salida del gen sigue ∝ a sus filas reducidas (0,29 GB por millón de filas medido en PEPC); (c) un gen es una tarea indivisible; (d) las etiquetas se parsean una vez por proceso del pool (forkserver); (e) marine toy sin validar; (f) no hay parámetro nuevo; (g) test_kernel_fuzz y test_discovery_kernel siguen sin servir; (h) las pruebas con Nextflow rotas antes (11 de test_wiring.py, 2 de test_core_batched.py, 4 de test_perm_replay_cli/runner con --caap_mode) siguen igual.
4. **Trampas de validación (siguen vigentes):** caché de ASR de la referencia a `<outdir>/scripts/asr_cache` con `cp -a`; lanzar pytest con el `Rscript` de micromamba en el PATH; `--continue-on-collection-errors` para la suite (17 errores de recolección previos); no lanzar dos corridas de Nextflow a la vez (14 GB de RAM).
