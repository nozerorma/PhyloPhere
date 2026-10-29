# WP2. Eliminar la reutilización de exports, desactivar POSENRICH y separar las opciones de limpieza

Estado: HECHO (2026-10-21), validado con PEPC c4 y el toy de cancer; marine toy validado después, con el código de WP4 (registro de WP4). Depende de: WP0. Independiente de WP1.
Conversación independiente. Antes de empezar, leer `docs/fases/README.md`.

## Objetivo

Tres cambios de borrado o reorganización sin efecto numérico: quitar la vía de reutilización de exports de perm-discovery y de nulos precalculados, desactivar POSENRICH con aviso, y sustituir `cleanup_policy` por dos opciones independientes.

## Decisiones aplicables (del usuario)

- La reutilización de exports no se usa y se elimina.
- POSENRICH se desactiva con aviso explícito cuando falta el nulo por posición; documento sobre si y cómo portarlo: `docs/POSENRICH_PORTE_EVALUACION.md` (ya redactado, revisar al cerrar la fase).
- `remove_intermediates` y `remove_work` pasan a ser opciones independientes. Se decide más adelante (WP8) si `remove_intermediates` conserva sentido una vez cambiado el nulo.

## Alcance (ficheros)

**Reutilización de exports** (todo lo que la toca):
- `subworkflows/CT/caas_permulation.nf`: `permDiscFiles`, `perm_disc_in`, `prune_exports`, la publicación de `perm_disc`, el canal `reuse_disc` y la segunda vía de `CAAS_CORE_BATCHED` (~163 a 192, 225 a 292, 345 a 374); el comentario de cabecera (~34).
- `main.nf` (~321 a 391, descubrimiento de exports y comprobación de `b_0`; `perm_reuse_ch` ~309 y ~386), `workflows/scoring.nf` (~183 a 213), `workflows/enrichment.nf` (~184).
- Parámetros `caas_perms_file` y `caas_pos_detail_file`: `conf/scoring.config:46,52`, `workflows/help.nf:343-344`, `README.md:442`; el mensaje de nulo obsoleto de `fcs_enrich.R:679-682`.
- GUI: `gui/models/precomputed.py:124,131-139` (`caas_perms_file`, `caas_pos_detail_file`, `perm_dir`), `gui/generation/templates/run_single.sh.j2:127-150, 243-259, 646, 652`, `sbatch_array.sh.j2` (PRECOMP_*), `gui/generation/report_registry.py:323`.
- Se conservan: `caas_cycle_roster_file`, `caas_pos_cycle_caas_file`, `caas_pos_sample_file`, `caas_pos_quantiles_file` (`conf/scoring.config:47-53`).
- Claves retiradas: mecanismo `_RETIRED_*` de `gui/models/serialization.py` (`migrate()` ~74 a 107) sin subir `SCHEMA_VERSION`.

**POSENRICH:** `subworkflows/ENRICHMENT/posenrich.nf` y `workflows/enrichment.nf` (omitir el módulo con aviso), informe 14.

**Limpieza:** `conf/common.config:63-65`, `workflows/help.nf:165`, `subworkflows/CT/caas_permulation.nf:48-53` (`shouldPruneIntermediates`) y usos en ~208, ~210 y ~292 a 297, `gui/models/runtime.py:51`, `gui/widgets/tabs/runtime_tab.py:43-100, 317`, `gui/generation/validate.py:66`, `gui/generation/templates/run_single.sh.j2:366, 513, 797-812`, `sbatch_array.sh.j2:72`, todas las plantillas JSON (valores actuales: `remove_work` en `carn`, `neoplasia_def`, `malignant_def`, `tier1_pepc_c4`, `tier1_pepc_c4_nocons`, `tier2_marine_toy`; `prune_intermediates` en `cancer_no_prune_multi_toy`; `none` en `LQ_fabio` y `log_LQ_fabio`; sin clave en `cancer_def`, `cancer_report`, `carn_toy`, `diet_reports`, `tier2_marine`).

## Cambios

1. Borrar la reutilización de exports y los parámetros y campos de GUI asociados; claves antiguas por `migrate()`.
2. POSENRICH omitido con aviso cuando no hay nulo por posición.
3. `cleanup_policy` → dos booleanos, `remove_intermediates` y `remove_work`, en todos los puntos anteriores. `migrate()`: `none` → (F, F), `prune_intermediates` → (T, F), `remove_work` → (T, T) (conserva el comportamiento actual de `remove_work`; los controles nuevos son independientes). `help.nf` y la GUI describen cada opción por separado. Hoy `shouldPruneIntermediates()` solo actúa en `CAAS_CORE_BATCHED` (poda de exports y `mv` de trozos); el resto de `rm` de las tareas es higiene sin condición.
4. Añadir pruebas de `migrate()` y de los scripts generados (hoy no existen para la limpieza).

## Interfaz hacia las fases siguientes

- `CAAS_CORE_BATCHED` con una sola vía (réplica del descubrimiento y del replay en la misma tarea). WP3 parte de aquí.

## Validación (criterios de paso)

- Es borrado puro: las tablas de marine toy y de PEPC idénticas a la referencia (`compare_null.py`, tolerancia 1e-12) y el DAG de Nextflow sin cambios (comparar `pipeline_dag` o la lista de procesos).
- Pruebas nuevas de `migrate()` (los tres valores antiguos y proyecto sin la clave) y de los scripts generados con cada combinación.
- La GUI abre y guarda una plantilla antigua sin perder campos.

## Fuera de alcance

Bloques de ciclos (WP3), formato de teselas (WP5), decisión final sobre `remove_intermediates` (WP8).

## Riesgos y puntos a verificar

- Quitar `caas_pos_detail_file` elimina la reconstrucción del merge dentro de SCORING (`scoring.nf:187-196, 207-213`): comprobar que ningún informe la usa.
- Validar `remove_work` en una corrida real no se hizo todavía (comparación de huellas y ahorro de disco): se hace en WP8 o aquí si hay una corrida disponible.

## Registro de la fase

```
- Conversación y fechas: WP2, 2026-10-21, en la misma conversación que WP0 (a petición del usuario), local, rama remove_caap.
- Commit de partida / commit de cierre: b8e8781 (docs de WP0 cerrado; código de la tubería de a92fd2d) / e2f9515.
- Qué se hizo (ficheros y cambios):
  * Reutilización de exports eliminada. subworkflows/CT/caas_permulation.nf: sin permDiscFiles, perm_disc_in, publicación de perm_disc, canal reuse_disc ni emit perm_discovery; CAAS_CORE_BATCHED con una sola vía (réplica y ASR en la misma tarea); la entrada de CAAS_CORE ya no lleva reuse_disc. main.nf: fuera el descubrimiento de exports (perm_reuse_ch, candidate_base_dirs, comprobación de b_0 de los exports) y la condición por caas_perms_file/caas_pos_detail_file; se conserva la vía de arranque desde un fichero de remuestreo (--resample_from, o el que está junto a --discovery_from o --background_input). workflows/scoring.nf: sin reconstrucción del nulo (CAAS_CORE_MERGE importado) ni importación de caas_perms.rds; el nulo es el vivo o el centinela NO_FILE. workflows/enrichment.nf: caas_perms_resolved sin parámetro. conf/scoring.config, workflows/help.nf, README.md: fuera caas_perms_file y caas_pos_detail_file. fcs_enrich.R: el aviso de nulo obsoleto remite a repetir CAAS_CORE_MERGE. GUI: gui/models/precomputed.py (sin caas_perms_file ni caas_pos_detail_file), run_single.sh.j2 (bloques PRECOMP y reuse_exploratory, claves de params.json), comentario de modules.py. Se conservan caas_cycle_roster_file, caas_pos_cycle_caas_file, caas_pos_sample_file y caas_pos_quantiles_file. gui/generation/report_registry.py:323 se deja: ese caas_perms_file es un parámetro del Rmd de informe que lee el caas_perms.rds publicado, no una vía de reutilización de la tubería.
  * POSENRICH se omite con aviso cuando no hay nulo por posición (workflows/enrichment.nf: posenrich_null_available = nulo vivo de CAAS_CORE_MERGE o --caas_pos_cycle_caas_file existente; sin él, log.warn y el módulo no se ejecuta; los consumidores posteriores ya toleraban posenrich_out nulo). Hoy el nulo por posición existe en las corridas normales, así que POSENRICH sigue ejecutándose en ellas: se desactivará de forma efectiva cuando WP5 retire perm_pos_cycle_caas.tsv.gz. docs/POSENRICH_PORTE_EVALUACION.md revisado (ver abajo).
  * Limpieza: cleanup_policy sustituido por remove_intermediates y remove_work (booleanos independientes) en conf/common.config, workflows/help.nf, caas_permulation.nf (removeIntermediatesOn(): poda de exports y mv de trozos), gui/models/runtime.py, gui/widgets/tabs/runtime_tab.py (dos casillas con descripción propia), gui/generation/validate.py (comprobación del combo retirada), run_single.sh.j2 (params.json con ambos valores; la retirada de work/ depende solo de remove_work), sbatch_array.sh.j2, README.md y las 9 plantillas de gui/templates que tenían clave (valores según migrate()). gui/models/serialization.py: migrate() convierte cleanup_policy (none a F,F; prune_intermediates a T,F; remove_work a T,T; las claves nuevas explícitas ganan; un valor desconocido avisa y deja ambas en F); no se sube SCHEMA_VERSION.
  * Pruebas: PhyloPhere_validation/wiring/test_cleanup_options.py (nuevo, 12 pruebas: migrate() con los tres valores, sin clave, claves explícitas y valor desconocido; scripts generados con las cuatro combinaciones; la retirada de work/ solo por remove_work; ausencia de los parámetros retirados en los ficheros versionados). Adaptados mini_core.nf, mini_core_chain.nf, test_wiring.py y test_core_batched.py a la entrada única (se retiran las dos pruebas de exports reutilizados).
- Comprobaciones ejecutadas (comando, conjunto de datos, tolerancia, resultado):
  * PEPC c4, con el código de WP2 (scripts renderizados con render_tier1_scripts.py desde validation/tier1/templates/tier1_pepc_c4.json, que aún lleva cleanup_policy y pasa por migrate()): compare_null.py contra golden_581b999/pepc_c4_tolerance/results/c4_complete, tolerancia 1e-12: 11 passed, 0 failed, delta 0.0 y 0 celdas no idénticas; procesos del trace idénticos (34 y 34, mismos nombres y recuentos). Se ejecutaron también los dos rasgos (c4_phenotypic no se compara).
  * Toy de cancer (45 genes, 100 ciclos), con la caché de ASR de la referencia: 16 passed, 0 failed, 5 skipped, delta 0.0 en gene_cycle_scores, perm_pos_*, caas_perms.rds (all.equal TRUE), b0/, scoring/position_scores.tsv, fcs_stats*.tsv y las tablas FCS (p.perm incluido); procesos idénticos (33 y 33).
  * Primera corrida del toy de cancer con caché de ASR vacía: 11 tablas distintas (delta máximo 0,088). Causa medida: codeml no es reproducible entre ejecuciones (ver WP0_linea_base.md, punto 5); no es un efecto del código de WP2. Registros en wp2_check/cancer_run1_fresh_asr*.
  * Suite wiring (10 ficheros sin Nextflow: test_gui_*, test_param_audit, test_retired_names, test_tier1_render): 55 failed, 98 passed en el código de partida (a92fd2d, checkout limpio) y 55 failed, 98 passed con WP2, el mismo conjunto (lista en golden_581b999/pytest_baseline/wiring_gui_failing_a92fd2d.txt). Los 12 tests nuevos de test_cleanup_options.py pasan.
- Mediciones (recurso antes y después): PEPC c4 y cancer toy sin cambio de procesos ni de tablas. Sin medición de disco de remove_intermediates / remove_work en corrida real (queda para WP8).
- Desviaciones del plan y decisiones tomadas (con su justificación):
  * No he editado validation/tier*/templates (plantillas y scripts de validation/ los gestiona el usuario): conservan cleanup_policy, que migrate() convierte al abrirlas. El usuario debe sincronizarlas si quiere las claves nuevas en esos ficheros.
  * wiring/test_wiring.py y los tests de Nextflow (mini_*.nf) se adaptaron pero no se ejecutaron (necesitan Nextflow y ~20 min); la lista de procesos de las dos corridas reales sí coincide.
  * POSENRICH: se interpreta "omitido con aviso cuando no hay nulo por posición" como condición del módulo, no como apagado incondicional; así WP2 sigue siendo borrado sin efecto numérico.
  * Marine toy: la comparación del plan (marine toy y PEPC) se hizo con PEPC y el toy de cancer; marine toy sigue pendiente (disco externo).
- Pendiente que pasa a otra fase: decisión sobre remove_intermediates (WP8); medir remove_work en una corrida real (WP8); POSENRICH deja de ejecutarse de forma efectiva al retirar perm_pos_cycle_caas.tsv.gz (WP5); sincronizar validation/tier*/templates y marine toy (usuario).
```
