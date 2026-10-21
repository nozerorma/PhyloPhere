# WP2. Eliminar la reutilización de exports, desactivar POSENRICH y separar las opciones de limpieza

Estado: pendiente. Depende de: WP0. Independiente de WP1.
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
- Conversación y fechas:
- Commit de partida / commit de cierre:
- Qué se hizo (ficheros y cambios):
- Comprobaciones ejecutadas (comando, conjunto de datos, tolerancia, resultado):
- Mediciones (recurso antes y después):
- Desviaciones del plan y decisiones tomadas (con su justificación):
- Pendiente que pasa a otra fase:
```
