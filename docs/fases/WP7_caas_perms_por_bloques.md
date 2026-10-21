# WP7. Fase 5: `caas_perms` por bloques y exportación opcional

Estado: pendiente. Depende de: WP5.
Conversación independiente. Antes de empezar, leer `docs/fases/README.md`.

## Objetivo

Sustituir el fichero único `caas_perms.rds` (matrices genes × N) por una matriz por bloque de ciclos, con el roster y la estampa, y dejar la reconstrucción del RDS único como exportación opcional.

## Decisión del usuario y razones

Bloques canónicos + exportación opcional (`--export-single`). Razones discutidas: (1) con bloques el FCS ya no necesita el objeto entero; (2) a 10.000 ciclos serían ~3 GB que un paso serie tendría que escribir y recargar, justo el problema que se quita; (3) `scoring_compute.R:520` solo lee el roster, que sale de `cycle_roster.txt`; (4) las vías de reutilización por parámetro desaparecen (WP2). Coste de no conservar el RDS único: se rompe cualquier script externo que lo lea (incluidos los scripts de control en `~/.cc` y `adhoc_scripts/`); por eso la exportación a demanda.

## Estado de partida

- Consumidores hoy: `fcs_enrich.R:651` (lee el RDS completo, usa solo `caas_corStat_byrank`), `scoring_compute.R:520` (roster). Sin lectores del eje ASR (WP5 deja de construirlo).
- Scripts de `adhoc_scripts/` que leen `caas_perms.rds` por la ruta del run: `fcs_calibrated_ranking_control.sbatch`, `fcs_rawmax_ranking_control.sbatch`, `fcs_null_as_observed.sbatch`, `fcs_scaling.sbatch`.
- Medido: `scoring_caas_perms.R` ocupa ~3,4 MB por ciclo (lee la tabla ancha entera con `read_tsv` y construye matrices densas).

## Alcance (ficheros)

`subworkflows/SCORING/local/src/scoring_caas_perms.R`, `subworkflows/CT/caas_permulation.nf` (salidas publicadas del merge), `fcs_enrich.R` (lectura), `scoring_compute.R` (roster), `gui/models/precomputed.py` (si queda algún campo), documentación, scripts de control de `adhoc_scripts/` que lean el RDS.

## Cambios

1. Formato canónico: `blocks/blk_k.matrices.rds` + `cycle_roster.txt` + estampa (`gene_stat`, y el estadístico base elegido en WP9).
2. `scoring_caas_perms.R` lee la tabla por trozos y rellena solo las matrices de CAAS (sin eje ASR).
3. `--export-single`: reconstruye un RDS único bajo demanda (lectura por bloques y escritura serie).
4. Actualizar los scripts de control y la documentación a la nueva ruta o a la exportación.

## Interfaz

Los consumidores leen por bloque; el roster es la fuente de N.

## Validación (criterios de paso)

- El RDS único exportado es idéntico (`all.equal`, tolerancia 1e-12, mismas filas y columnas) al `caas_perms.rds` de referencia en marine toy y PEPC.
- FCS y SCORING con bloques dan los mismos resultados que con el RDS único.
- Medir RSS de `scoring_caas_perms.R` por bloque: debe depender de B y no de N.

## Fuera de alcance

La partición por bloques del FCS (WP6), decisiones de limpieza (WP8).

## Riesgos y puntos a verificar

- Roster y estampa: un bloque que falta o un ciclo sin detección debe seguir contando en N.

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
