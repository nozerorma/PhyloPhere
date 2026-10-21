# WP0. Línea base congelada

Estado: pendiente (preparación hecha). Depende de: commit del usuario (hecho: `581b999`, rama `remove_caap`).
Conversación independiente. Antes de empezar, leer `docs/fases/README.md`.

## Objetivo

Dejar referencias verificables del comportamiento actual (nulo, scoring, FCS) contra las que se comparará cada fase posterior, y la suite de pruebas existente ejecutada y registrada.

## Estado de partida (2026-10-21)

- Commit `581b999` del usuario. Checkout del cluster en `d246e82` con una entrada sin confirmar (`~/scratch/0.Phylophere/.tmp/golden/checkout_revision.txt`). Entre `d246e82` y `581b999` los ficheros de la cadena del nulo no cambian (`caas_permulation.nf`, `gene_wrapper.py`, `reaggregate_perm_scores.py`, `scoring_caas_perms.R`, `fcs_enrich.R`, `fcs.nf`, `rer_gene_lists.nf`). Cambian `permulations.R`, `lean_contrast_selector.R` y `ct_resample.nf` (solo la opción `pair_center_gate`, que con `none` es un no-op: el bloque se ejecuta únicamente si el valor es distinto de `none`) y `scoring_compute.R` (solo la lectura de la acumulación exacta). Por tanto las salidas del cluster del 19/10 valen como referencia de `581b999`.
- Ya hecho: copia congelada del código del FCS en `PhyloPhere_validation/unification/legacy_fcs/` (`fcs_enrich.R`, `fcs_compute.R`, `percentile_flags.R`, `SHA256SUMS`, `FROZEN_FROM.txt`), idéntica al commit (verificado); huellas SHA-256 de las salidas de neoplasia y maligno completo (incluido `strict/`) en `~/scratch/0.Phylophere/.tmp/golden/{neoplasia,malignant}_manifest.sha256` (trabajo 7568379; los directorios de trozos por gen se listan en `*_shards.txt` y no se hashean).

## Tareas

1. **Referencias locales** en `~/.cc/golden_581b999/`:
   - Marine toy (121 genes, 100 ciclos): `validation/tier2/output/marine_toy/results/{marine,control}_complete` y `strict`. Generadas el 2026-10-20 09:33 con HEAD `d246e82` (informe `validation/tier2/reports/marine_toy_resources.md`; si el árbol tenía cambios no consta). Candidata a referencia; confirmar.
   - PEPC tier 1 (1 gen, 1000 ciclos, 100 hipótesis): `validation/tier1/output/pepc/results/c4_complete`, generada el 2026-10-15. Puede ser anterior a cambios del nulo (comparar con `git log` por fecha); regenerar si no se puede demostrar equivalencia (`validation/tier1/input/pepc/scripts/run_tier1_pepc_*_complete.sh`).
   - `validation/unification/cancer_b0_toy` (45 genes, 100 ciclos): de código anterior a la unificación; regenerar con `gui/templates/cancer_no_prune_multi_toy.json` si se quiere usar.
   - Para cada una, anotar el comando y el commit con que se generó. Decidir en la conversación si se regenera en local o en el cluster.
2. **Huellas del cluster:** traer a local una copia de `~/scratch/0.Phylophere/.tmp/golden/` y guardarla en `~/.cc/golden_581b999/cluster_manifest/`. Si hace falta una copia física de salidas (por ejemplo `fcs_stats.tsv`, `caas_perms.rds`, un `.continuous.perms.rds` de RER y las tablas FCS antiguas) hacerla por `sbatch` hacia `~/scratch/0.Phylophere/.tmp/golden/` y anotar el tamaño.
3. **Suite existente:** `cd /home/miguel/IBE-UPF/PhD/PhyloPhere_validation && RUN_SLOW=1 TMPDIR=~/.cc/tmp python3 -m pytest -q unification resources` (codeml y R en el `PATH`; el submódulo `wiring/` necesita Nextflow y tarda ~20 min). Anotar qué pasa y qué no, para distinguir luego los fallos nuevos.
4. **Ampliar `compare_null.py`** (`PhyloPhere_validation/unification/`; tablas y `--extra`) a `scoring/position_scores.tsv` y a las tablas FCS, con su prueba (`test_compare_null.py`).
5. **Confirmar la copia del FCS:** `cd PhyloPhere_validation/unification/legacy_fcs && sha256sum -c SHA256SUMS`.
6. **Dejar escrito el índice** de referencias (qué, dónde, commit, comando) al final de este fichero.

## Interfaz hacia las fases siguientes

- `~/.cc/golden_581b999/` con índice; `legacy_fcs/` (código de FCS de referencia); `compare_null.py` extendido; lista de pruebas de la suite que pasan o fallan hoy.

## Validación (criterios de cierre)

- Referencias con origen documentado (commit y comando) y comprobables por huella.
- Suite ejecutada y su resultado anotado.
- `compare_null.py` extendido con prueba propia que pasa.

## Fuera de alcance

Cambios de código del pipeline. Cualquier cómputo del cluster solo por `sbatch`.

## Riesgos y puntos a verificar

- Ya se detectó que `test_null_consumers`, `test_scores` y `test_compare_null` asumen `gene_cycle_scores.tsv`: se portan en WP5.
- Si una referencia no puede atarse a un commit, no sirve como referencia (regenerarla).

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
