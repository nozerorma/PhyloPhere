# WP0. Línea base congelada

Estado: CASI CERRADA (2026-10-21). Hecho: huellas del cluster, suite, `compare_null.py`, PEPC y cancer_b0_toy regenerados; pendiente solo marine toy (disco externo). Depende de: commit del usuario (hecho: `581b999`, rama `remove_caap`).
Conversación independiente. Antes de empezar, leer `docs/fases/README.md`.

## Objetivo

Dejar referencias verificables del comportamiento actual (nulo, scoring, FCS) contra las que se comparará cada fase posterior, y la suite de pruebas existente ejecutada y registrada.

## Estado de partida (2026-10-21)

- Commit `581b999` del usuario. Checkout del cluster en `d246e82` con una entrada sin confirmar (`~/scratch/0.Phylophere/.tmp/golden/checkout_revision.txt`). Entre `d246e82` y `581b999` los ficheros de la cadena del nulo no cambian (`caas_permulation.nf`, `gene_wrapper.py`, `reaggregate_perm_scores.py`, `scoring_caas_perms.R`, `fcs_enrich.R`, `fcs.nf`, `rer_gene_lists.nf`). Cambian `permulations.R`, `lean_contrast_selector.R` y `ct_resample.nf` (solo la opción `pair_center_gate`, que con `none` es un no-op: el bloque se ejecuta únicamente si el valor es distinto de `none`) y `scoring_compute.R` (solo la lectura de la acumulación exacta). Por tanto las salidas del cluster del 19/10 valen como referencia de `581b999`.
- Ya hecho: copia congelada del código del FCS en `PhyloPhere_validation/unification/legacy_fcs/` (`fcs_enrich.R`, `fcs_compute.R`, `percentile_flags.R`, `SHA256SUMS`, `FROZEN_FROM.txt`), idéntica al commit (verificado); huellas SHA-256 de las salidas de neoplasia y maligno completo (incluido `strict/`) en `~/scratch/0.Phylophere/.tmp/golden/{neoplasia,malignant}_manifest.sha256` (trabajo 7568379; los directorios de trozos por gen se listan en `*_shards.txt` y no se hashean).

## Tareas

1. **Referencias locales** en `/home/miguel/IBE-UPF/PhD/PhyloPhere_validation/golden_581b999/`:
   - Marine toy (121 genes, 100 ciclos): `validation/tier2/output/marine_toy/results/{marine,control}_complete` y `strict`. Generadas el 2026-10-20 09:33 con HEAD `d246e82` (informe `validation/tier2/reports/marine_toy_resources.md`; si el árbol tenía cambios no consta). Candidata a referencia; confirmar.
   - PEPC tier 1 (1 gen, 1000 ciclos, 100 hipótesis): `validation/tier1/output/pepc/results/c4_complete`, generada el 2026-10-15. Puede ser anterior a cambios del nulo (comparar con `git log` por fecha); regenerar si no se puede demostrar equivalencia (`validation/tier1/input/pepc/scripts/run_tier1_pepc_*_complete.sh`).
   - `validation/unification/cancer_b0_toy` (45 genes, 100 ciclos): de código anterior a la unificación; regenerar con `gui/templates/cancer_no_prune_multi_toy.json` si se quiere usar.
   - Para cada una, anotar el comando y el commit con que se generó. Decidir en la conversación si se regenera en local o en el cluster.
2. **Huellas del cluster:** traer a local una copia de `~/scratch/0.Phylophere/.tmp/golden/` y guardarla en `/home/miguel/IBE-UPF/PhD/PhyloPhere_validation/golden_581b999/cluster_manifest/`. Si hace falta una copia física de salidas (por ejemplo `fcs_stats.tsv`, `caas_perms.rds`, un `.continuous.perms.rds` de RER y las tablas FCS antiguas) hacerla por `sbatch` hacia `~/scratch/0.Phylophere/.tmp/golden/` y anotar el tamaño.
3. **Suite existente:** `cd /home/miguel/IBE-UPF/PhD/PhyloPhere_validation && RUN_SLOW=1 TMPDIR=/home/miguel/IBE-UPF/PhD/PhyloPhere_validation/.tmp python3 -m pytest -q unification resources` (codeml y R en el `PATH`; el submódulo `wiring/` necesita Nextflow y tarda ~20 min). Anotar qué pasa y qué no, para distinguir luego los fallos nuevos.
4. **Ampliar `compare_null.py`** (`PhyloPhere_validation/unification/`; tablas y `--extra`) a `scoring/position_scores.tsv` y a las tablas FCS, con su prueba (`test_compare_null.py`).
5. **Confirmar la copia del FCS:** `cd PhyloPhere_validation/unification/legacy_fcs && sha256sum -c SHA256SUMS`.
6. **Dejar escrito el índice** de referencias (qué, dónde, commit, comando) al final de este fichero.

## Interfaz hacia las fases siguientes

- `/home/miguel/IBE-UPF/PhD/PhyloPhere_validation/golden_581b999/` con índice; `legacy_fcs/` (código de FCS de referencia); `compare_null.py` extendido; lista de pruebas de la suite que pasan o fallan hoy.

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
- Conversación y fechas: WP0, 2026-10-21 (15:40 a 16:30 CEST), local, rama remove_caap.
- Commit de partida / commit de cierre: 581b999 (HEAD e5b924b solo añade docs; el código de la tubería es idéntico) / a92fd2d. Las referencias nuevas se generaron con 581b999 MÁS la corrección de una línea de accum_gene_lists.nf (commit a92fd2d); el commit de cierre a92fd2d la incluye (hecho por Claude a petición del usuario; solo la corrección y docs/fases, sin validation/).
- Qué se hizo (ficheros y cambios):
  * Huellas del cluster copiadas a PhyloPhere_validation/golden_581b999/cluster_manifest/ (51 KB, scp de 7 ficheros; checkout del cluster d246e82, 1 entrada sin confirmar).
  * legacy_fcs verificado: sha256sum -c SHA256SUMS, 3 OK.
  * unification/compare_null.py: nuevas listas SCORING_TABLES (scoring/position_scores.tsv) y FCS_TABLES (scoring/fcs_stats{,_rer,_accum}.tsv y las tres tablas fcs_results de fcs/ y de scoring/rer/), comparadas siempre (se omiten si faltan en ambos); --extra se mantiene. Prueba nueva test_position_scores_and_fcs_tables_are_compared_by_default; test_a_stored_run_equals_itself... copia ahora también esas tablas.
  * Suite ejecutada y su línea base guardada en golden_581b999/pytest_baseline/.
  * Rutas ~/.cc/golden_* y TMPDIR corregidas en README.md, WP0 y PLAN_APROBADO (partición IBE).
  * Único cambio de código (autorizado por el usuario): subworkflows/CT_ACCUMULATION/accum_gene_lists.nf línea 52, '^PValue$' pasa a '^PValue\$' (escape del $ dentro del bloque de comillas triples).
  * PEPC regenerado en golden_581b999/pepc_c4_tolerance/ (c4_complete y c4_phenotypic_complete; scripts renderizados con render_tier1_scripts.py, copiados a la raíz del repo para ejecutarlos y retirados después): ambas corridas exit 0 (34 y 35 tareas, 0.7 CPU-h en total, 690 MB).
- Comprobaciones ejecutadas (comando, conjunto de datos, tolerancia, resultado):
  * cd PhyloPhere_validation && RUN_SLOW=1 TMPDIR=$PWD/.tmp python -m pytest -rA -q --continue-on-collection-errors unification resources (entorno phylophere, codeml y R en PATH; wiring/ no ejecutado): 290 passed, 252 failed, 17 errors, 273 s. Listas exactas: pytest_baseline/{passing,failing}_tests.txt; log completo pytest_full.log.
    Causas medidas (primer error de cada fallo): NameError SCHEMES/US (test_kernel_fuzz, test_b0_discovery, test_discovery_kernel: tests escritos para el código con CAAP/GS); caas_id(caap_group=...) y campo caap_group en tablas; GeneUnit con 5 argumentos (test_postproc); process_all_genes_perms sin chunk_threshold; "tax_id shared by several species" (Eleocharis_baldwinii_*; test_asr_loader, test_asr_concurrency, test_species_mapping, test_driver...; curación de especies); pool_domains sin pares de diseño (test_modal_tie, test_pool_domains_order); 2 errores de recolección por importar la acumulación eliminada (randomize.py: test_accum_null_counting, test_accum_occupancy_universe). Leído en los mensajes de error, no analizado test a test: los fallos son de interfaz, no de valores numéricos. test_compare_null.py: 10 passed.
  * compare_null.py, salida PEPC del 15/10 (HEAD 9ae9bd2, validation/tier1/output/pepc/results/c4_complete) contra la regenerada (c4_complete): 11 passed, 0 failed, 10 skipped (tablas ausentes en ambas), tolerancia 1e-12, max_abs_delta 0.0 y 0 celdas no idénticas bit a bit en gene_cycle_scores (1000 filas), perm_pos_sample, perm_pos_quantiles, perm_pos_cycle_caas (40.693), caas_perms.rds (all.equal TRUE), la rebanada b0/, scoring/position_scores.tsv (36 filas) y scoring/fcs_stats.tsv. Es una equivalencia medida del nulo y del scoring de PEPC c4 entre 9ae9bd2 y 581b999 + corrección (no cubre c4_phenotypic, strict ni FADE).
  * compare_null.py contra marine_complete consigo mismo: 19 passed, 2 skipped, delta 0 (solo comprueba que las tablas se leen y comparan; no es una prueba de equivalencia).
- Mediciones (recurso antes y después): suite 273 s en local. Disco: golden_581b999 < 30 MB.
- Desviaciones del plan y decisiones tomadas (con su justificación):
  * Los scripts run_tier1_pepc_* del árbol eran anteriores a los parámetros actuales: se renderizaron de nuevo. El primer lanzamiento falló en la compilación de Nextflow 25.10.4: accum_gene_lists.nf línea 52, grep('^PValue$', ...) dentro de un bloque de comillas triples (el $' se interpreta como inicio de GString). La línea entró en 581b999 (en d246e82 era 'PValueEmpirical'). El módulo se importa aunque accumulation esté desactivada, así que con 581b999 sin corrección ningún run local arrancaba. Tras la corrección el pipeline compila y corre; búsqueda de otro $' o $" sin escapar en los .nf versionados: ninguno. Log del fallo en pepc_c4_tolerance/failed_attempt/.
  * La salida PEPC del 15/10 no se usa como referencia por sí misma (HEAD 9ae9bd2), pero la comparación anterior demuestra que coincide con la regenerada.
  * Marine toy: pendiente por decisión del usuario (el disco externo con los alineamientos toy, /media/miguel/external/phylophere_tier2/input/alignments/toy, no está disponible). Las salidas del 20/10 (HEAD d246e82, relanzadas con -resume, estado del árbol no registrado) no se copian como referencia hasta confirmarlo.
  * cancer_b0_toy: regenerado en local con a92fd2d (ver índice). Parámetros distintos de los del cluster: toy_n 45, AMI y STRING apagados (el fondo de 31 genes no tiene aristas de STRING a 700) y VEP y POSENRICH apagados. Lanzamiento con -resume tras el fallo de DOMINO. Que los alineamientos phylip-relaxed sean los mismos que los BMGE en fasta del cluster lo afirma el usuario; no lo he verificado.
- Pendiente que pasa a otra fase: ver "Pendiente de WP0" y "Qué comprobar primero" abajo.
```

## Índice de referencias

Directorio: `/home/miguel/IBE-UPF/PhD/PhyloPhere_validation/golden_581b999/` (fuera del repositorio).

| Qué | Dónde | Commit / origen | Estado |
|---|---|---|---|
| Huellas SHA-256 de neoplasia y maligno completo (con `strict/`) y lista de directorios de trozos | `cluster_manifest/{neoplasia,malignant}_manifest.sha256`, `*_shards.txt`, `checkout_revision.txt`, `slurm-7568379.*` | Cluster, checkout `d246e82` (+1 entrada sin confirmar); trabajo 7568379, 2026-10-21 15:32; cadena del nulo idéntica a `581b999` | Hecho |
| Código del FCS de referencia | `../unification/legacy_fcs/` (`fcs_enrich.R`, `fcs_compute.R`, `percentile_flags.R`, `SHA256SUMS`) | `581b999`; `sha256sum -c SHA256SUMS` OK | Hecho |
| Línea base de la suite | `pytest_baseline/{passing_tests,failing_tests}.txt`, `pytest_full.log` | `581b999`; comando en el registro | Hecho |
| WP9-A1 (estadístico génico) | `a1_gene_statistic/` | `581b999` | Hecho (preparación) |
| PEPC tier 1, c4 y c4_phenotypic, tolerante, 1000 ciclos | `pepc_c4_tolerance/results/{c4,c4_phenotypic}_complete` (+ scripts renderizados, `PROVENANCE.txt`, `run.log`, `failed_attempt/`) | `581b999` + corrección de `accum_gene_lists.nf` (commit de cierre), 2026-10-21 | Hecho; igual bit a bit a la salida del 15/10 en c4 |
| Marine toy (121 genes, 100 ciclos, 4 corridas) | `validation/tier2/output/marine_toy/results/` (no copiado) | `d246e82`, 2026-10-20, relanzamientos con `-resume`; lanzadores en la raíz del repo (`run_tier2_marine_toy_{local,single}_{complete,strict}.sh`, ignorados por git) | PENDIENTE (falta el disco externo) |
| cancer_b0_toy local (template `gui/templates/cancer_no_prune_multi_toy.json` con rutas locales, `toy_n` 45 = 31 genes tras filtros, 100 ciclos, VEP, POSENRICH, AMI y STRING apagados, alineamientos `Primate_alignments` en `phylip-relaxed`) | `cancer_b0_toy/results/neoplasia_prevalence_toy_complete` (+ `build_project.py`, `project.json`, scripts, `PROVENANCE.txt`, `run.log`, `run_attempt1_domino_failed.log`) | `a92fd2d`, 2026-10-21; el segundo lanzamiento usó `-resume` sobre el primero, que solo falló en `DOMINO_BUILD_NETWORK` | Hecho (exit 0; 272 MB). Sin referencia previa comparable atada a commit; la carpeta `unification/cancer_b0_toy` es anterior a la unificación |

Comparador: `PhyloPhere_validation/unification/compare_null.py` (tablas del nulo, `b0/`, `caas_perms.rds`, `scoring/position_scores.tsv`, tablas FCS; `--extra` para más), tolerancia 1e-12.

## Pendiente de WP0

1. Marine toy: con el disco montado, regenerar las cuatro corridas (acumulado medido el 20/10: 51 min lenient marine, 1 h 16 min lenient control, 10 min cada strict) o aceptar las salidas del 20/10 como provisionales.

4. Hallazgo: DOMINO debe saltarse el paso (skip silencioso) sin aristas suficientes (detalle en `README.md`, sección de hallazgos pendientes). En el toy de cancer se evitó con `scoring_ami` y `scoring_string` apagados.

## Qué comprobar primero en WP1 y WP2 (pueden ir en paralelo)

- La tubería arranca con Nextflow 25.10.4 solo con la corrección de `accum_gene_lists.nf` (incluida en el commit de cierre de WP0): partir de ese commit. Los lanzadores generados deben estar en la raíz del repo (el runner single toma REPO_DIR de su ubicación).
- Comparar con `compare_null.py` contra `pepc_c4_tolerance/results/c4_complete` (tolerancia 1e-12) tras cualquier cambio del nulo.
- WP1 (FCS/RER, `fcs_enrich.R`, `fcs_compute.R`): comparar contra `legacy_fcs/` (código) y, para números, contra las huellas del cluster; en local aún no hay salida FCS de referencia atada a `581b999`.
- WP2 (exports, POSENRICH, limpieza): los tests que tocan estas partes ya fallan hoy por la interfaz (lista en `pytest_baseline/failing_tests.txt`); distinguir los fallos nuevos con esa lista, no con verde/rojo.
