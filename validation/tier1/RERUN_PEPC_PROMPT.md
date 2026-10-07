# Prompt: repetir Tier 1 PEPC con el nulo construido dentro del propio run

Copiar el bloque siguiente en una sesión nueva de Claude Code abierta en `/home/miguel/IBE-UPF/PhD/PhyloPhere`.

````text
Trabajas en el repositorio PhyloPhere (rama scoring_v2, a partir del commit 8faee85). Lee antes `git log 5bafadb..HEAD`,
`docs/scoring_v2_p_emp.md` y `validation/README.md` (sección "Tier 1 procedure").

OBJETIVO
Rehacer por completo el caso Tier 1 PEPC (rasgos `c4` genotípico y `c4_phenotypic`) con el pipeline actual. El nulo de
permulación de cada run (caas_perms.rds, perm_pos_cycle_caas.tsv.gz, gene_cycle_scores.tsv) debe construirse dentro del propio
run, a partir de sus propias permulaciones: ningún nulo, descubrimiento, caché de ASR ni resultado de FADE importado de otro run.
Después, calcular las tablas de posiciones con las columnas nuevas y comprobar la calibración del p empírico.

QUE HA CAMBIADO EN EL PIPELINE (resumen)
- caas_score_aggregation = cumulative por defecto (suma de asr_path_score sobre los 5 esquemas / 5; el esquema que no detecta
  cuenta 0). perm_pos_cycle_caas.tsv.gz lleva la columna score_aggregation y scoring_compute.R se detiene si el nulo y el
  observado usan reglas distintas.
- p.adj_sam y el p de gen poolado (scoring_gene_perm_pooled) ya no existen.
- p.emp = 1 cuando el score observado es 0. gene_caas_pperm_adj es BH sobre todo el universo del nulo.
- Columnas nuevas de position_scores.tsv: p.emp_fact, p.adj_bh_fact (p factorizado y su BH), amino_encoded, n_conserved_pairs;
  de gene_scores.tsv: gene_caas_pperm_fact{,_top,_bottom} y gene_caas_pperm_fact_adj{,_top,_bottom}.
- fcs_caas_score = raw por defecto (fact reordena los rankings CAAS del FCS por -log10 del p de gen factorizado).

REGLAS
- No hagas commit salvo que se pida. No uses /tmp: usa el directorio scratchpad de la sesión para intermedios.
- No toques `validation/tier1/previous_work/` ni `PhyloPhere_validation` salvo lo indicado en el paso 9.
- Estilo de los informes: descriptivo del estado actual, sin historia ("antes", "ahora", fechas, commits), sin rayas largas.
- Si algo falla, diagnostica con los logs (run.log, .nextflow.log, work/) y cuéntalo; no parchees el pipeline en silencio.

PASOS
1. Estado. `git status --short` (solo debe aparecer gui/templates/carn.json sin seguimiento) y `git log -1 --format=%h`.

2. Guardar lo mínimo del run anterior y borrarlo. Si existe `validation/tier1/output/pepc`:
   a) copia a scratchpad los position_scores.tsv, gene_scores.tsv y fcs_stats.tsv de results/c4_complete/scoring y
      results/c4_phenotypic_complete/scoring (servirán para comparar antes y después);
   b) muestra `du -sh validation/tier1/output/pepc` (unos 427 MB: results, work, asr_cache y los scripts generados);
   c) borra solo ese directorio: `rm -rf validation/tier1/output/pepc`.
   No borres nada más (previous_work, reports, templates, input).

3. Fixture. Comprueba con `ls -a validation/tier1/input` que existen `pepc/` (genotípico) y `.pepc_phenotypic/`, y que están los
   ficheros que lista `validation/tier1/input/pepc/README.md` (align, align_cds, tree.nwk, my_traits.tsv, gene_ensembl.tsv,
   taxid.tsv...). Si falta algo, reconstrúyelo con los scripts de `validation/tier1/input/pepc/scripts/` (build.py, build_cds.py).

4. Renderizar los scripts (directorios nuevos de resultados, work y caché de ASR):
   python3 validation/tier1/scripts/render_tier1_scripts.py --template gui/templates/tier1_pepc_c4.json \
       --outdir validation/tier1/output/pepc --evidence-top-n 30
   Revisa el JSON de parámetros del script generado y comprueba: caas_score_aggregation = cumulative, ct_tool_discovery,
   ct_tool_resample y caas_permulation_enrichment = true, caas_full_perms = 1000, y que NO están definidos caas_perms_file,
   caas_pos_detail_file, caas_pos_cycle_caas_file, discovery_from, fade_json_dir_top, fade_json_dir_bottom ni ningún
   *_perms_file/*_null_file. Si alguno está, páralo y avisa.

5. Ejecutar (entorno `phylophere`; unas 1.1 horas de CPU en el run anterior):
   cd validation/tier1/output/pepc
   ( date '+%F %T' > started.at; RESUME=0 bash run_tier1_pepc_local_complete.sh > run.log 2>&1; echo $? > exit.code ) &
   Espera a que termine. Éxito: exit.code = 0 y la ejecución de ambos rasgos completa (c4_complete y c4_phenotypic_complete en results/).

6. Comprobar que el nulo es el del propio run. Para cada rasgo, en results/<rasgo>_complete/caas_permulation/:
   - caas_perms.rds (genes x 1000 ciclos), perm_pos_cycle_caas.tsv.gz y gene_cycle_scores.tsv existen y sus fechas caen dentro
     de la ventana del run (started.at ... fin);
   - perm_pos_cycle_caas.tsv.gz tiene la columna score_aggregation con el valor cumulative y 1000 ciclos distintos entre
     "b_1".."b_1000" (más los que no dejaron filas, que cuentan en N);
   - scoring/position_scores.tsv contiene p.emp, p.adj_bh, p.emp_fact, p.adj_bh_fact, amino_encoded, n_conserved_pairs y NO
     contiene p.adj_sam; scoring/gene_scores.tsv contiene gene_caas_pperm_fact* y no contiene *pooled*.

7. Análisis de posiciones (escribe un script reproducible en validation/tier1/scripts/, p. ej. pepc_pvalue_tables.py):
   a) Conjunto de verdad: validation/truthsets/tier1/pepc_c4.sites.tsv (10 sitios, numeración de maíz PEPC1 = Position + 1).
      Para cada rasgo y sitio: CAAS_score, rango por score, p.emp, p.adj_bh, p.emp_fact, p.adj_bh_fact. Cuántos sitios de
      verdad quedan a p.adj_bh <= 0.05 y a p.adj_bh_fact <= 0.05; cuántas posiciones que no son de verdad quedan también
      por debajo, con sus valores.
   b) Distribución: número de posiciones con p.emp y p.emp_fact en los tramos <=0.001, <=0.01, <=0.05, y en el suelo
      (p.emp = 1/(N+1)).
   c) Comparación con el run anterior (copias del paso 2): puestos de las posiciones, scores y p. Aclara que la agregación
      cambió de media a acumulativo y por eso los scores no son comparables uno a uno.

8. Calibración con el nulo propio. Para cada rasgo, tomando cada ciclo nulo como observado frente a los demás
   (perm_pos_cycle_caas.tsv.gz, estadístico = máximo sobre lados de caas_score):
   - p.emp pseudo-observado = (ciclos que detectan la posición con score >= el suyo, sin contarse, + 1) / N;
   - p factorizado pseudo-observado: la regla de scoring_compute.R (helpers .fact_fit, .fact_assign, .fact_p) con la clase
     de propensión de la posición y la reserva de la clase sin el ciclo evaluado;
   - la tasa de pares (posición, ciclo) con p <= alpha (alpha = 0.001, 0.01, 0.05) sobre toda la familia (posiciones
     detectadas por el nulo o el observado; las no detectadas cuentan p = 1) debe ser <= alpha, también por clases de
     propensión (nd <= 5, 6-20, 21-100, > 100). Reporta cualquier violación.

9. Informes y herramientas (descriptivo, sin historia):
   - reescribe validation/tier1/reports/pepc_results.md y pepc_genotypic_vs_phenotypic.md con los números nuevos: quita todo
     lo de p.adj_sam (columnas, conteos, calibración) e incluye p.emp_fact y p.adj_bh_fact;
   - actualiza validation/tier1/scripts/compare_pepc_runs.py (quita la columna sam, añade p.emp_fact) y
     /home/miguel/IBE-UPF/PhD/PhyloPhere_validation/wiring/test_compare_pepc_runs.py; ejecuta esa prueba y la suite
     `unification` de PhyloPhere_validation (PATH con el entorno phylophere; README de esa carpeta);
   - si cambia el uso de recursos, regenera las tablas con validation/tier1/scripts/pepc_resources_tables.py.

10. OPCIONAL, solo si queda tiempo: repetir el run con caas_score_aggregation = mean en validation/tier1/output/pepc_mean
    (copia el template, fija modules.scoring.caas_score_aggregation = "mean", renderiza y ejecuta como arriba) para aislar el
    efecto de la agregación sobre los 10 sitios de verdad.

ENTREGA
Un resumen corto con: estado del run (éxito o fallo y dónde), tabla de los 10 sitios por rasgo (score, p.emp, p.emp_fact y sus
BH), recuentos de falsos positivos, resultado de la calibración (tasas frente a alpha, por clases), ficheros creados o
modificados, y cualquier cosa que no cuadre. Sin commit.
````
