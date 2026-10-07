# Prompt: analizar el toy de Carn con los p empíricos nuevos

```text
Analiza el run toy de Carn en correfoc (repositorio PhyloPhere, rama scoring_v2, desde el commit 8faee85; lee
docs/scoring_v2_p_emp.md).

DONDE ESTA: /homes/users/mramon/scratch/3.Work_dirs_carn_toy2/Carn_toy_complete es el directorio de TRABAJO de Nextflow
(solo subcarpetas hexadecimales: 02, 09, 0d...), no los resultados. Los resultados (scoring/, caas_permulation/, fcs/...)
están en el outdir del run: búscalo en el .nextflow.log o el script de lanzamiento (params.outdir / results_dir), o en
~/scratch/2.Primates/2.Primates_results/CAAS_RESULTS/carn/ (el run completo estaba en Carn_complete; el toy debería ser
Carn_toy_complete). Si hay varios candidatos, elige el más reciente y dime cuál. Si no lo encuentras, pregunta; no analices el
directorio de trabajo.

REGLAS: nada de cómputo en el nodo de entrada (solo ls, head, tail, cat de logs y scp de tablas pequeñas); lo pesado, en
local con las tablas copiadas al scratchpad (nunca /tmp a pelo) o con sbatch. Si ssh se cuelga: ssh -o ControlPath=none
-o ControlMaster=no correfoc. Sin commit. Estilo descriptivo, sin rayas largas.

CONTEXTO (lo que cambió y hay que tener en cuenta)
- caas_score_aggregation = cumulative por defecto; el nulo trae la columna score_aggregation y el scoring se detiene si no
  coincide con el observado.
- p.emp se mantiene (con p = 1 si el score es 0). Nuevos: p.emp_fact y p.adj_bh_fact (posición), gene_caas_pperm_fact{,_top,
  _bottom} y _fact_adj* (gen). p factorizado = (ciclos que puntúan la unidad + 1)/(N + 1) x fracción de detecciones de su clase
  (propensión, y tamaño en genes) con score >= el observado; no está limitado por 1/(N + 1).
- Ya no existen p.adj_sam ni el p de gen poolado. gene_caas_pperm_adj es BH sobre todo el universo del nulo.
- En el toy N (permulaciones) es pequeño: el suelo de p.emp es 1/(N + 1) y las clases tienen reservas pequeñas (una clase de
  gen con menos de 200 detecciones vuelve a la clase de propensión). Lee N de ncol(caas_perms.rds$caas_corStat_byrank$global)
  y juzga todo con ese N, no con el de producción.

QUE HACER
1. Salud del run: exit code, errores en .nextflow.log o logs; que position_scores.tsv tenga p.emp, p.adj_bh, p.emp_fact,
   p.adj_bh_fact, amino_encoded, n_conserved_pairs y no p.adj_sam; que gene_scores.tsv tenga gene_caas_pperm_fact*; que
   perm_pos_cycle_caas.tsv.gz tenga score_aggregation = cumulative. Si fcs_caas_score = fact, que exista
   scoring/caas_perms_fcs.rds. Di qué valor tuvieron caas_score_aggregation y fcs_caas_score.
2. Distribución de p: por tramos (suelo = 1/(N+1), <=1e-3, 1e-2, 5e-2, = 1) para p.emp y p.emp_fact, y recuento a BH <= 0.05 y
   <= 0.25 de p.adj_bh y p.adj_bh_fact. Cuántas posiciones tienen score 0 y p = 1.
3. Calibración con el propio nulo (cada ciclo nulo como observado frente a los demás, sin contarse): tasa de pares
   (posición, ciclo) con p <= alpha (1e-3, 1e-2, 5e-2) sobre toda la familia, para p.emp y p.emp_fact, por clase de propensión
   (ciclos que puntúan la posición: <=5, 6-20, 21-100, >100). Debe ser <= alpha. Lo mismo para genes (global, top, bottom) por
   propensión y por tamaño del gen (1, 2-5, 6-20, >20 posiciones). Reporta violaciones y el tamaño de las reservas.
4. Observado frente a nulo: número de posiciones con p <= alpha en el observado frente a lo que da un ciclo nulo (media y
   percentil 95). En controles negativos de PEPC el factorizado heredaba un exceso de cola de p.emp (unas 2 veces a alpha 0.01);
   mira si aquí aparece algo parecido o es señal.
5. Las mejores posiciones y genes por p_fact: score, n_schemes, nd (ciclos que las puntúan), n_conserved_pairs, residuos por
   especie, amino_encoded. Distingue las que son de propensión alta (columnas variables, tipo MHC) de las de propensión baja
   y di cuáles te parecen convergencia real y cuáles no.
6. Genes: gene_caas_pperm frente a gene_caas_pperm_fact (correlación de rangos, mínimos, nº con p <= 1e-4, BH sobre el universo);
   confirma que ningún gen sale a BH <= 0.05 si es lo que se espera con este N, y cuál es el mejor.
7. FCS: si fcs_caas_score = raw, compara offline el ranking crudo (gene_caas_score) con -log10(gene_caas_pperm_fact): rho de
   Spearman, empates en el top, genes que más se mueven. Si es fact, comprueba que score_global/top/bottom de fcs_stats.tsv
   son -log10 de esos p y mira las diferencias en los términos enriquecidos frente al run de Carn anterior si lo hay.

ENTREGA: resumen corto con tablas, lo que cuadra, lo que no, y qué decidirías cambiar. Indica siempre N y el número de
unidades de cada reserva cuando una cifra dependa de ellos.
```
