# WP1. Estadísticos del FCS y de RER, informes y reducción de redundancia

Estado: HECHO (2026-10-21), validado con el toy de cancer y las tablas reales; pendientes: medir RER_CONT a 10.000 permulaciones (cluster) y la diferencia de recuentos 1, 0, 6, 6 (ver registro). Depende de: WP0. Independiente de WP2 y WP3 (ficheros disjuntos).
Conversación independiente. Antes de empezar, leer `docs/fases/README.md`.

## Objetivo

Cambiar la capa de significación del FCS y el ranking de RER, actualizar los informes y añadir una reducción de redundancia a posteriori. No toca el nulo de CAAS.

## Decisiones aplicables (del usuario)

- `fcs_min_genes` 10 y `fcs_max_genes` sin límite (0), en configuración, GUI, plantillas, ayuda, README y valores por defecto de R.
- Criterio principal: BH sobre el p de permulación del Wilcoxon por ranking y base de datos, q < 0,15 (`fcs_fdr`). Segundo umbral: `fcs_pperm_thr` 0,025 (Redlich et al. 2024: p ajustado < 0,15 y p de permulación < 0,025, justificado porque el p de permulación corrige la no independencia de los rangos de genes y de la filogenia; ellos reconocen que parte del refinamiento viene de ese filtro). Voto bilateral en rankings firmados. Se conserva el q analítico como columna.
- Lachenbruch fuera del FCS. POSENRICH no se toca aquí (tiene su propia copia del test y cae a `fcs_fdr` si desaparece `fcs_fdr_lachenbruch`).
- RER: solo ranking global y firmado; estadístico paramétrico `sign(Rho) · -log10(P)` en observado y nulo (como `getStat`); permulaciones a 10.000 (medir antes `RER_CONT`).
- El FCS se presenta como resultado.
- Informe 15: el peso por evidencia pasa a ser el número de clústeres significativos distintos que contienen el gen (tope 2); `evidence_score` pasa a −log10 del menor `q.perm`.
- Redundancia: solo a posteriori (no modifica ningún q). La inferencia principal sigue sobre las colecciones completas por base de datos.

## Alcance (ficheros)

`subworkflows/ENRICHMENT/local/src/fcs_enrich.R` (funciones de Lachenbruch ~233 a 281; `fcs_classify_evidence` ~564; esquema de columnas ~538 y ~635 a 644; `fcs_empirical_corstat` ~593; `fcs_run_all` ~620 a 808), `fcs_compute.R`, `subworkflows/ENRICHMENT/fcs.nf`, `conf/enrichment.config`, `conf/rerconverge.config`, `subworkflows/RERCONVERGE/rer_gene_lists.nf` (líneas ~72 a 86), `continuous_rer.R`, `scoring_compute.R` (~1176 a 1186, escribe `fcs_stats_rer.tsv` desde el p de permulación: hace falta un `rer_P` paramétrico aguas arriba), `12.FCS_general_report.Rmd` (es el mayor consumidor de las etiquetas y de `lach_*`), `15.Comparison_report.Rmd` (`.evidence_weight_of` ~139 a 146; ~434 a 460; ~688 a 725; ~1046 a 1105), pestaña de enriquecimiento de la GUI (`enrichment_tab.py`, `modules.py`), `gui/i18n.py`, `gui/templates/*.json`, `validation/tier1/templates`, `validation/tier2/templates`, `README.md`, `workflows/help.nf`, `run_single.sh.j2` (valores por defecto 5 y 500) y `sbatch_array.sh.j2`.

## Cambios

1. Parámetros: `fcs_min_genes` 10, `fcs_max_genes` 0 en todas partes (hoy hay 5 y 1000 en config y GUI, 5 y 500 en `run_single.sh.j2`, y 500 por defecto en varias funciones de R; `malignant_def.json` tiene mínimo 2).
2. Lachenbruch: eliminar funciones, columnas `lach_*`, `sig_*` y `evidence_*`, y el parámetro `fcs_fdr_lachenbruch` (GUI, i18n, plantillas).
3. Clasificación nueva `fcs_classify_significance`: BH sobre `p.perm` → `q.perm`; `sig_perm = q.perm < fcs_fdr`; `sig_top = sig_perm & p.perm < fcs_pperm_thr`; dirección = `sign(stat)` en rankings firmados.
4. RER: ranking global únicamente (se eliminan acelerado y decelerado, `rer_gene_lists.nf:81-82` y las pestañas RER de top y bottom del informe 15); `score_global = sign(Rho) · -log10(P)`; el nulo es `corStat` tal cual (se elimina `fcs_empirical_corstat` y el camino de `corRho`); definir el leading edge de un ranking firmado (hoy `12.Rmd:318-322` descarta `drv_score <= 0`). Comprobar que `continuous_rer.R:165` usa el mismo método de correlación que `getPermsContinuous` (`method = "p"`). El RER binario comparte `rer_perm_batches`.
5. Informe 15: sustituir los pesos Supported = 2 y Exploratory = 1.
6. Reducción de redundancia a posteriori: nuevo `fcs_cluster.R` llamado en el informe 12 tras `sig`. Asignación voraz por significancia: el término con menor `q.perm` es representante; cada término siguiente entra en el primer representante con similitud de leading edge ≥ umbral, o crea uno nuevo (determinista, sin encadenamiento). Salida `cluster_id` y `is_representative`. Similitud por defecto: Jaccard (el coeficiente de solapamiento colapsa en exceso las bases jerárquicas: Reactome 1236 → 26 y HPO 3113 → 1 con asignación voraz por tamaño; con Jaccard ≥ 0,5, 657 y 1858). Parámetros expuestos en la GUI. Calibrar el umbral con las tablas existentes (recomendaciones externas: EnrichmentMap, solapamiento 0,5 o Jaccard 0,25, son de su documentación).
7. Comprobación barata, sin recomputar: recalcular BH sobre las tablas FCS existentes quitando copias exactas y casi copias (Jaccard ≥ 0,9 de la composición completa, independiente de los datos) y anotar cuánto cambia el número de significativos. No entra por defecto.

## Interfaz hacia las fases siguientes

- Columnas de salida: `q.perm`, `sig_perm`, `sig_top`, `cluster_id`, `is_representative`, q analítico. WP6 reutiliza la clasificación como cierre (`fcs_close()`).

## Validación (criterios de paso)

- `p.perm` de CAAS igual al de referencia hasta 1e-14 relativo (el nulo no cambia): comparar con las tablas congeladas y con `legacy_fcs/` (WP0).
- Clasificación igual a una recomputación independiente con `p.adjust` desde la tabla antigua; recuentos de referencia 1, 0, 6 y 6 conjuntos (neoplasia estricto, maligno estricto, neoplasia completo, maligno completo) con la puerta de 0,025 apagada.
- Estadístico observado de RER frente a `RERconverge::getStat` del commit 2bd328f hasta 1e-12; identidad numérica `corStat == sign(corRho) · -log10(corP)` en el objeto de permulaciones.
- Agrupación: conjuntos sintéticos de solapamiento conocido, estabilidad al orden de filas, partición determinista.
- Informes 12 y 15 renderizan sin etiquetas ni columnas `lach_*`; POSENRICH sigue funcionando con `fcs_fdr`.

## Fuera de alcance

El nulo de CAAS, sus bloques y `GENE_P` (WP5 y WP6); POSENRICH (WP2).

## Riesgos y puntos a verificar

- El informe 15 depende de las etiquetas para ponderar: la nueva regla cambia el score compuesto de genes; revisar su efecto en la salida.
- Subir las permulaciones de RER a 10.000: `combinePermData` puede duplicar memoria de forma transitoria; medir `RER_CONT` antes (hoy 8 CPU y 32 GB × intento, tarea serie con `errorStrategy 'ignore'`).

## Registro de la fase

```
- Conversación y fechas: WP1, 2026-10-21, en la misma conversación que WP0 y WP2 (a petición del usuario), local, rama remove_caap.
- Commit de partida / commit de cierre: e2f9515 (WP2) / b0786dc.
- Qué se hizo (ficheros y cambios):
  * Parámetros: fcs_min_genes 10 y fcs_max_genes 0 en conf/enrichment.config, workflows/help.nf, README.md, gui/models/modules.py, run_single.sh.j2 (valores por defecto del JSON de parámetros), defaults de las funciones de R y las 14 plantillas de gui/templates (malignant_def tenía mínimo 2).
  * Lachenbruch fuera del FCS. fcs_enrich.R: eliminadas fcs_lach_prepare, fcs_lach_index, fcs_lach_stats, fcs_run_lachenbruch, fcs_compute_lach_p_perm, fcs_classify_evidence y las columnas lach_*, sig_wilcoxon, sig_lachenbruch, evidence_count y evidence_label. Parámetros fcs_fdr_wilcoxon y fcs_fdr_lachenbruch retirados (config, fcs.nf, GUI, run_single.sh.j2, plantillas; migrate() los descarta: _RETIRED_FIELDS["enrichment"]). POSENRICH no se toca: sigue cayendo a fcs_fdr (params.fcs_fdr_lachenbruch ?: params.fcs_fdr).
  * Clasificación nueva fcs_classify_significance (fcs_enrich.R): q.perm = BH sobre p.perm por (ranking, base de datos); sig_perm = q.perm < fcs_fdr; sig_top = sig_perm y p.perm < fcs_pperm_thr; direction = signo de stat en los significativos. El q analítico (p.adj) se conserva como columna y no decide. Sin nulo en un ranking: sin p.perm y sin significativos. Columnas de salida: ranking, database, pathway, stat, pval, p.adj, p.perm, q.perm, num.genes, gene.vals, sig_perm, sig_top, direction (+ description, cluster_id, is_representative y sig en el informe 12).
  * RER: solo el ranking global y firmado. rer_gene_lists.nf: score_global = sign(Rho) * -log10(P) paramétrico, sin score_accelerating ni score_decelerating (las listas de genes y los flags flag_rer_acc/decc se mantienen). scoring_compute.R: nueva columna rer_P (P paramétrico) y fcs_stats_rer.tsv con score_global de esa columna. fcs_enrich.R: el nulo de RER es corStat tal cual (eliminados fcs_empirical_corstat y el camino de corRho); fcs_resolve_corstat_rk usa corStat para "global".
  * Reducción de redundancia a posteriori: nuevo subworkflows/ENRICHMENT/local/src/fcs_cluster.R (fcs_set_similarity, fcs_cluster_sets, fcs_add_clusters): asignación voraz por q.perm dentro de cada ranking, similitud de borde principal Jaccard (u overlap), umbral 0,5; columnas cluster_id (<ranking>_C<k>) e is_representative. Parámetros fcs_cluster_similarity y fcs_cluster_threshold en conf/enrichment.config, fcs.nf, help.nf, README.md, GUI (campos en la pestaña de enriquecimiento), run_single.sh.j2 y sbatch_array.sh.j2.
  * Informe 12 reescrito (una prueba, orden por q.perm, tablas de significativos con Group y Rep, gráfico de puntos por ranking solo con representantes, borde principal con signo: en un ranking firmado, los miembros cuyo score tiene el signo del enriquecimiento del conjunto; en uno de magnitud, los de score > 0). Informe 15: peso por gen = número de clústeres significativos distintos que contienen el gen (tope 2; columnas CAAS clusters y RER clusters), evidence_score = -log10 del menor q.perm; sin rankings acelerado y decelerado (RER solo contribuye a Global); sin etiquetas Supported y Exploratory.
  * Pruebas: PhyloPhere_validation/unification/test_fcs_no_null.py (reescrita), test_fcs_cluster.py (nueva, 7 pruebas), test_rer_statistic.py (nueva, 3 pruebas), test_enrichment_tests_definition.py (las tres de equivalencia con el Lachenbruch de POSENRICH usan ahora legacy_fcs/fcs_enrich.R); wiring/test_gui_params_json.py adaptada a los campos nuevos.
- Comprobaciones ejecutadas (comando, conjunto de datos, tolerancia, resultado):
  * p.perm de CAAS igual a la referencia (toy de cancer, 183 filas, fcs_compute.R con --num-g 5 --max-g 1000 sobre el caas_perms.rds de golden_581b999/cancer_b0_toy): diferencia relativa máxima de stat, pval, p.adj, p.perm y num.genes entre 0 y 2,5e-16 (criterio 1e-14). Repetido dentro del pipeline completo (-resume, mismos parámetros 5/1000): diferencia 0.0 en las cinco columnas.
  * q.perm igual a un p.adjust independiente (R y Python), diferencia máxima 2,2e-16, sobre la tabla del toy y sobre las cuatro tablas reales de ~/.cc/caas_fcs.
  * Estadístico de RER: corStat == sign(corRho) * -log10(corP) en el objeto de permulaciones del toy marino (118.935 celdas, mismo patrón de NA, diferencia máxima 0); estadístico observado igual a RERconverge::getStat (diferencia 0 sobre el objeto y 1,8e-15 sobre el resumen). En la corrida del pipeline, fcs_stats_rer.tsv frente a sign(Rho) * -log10(P) del resumen: 3,2e-13 (precisión del TSV).
  * Pipeline completo con RER (toy de cancer, 2 lotes x 50 permulaciones, caché de ASR de la referencia), exit 0: informes 12 (CAAS y RER) y 15 sin errores; RER con un solo ranking (global, 14.331 conjuntos, 6.892 con stat negativo, p.perm sin NA, mínimo q.perm 0,2475 por el suelo de 100 ciclos: ninguno significativo); borde principal de RER: 561.154 filas y el signo del score coincide con el del enriquecimiento en todas (ninguna con score 0).
  * Recuentos sobre las cuatro tablas reales (BH sobre p.perm por ranking y base de datos, q < 0,15, puerta de 0,025 apagada): neoplasia estricto 0, maligno estricto 0, neoplasia completo 6, maligno completo 3; con la puerta de 0,025 encendida: 0, 0, 6, 3. NO reproducen los recuentos de referencia del plan (1, 0, 6, 6): esos salieron del script de control ctrl4 (~/.cc/enrich/ctrl_calib4.R), que recalcula el p de permulación sobre los ciclos, y no de las tablas publicadas; coinciden solo en neoplasia completo. Diferencia sin resolver (ver pendiente).
  * Comprobación barata de redundancia (golden_581b999/wp1_checks/dedupe_bh_check.py): quitando copias exactas y casi copias (Jaccard >= 0,9 de la composición completa, selección independiente de los datos) desaparecen el 10,5 a 10,8 % de las filas de las tablas (disgenet 31,6 %, HPO 22,3 %, OMIM 25,6 %, Reactome 11,9 %) y los significativos no cambian (0, 0, 6, 3). No entra por defecto.
  * Suites: unification + resources: 252 failed, 300 passed, 17 errors, el mismo conjunto de fallos que la línea base de WP0 (0 nuevos; 290 a 300 passed por las pruebas nuevas). wiring (10 ficheros sin Nextflow): 57 failed con dos fallos nuevos (pruebas que usaban fcs_fdr_wilcoxon y fcs_fdr_lachenbruch), adaptadas; tras adaptarlas, los fallos son los 55 de la línea base de WP2.
  * Observación de código (leída, no medida): en la instalación de RERconverge, correlateWithContinuousPhenotype y getPermsContinuous usan la misma correlación (getAllCor con method = "p" y min.pos = 0), pero el nulo (getNullCor) llama a getAllCor con min.sp por defecto (10) mientras el observado recibe rer_minsp (15 por defecto): puede cambiar qué genes tienen estadístico en el nulo y en el observado, no su valor.
- Mediciones (recurso antes y después): sin medición de RER_CONT a 10.000 permulaciones (ver desviaciones). Toy de cancer con RER y 100 permulaciones: corrida completa de unos 35 min en local.
- Desviaciones del plan y decisiones tomadas (con su justificación):
  * Permulaciones de RER: el valor por defecto (rer_perm_batches 10 x rer_perms_per_batch 100 = 1.000) NO se sube a 10.000: el plan exige medir antes RER_CONT (combinePermData puede duplicar memoria) y esa medición necesita el cluster (solo sbatch). La perilla ya existe en config, GUI y plantillas.
  * Umbral de agrupación por defecto 0,5 con Jaccard, el del plan (Reactome 1236 a 657 y HPO 3113 a 1858); no se recalibró aquí porque las tablas locales no llevan el borde principal.
  * Los flags flag_rer_acc y flag_rer_decc se siguen escribiendo en el fcs_stats.tsv de RER_GENE_LISTS (alimentan las columnas RER_dir del informe 15); fcs_stats_rer.tsv (el ranking del FCS) solo lleva gene y score_global. Las pestañas top y bottom de RER del informe 15 desaparecen porque RER ya no tiene esos rankings.
  * Los toys de 31 genes no tienen conjuntos de al menos 10 genes: con los nuevos valores por defecto la tabla de FCS sale vacía (probado, el informe 12 la maneja); la validación numérica se hizo con mínimo 5 y máximo 1000 como en la referencia.
  * validation/tier*/templates no se editaron (los gestiona el usuario): conservan fcs_fdr_wilcoxon, fcs_fdr_lachenbruch y los mínimos antiguos; migrate() descarta los campos retirados.
  * No se ejecutaron los tests de wiring que necesitan Nextflow.
- Pendiente que pasa a otra fase: (1) medir RER_CONT a 10.000 permulaciones por sbatch y decidir el valor por defecto; (2) resolver la diferencia de recuentos 1, 0, 6, 6 frente a 0, 0, 6, 3 (¿tablas o control?); (3) WP6 reutiliza fcs_classify_significance como cierre (fcs_close()) y las columnas q.perm, sig_perm, sig_top, cluster_id, is_representative y p.adj; (4) sincronizar validation/tier*/templates; (5) POSENRICH conserva su Lachenbruch (WP2/WP5).
```

## Adenda (2026-10-21): regla analítica de Redlich como columna

A petición del usuario, junto a la regla de q.perm se informa la regla publicada de Redlich et al. 2024 (leída en su suplemento, líneas 27 a 38 y 420 a 427, y los pies de las Fig. S5 y S7): BH del p analítico del test de enriquecimiento (p.adj) por debajo de un nivel, y p de permulación bruto (sin corregir) por debajo de 0,025 (fcs_pperm_thr). Se añade la columna redlich_level = el menor de los niveles de fcs_analytic_levels (por defecto 0.05,0.1,0.15; campo en la GUI) que cumple el conjunto, vacía si no cumple ninguno. No interviene en sig_perm, sig_top ni en los grupos de redundancia. El informe 12 muestra una tabla de recuentos por regla y nivel (acumulados) y una sección con los conjuntos de la regla analítica; se requirió p.perm < 0,025 sin exigir stat > 0 (en un ranking de magnitud un conjunto empobrecido tiene p.perm >= 0,5; en uno firmado se conservan los dos lados).

Comprobado: sobre las cuatro tablas reales, fcs_classify_significance da el mismo recuento que un cálculo independiente (p.adj < nivel y p.perm < 0,025 y stat > 0): 192, 248 y 299 (neoplasia estricto), 155, 219 y 255 (maligno estricto), 67, 97 y 137 (neoplasia completo), 80, 111 y 144 (maligno completo), a 0,05, 0,10 y 0,15. Con q.perm < 0,05 y < 0,10 la regla nueva da 0 conjuntos en las cuatro (mínimo q.perm 0,108 a 0,344). Pruebas nuevas en test_fcs_no_null.py (nivel mínimo, niveles como parámetro, independencia de sig_perm, ranking firmado); sin fallos nuevos en unification (302 passed, mismo conjunto de 252 fallos y 17 errores) ni en wiring (55 fallos de la línea base; 110 passed).

Lachenbruch (pregunta abierta del usuario, sin cambios de código): medido sobre las tablas reales, dentro de los 9 conjuntos de la regla nueva (6 de neoplasia completo y 3 de maligno completo) lo pasa la puerta antigua (lach_p.adj < 0,15 y lach_p.perm < 0,025) en 1 y 2 respectivamente; frente al Wilcoxon antiguo (299, 255, 137, 144) exigir también Lachenbruch dejaba 73, 19, 44 y 68. Su p de permulación está inflado en los estrictos (proporción con p < 0,05: 0,083 y 0,073 frente a 0,035 y 0,037 del Wilcoxon). Opciones registradas: devolverlo como anotación sin voto, o como filtro opcional apagado por defecto; la decisión es del usuario.

## Adenda (2026-10-21): medición de RER_CONT (permulaciones de RER)

Qué se midió: RER_CONT (`continuous_rer.R`, el del commit `d246e82` del cluster, sin cambios) sobre las entradas reales de neoplasia (RER de 16.109 genes, `rer_minsp` 15, `ha_logit`), 4 CPU y 32 GB pedidos, en 14 tareas independientes de un array de SLURM (trabajo 7571897): los primeros G genes de un orden aleatorio fijo (semilla 1998, subconjuntos anidados) con N permulaciones en lotes de 100 (y una con lotes de 1.000). Cada tarea registra tiempo, memoria máxima (VmHWM), CPU y tamaño del objeto. Datos y scripts: `PhyloPhere_validation/golden_581b999/wp1_checks/{rer_cont_scaling.tsv,rer_cont_scaling_raw.txt,rer_cont_rally.sh,rally_params.tsv}`.

Resultados medidos (tiempo, pico de memoria y objeto de permulaciones):
- 16.109 genes: 500 permulaciones 1.102 s, 1,07 GB, 184 MB; 1.000: 2.153 s, 1,58 GB, 368 MB; 2.000: 4.255 s, 2,46 GB, 737 MB.
- 4.000 genes, 10.000 permulaciones: lotes de 100 5.477 s y 3,29 GB; lotes de 1.000 5.467 s y 2,80 GB (objeto de 914 MB).
- 1.000 genes, 10.000 permulaciones: 1.682 s y 1,16 GB. 200 genes, 10.000: 650 s y 649 MB.
- El tiempo de CPU es igual al de reloj en todas las tareas (un solo hilo efectivo). El tamaño del objeto es 24·G·N bytes (corP, corRho y corStat).

Lectura (ajuste sobre las tareas de 1.000 permulaciones con G de 1.000 en adelante): T ≈ 60 s + 0,13 s por gen por cada 1.000 permulaciones. Predice 5.256 s para 4.000 genes y 10.000 permulaciones frente a 5.467 y 5.477 medidos (4 % por debajo). El reparto en lotes no cambia el tiempo (5.467 frente a 5.477 s) y cambia la memoria: el pico es unos 2,55 bytes por byte de objeto con lotes de 1.000 y unos 3,08 con lotes de 100, además de ~470 MB de base.

Extrapolación (NO medida) a 16.109 genes y 10.000 permulaciones: unas 5,8 a 6,1 h (lineal en N, validado hasta 10.000 con 4.000 genes) y un pico de unos 9,9 GB con lotes de 1.000 o 11,8 GB con lotes de 100 (objeto de 3,7 GB). Cabe en los 32 GB de RER_CONT; el tiempo cabe en 8 h del primer intento con poco margen. No se midió la esquina completa (6 h) ni RER_BIN (el rasgo binario comparte `rer_perm_batches` y usa otra permulación; en el toy marino costó 16 min por 1.000 permulaciones con 121 genes).

Pendiente: decidir el valor por defecto (la medición apoya 10.000 permulaciones con 10 lotes de 1.000) y, si se quiere, ajustar los recursos de RER_CONT (1 CPU, más tiempo); `rer_perm_batches` y `rer_perms_per_batch` son campos de la GUI y de la configuración.

