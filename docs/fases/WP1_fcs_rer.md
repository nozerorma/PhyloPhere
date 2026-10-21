# WP1. Estadísticos del FCS y de RER, informes y reducción de redundancia

Estado: pendiente. Depende de: WP0. Independiente de WP2 y WP3 (ficheros disjuntos).
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
- Conversación y fechas:
- Commit de partida / commit de cierre:
- Qué se hizo (ficheros y cambios):
- Comprobaciones ejecutadas (comando, conjunto de datos, tolerancia, resultado):
- Mediciones (recurso antes y después):
- Desviaciones del plan y decisiones tomadas (con su justificación):
- Pendiente que pasa a otra fase:
```
