# Revisión del FCS frente a RERconverge

## Contrato vigente tras WP8 (2026-10-23)

Esta sección sustituye las propuestas y descripciones históricas de implementación de este documento. WP5, WP6 y WP7 están implementados y validados de forma dirigida; no quedan por construir sus bloques ni su FCS productivo. La aceptación integrada y las medidas reales a escala siguen pendientes.

- Nulo CAAS canónico: `caas_permulation/blocks/blk_*/matrices.rds`, `cycle_roster.txt`, `null_metadata.json` y contadores/tablas WP5. Cada bloque conserva genes, ciclos, unidades, `gene_stat`, `stat` y `parts` (max/fmax/n). Completo y strict usan el mismo contrato. No hay tabla pública de cada posición × ciclo.
- `caas_export_single=false`: el camino normal no ensambla `caas_perms.rds`. La exportación explícita reconstruye el objeto histórico con coste global deliberado. Se conserva la importación del RDS histórico; el almacén canónico completo permite reutilizar el nulo sin él.
- Ranking CAAS productivo: −log10(p de permulación del gen), con observado y leave-one-out en `gene_p.R`. `gene_stat` sigue parametrizado; WP9 decide entre máximo, F(max) y F(max)^n. El ranking crudo es un control, no el fallback productivo.
- Rutas: `stat = AUC − 0,5`, cola superior CAAS y bilateral RER; empates inclusivos, `(k+1)/(N_válidas+1)`. `sig_perm` es BH sobre `p.perm`, por ranking/base. `redlich_level` es independiente y `sig_top` está retirado. El observado interno usa RDS para conservar empates y recuentos; el TSV es público.
- RER continuo/binario genera bloques con semillas por lote. La disponibilidad nula se evalúa por ciclo y no se restringe a genes disponibles en el observado. Cambiar tamaño de lote cambia las simulaciones; conservar ajustes, semilla y partición al comparar.
- No se afirma memoria constante del pipeline: GENE_P retiene todos los ciclos de un lote de genes; el bosque RER y el universo de posiciones tienen costes separados. Las medidas sintéticas de WP6 no sustituyen recursos reales.
- `remove_intermediates` retirado: era inerte tras WP5. La higiene de tareas es incondicional; ningún bloque, roster o índice público se poda. Los proyectos antiguos conservan solo su elección de `remove_work`. Este último pertenece al launcher: tras éxito elimina work, directorio de ejecución y NXF_HOME privado, con guardas de rutas; no actúa por invocar Nextflow directamente.
- POSENRICH se omite con aviso en el flujo productivo tras WP5. Decisión del usuario en WP8: preparar una fase propia de AUC por posición (opción A), sin activarla aquí; `docs/fases/WP10_posenrich_auc.md`.

Controles: máximo crudo frente a ranking calibrado, nulo como observado, tres estadísticos base y p/BH por gen se distinguen de las pruebas por posición (`p.emp`, `p.emp_fact`) y del AUC por conjunto. Los cuatro `adhoc_scripts/fcs*.sbatch` son oráculos de memoria completa; `tests/test_wp6_blocks.R`, `tests/test_fcs_auc_contract.R` y `tests/test_wp7.py` cubren contratos dirigidos. Los seis RDS PEPC/marine/cáncer completo/strict y los 399 resultados marinos de WP7 están registrados en su fase, no se han vuelto a ejecutar en WP8.

La sensibilidad de agrupación de WP8 usa una tabla histórica con leading edges, no recalcula inferencia productiva ni calibra potencia: Jaccard 0,25/0,5/0,75 conserva los tres grupos BH de maligno completo; overlap da 2/2/3. Se conserva Jaccard 0,5 y la agrupación solo de presentación. Criterios y comandos para aceptación: `docs/fases/WP8_aceptacion_integrada.md`.


Fecha de la revisión: 20/10/2026.
Estado: revisión histórica del 20/10; implementación vigente y pendientes en el contrato anterior.

## 1. Criterio de trabajo

RERconverge no es un método propio. Cada punto en que el FCS de PhyloPhere se aparta de él se resuelve de una de dos formas:
se alinea con RERconverge, o se conserva con una justificación escrita (fuente publicada o medida).
No se toman decisiones ad hoc. Línea base: commit `2bd328f7530b4aca9b48c0b3997875c9b77a7026` (`environment/install_env.sh`).

## 2. Fuentes verificadas

- Código del paquete instalado en correfoc (0.3.0), deserializado desde su `.rdb`: `fastwilcoxGMT`, `simpleAUCgenesRanks`, `getStat`,
  `getEnrichPerms`, `permpvalenrich`, `permpvalcor`, `getPermsContinuous`. Coincide con el commit fijado.
- Walkthrough de permulaciones de Saputra et al. 2021 (suplemento), pp. 11 a 13:
  `validation/tier2/references/marine/saputra2021_msab068/supplementary_data/Supplementary Walkthrough.pdf`.
- Artículo de Saputra et al. 2021, sección "Empirical P Values for Pathway Enrichment" (doi 10.1093/molbev/msab068).
- Suplemento de Redlich et al. 2024 (doi 10.1093/molbev/msae210), texto de la Fig. S7 (origen de los umbrales 0,15 y 0,025) y de la Fig. S5.
- Comprobaciones de código propias (sintéticas, 400 genes): el p de nuestro nulo vectorizado, el del lado observado y `stats::wilcox.test`
  coinciden con diferencia relativa máxima de 1e-14; los dos ayudantes del p de permulación coinciden con su definición;
  la puntuación empírica de cada columna nula coincide con `permpvalcor` con la columna excluida (diferencia máxima 0,003 en unidades de -log10 p).

## 3. Resumen de puntos

| ID | Punto | Estado | Aplica a CAAS |
|---|---|---|---|
| P1 | Tamaño mínimo de conjunto = 10 | Decidido | Sí |
| P2 | Tamaño máximo de conjunto | Propuesto: quitar | Sí |
| P3 | Voto bilateral en rankings firmados | Decidido | No |
| P4 | Criterio de significación: BH sobre el p de permulación | Decidido en principio, regla unificada pendiente | Sí |
| P5 | Umbrales 0,15 y 0,025 | 0,15 aceptado, 0,025 pendiente | Sí |
| P6 | Lachenbruch | Pendiente | Sí (se conserva) |
| P7 | Votos y etiqueta "Supported" | Pendiente | Sí |
| P8 | Rankings acelerado y decelerado en RER | Propuesto: solo global | No |
| P9 | Ranking por gen en RER: p de permulación o P paramétrico | Pendiente | No |
| P10 | Número de permulaciones | Pendiente | Sí |
| P11 | Corrección de empates en el Wilcoxon | Se conserva | Sí |
| P12 | Pseudo-cuenta +1 en el p de permulación | Se conserva | Sí |
| P13 | Nulo vectorizado | Verificado | Sí |
| P14 | Parámetros distintos entre rasgos | Alinear | Sí |
| P15 | Calibración del nulo de CAAS | Abierto (criterio estricto del nulo verificado; falta el control de nivel) | Solo CAAS |

## 4. Detalle por punto

### P1. Tamaño mínimo de conjunto

- RERconverge: `num.g = 10` por defecto en `fastwilcoxGMT`, y llamada explícita `num.g=10` en el walkthrough.
- PhyloPhere: `conf/enrichment.config:18` (5), plantillas `gui/templates/malignant_def.json:202` (2) y `neoplasia_def.json:202` (5);
  `fcs_enrich.R:153` usa 10 por defecto.
- Medido: el 51 % de los conjuntos probados en maligno y el 25 % en neoplasia tienen menos de 10 genes. Con mínimo 10, los votos de Wilcoxon
  de RER maligno decelerado pasan de 35 (mínimo 2) a 66, porque la familia de BH pasa de 20.746 a 10.125 conjuntos.
- Propuesta: 10 en configuración, plantillas y GUI, igual en los dos rasgos.
- CAAS: los "Supported" cambian poco (global de maligno no estricto 19 a 19; neoplasia estricto 32 a 34 en global y 38 a 49 en top).

### P2. Tamaño máximo de conjunto

- RERconverge: no existe en `fastwilcoxGMT`.
- PhyloPhere: 1000 en `conf/enrichment.config:19` y plantillas, 500 por defecto en `fcs_enrich.R:153`. No hay justificación publicada.
- Medido (16 GMT, universo de RER): no hay conjuntos por encima de 500 en 10 de 16 bases. Cambia en `TCGA.PanCancer.OncoVar`
  (8 de 9 conjuntos superan 500; el universo anotado pasa de 15.709 a 435 genes), `c8` (42), factores de transcripción (36),
  HPO (194), Reactome (22) y DisGeNET (15).
- Propuesta: sin límite por defecto (valor 0), con `OncoVar` tratado aparte.
- CAAS: mismos GMT. La medida usó el universo de RER (16.109 genes); el de CAAS es de 11.530 a 12.473 genes.

### P3. Bilateralidad en rankings firmados

- RERconverge: contraste bilateral por defecto sobre un único ranking firmado; la dirección es el signo del estadístico.
- PhyloPhere: `fcs_enrich.R:567` exige `stat > 0` en todos los rankings. En el ranking firmado global se descarta el sentido negativo:
  6 conjuntos de neoplasia y 16 de maligno cumplen BH < 0,15 y p de permulación < 0,025 y quedan marcados como no significativos
  (por ejemplo HALLMARK_INFLAMMATORY_RESPONSE, BH 2,9e-4).
- Propuesta: voto bilateral en rankings firmados, con la dirección tomada del signo de `stat`.
- CAAS: no aplica; los tres rankings (global, top, bottom) son no negativos y unilaterales.

### P4. Criterio de significación: BH sobre el p de permulación

Tres p-valores intervienen:

1. p del gen, que ordena el ranking. RERconverge: P paramétrico. PhyloPhere (RER): p de permulación del gen.
2. p analítico del conjunto: Wilcoxon con aproximación normal sobre esos rangos.
3. p de permulación del conjunto: fracción de fenotipos nulos con un estadístico del conjunto igual de extremo.

La lista que se ajusta es la de conjuntos de una base de datos.

- RERconverge: el walkthrough (pp. 12 y 13) aplica `p.adjust(permpval, "BH")` por lista de anotación y lo guarda como `permpvaladj`;
  indica que el p de permulación se puede tratar como un p normal, incluida la corrección múltiple. El p de la Tabla 1 de Saputra va ajustado.
- PhyloPhere: `fcs_enrich.R:564` aplica BH al p analítico y usa el p de permulación como puerta fija sin ajustar (< 0,025).
- El walkthrough indica que el p paramétrico de enriquecimiento tiende a estar subestimado, de modo que el BH recae sobre el p peor calibrado.
- Medido en RER, ranking acelerado de neoplasia (conjuntos de al menos 10 genes): BH analítico < 0,15, 45 conjuntos; con p de permulación < 0,025,
  38; voto de Lachenbruch, 27; "Supported", 23.
  De los 29 "Supported" de neoplasia, 16 cumplen BH(p de permulación) < 0,15 (25 con < 0,25, ninguno con < 0,05).
  En maligno, 2 de 27 (8 con < 0,25).
  HALLMARK_PI3K_AKT_MTOR_SIGNALING en neoplasia: p analítico 0,00154, BH sobre 50 Hallmark 0,0256; p de permulación 0,00999 = (9+1)/1001;
  BH sobre el p de permulación 0,237.
- Diferencia entre RER y CAAS (conjuntos de al menos 10 genes, todos los rankings): conjuntos con p de permulación del Wilcoxon <= 0,001,
  observados frente a esperados por azar: RER neoplasia 73 frente a 30; RER maligno 71 frente a 30; CAAS neoplasia no estricto 27 frente a 28;
  CAAS maligno no estricto 32 frente a 27; CAAS neoplasia estricto 7 frente a 28; CAAS maligno estricto 0 frente a 27.
  El BH sobre el p de permulación solo puede pasar conjuntos si hay un exceso de p pequeños. En RER lo hay; en CAAS no estricto la cola es igual a la
  esperada por azar; en CAAS estricto es menor.
- La resolución (p mínimo 1/1001) pesa en las bases grandes: con unos 3000 conjuntos (HPO) harían falta unos 20 en el suelo. Ningún conjunto de HPO
  ni de DisGeNET aporta en RER.
- Propuesta: BH sobre el p de permulación por ranking y base de datos, conjuntos de al menos 10 genes, q < 0,15, igual en RER y en CAAS.
  El q analítico se sigue mostrando.
- CAAS: los "Supported" de maligno no estricto, maligno estricto, neoplasia no estricto y neoplasia estricto (68, 19, 44 y 73) pasan a 2, 0, 1 y 1;
  con cualquier voto, a 12, 7, 18 y 102. Solo 2 conjuntos de los "Supported" actuales cumplen también la regla nueva (ver P15 para el caso estricto).
- Aplica también a POSENRICH: `posenrich_enrich.py:117` y `:657` usan la misma puerta (`lach_p.adj` y `lach_p.perm < pperm_thr`).

### P5. Umbrales 0,15 y 0,025

- Origen: suplemento de Redlich et al. 2024, texto de la Fig. S7, dentro de la comparación de métodos de dieta (rutas de MGI; categórico y binario
  frente a phylANOVA). El 0,15 se elige porque el análisis tiene poca potencia (solo dos rutas bajaban de 0,05 con phylANOVA). El 0,025 se añade
  "para obtener un número razonable de rutas de cabeza"; el texto lo justifica porque el p de permulación corrige la no independencia de los rangos
  de genes y la filogenia. La Fig. S5 del mismo suplemento usa p ajustado <= 0,10 y p de permulación <= 0,05, con 10.000 permulaciones.
  No hay otra fuente de la regla en el texto principal ni en la viñeta de permulaciones categóricas.
- PhyloPhere: `conf/enrichment.config:20` a `:23`; `README.md:464`.
- Propuesta: se conserva 0,15. La puerta de 0,025 deja de ser criterio de significación; si se mantiene como columna, se documenta como filtro
  descriptivo de ese suplemento.
- CAAS: igual.

### P6. Lachenbruch

- RERconverge: no existe.
- PhyloPhere: se aplica a rankings con ceros: CAAS (global, top, bottom) y RER acelerado y decelerado. No corre en el ranking firmado global.
- RER: el cero del ranking acelerado mezcla genes con rho de signo contrario y genes sin evidencia; la prevalencia pasa a ser casi un test de signo
  (en PI3K 66 % de rho positivo frente al 51 % del fondo, y el 83 % del chi cuadrado). El test sigue siendo válido, pero no mide prevalencia de señal.
- CAAS: el cero significa ausencia de detección, que es lo que asume el test.
- Como filtro de concordancia aporta: en RER neoplasia acelerado, 13 de 23 "Supported" cumplen BH(p de permulación) < 0,15, frente a 4 de 15 con solo
  Wilcoxon y 0 de 7 de los que fallan la puerta de 0,025. En maligno apenas discrimina (2 de 26).
- Propuesta: no usarlo en RER si se adopta P8; conservarlo en CAAS.

### P7. Votos y etiqueta

- PhyloPhere: `fcs_enrich.R:564` a `:590`; "Supported" = dos votos, "Exploratory" = uno.
- Los dos votos leen los mismos genes y no son confirmaciones independientes. La etiqueta no es una tasa de error controlada.
- Propuesta pendiente: el criterio formal pasa a ser el q de permulación (P4); la concordancia entre tests se mantiene como columna, con un nombre que no
  sugiera independencia.

### P8. Rankings acelerado y decelerado en RER

- RERconverge: un único ranking firmado; la dirección es el signo del estadístico.
- PhyloPhere: `rer_gene_lists.nf:81` a `:82` define `score_accelerating` y `score_decelerating` (puntuación solo si rho > 0 o < 0, cero en el resto).
  Su origen es el contraste unilateral del ómnibus categórico de Redlich (commit `ef00fa8`), donde el estadístico es no negativo por construcción.
- Medido (PI3K, universo de 16 mil genes): el p varía unas 10 veces según la combinación de decisiones: 0,0068 con el test de RERconverge, 0,0037 al
  quedarse solo con rho positivo y test unilateral, 0,0022 con corrección de empates, 0,0014 con el p de permulación como puntuación, 0,0007 con ambas.
- Propuesta: RER solo con el ranking global (Wilcoxon, bilateral, BH sobre el p de permulación). Acelerado y decelerado quedan como partición descriptiva
  del borde de avance. Con el criterio fiel (global, al menos 10 genes) salen 1 conjunto con q < 0,15 y 21 con q < 0,25 en neoplasia, y 3 y 11 en maligno.
- CAAS: no aplica; top y bottom se definen por lado, no por signo de un coeficiente.
- Impacto aguas abajo: informe 15 (`rer_dir`, evidencia RER), informe 12.

### P9. Ranking por gen en RER

- RERconverge, código: `getStat = sign(Rho) * -log10(P paramétrico)` y el mismo estadístico en cada permulación (`getPermsContinuous`).
  Suplemento de Redlich: el nulo se calcula sobre "null sets of gene correlation statistics", con 10.000 permulaciones.
  Texto de Saputra: el nulo puntúa cada gen con el log de su p nulo empírico por el signo; no describe el lado observado.
- PhyloPhere: `rer_gene_lists.nf:81` usa el p de permulación del gen; `fcs_enrich.R:593` calcula la puntuación empírica de cada columna nula.
- Medido: el suelo del p de permulación de gen es 0,002; el 0,6 % de los genes está empatado en él.
- Argumento para alinear con el paquete: el p de permulación del conjunto ya corrige la calibración, de modo que la transformación a p empírico por gen
  es redundante para la validez y pierde resolución.
- Propuesta: paramétrico en ambos lados, salvo que una sensibilidad sobre pocos conjuntos muestre otra cosa. La sensibilidad requiere el RDS de
  permulaciones en correfoc (trabajo de cómputo, no en el nodo de entrada).
- CAAS: no aplica; la puntuación es propia y el nulo está en su escala.

### P10. Número de permulaciones

- RERconverge: al menos 1000, idealmente tantas como sea posible; el p mínimo es 1/numperms (walkthrough). Redlich usa 10.000.
- PhyloPhere: `rer_perm_batches` 10 por `rer_perms_per_batch` 100 en las plantillas (1000).
- Propuesta pendiente: subir las de RER si el cómputo lo permite. Para CAAS no parece resolver la cola vacía de P4.

### P11. Corrección de empates

- RERconverge: `simpleAUCgenesRanks` usa la varianza sin término de empates.
- PhyloPhere: `fcs_tie_term` y `fcs_rank_sum_sd` (`fcs_enrich.R:130` a `:139`).
- Justificación: fórmula estándar de Mann-Whitney con empates; con rankings de CAAS (92 % de ceros) la varianza sin el término es unas 5,75 veces mayor
  (medida de la nota del 08/10/2026). En RER acelerado el efecto es menor (en PI3K el p pasa de 0,0037 a 0,0022). El resultado coincide con
  `wilcox.test(correct = FALSE)` hasta 1e-14. Se conserva.

### P12. Pseudo-cuenta +1

- RERconverge: `permpvalenrich` cuenta `abs(nulo) > abs(obs)` sin +1; `permpvalcor` (genes) usa `(num + 1) / (denom + 1)`.
- PhyloPhere: `(1 + #nulo extremo) / (N + 1)` (`fcs_enrich.R:405` a `:424` y `:459` a `:474`). Se conserva: evita p = 0 y coincide con la
  versión de genes de RERconverge.

### P13. Nulo vectorizado

- El nulo se calcula con una implementación vectorizada propia (`fcs_enrich.R:355`), no con `fastwilcoxGMTall` por permulación.
- Verificado: AUC (prueba existente), p-valor frente a `wilcox.test`, ayudantes del p de permulación y escala empírica frente a `permpvalcor`.
- Límite: datos sintéticos de 400 genes y 25 columnas; no se ha probado con 16 mil genes y 1000 columnas.

### P14. Parámetros distintos entre rasgos

- `malignant_def.json:202` usa `fcs_min_genes` 2 y `neoplasia_def.json:202` usa 5. Los dos rasgos no son comparables entre sí hasta alinear el valor (P1).

## 5. CAAS: calibración del nulo (P15)

Hipótesis de trabajo: que en CAAS haya pocos conjuntos significativos tras la permulación, frente a los observados con el p analítico, puede indicar
que alguna parte de la base del FCS de CAAS no se está haciendo como corresponde.

Observaciones:

- Calibración del p de permulación entre conjuntos de al menos 10 genes (esperado bajo un nulo uniforme: 5 % con p <= 0,05 y 25 % con p <= 0,25).
  RER: cerca de lo esperado (Lachenbruch 5,0 a 5,2 % con p <= 0,05). CAAS no estricto: cerca de lo esperado.
  CAAS estricto: el Wilcoxon queda por debajo (2,2 a 4,7 % con p <= 0,05) y el Lachenbruch por encima (7,6 a 8,6 %, y 35 a 40 % con p <= 0,25).
- Con BH sobre el p de permulación solo, neoplasia estricto da 102 conjuntos con q de Lachenbruch < 0,15 (WikiPathways, quinasas, KEGG), sin respaldo del p analítico.
- Los temas recurrentes del FCS de CAAS (complemento y coagulación, reguladores de proteasas, galactosa, lumen vacuolar, unión a cofactores, colágeno y matriz)
  tienen q de permulación mayoritariamente entre 0,2 y 0,8. Mejores casos: galactosa en neoplasia estricto (q de Lachenbruch 0,06), complemento y coagulación
  en neoplasia estricto (0,10), lumen vacuolar en las corridas no estrictas (0,11 a 0,17). Ninguno cumple ambas pruebas.
  El único conjunto con q < 0,15 en las dos pruebas es HALLMARK_IL2_STAT5_SIGNALING en neoplasia no estricto (global).

Causas posibles ya identificadas en el pipeline:

- Desajuste de nivel entre observado y ciclos nulos. En Carn, el observado tenía 1135 genes puntuados frente a una mediana de 190 por ciclo nulo, y el
  control "ciclo nulo como observado" dio tasas de falso positivo de la puerta de permulación del 19 % con el estadístico AUC (nominal 2,5 %) para ciclos
  de nivel parecido; ese control es anterior a la corrección de empates. No se ha medido en las corridas de cáncer.
- Asimetría de pares entre observado y nulo (selector del nulo se detiene en K, el observado sigue hasta que Dunn para). `perm_match_pss` lo corrige para
  rasgos continuos, pero no se aplica a rasgos de conteo (`conf/ct.config:19`), que son los de estas corridas. En el modelo de cáncer de juguete la asimetría
  medida era pequeña y el hueco, borderline.
- Criterio estricto del nulo (verificado el 20/10/2026): las salidas estrictas (`caas_permulation`, `ct_disambiguation`, `meta_caas`) están en
  `<corrida>_complete/strict/` y salen del mismo replay que las no estrictas. La definición es la misma en el observado y en el nulo:
  - Mismo descubrimiento: `check_pattern` (`caas_id.py`) con `--max_conserved = int(n_pares * (1 - min_divergent_fraction))` y `--min_non_overlapping`
    (`caas_permulation.nf:244` a `:245`). El observado es la rebanada `b_0` del mismo replay.
  - Mismo criterio de fila estricta: el observado usa `is_strict_row` (`observed.py:32`: `conserved_pair` vacío o `0:`); el nulo usa `conserved = bool(conserved_pair)`
    (`disambiguate_single.py:1005`) con el mismo recorte del prefijo `{n}:` (`gene_wrapper.py`).
  - Mismo agrupador: `pool_labelings(strict=True)` (`driver.py:278`) deja fuera las hipótesis con pares conservados y las cuenta como hipótesis que no detectan,
    con el mismo diseño que el observado.
  - Comprobación con datos: las puntuaciones génicas de `b_0` en el estricto (`strict/caas_permulation/b0/gene_cycle_scores.tsv`) coinciden con
    `gene_scores.tsv` del observado estricto en las cuatro corridas: mismos conjuntos de genes puntuados (neoplasia 2.480, 1.224 y 1.658 en global, top y bottom;
    maligno 2.212, 726 y 1.781) y diferencia máxima 1,1e-16. El caso no estricto da el mismo resultado.
  - Límite: la comparación con datos cubre el ciclo `b_0` (etiquetado observado). Para los ciclos permulados la igualdad se apoya en que ejecutan el mismo código
    con `strict=True`; no hay un análogo observado con el que compararlos.
- Nivel de detección del observado frente a los 1000 ciclos nulos (posiciones detectadas, `perm_pos_quantiles.tsv`, `n_total`):

  | Corrida | Observado | Mediana nula | Cociente | Ciclos nulos con al menos el observado |
  |---|---|---|---|---|
  | Neoplasia no estricto | 81.355 | 55.954 | 1,45 | 116 de 1000 |
  | Neoplasia estricto | 5.004 | 2.252 | 2,22 | 70 de 1000 |
  | Maligno no estricto | 69.449 | 45.660 | 1,52 | 128 de 1000 |
  | Maligno estricto | 4.393 | 1.542 | 2,85 | 32 de 1000 |

  El desajuste de nivel es moderado en las corridas no estrictas y mayor en las estrictas. Coincide en orden con la calibración del Lachenbruch de permulación
  (inflado en estricto, cercano a lo esperado en no estricto). Es una asociación, no una demostración: un nivel observado alto puede ser señal real
  y no desajuste del nulo.

Lectura: hay dos explicaciones compatibles con los datos. Una, que no haya señal específica de ruta más allá de lo esperable (apoyada por el p de permulación
casi uniforme en no estricto). Otra, que el nulo no sea comparable con el observado (apoyada por la contradicción entre Wilcoxon y Lachenbruch en estricto y por el
desajuste de nivel, documentado en Carn y medido ahora en las corridas de cáncer, mayor en estricto). La definición del criterio estricto en el nulo queda descartada
como causa: coincide con la del observado. Lo que distingue las dos explicaciones es el control "ciclo nulo como observado" con la regla vigente, sobre las corridas de cáncer.

Control "ciclo nulo como observado", neoplasia estricto (20/10/2026; `adhoc_scripts/fcs_null_as_observed.sbatch`, 20 ciclos, 5 por estrato de nivel,
código del FCS sin cambios, nulo de los otros 999 ciclos; conjuntos de al menos 10 genes). Resultados:

- El Lachenbruch de permulación depende del nivel del ciclo (Spearman 0,85 entre el nivel y la proporción de conjuntos con p <= 0,05). Con p <= 0,05, de 3,3 a 4,2 %
  en el estrato bajo, de 4,5 a 5,6 % en el medio, de 7,0 a 7,7 % en el alto y de 6,9 a 9,7 % en el parecido al observado (con p <= 0,25, de 32 a 39 %).
  El observado (8,3 a 8,6 % y 35 a 40 %) coincide con lo que dan los ciclos nulos de su mismo nivel: la inflación del Lachenbruch en el estricto se explica por el
  nivel, sin señal.
- El Wilcoxon de permulación es casi invariante al nivel y algo conservador: 3,6 a 5,2 % con p <= 0,05 en todos los estratos (media 4,2 a 4,4 %; Spearman -0,54).
  El déficit del observado (2,2 a 4,7 %) se reproduce en los ciclos nulos de su nivel (3,6 a 3,8 %).
- Conjuntos "Supported" con la regla vigente por ciclo nulo (sin señal): 6,8 (4 a 10) en el estrato bajo, 7,4 (2 a 13) en el medio, 28,2 (8 a 52) en el alto y
  35,4 (23 a 47) en el parecido al observado. El observado da 73 con mínimo de conjunto 5 y 86 con mínimo 10 (mismo criterio que el control). Con el nivel del
  observado, por tanto, la regla vigente da 35 "Supported" sin señal, del orden del 40 % de los observados.
- Regla propuesta (BH sobre el p de permulación, q < 0,15, por ranking y base): con el Wilcoxon, de 2,2 a 6,0 conjuntos por ciclo nulo (0 a 19), compatible con
  el control del FDR en unas 51 familias; con el Lachenbruch, de 2,0 (estrato bajo) a 82,2 (alto; hasta 300) y 31,6 (parecido al observado; 7 a 69), es decir, no
  controla el FDR cuando el nivel del observado supera al de los ciclos nulos. El observado da 1 conjunto con el Wilcoxon y 102 con el Lachenbruch.
- Límite: 5 ciclos por estrato en neoplasia estricto; la dispersión entre ciclos del mismo estrato es grande.

Las cuatro corridas (20/10/2026; 70 ciclos nulos en total: 20 en cada estricto y 15 en cada no estricto; estratos definidos por el cociente entre el nivel del
ciclo y el del observado de su corrida; conjuntos de al menos 10 genes; mismas reglas aplicadas al observado):

| Corrida | "Supported" del observado | "Supported" por ciclo nulo con nivel >= 0,8 del observado | Ciclos | Proporción |
|---|---|---|---|---|
| Neoplasia estricto | 86 | 35,4 | 5 | 41 % |
| Maligno estricto | 19 | 47,6 | 5 | 251 % |
| Neoplasia no estricto | 45 | 36,7 | 13 | 82 % |
| Maligno no estricto | 72 | 43,2 | 10 | 60 % |

- Wilcoxon de permulación: de 3,7 a 5,3 % con p <= 0,05 en todos los niveles; el observado, de 3,2 a 4,9 %.
- Lachenbruch de permulación: en los estrictos depende del nivel (Spearman 0,85 y 0,88; de 8,2 a 11,5 % con p <= 0,05 a nivel del observado, que da de 8,0 a 8,4 %);
  en los no estrictos queda cerca del nominal (de 4,4 a 5,5 %; observado de 5,5 a 6,0 %).
- Regla propuesta (BH sobre el p de permulación) con ciclos de nivel >= 0,8 del observado: el Wilcoxon da de 1,0 a 4,5 conjuntos por ciclo, compatible con el control del FDR;
  el Lachenbruch da 31,6 y 57,0 en los estrictos y 8,9 y 13,5 en los no estrictos. El observado da con el Wilcoxon 1, 0, 6 y 6 conjuntos (neoplasia estricto, maligno estricto,
  neoplasia no estricto, maligno no estricto), frente a 2,2, 1,0, 3,9 y 4,5 por ciclo nulo.
- Lectura: con un criterio calibrado (Wilcoxon con BH sobre el p de permulación), ninguna corrida muestra enriquecimiento de conjuntos por encima del de ciclos nulos del mismo
  nivel. Los "Supported" de la regla vigente no son evidencia de enriquecimiento; es información de priorización. Esto no demuestra ausencia de señal: la potencia con 1000 ciclos
  y familias grandes es limitada.
- Límites: una semilla; entre 1 y 13 ciclos por nivel; en maligno estricto los nulos a nivel del observado van de 30 a 71 "Supported" por ciclo. Un ciclo de nivel casi nulo
  (0,03 del observado) da 601 conjuntos con el Wilcoxon: el test es inestable con muy pocos genes puntuados.
- Memoria máxima de las tareas del clúster: 5,3 GB.

Alcance del control: evalúa la calibración del estadístico del FCS (el p de permulación del conjunto) cuando el observado y los ciclos nulos difieren en nivel.
No cuestiona la generación del nulo para rasgos de conteo, que se considera adecuada para su pregunta (cuántas posiciones convergentes aparecen por azar con
etiquetados aleatorios). Si el control muestra descalibración, la corrección se plantea en el uso del p de permulación dentro de la regla de decisión
del FCS, no en la forma de generar los ciclos.

Construcción de los scores y p de CAAS (lectura del código, 20/10/2026):

- Score de posición: `CAAS_score` = `asr_path_score` de la fila, en (0, 1] (media sobre los dominios de la media ponderada por PSS de las especies de cada dominio). Mismo valor en el observado y en el nulo.
- Score de gen: `size_adj_max = F(máximo)^n` (`CT_DISAMBIGUATION/local/src/core/scores.py`). F es la distribución empírica de los scores de posición del mismo etiquetado y dirección;
  n es el número de posiciones puntuadas del gen en esa dirección. Es el cuantil del máximo de n posiciones sacadas del pool del propio ciclo, por lo que dentro de un ciclo es casi
  uniforme entre los genes puntuados. No lleva nivel ni propensión; los genes sin posición puntuada no tienen valor.
- `p.emp` de posición: (k+1)/(N+1), con k los ciclos nulos que detectan la posición y tienen un score del lado máximo mayor o igual al observado (los que no la detectan cuentan como menos infinito);
  N sale del roster de ciclos; la familia del BH es el universo posicional (posiciones probadas con detección en el nulo más las observadas, con p = 1 para las no detectadas).
  `p.emp_fact` factoriza P(detectar) por P(score | clase de propensión), con 20 clases.
- El p de permulación de gen que se retiró (`gene_caas_pperm*`) combinaba el score con el número de descubrimientos por gen. El p de gen del control (`fcs_calibrated_ranking_control.sbatch`) es
  p = (1 + ciclos nulos con score del gen >= el observado) / (N + 1), con p = 1 si el gen no se puntúa; el ranking es -log10 p, como en el FCS de RER.

Nivel de gen, neoplasia estricto (ranking global): el 27 % de los genes no se puntúa nunca en 1000 ciclos nulos; AUC de "gen puntuado frente a su propensión" 0,868 en el observado y 0,864 en los
nulos (0,79 a 0,90); entre los genes puntuados, Spearman del score con la propensión 0,05 (observado) y 0,03 (nulos); genes puntuados 2.480 frente a una mediana nula de 1.090 (25 de 1000 ciclos
con al menos el observado). El observado puntúa más genes, pero del mismo tipo que un ciclo nulo.

Variabilidad del estadístico de conjunto (z del p analítico, ranking global): el Wilcoxon da una desviación típica de 1,2 a 1,5 en los ciclos nulos y de 1,3 a 1,6 en el observado; el Lachenbruch, de 1,9 a 3,6
en los estrictos y de 1,2 a 1,3 en los no estrictos. El z del observado y el z medio de los nulos se correlacionan entre conjuntos (Spearman de 0,63 a 0,80 en el Wilcoxon): una parte fija del conjunto
(propensión media de sus genes, desviación típica de 0,8 a 1,3) más un ruido propio de cada ciclo (0,8 a 0,9). Entre el 10 y el 20 % de los conjuntos que llama el Wilcoxon en el observado aparecen en
al menos la mitad de los ciclos nulos.

Pendiente, ver la sección 7 y la nota de traspaso de la memoria: control posicional de `p.emp` (cáncer y carnivoría), prueba de `^n` y del máximo de posición en el score de gen, control del ranking
calibrado por gen (trabajos 7566281 y 7566282 en correfoc).

Hasta ese control, los resultados del FCS de CAAS se tratan como priorización con concordancia entre tests, no como resultado permulado.

## 6. Cambios que implicaría (si se aprueban)

- `conf/enrichment.config`: `fcs_min_genes` 10; `fcs_max_genes` 0 (P1, P2).
- `gui/templates/malignant_def.json` y `neoplasia_def.json`: mismos valores; revisar los campos de la GUI (P1, P2, P14).
- `subworkflows/ENRICHMENT/local/src/fcs_enrich.R`: `fcs_classify_evidence` (q de permulación como criterio, voto bilateral en rankings firmados); valores por
  defecto de `fcs_run_ranking` (P2, P3, P4).
- `subworkflows/ENRICHMENT/local/src/posenrich_enrich.py`: misma regla de significación (P4).
- `subworkflows/ENRICHMENT/local/12.FCS_general_report.Rmd` y `15.Comparison_report.Rmd`: etiquetas, columnas de q de permulación y parámetro `pperm_thr` (P4, P5, P7, P8).
- `subworkflows/RERCONVERGE/rer_gene_lists.nf`: ranking por gen y rankings dirigidos (P8, P9).
- `README.md`: descripción de `fcs_pperm_thr` y de la regla de significación (P5).

## 7. Decisiones pendientes

1. Regla unificada de P4 para RER y CAAS, con CAAS estricto marcado como no concluyente.
2. Ranking por gen en RER (P9) y su sensibilidad.
3. RER solo con el ranking global (P8) y sus consecuencias en los informes 12 y 15.
4. Número de permulaciones de RER (P10).
5. Tamaño máximo de conjunto (P2).
6. Control "ciclo nulo como observado" en las corridas de cáncer (P15): se prepara como `sbatch`; no se ejecuta en el nodo de entrada.

## 8. Reproducibilidad

Los datos y scripts de esta revisión están fuera del repositorio:

- `~/.cc/rer_installed/`: copia de los tres archivos del paquete, funciones deserializadas (`dump/`), `check_null_p.R`, `check_empirical_scale.R`.
- `PhyloPhere_validation/unification/check_fcs_null_equivalence.R`: equivalencia del AUC del nulo.
- `~/.cc/pi3k/`: tablas del FCS de RER y filas de borde de avance de PI3K.
- `~/.cc/caas_fcs/`: tablas del FCS de CAAS de las cuatro corridas.
- `~/.cc/msae210/`: texto del artículo y del suplemento de Redlich.
- `~/.cc/enrich/ctrl_calib4.R` y `ctrl_sum.R`: control "ciclo nulo como observado" de Carn.

## 9. Implementación (WP1, 2026-10-21)

Aplicado en `581b999` + commits de WP0 a WP2: significación por BH sobre el p de permulación (`q.perm`, `sig_perm`; `sig_top` retirado después), Lachenbruch fuera del FCS, RER solo con el ranking global firmado `sign(Rho) * -log10(P)` y su nulo `corStat`, reducción de redundancia a posteriori (`fcs_cluster.R`) e informes 12 y 15 reescritos. Registro completo, comprobaciones y diferencias abiertas: `docs/fases/WP1_fcs_rer.md`. Lo no aplicado: subir las permulaciones de RER por defecto (P10; falta medir `RER_CONT`).

Medido al aplicarlo: con la regla nueva, las cuatro tablas reales dan 0, 0, 6 y 3 conjuntos significativos (neoplasia estricto, maligno estricto, neoplasia completo, maligno completo), frente a 1, 0, 6 y 6 de la sección 5, que salieron del script de control; quitar copias y casi copias (Jaccard >= 0,9) retira el 10,5 a 10,8 % de las filas sin cambiar esos recuentos.

Medición de RER_CONT (P10, 2026-10-21): el coste de las permulaciones de RER es lineal en N (unos 2,1 a 2,2 s por permulación con 16.109 genes; 4.000 genes con 10.000 permulaciones: 5.470 s y 2,8 a 3,3 GB según el tamaño de lote), la memoria crece con el objeto (24·G·N bytes) y el proceso usa un solo hilo; detalle y extrapolación señalada en `docs/fases/WP1_fcs_rer.md`.



## WP9: control dirigido del estadístico base (2026-10-23)

Se comparan las seis variantes en marino completo/strict (bloques, 100 ciclos, tres rankings) y neoplasia estricto histórico (1000, global, Hallmark+KEGG). El control usa gene_p productivo y AUC centrado con comparación exacta; añade mitades disjuntas porque la tasa marginal LOO es casi una identidad de rangos. F(max)^n reduce ρ(score,n) de 0,225 a −0,020 en neoplasia, con tasas holdout de rutas de nivel observado ~8,25 % para las seis variantes: cambiar exponente no elimina la asociación con nivel. No es una estimación de FDR. Se recomienda mantener provisionalmente el default, sin cerrarlo antes de A2. Resultados, universos, medidas y límites: `docs/fases/WP9_resultados_dirigidos_2026-10-23.md`; comandos y pendientes en `WP9_estadistico_gen.md`. Sin cambios del cálculo productivo ni corrida 10k ejecutada.
