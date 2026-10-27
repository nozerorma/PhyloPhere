# POSENRICH: diseño de preguntas, universo y calibración

Estado: diseño acordado, sin implementar. Sustituye el contrato inferencial de `POSENRICH_MODULAR_CONTRATO.md` y `POSENRICH_MODULAR_REVISION.md`; las componentes descriptivas existentes (registro, caracterización, proximidad) se conservan y se adaptan según la sección 8. Los puntos abiertos están en la sección 10.

Las cifras citadas se midieron en el run `malignant_def/malignant_prevalence_complete` (10.000 ciclos, correfoc). Son comprobaciones de las entradas, no resultados de enriquecimiento.

## 1. Qué se pregunta

Una posición es **positiva** si tiene score CAAS canónico > 0. Todas las preguntas se refieren a posiciones positivas y a su relación con conjuntos de posiciones (colecciones).

| Nivel | Pregunta | Región | Inferencia |
|---|---|---|---|
| Caracterización | ¿Cómo se reparten las positivas entre las capas? | capas y términos | ninguna (tablas con denominadores) |
| 1, no acotado | ¿Un gen acumula más positivas que su cuota de U? | el gen entero | hipergeométrico exacto (sección 5) |
| 2, acotado | ¿Las positivas se concentran en un término dentro de cada gen, más allá de lo que el gen ya acumula? | término ∩ gen | estadístico agregado calibrado con el nulo (sección 6); tabla por gen descriptiva |
| Proximidad | ¿Las detecciones están más cerca de los sitios o regiones de su proteína? | misma proteína | descriptiva hasta validar el control |

El nivel 2 condiciona en el nivel 1: se fija el número de positivas de cada gen y se mira solo su reparto dentro del gen. Así los dos niveles no cuentan dos veces la misma señal. El contraste global de un conjunto contra U mezcla ambas fuentes de señal y no se usa como test. Sus recuentos y proporciones se publican como descripción.

## 2. Universo U

U es el universo de oportunidad de una detección: las posiciones probadas por CAAStools (`caastools/background.output`) en los genes del fondo limpio, que puntúan > 0 en algún ciclo del nulo filtrado (`caas_permulation/null_positions.parquet`), más las positivas observadas. Es el universo que ya usan SCORING y Accumulation, y debe construirse con un único constructor compartido.

Medido: el fondo probado tiene 6,5 M de posiciones y U tiene 178.365 (2,7 %). Las 27.080 positivas observadas están todas en la intersección entre el nulo y lo probado, es decir, U no recorta la señal observada. La tasa base de positivas en U es 0,152.

La detectabilidad por posición (`nd/N`) es muy heterogénea dentro de U: mediana 0,063, percentil 5 en 0,0005, percentil 95 en 0,25. Un contraste que trate las posiciones de U como equivalentes es un nulo débil, y es la razón de calibrar con los ciclos.

Posiciones fuera de U no se tratan como negativas. No son evaluables para frecuencia, pero sí son blancos válidos para proximidad.

## 3. Entradas

- `caas_permulation/blocks/blk_*/positions.parquet`: (gene, cycle, position, `all`) por posición puntuada y ciclo. `all` es el mejor lado (`fmax(top, bottom)`, `null_merge.py`), igual que el score global del observado. Solo hay score global: el nulo no calibra por lado.
- `caas_permulation/null_positions.parquet` y `blocks/*/detections.parquet`: `nd` y primer ciclo por posición.
- `scoring/position_scores.tsv`: scores y lado del observado. Solo 35 de 27.080 posiciones (0,13 %) puntúan en ambos lados, de modo que el lado se conserva en el observado y se puede describir la composición top/bottom por término.
- `caastools/background.output`, `universes/cleaned_background_main_*.txt`, `filtered_discovery.tsv`.
- Colecciones: ver sección 7.

## 4. Capas

Se conservan: `pfam_domains`, `pfam_clans`, `ucr_core`, `ucr_flank`, `fubar_positive`, `fubar_purifying`, `fade_top_sig`, `fade_bottom_sig` y las listas de sitios (`custom`). VEP, COSMIC y PrimateAI se mantienen como caracterización cruda.

Se retiran como colecciones de test: ortogrupo, bin genómico, COSMIC/PAI3D por ortogrupo y capas de selección por ortogrupo. Un conjunto por ortogrupo es un gen, y ese nivel es Accumulation. Las ventanas de 1 Mb solo se usan como atributo de posición para la proximidad genómica.

Una colección tiene valor para el test cuando su vocabulario se repite entre genes (dominio, clan, región conservada, sitio). Así se puede agrupar entre genes y estratificar por gen.

## 5. Nivel 1: Accumulation

Se mueve a POSENRICH. Para cada dirección (all, top, bottom): sea M = |U|, K el número de posiciones positivas distintas de la dirección y n_g las posiciones de U del gen g. El recuento observado del gen sigue un hipergeométrico `Hyper(M, n_g, K)`, y el p exacto de cola superior se ajusta con BH sobre los genes probados. La implementación actual (`accumulation_test.py`) se conserva. Sus entradas no incluyen scores, así que no depende de SCORING.

Limitación: el hipergeométrico trata las posiciones de U como equivalentes. El nivel 1 no se calibra con el nulo en esta versión. Una variante calibrada contra los recuentos por (gen, ciclo) del nulo queda anotada como posible mejora posterior; no forma parte del diseño actual.

## 6. Nivel 2: concentración dentro del gen

**Estadístico (esquema).** Un solo estadístico graduado por score. Cada posición i de U aporta un peso w_i: su score global (`all`) si es positiva y 0 si no se detecta. Para un término S y cada gen g con posiciones dentro y fuera de S, sea W_g la suma de pesos del gen, n_g sus posiciones en U y m_g las que están en S. El residuo del gen es

    r_g = (suma de w_i sobre S ∩ g) − W_g · m_g / n_g

y el estadístico principal es la suma de r_g / s_g sobre genes, donde s_g es la desviación típica de r_g en los ciclos del nulo. La suma simple de r_g se calcula como sensibilidad. Se calcula igual sobre el observado y sobre cada ciclo del nulo, con los márgenes (W_g) de ese ciclo.

Comprobación previa a publicar p/q: con el control "ciclo nulo como observado" se compara la tasa de falsos positivos de las dos versiones. Si la normalizada no mantiene el nivel nominal (por ruido en genes con pocas posiciones), se revisa antes de ver resultados del observado.

**Exclusión de genes.** Un gen entra en el estadístico solo si tiene posiciones dentro y fuera del término y una s_g estimable (no nula, con peso distinto de cero en una fracción suficiente de ciclos; el umbral se fija antes de ver resultados del observado). Los demás se excluyen de forma explícita: cada término publica la lista de genes excluidos con su motivo (sin posiciones fuera del término, sin positivas, s_g nula o inestable). Nunca se imputa un residuo.

**Calibración.** p = (1 + nº de ciclos con estadístico ≥ observado) / (N + 1), con N = 10.000 y suelo 1e-4. Se calcula por bloque (500 ciclos, 9,3 M de filas por bloque) y los recuentos de excedencias se suman. Universo: U. BH por familia (colección y ranking).

**Qué cubre el nulo sin parámetros nuevos:** la detectabilidad por posición, la co-ocurrencia de detecciones dentro de un gen y el filtro de genes del nulo.

**Informatividad.** Solo aportan los genes con posiciones dentro y fuera del término y con al menos una positiva. Medido: la mediana de un gen es 10 posiciones de U con 2 positivas; 5379 de 10.797 genes tienen al menos 2 positivas. La tabla por gen es descriptiva casi siempre, y el agregado entre genes es la unidad de inferencia. Los términos que cubren casi todo un gen (por ejemplo una lista de sitios definida por gen) son equivalentes al nivel 1 y se marcan como tales.

**Frecuencia e intensidad.** El estadístico graduado responde a la vez si el término concentra detecciones y si concentra detecciones fuertes: con pesos binarios es un recuento, y con el score es el estadístico graduado. Hay un solo p por término y la familia de BH no se duplica. El p no separa frecuencia de intensidad. Esa separación se hace con dos tamaños de efecto descriptivos, sin p:

- **Frecuencia:** proporción de positivas entre las posiciones de U del término, frente a la del resto de sus genes.
- **Intensidad:** diferencia de score medio entre las positivas del término y las positivas del resto de sus genes.

**Ranking.** El peso es el score crudo (`all`), que tiene equivalente en cada ciclo del nulo. `p.emp` no se usa en el estadístico (no se guarda por ciclo) y se mantiene en las vistas descriptivas.

## 7. Colecciones y registro

- El universo de test es U. La geometría de las anotaciones no se recorta a U: los sitios y regiones fuera de U siguen siendo blancos para proximidad.
- Valores únicos por posición (`fubar_tested`, ventana de 1 Mb, región UCR, llamada FUBAR) son columnas de `positions.tsv`, no filas de membresía.
- Las regiones (Pfam, UCR) se guardan como intervalos por (gen, término) y se expanden solo dentro de U para localización. Los sitios se guardan por posición.
- Medido: `members.tsv` tiene 3,7 GB, dominado por Pfam (6,6 M de filas), `fubar_tested` (6,4 M), ventanas (6,4 M), UCR core (5,5 M), FUBAR purificador (3,0 M) y clanes (2,5 M). Los sitios son unas pocas miles de filas.
- Ventanas genómicas: el tamaño es el de `scoring_window_size_bp` (`conf/scoring.config`, 1.000.000 por defecto), no un valor fijo en el constructor, y el identificador de colección no lleva el tamaño. Solo son un atributo de posición; la proximidad genómica mide distancias en nt.
- Lista de genes de U: el fondo limpio de CT_POSTPROC (`postproc/cleaned_backgrounds/cleaned_background_main.txt`, 11.530 genes en el run). Es el fondo de CAAStools (`caastools/background_genes.output`, 11.710 genes) sin los genes que el filtrado de genes retira del observado (`postproc/gene_filtering/removed_genes_summary.tsv`, categorías Extreme y Dubious; 180 en el run), que por construcción no pueden tener positivas observadas. Con esta lista |U| = 165.468 en el run; con la del fondo de CAAStools, 178.365 (12.897 posiciones más del nulo y ninguna positiva observada). En una ejecución integrada, el canal `pp_cleaned_bg` de `main.nf` alimenta a todos los consumidores. La ejecución de SCORING del run que escribió `null_universe_summary.tsv` usó la lista de 11.710 genes; qué parámetro o ruta la trajo no se ha verificado. `posenrich_universe.py` recibe la lista como entrada explícita y escribe su huella y |U| en el manifiesto. SCORING registra la huella de su lista en `null_universe_summary.tsv` y el registro compara ambas huellas, con aviso si difieren.

## 8. Cambios en las componentes existentes

- **Registro:** cambios de la sección 7. Con un único constructor de U compartido con Accumulation y SCORING.
- **Caracterización:** sin cambios de contrato. Pasa a leer los atributos desde columnas.
- **Localización y asociación global:** el control de redistribución por proteína y el AUC global descriptivo se sustituyen por el estadístico de la sección 6, calibrado con los ciclos. Las salidas con `p` retenido desaparecen.
- **Proximidad:** se conserva. El control de distancias sigue sin validar.
- **Informes:** el informe modular consume las nuevas tablas. El informe histórico 14 sigue omitido.

## 9. Accumulation fuera de SCORING

Decisiones: Accumulation se mueve a POSENRICH y se elimina de SCORING; no se integra en los informes existentes; se elimina `flag_accum` de todas partes; `fcs_stats_accum` se retira.

Qué toca (inventario, por verificar al implementar):

- `main.nf`: llamada a CT_ACCUMULATION (l.287, 520-556), canal a SCORING (l.727, 766), `fcs_stats_accum` a ENRICHMENT (l.845, 907); parámetros `ct_accumulation`, `accumulation_background_input`, `accumulation_caas_input`.
- `workflows/scoring.nf`: entrada `accum_ch`, `scoring_accum_dir`.
- `scoring_compute.nf` y `scoring_compute.R`: sección 4b, columnas `accum_p`, `accum_fdr`, `accum_significant` (all, top, bottom), banderas `flag_accum*`, `accum_available*`, correlaciones, `fcs_stats_accum.tsv`, parámetro `accumulation_fdr`.
- Informe 11 (`11.Scoring_report.Rmd`, `report_evidence.R`, `scoring_report.nf`).
- ENRICHMENT: `scoring_enrichment.nf`, `workflows/enrichment.nf` (`fcs_stats_accum`), informes 12, 13, 14 y 15.
- `rer_gene_lists.nf` y `accum_gene_lists.nf`.
- `genomic_window_accum_gene_caas.tsv` es una tabla para un gráfico por ventanas del informe 11, no un test, y queda fuera.
- GUI, plantillas, ayuda y tests de las partes anteriores.

Se hace como fase propia, con tests, después de cerrar el diseño. La integración de POSENRICH en el informe comparativo se decide cuando POSENRICH esté construido.

## 10. Puntos abiertos

1. **Umbral de estimabilidad de s_g.** Fijar, antes de ver resultados del observado, la fracción mínima de ciclos con peso distinto de cero en el gen para que s_g se considere estimable (sección 6, exclusión de genes).
2. **Lado.** El nulo guarda solo `all`. Una bandera de lado de un byte permitiría calibrar por dirección. No se hace ahora; el lado se describe en el observado.
3. **Nivel 1 calibrado.** Anotado como mejora posible (sección 5). No forma parte del diseño actual.
4. **Alineación de SCORING con la lista de U.** Decidido: la lista es el fondo limpio de CT_POSTPROC. Falta implementar el registro de la huella en SCORING y la comparación en el registro (fase P1b del plan), y averiguar por qué la ejecución de SCORING del run recibió la lista de 11.710 genes.
5. **Validación.** Ciclo nulo como observado: tasa de falsos positivos por término y colección, y comparación de la versión normalizada con la suma simple, antes de publicar p/q.

## 11. PP2A

Ensayo propio, fuera de POSENRICH y de la GUI de PhyloPhere (sección PP2A de la pestaña Enrichment y `posenrich_custom_file` de las plantillas). Pregunta principal: proximidad de las detecciones a los sitios de su proteína. Pregunta secundaria: solapamiento exacto.

Medido con los sitios mapeados (`20261024/ppga/pp2a`): de 5284 posiciones únicas de sitios, 3815 están en genes del universo, 3772 son posiciones probadas y 76 caen en U (2,0 %; 22 positivas). El fondo medido de Sapkota tiene 22.383 posiciones, 254 en U. El solapamiento exacto tiene poca potencia por esto. La proximidad usa los 3772 sitios como blancos y no está limitada por U. Los sitios mapeados quedan simétricamente por debajo en U respecto a sus vecinos a ±1 y ±2 residuos (76 frente a 90-106), compatible con que las columnas S/T/Y se detecten menos. No se ha comprobado con el residuo del alineamiento.

La proyección a la referencia productiva (accession, isoforma, péptido único, MAP) la hace `bin/map_pp2a_to_caas.py`; los registros que no validan quedan como pendientes, nunca como negativos.

Pendiente de diseñar: el control de proximidad que respete proteína, densidad de sitios y cobertura medida, y el fondo por estudio.

## 12. Fuera de alcance

Lachenbruch y el hipergeométrico sin condicionar en el gen como tests; conjuntos por ortogrupo; calibración por lado; validación inferencial a escala; FCS de posiciones.
