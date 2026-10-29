# POSENRICH: diseño de preguntas, universo y calibración

Estado: implementado (universo, registro, caracterización, niveles 1 y 2, solapamientos y proximidad, informe 14). Falta la validación del nivel 2 a escala con 10.000 ciclos (sección 10). Es el documento vigente de POSENRICH; las fuentes de posiciones y la preparación de PP2A están en `POSITION_COLLECTIONS.md`.

Las cifras citadas se midieron en el run `malignant_def/malignant_prevalence_complete` (10.000 ciclos, correfoc). Son comprobaciones de las entradas, no resultados de enriquecimiento.

## 1. Qué se pregunta

Una posición es **positiva** si tiene score CAAS canónico > 0. Todas las preguntas se refieren a posiciones positivas y a su relación con conjuntos de posiciones (colecciones).

| Nivel | Pregunta | Región | Inferencia |
|---|---|---|---|
| Caracterización | ¿Cómo se reparten las positivas entre las capas? | capas y términos | ninguna (tablas con denominadores) |
| 1, no acotado | ¿Un gen se queda con una proporción de las detecciones mayor que en ciclos al azar? | el gen entero | proporción del gen calibrada con el nulo (sección 5) |
| 2, acotado | ¿Dentro de los genes que llevan un término, las positivas caen en el término más que en el resto del mismo gen? | término ∩ gen | estadístico agregado calibrado con el nulo (sección 6); tabla por gen descriptiva |
| Solapamientos y proximidad | ¿Qué detecciones coinciden con un sitio o región de su proteína, o están cerca? | misma proteína | listas descriptivas, sin p |

El nivel 2 condiciona en el gen: se fija el número de positivas de cada gen (en el observado y en cada ciclo) y se mira solo su reparto dentro del gen. Así los dos niveles no cuentan dos veces la misma señal. El contraste global de un conjunto contra U mezcla ambas fuentes de señal y no se usa como test. Sus recuentos y proporciones se publican como descripción.

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

Reglas de las entradas. El score de posición es el score CAAS canónico (`core.scores.position_score`): solo es positivo si el valor es positivo; cero, negativo o ausente no son un score. `position_scores.tsv` tiene una fila por (Gene, Position, side) con score > 0, y el `p.emp` es compartido entre lados (máximo): las filas top y bottom no son dos pruebas direccionales. `filtered_discovery.tsv` tiene una fila por (Gene, Position, side) evaluada por CT, con `asr_path_score` ≤ 0 cuando se evaluó sin detección positiva. Una posición sin score no se convierte en cero ni en negativo si no se demuestra que se evaluó. La posición es cero-based en el alineamiento recortado de CAAS, con una unidad por posición biológica, sin multiplicar por hipótesis, transcritos ni lados.

## 4. Capas

Se conservan: `pfam_domains`, `pfam_clans`, `ucr_core`, `ucr_flank`, `fubar_positive`, `fubar_purifying`, `fade_top_sig`, `fade_bottom_sig` y las listas de sitios (`custom`). VEP, COSMIC y PrimateAI se mantienen como caracterización cruda.

Se retiran como colecciones de test: ortogrupo, bin genómico, COSMIC/PAI3D por ortogrupo y capas de selección por ortogrupo. Un conjunto por ortogrupo es un gen, y ese nivel es Accumulation. Las ventanas de 1 Mb solo se usan como atributo de posición para la proximidad genómica.

Una colección tiene valor para el test cuando su vocabulario se repite entre genes (dominio, clan, región conservada, sitio). Así se puede agrupar entre genes y estratificar por gen.

## 5. Nivel 1: Accumulation

Se calcula en POSENRICH (`accumulation_test.py`). Sea K el número de posiciones positivas observadas en U y c_g las del gen g. La proporción del gen es c_g / K. En cada ciclo b del nulo se calcula la misma proporción con el total de ese ciclo, c_g(b) / K_b, y

    p_g = (1 + nº de ciclos con proporción ≥ la observada) / (N + 1)

con suelo 1/(N + 1). Solo se prueban los genes con al menos una positiva observada; q es BH sobre ellos. Los ciclos sin ninguna positiva en U no definen proporción y se descartan (se cuentan en el manifiesto).

Motivo de comparar proporciones y no recuentos: el observado tiene más positivas que un ciclo típico (en `malignant_prevalence_complete`, 27.080 frente a 15.142 de media por ciclo dentro de U), y el recuento bruto trasladaría ese exceso global a todos los genes. El nulo recoge lo que un hipergeométrico ignora: la detectabilidad desigual de las posiciones y la co-ocurrencia de detecciones dentro de un gen.

Accumulation y el motor del nivel 2 leen la misma tabla `universe.tsv.gz` de `posenrich_universe.py` (una fila por posición de U con `nd`, `observed_positive`, `observed_top` y `observed_bottom`) y los mismos bloques del nulo. Las posiciones observadas son las filas con score positivo de `filtered_discovery.tsv`, de modo que el nivel 1 no depende de SCORING. Los genes de la lista sin posiciones en U no tienen fila.

El nulo guarda solo el mejor lado (`all`), de modo que el test se aplica a detecciones en cualquier lado. Los recuentos observados en los lados top y bottom se publican como columnas descriptivas.

Salidas: `accumulation/accumulation_results.tsv` (`gene`, `n_positions`, `positives`, `positives_top`, `positives_bottom`, `share`, `null_share_mean`, `null_share_sd`, `share_ratio`, `p`, `q`, `n_null_cycles`) y `accumulation_manifest.json` (con la distribución de positivas por ciclo del nulo).

## 6. Nivel 2: concentración dentro del gen

**Estadístico.** Para un término S y un gen g con posiciones dentro y fuera de S, sea n_g sus posiciones en U, m_g las que están en S y k_g sus positivas. Si las positivas de g se repartieran por igual entre sus posiciones, caerían k_g m_g / n_g en S. El residuo del gen es

    r_g = (positivas de g en S) − k_g m_g / n_g

y el estadístico del término es la suma de r_g sobre sus genes elegibles. Una posición pesa lo mismo (es positiva o no); el score solo interviene en la intensidad descriptiva. Se calcula igual sobre el observado y sobre cada ciclo del nulo, con los k_g de ese ciclo.

**Elegibilidad, fijada por la anotación y por U, nunca por las positivas observadas.** Un gen es elegible para un término si tiene posiciones dentro y fuera del término. Un término se prueba si tiene miembros en al menos `--min-genes` genes de U (2 por defecto): un término en un solo gen es un test dentro de ese gen, no una concentración entre genes. Los demás términos se listan con su estado (`below_min_genes`, `not_identifiable`) y sin p.

**Calibración.** p = (1 + nº de ciclos con estadístico ≥ observado) / (N + 1), con suelo 1/(N + 1). BH por colección. Se calcula por bloque (500 ciclos) y los recuentos de excedencias se suman, de modo que el resultado no depende de cómo se reparten los ciclos (comprobado con 1, 7 y 42 bloques en `tests/test_posenrich_within_gene.py`). Universo: U.

**Implementación.** `posenrich_within_gene.py`. Entradas: `universe.tsv.gz`, la tabla observada con `CAAS_score`, el directorio del registro (`--registry-dir`) y los bloques `blk_*` (`positions.parquet` y `cycles.txt`). Salidas: `within_gene_results.tsv` (un término por fila: `statistic`, `expected_positives_members`, `p`, `q`, tamaños de efecto, recuento de genes elegibles y excluidos, `status`), `within_gene_by_gene.tsv` (un par término-gen por fila, con su residuo y el motivo de exclusión) y `within_gene_manifest.json`. `--holdout K` reserva K ciclos del nulo como observado y escribe `holdout_fpr.tsv`.

**Qué cubre el nulo sin parámetros nuevos:** la detectabilidad por posición, la co-ocurrencia de detecciones dentro de un gen y el filtro de genes del nulo.

**Frecuencia e intensidad.** El p no separa frecuencia de intensidad. Esa separación se hace con dos tamaños de efecto descriptivos, sin p:

- **Frecuencia:** proporción de positivas entre las posiciones de U del término, frente a la del resto de sus genes elegibles.
- **Intensidad:** diferencia de score medio entre las positivas del término y las positivas del resto de sus genes.

**Ranking.** El test es sobre posiciones positivas; `p.emp` no se usa en el estadístico y se mantiene en las vistas descriptivas.

## 7. Colecciones y registro

- El universo de test es U. La geometría de las anotaciones no se recorta a U: los sitios y regiones fuera de U siguen siendo blancos para proximidad.
- El registro guarda las pertenencias según lo que son, no como una fila por término y posición (`posenrich_members.py` las expande cuando un consumidor las pide):
  - **Atributos** (un valor por posición): columnas de `positions.tsv`: `ucr_core`, `ucr_flank`, `ucr_region`, `ucr_variability`, `fubar_call` (positive, purifying, neutral o vacío si no se probó), `fade_top_sig`, `fade_bottom_sig` y `window`. Las colecciones `ucr_core`, `ucr_flank`, `fubar_positive`, `fubar_purifying`, `fubar_tested`, `fade_top_sig`, `fade_bottom_sig` y `genomic_windows` se derivan de ellas.
  - **Regiones** (Pfam: dominios y clanes): `regions.tsv`, un registro por tramo de columnas CAAS consecutivas de un término en un gen (inicio y fin inclusivos). Los huecos de la referencia parten el tramo, de modo que la expansión da exactamente las posiciones que daba la construcción por posición.
  - **Sitios** (listas personalizadas): `sites.tsv`, una fila por término, gen y posición.
- El registro no guarda la tabla expandida de pertenencias (una fila por término y posición). Las colecciones de cobertura (`fubar_tested` y `genomic_windows`) no cuentan como anotación en `n_positions_annotated`.
- Tamaño de ventana: `--window_size_bp`, con el valor de `scoring_window_size_bp`; la colección se llama `genomic_windows`, sin el tamaño. El nombre de cada ventana cuenta su inicio en Mb, kb o bp (1 Mb: `chr1_12M`; 100 kb: `chr1_12300k`).
- Motivo de la forma compacta: la tabla expandida del run ocuparía 3,7 GB, dominada por Pfam (6,6 M de filas), `fubar_tested` (6,4 M), ventanas (6,4 M), UCR core (5,5 M), FUBAR purificador (3,0 M) y clanes (2,5 M). Los sitios son unas pocas miles de filas.
- Lista de genes de U: el fondo limpio de CT_POSTPROC (`postproc/cleaned_backgrounds/cleaned_background_main.txt`). Es el fondo de CAAStools (`caastools/background_genes.output`) sin los genes que el filtrado de genes retira del observado (`postproc/gene_filtering/removed_genes_summary.tsv`, categorías Extreme y Dubious), que por construcción no pueden tener positivas observadas. `main.nf` la resuelve desde `pp_cleaned_bg` (CT_POSTPROC en vivo) o, en una ejecución con resultados precomputados, desde `scoring_background_input` o `accumulation_background_input`; `--background_input` es la lista cruda de CAAStools y nunca se toma como fondo limpio. `posenrich_universe.py` recibe la lista como entrada explícita y escribe su huella y |U| en el manifiesto. SCORING registra la huella de su lista en `null_universe_summary.tsv` y el registro compara ambas huellas, con aviso si difieren.

## 8. Componentes descriptivas

- **Registro:** sección 7, con el constructor de U compartido con Accumulation y SCORING. Las fuentes se leen con `position_sources.py`.
- **Caracterización:** lee los atributos desde las columnas del registro y no hace inferencia.
- **Solapamientos y proximidad:** proteica (residuos de referencia) y genómica (nucleótidos), separadas. Sin p. Las listas de solapamientos (identidad exacta de la posición con un miembro del término) y de posiciones cercanas (dentro del radio mayor, con el miembro más próximo) se publican por detección y por gen (`proximity_protein_hits.tsv`, `proximity_protein_by_gene.tsv`).
- **Informe 14:** consume las tablas de todas las componentes (sección 13).

## 9. Accumulation y SCORING

Accumulation es el nivel 1 de POSENRICH. SCORING no calcula ni lee Accumulation: no tiene columnas `accum_*` ni banderas `flag_accum*`, y no existe `fcs_stats_accum`. ENRICHMENT y los informes 11, 12 y 13 no lo consumen. `genomic_window_accum_gene_caas.tsv` es una tabla del informe 11 para un gráfico por ventanas, no un test.

## 10. Puntos abiertos

1. **Validación del nivel 2 a escala.** Ciclo nulo como observado con los 10.000 ciclos y varios conjuntos de términos: tasa de falsos positivos por término y colección antes de publicar p y q como resultado. Se lanza por sbatch.
2. **Lado.** El nulo guarda solo `all`. Una bandera de lado de un byte permitiría calibrar por dirección; no se hace, y el lado se describe en el observado.
3. **Control de proximidad.** Pendiente de diseñar (sección 11).

## 11. PP2A

El ensayo PP2A propio no está construido. Mientras tanto, la lista de sitios entra como conjunto personalizado (`posenrich_custom_file`, campo "Custom position sets" de la pestaña POSENRICH, sin mención a PP2A), que POSENRICH caracteriza y prueba como cualquier término; el ensayo propio sustituirá ese uso. Pregunta principal: proximidad de las detecciones a los sitios de su proteína. Pregunta secundaria: solapamiento exacto.

Medido con los sitios mapeados (`20261024/ppga/pp2a`): de 5284 posiciones únicas de sitios, 3815 están en genes del universo, 3772 son posiciones probadas y 76 caen en U (2,0 %; 22 positivas). El fondo medido de Sapkota tiene 22.383 posiciones, 254 en U. El solapamiento exacto tiene poca potencia por esto. La proximidad usa los 3772 sitios como blancos y no está limitada por U. Los sitios mapeados quedan simétricamente por debajo en U respecto a sus vecinos a ±1 y ±2 residuos (76 frente a 90-106), compatible con que las columnas S/T/Y se detecten menos. No se ha comprobado con el residuo del alineamiento.

La proyección a la referencia productiva (accession, isoforma, péptido único, MAP) la hace `bin/map_pp2a_to_caas.py`; los registros que no validan quedan como pendientes, nunca como negativos.

Pendiente de diseñar: el control de proximidad. Debe contrastar si las detecciones están más cerca de los sitios de su misma proteína de lo esperado bajo un control que respete proteína, elegibilidad, cobertura medida y densidad de sitios, con radios predefinidos y una sensibilidad que excluya las coincidencias exactas. Incluye también el fondo por estudio: el fondo medido de Sapkota no es el de Kauko ni el de Stewart, ni equivale a todas las S/T/Y. La distancia proteica en residuos no se mezcla con la genómica ni se presenta como cercanía 3D.

## 12. Fuera de alcance

Lachenbruch y el hipergeométrico sin condicionar en el gen como tests; conjuntos por ortogrupo; calibración por lado; validación inferencial a escala; FCS de posiciones.

## 13. Módulo POSENRICH y estructura del código

POSENRICH es un módulo propio, separado de ENRICHMENT: `workflows/posenrich.nf`, `subworkflows/POSENRICH/` y `conf/posenrich.config`. Contiene el universo, el registro, la caracterización, el nivel 1 (Accumulation), el nivel 2 (`posenrich_within_gene.py`), los solapamientos y la proximidad, la generación de las tablas Pfam y UCR y el informe 14. Su interruptor es `--posenrich`.

ENRICHMENT conserva FCS, AMI/DOMINO y los informes 12, 13 y 15. El orden de datos es SCORING, POSENRICH, ENRICHMENT. El informe 15 recibe la caracterización de POSENRICH como entrada opcional; sin ella omite la sección de posiciones interesantes.

**Informe 14** (`14.POSENRICH_report.Rmd`). Secciones: universo y cobertura, caracterización cruda de las positivas por capa, nivel 1 por dirección, nivel 2 (términos con p, q, tamaños de efecto y genes excluidos con su motivo), proximidad, métodos y límites. La GUI tiene una sola especificación de informe regenerable para POSENRICH.

**Código** (`subworkflows/POSENRICH/local/src/`):

| Script | Función |
|---|---|
| `position_sources.py` | Módulo de lectura de fuentes (Pfam, UCR, FUBAR, FADE, ventanas, listas personalizadas) con proyección a columnas CAAS a través del MAP. No es un programa ni escribe archivos |
| `posenrich_registry.py` | Registro: posiciones con atributos como columnas, regiones, sitios, colecciones, cobertura y manifiesto |
| `posenrich_members.py` | Expande atributos, regiones y sitios a pertenencias cuando un consumidor las pide |
| `posenrich_universe.py` | Universo U (`universe.tsv.gz`) y su manifiesto |
| `accumulation_test.py` | Nivel 1 (proporción por gen calibrada con el nulo) |
| `posenrich_within_gene.py` | Nivel 2 |
| `posenrich_characterize.py` | Caracterización descriptiva |
| `posenrich_proximity.py` | Proximidad proteica y genómica, listas de solapamientos y de posiciones cercanas |

Las colecciones por ortogrupo, COSMIC y PAI3D por ortogrupo y las exportaciones a GMT no forman parte del módulo (sección 4). Las pruebas están en `tests/test_posenrich_*` y `tests/test_posenrich_report.R`; el smoke de Nextflow es `tests/test_posenrich_smoke.py`.

## 14. Validación a escala

`adhoc_scripts/posenrich_scale_check.sbatch` ejecuta los scripts del módulo directamente sobre los ficheros de un run terminado (sin Nextflow) y escribe todo bajo `.tmp/posenrich_check/job_<id>/`; el directorio del run solo se lee. Mide tiempo y memoria máxima de cada paso, tamaño de las tablas del registro, la tasa de falsos positivos del nivel 2 con los 10.000 ciclos (`--holdout`), el nivel 1 (genes probados, q < 0,1 y concordancia de los recuentos con la Accumulation publicada por el run), y renderiza el informe 14.

**Cómo leer el resultado.**

- **Universo.** |U| y positivas observadas coinciden con las del run. Si la lista de genes del registro y la de `null_universe_summary.tsv` difieren, `gene_list_check` del registro da `mismatch`: es la detección funcionando, no un fallo, y apunta a que SCORING usó otra lista.
- **Nivel 1.** `n_positions` y `positives` iguales a `PoolPositions` y `ActualCount` de la Accumulation publicada en todos los genes presentes en ambas tablas. Los p no son comparables: la publicada es hipergeométrica.
- **Registro.** Compara el tamaño de `positions.tsv`, `regions.tsv` y `sites.tsv` con el de `positions.tsv` y `members.tsv` del registro anterior del run (611 MB y 3,7 GB). Si la memoria máxima del paso `registry` o `characterization` se acerca al límite del trabajo, esos procesos necesitan optimización antes de usarlos en otros runs.
- **Nivel 2, falsos positivos.** Los ciclos reservados son nulos por construcción, así que la tasa por colección a α = 0,05 y 0,01 debe ser compatible con α. Los términos comparten posiciones y cada ciclo desplaza todos sus términos, de modo que el error de la tasa lo marca el número de ciclos reservados (1.000: desviación típica de unos 0,007 a α = 0,05 y 0,003 a α = 0,01, que subestima el error real). Como criterio heurístico, no demostrado: una tasa entre α − 0,015 y α + 0,015 (a α = 0,05) y entre 0 y 0,015 (a α = 0,01) es compatible con el nivel nominal; una tasa claramente superior (por encima de 0,07 a α = 0,05 o de 0,02 a α = 0,01) indica un estadístico anticonservador y sus p y q no deben interpretarse; una tasa inferior es conservadora y aceptable.
- **Nivel 2, observado.** Número de términos probados, por debajo del mínimo de genes y no identificables, genes elegibles y excluidos por término, y el suelo de p (1/(N+1) con N = ciclos de calibración). Un término con q < 0,1 solo es una señal de priorización si la validación anterior es satisfactoria.
- **Informe 14.** Renderiza sin errores, con el panel de falsos positivos y las seis secciones.

