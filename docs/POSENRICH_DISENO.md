# POSENRICH: diseño de preguntas, universo y calibración

Estado: implementado (universo, registro, caracterización, niveles 1 y 2, proximidad e informe 14). Falta la validación del nivel 2 a escala con 10.000 ciclos (sección 10). Es el documento vigente de POSENRICH; las fuentes de posiciones y la preparación de PP2A están en `POSITION_COLLECTIONS.md`.

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

Reglas de las entradas. El score de posición es el score CAAS canónico (`core.scores.position_score`): solo es positivo si el valor es positivo; cero, negativo o ausente no son un score. `position_scores.tsv` tiene una fila por (Gene, Position, side) con score > 0, y el `p.emp` es compartido entre lados (máximo): las filas top y bottom no son dos pruebas direccionales. `filtered_discovery.tsv` tiene una fila por (Gene, Position, side) evaluada por CT, con `asr_path_score` ≤ 0 cuando se evaluó sin detección positiva. Una posición sin score no se convierte en cero ni en negativo si no se demuestra que se evaluó. La posición es cero-based en el alineamiento recortado de CAAS, con una unidad por posición biológica, sin multiplicar por hipótesis, transcritos ni lados.

## 4. Capas

Se conservan: `pfam_domains`, `pfam_clans`, `ucr_core`, `ucr_flank`, `fubar_positive`, `fubar_purifying`, `fade_top_sig`, `fade_bottom_sig` y las listas de sitios (`custom`). VEP, COSMIC y PrimateAI se mantienen como caracterización cruda.

Se retiran como colecciones de test: ortogrupo, bin genómico, COSMIC/PAI3D por ortogrupo y capas de selección por ortogrupo. Un conjunto por ortogrupo es un gen, y ese nivel es Accumulation. Las ventanas de 1 Mb solo se usan como atributo de posición para la proximidad genómica.

Una colección tiene valor para el test cuando su vocabulario se repite entre genes (dominio, clan, región conservada, sitio). Así se puede agrupar entre genes y estratificar por gen.

## 5. Nivel 1: Accumulation

Se mueve a POSENRICH. Para cada dirección (all, top, bottom): sea M = |U|, K el número de posiciones positivas distintas de la dirección y n_g las posiciones de U del gen g. El recuento observado del gen sigue un hipergeométrico `Hyper(M, n_g, K)`, y el p exacto de cola superior se ajusta con BH sobre los genes probados (los que tienen al menos una positiva en la dirección).

Accumulation y el motor del nivel 2 leen la misma tabla `universe.tsv.gz` de `posenrich_universe.py`, una fila por posición de U con `nd`, `observed_positive`, `observed_top` y `observed_bottom` (0/1). K sale de la suma de la columna de la dirección (`observed_positive` para `all`, de modo que una posición detectada en ambos lados cuenta una vez) y n_g del recuento de filas por gen. Accumulation no construye universo propio ni lee el fondo, las posiciones probadas ni el nulo. Los genes de la lista sin posiciones en U (n_g = 0) no tienen fila en la salida y no se prueban. Las posiciones observadas son las filas con score positivo de `filtered_discovery.tsv`, de modo que el nivel 1 no depende de SCORING; con `position_scores.tsv` como observadas se obtiene la misma tabla (comprobado en el run `malignant_prevalence_complete`). Contra las salidas que daba el constructor propio de Accumulation (11.530 genes), el recuento, `PoolPositions` y el p coinciden bit a bit en los 10.619 genes con n_g > 0; los 911 restantes tenían n_g = 0 y ninguna positiva (`tests/test_posenrich_universe.py`). La salida conserva las columnas `Gene`, `ActualCount`, `PoolPositions`, `MeanNull`, `SD` y `PValue`.

Limitación: el hipergeométrico trata las posiciones de U como equivalentes. El nivel 1 no se calibra con el nulo en esta versión. Una variante calibrada contra los recuentos por (gen, ciclo) del nulo queda anotada como posible mejora posterior; no forma parte del diseño actual.

## 6. Nivel 2: concentración dentro del gen

**Estadístico (esquema).** Un solo estadístico graduado por score. Cada posición i de U aporta un peso w_i: su score global (`all`) si es positiva y 0 si no se detecta. Para un término S y cada gen g con posiciones dentro y fuera de S, sea W_g la suma de pesos del gen, n_g sus posiciones en U y m_g las que están en S. El residuo del gen es

    r_g = (suma de w_i sobre S ∩ g) − W_g · m_g / n_g

y el estadístico principal es la suma de r_g / s_g sobre genes, donde s_g es la desviación típica de r_g en los ciclos del nulo. La suma simple de r_g se calcula como sensibilidad. Se calcula igual sobre el observado y sobre cada ciclo del nulo, con los márgenes (W_g) de ese ciclo.

Comprobación previa a publicar p/q: con el control "ciclo nulo como observado" se compara la tasa de falsos positivos de las dos versiones. Si la normalizada no mantiene el nivel nominal (por ruido en genes con pocas posiciones), se revisa antes de ver resultados del observado.

**Exclusión de genes.** Un gen entra en el estadístico solo si tiene posiciones dentro y fuera del término y una s_g estimable (no nula, con peso distinto de cero en una fracción suficiente de ciclos; el umbral se fija antes de ver resultados del observado). Los demás se excluyen de forma explícita: cada término publica la lista de genes excluidos con su motivo (sin posiciones fuera del término, sin positivas, s_g nula o inestable). Nunca se imputa un residuo.

**Calibración.** p = (1 + nº de ciclos con estadístico ≥ observado) / (N + 1), con N = 10.000 y suelo 1e-4. Cada ciclo del nulo se normaliza con la s_g calculada sin ese ciclo (a partir de los mismos contadores de suma y suma de cuadrados), igual que el observado, que no entra en su s_g. Con la s_g de todos los ciclos, un ciclo extremo infla su propia dispersión y su estadístico queda subestimado frente al observado; el efecto es mayor en los pares con pocos ciclos no nulos. Se calcula por bloque (500 ciclos, 9,3 M de filas por bloque) y los recuentos de excedencias se suman. Universo: U. BH por familia (colección y ranking).

**Implementación.** `posenrich_within_gene.py`. Lee los bloques del nulo dos veces y solo cruzan bloques contadores aditivos, de modo que el resultado no depende de cómo se reparten los ciclos (comprobado con 1, 7 y 42 bloques). Entradas: `universe.tsv.gz` de `posenrich_universe.py`, la tabla observada con `CAAS_score`, el directorio del registro (`--registry-dir`) y los bloques `blk_*` (`positions.parquet` y `cycles.txt`). Salidas: `within_gene_results.tsv` (un término por fila, con p y q normalizados y simples, tamaños de efecto y recuento de genes elegibles y excluidos), `within_gene_by_gene.tsv` (un par término-gen por fila, con su residuo, s_g, ciclos no nulos y motivo de exclusión), `within_gene_manifest.json`. `--holdout K` reserva K ciclos del nulo como observado y escribe `holdout_fpr.tsv`. El umbral de estabilidad de s_g es `--min-nonzero-cycles` (30 por defecto): con menos valores no nulos la desviación de una variable con muchos ceros es inestable. Es una regla heurística, no demostrada.

**Validación medida.** Oráculo denso independiente (tres semillas), partición de bloques y ciclos reservados sobre datos sintéticos, en `tests/test_posenrich_within_gene.py`. Con un bloque real de 500 ciclos del run (U de 165.468 posiciones) y 60 términos aleatorios sin señal, 400 ciclos de calibración y 100 reservados: tasa de falsos positivos a α = 0,05 de 0,055 (normalizado) y 0,046 (simple), y a α = 0,01 de 0,010 y 0,007. Con 250 reservados y dos selecciones distintas: 0,043 y 0,027 (normalizado), 0,044 y 0,049 (simple) a α = 0,05. Los pares término-ciclo no son independientes (los términos comparten posiciones y cada ciclo desplaza todos sus términos), de modo que el error de estas tasas lo fija el número de ciclos reservados y no el de pares. Una validación concluyente necesita los 10.000 ciclos y varios conjuntos de términos, en el cluster por sbatch. Tiempo y memoria del bloque real con 474.450 pares: 24 s y 1,5 GB.

**Qué cubre el nulo sin parámetros nuevos:** la detectabilidad por posición, la co-ocurrencia de detecciones dentro de un gen y el filtro de genes del nulo.

**Informatividad.** Solo aportan los genes con posiciones dentro y fuera del término y con al menos una positiva. Medido: la mediana de un gen es 10 posiciones de U con 2 positivas; 5379 de 10.797 genes tienen al menos 2 positivas. La tabla por gen es descriptiva casi siempre, y el agregado entre genes es la unidad de inferencia. Los términos que cubren casi todo un gen (por ejemplo una lista de sitios definida por gen) son equivalentes al nivel 1 y se marcan como tales.

**Frecuencia e intensidad.** El estadístico graduado responde a la vez si el término concentra detecciones y si concentra detecciones fuertes: con pesos binarios es un recuento, y con el score es el estadístico graduado. Hay un solo p por término y la familia de BH no se duplica. El p no separa frecuencia de intensidad. Esa separación se hace con dos tamaños de efecto descriptivos, sin p:

- **Frecuencia:** proporción de positivas entre las posiciones de U del término, frente a la del resto de sus genes.
- **Intensidad:** diferencia de score medio entre las positivas del término y las positivas del resto de sus genes.

**Ranking.** El peso es el score crudo (`all`), que tiene equivalente en cada ciclo del nulo. `p.emp` no se usa en el estadístico (no se guarda por ciclo) y se mantiene en las vistas descriptivas.

## 7. Colecciones y registro

- El universo de test es U. La geometría de las anotaciones no se recorta a U: los sitios y regiones fuera de U siguen siendo blancos para proximidad.
- El registro guarda las pertenencias según lo que son, no como una fila por término y posición (`posenrich_members.py` las expande cuando un consumidor las pide):
  - **Atributos** (un valor por posición): columnas de `positions.tsv`: `ucr_core`, `ucr_flank`, `ucr_region`, `ucr_variability`, `fubar_call` (positive, purifying, neutral o vacío si no se probó), `fade_top_sig`, `fade_bottom_sig` y `window`. Las colecciones `ucr_core`, `ucr_flank`, `fubar_positive`, `fubar_purifying`, `fubar_tested`, `fade_top_sig`, `fade_bottom_sig` y `genomic_windows` se derivan de ellas.
  - **Regiones** (Pfam: dominios y clanes): `regions.tsv`, un registro por tramo de columnas CAAS consecutivas de un término en un gen (inicio y fin inclusivos). Los huecos de la referencia parten el tramo, de modo que la expansión da exactamente las posiciones que daba la construcción por posición.
  - **Sitios** (listas personalizadas): `sites.tsv`, una fila por término, gen y posición.
- El registro no guarda la tabla expandida de pertenencias (una fila por término y posición). Las colecciones de cobertura (`fubar_tested` y `genomic_windows`) no cuentan como anotación en `n_positions_annotated`.
- Tamaño de ventana: `--window_size_bp`, con el valor de `scoring_window_size_bp`; la colección se llama `genomic_windows`, sin el tamaño. El nombre de cada ventana cuenta su inicio en Mb, kb o bp (1 Mb: `chr1_12M`; 100 kb: `chr1_12300k`).
- Motivo de la forma compacta: la tabla expandida del run ocuparía 3,7 GB, dominada por Pfam (6,6 M de filas), `fubar_tested` (6,4 M), ventanas (6,4 M), UCR core (5,5 M), FUBAR purificador (3,0 M) y clanes (2,5 M). Los sitios son unas pocas miles de filas.
- Lista de genes de U: el fondo limpio de CT_POSTPROC (`postproc/cleaned_backgrounds/cleaned_background_main.txt`, 11.530 genes en el run). Es el fondo de CAAStools (`caastools/background_genes.output`, 11.710 genes) sin los genes que el filtrado de genes retira del observado (`postproc/gene_filtering/removed_genes_summary.tsv`, categorías Extreme y Dubious; 180 en el run), que por construcción no pueden tener positivas observadas. Con esta lista |U| = 165.468 en el run; con la del fondo de CAAStools, 178.365 (12.897 posiciones más del nulo y ninguna positiva observada). En una ejecución integrada, el canal `pp_cleaned_bg` de `main.nf` alimenta a todos los consumidores. La ejecución de SCORING del run que escribió `null_universe_summary.tsv` usó la lista de 11.710 genes; qué parámetro o ruta la trajo no se ha verificado. `posenrich_universe.py` recibe la lista como entrada explícita y escribe su huella y |U| en el manifiesto. SCORING registra la huella de su lista en `null_universe_summary.tsv` y el registro compara ambas huellas, con aviso si difieren.

## 8. Componentes descriptivas

- **Registro:** sección 7, con el constructor de U compartido con Accumulation y SCORING. Las fuentes se leen con `position_sources.py`.
- **Caracterización:** lee los atributos desde las columnas del registro y no hace inferencia.
- **Proximidad:** proteica (residuos de referencia) y genómica (nucleótidos), separadas. El control de distancias no está validado.
- **Informe 14:** consume las tablas de todas las componentes (sección 13).

## 9. Accumulation y SCORING

Accumulation es el nivel 1 de POSENRICH. SCORING no calcula ni lee Accumulation: no tiene columnas `accum_*` ni banderas `flag_accum*`, y no existe `fcs_stats_accum`. ENRICHMENT y los informes 11, 12 y 13 no lo consumen. `genomic_window_accum_gene_caas.tsv` es una tabla del informe 11 para un gráfico por ventanas, no un test.

## 10. Puntos abiertos

1. **Validación del nivel 2 a escala.** Ciclo nulo como observado con los 10.000 ciclos y varios conjuntos de términos: tasa de falsos positivos por término y colección, y comparación de la versión normalizada con la suma simple, antes de publicar p y q como resultado. Se lanza por sbatch.
2. **Umbral de estimabilidad de s_g.** El valor por defecto (`--min-nonzero-cycles` 30) es una regla heurística no demostrada. Debe fijarse antes de ver resultados del observado.
3. **Lado.** El nulo guarda solo `all`. Una bandera de lado de un byte permitiría calibrar por dirección; no se hace, y el lado se describe en el observado.
4. **Nivel 1 calibrado.** Anotado como mejora posible (sección 5); no forma parte del diseño.
5. **Lista de genes de SCORING.** La ejecución de SCORING del run que escribió `null_universe_summary.tsv` usó la lista de 11.710 genes en vez de la del fondo limpio (11.530). No se ha verificado qué parámetro o ruta la trajo. El registro compara las huellas de ambas listas y avisa si difieren.
6. **Control de proximidad.** Pendiente de diseñar (sección 11).

## 11. PP2A

Ensayo propio, fuera de POSENRICH y de la GUI de PhyloPhere (sección PP2A de la pestaña Enrichment y `posenrich_custom_file` de las plantillas). Pregunta principal: proximidad de las detecciones a los sitios de su proteína. Pregunta secundaria: solapamiento exacto.

Medido con los sitios mapeados (`20261024/ppga/pp2a`): de 5284 posiciones únicas de sitios, 3815 están en genes del universo, 3772 son posiciones probadas y 76 caen en U (2,0 %; 22 positivas). El fondo medido de Sapkota tiene 22.383 posiciones, 254 en U. El solapamiento exacto tiene poca potencia por esto. La proximidad usa los 3772 sitios como blancos y no está limitada por U. Los sitios mapeados quedan simétricamente por debajo en U respecto a sus vecinos a ±1 y ±2 residuos (76 frente a 90-106), compatible con que las columnas S/T/Y se detecten menos. No se ha comprobado con el residuo del alineamiento.

La proyección a la referencia productiva (accession, isoforma, péptido único, MAP) la hace `bin/map_pp2a_to_caas.py`; los registros que no validan quedan como pendientes, nunca como negativos.

Pendiente de diseñar: el control de proximidad. Debe contrastar si las detecciones están más cerca de los sitios de su misma proteína de lo esperado bajo un control que respete proteína, elegibilidad, cobertura medida y densidad de sitios, con radios predefinidos y una sensibilidad que excluya las coincidencias exactas. Incluye también el fondo por estudio: el fondo medido de Sapkota no es el de Kauko ni el de Stewart, ni equivale a todas las S/T/Y. La distancia proteica en residuos no se mezcla con la genómica ni se presenta como cercanía 3D.

## 12. Fuera de alcance

Lachenbruch y el hipergeométrico sin condicionar en el gen como tests; conjuntos por ortogrupo; calibración por lado; validación inferencial a escala; FCS de posiciones.

## 13. Módulo POSENRICH y estructura del código

POSENRICH es un módulo propio, separado de ENRICHMENT: `workflows/posenrich.nf`, `subworkflows/POSENRICH/` y `conf/posenrich.config`. Contiene el universo, el registro, la caracterización, el nivel 1 (Accumulation), el nivel 2 (`posenrich_within_gene.py`), la proximidad, la generación de las tablas Pfam y UCR y el informe 14. Su interruptor es `--posenrich`.

ENRICHMENT conserva FCS, AMI/DOMINO y los informes 12, 13 y 15. El orden de datos es SCORING, POSENRICH, ENRICHMENT. El informe 15 recibe la caracterización de POSENRICH como entrada opcional; sin ella omite la sección de posiciones interesantes.

**Informe 14** (`14.POSENRICH_report.Rmd`). Secciones: universo y cobertura, caracterización cruda de las positivas por capa, nivel 1 por dirección, nivel 2 (términos con p, q, tamaños de efecto y genes excluidos con su motivo), proximidad, métodos y límites. La GUI tiene una sola especificación de informe regenerable para POSENRICH.

**Código** (`subworkflows/POSENRICH/local/src/`):

| Script | Función |
|---|---|
| `position_sources.py` | Módulo de lectura de fuentes (Pfam, UCR, FUBAR, FADE, ventanas, listas personalizadas) con proyección a columnas CAAS a través del MAP. No es un programa ni escribe archivos |
| `posenrich_registry.py` | Registro: posiciones con atributos como columnas, regiones, sitios, colecciones, cobertura y manifiesto |
| `posenrich_members.py` | Expande atributos, regiones y sitios a pertenencias cuando un consumidor las pide |
| `posenrich_universe.py` | Universo U (`universe.tsv.gz`) y su manifiesto |
| `accumulation_test.py` | Nivel 1 |
| `posenrich_within_gene.py` | Nivel 2 |
| `posenrich_characterize.py` | Caracterización descriptiva |
| `posenrich_proximity.py` | Proximidad proteica y genómica |

Las colecciones por ortogrupo, COSMIC y PAI3D por ortogrupo y las exportaciones a GMT no forman parte del módulo (sección 4). Las pruebas están en `tests/test_posenrich_*` y `tests/test_posenrich_report.R`; el smoke de Nextflow es `tests/test_posenrich_smoke.py`.

## 14. Validación a escala

`adhoc_scripts/posenrich_scale_check.sbatch` ejecuta los scripts del módulo directamente sobre los ficheros de un run terminado (sin Nextflow) y escribe todo bajo `.tmp/posenrich_check/job_<id>/`; el directorio del run solo se lee. Mide tiempo y memoria máxima de cada paso, tamaño de las tablas del registro, la tasa de falsos positivos del nivel 2 con los 10.000 ciclos (`--holdout`), la concordancia del nivel 1 con la Accumulation publicada por el run, y renderiza el informe 14.

**Cómo leer el resultado.**

- **Universo.** |U| y positivas observadas coinciden con las del run (en `malignant_prevalence_complete`: 165.468 y 27.080 con la lista de 11.530 genes). Si la lista de genes del registro y la de `null_universe_summary.tsv` difieren, `gene_list_check` del registro da `mismatch`: es la detección funcionando, no un fallo, y apunta a que SCORING usó otra lista.
- **Nivel 1.** `PoolPositions` y `ActualCount` iguales y diferencia de p nula en todos los genes presentes en ambas tablas.
- **Registro.** Compara el tamaño de `positions.tsv`, `regions.tsv` y `sites.tsv` con el de `positions.tsv` y `members.tsv` del registro anterior del run (611 MB y 3,7 GB). Si la memoria máxima del paso `registry` o `characterization` se acerca al límite del trabajo, esos procesos necesitan optimización antes de usarlos en otros runs.
- **Nivel 2, falsos positivos.** Los ciclos reservados son nulos por construcción, así que la tasa por colección y estadístico a α = 0,05 y 0,01 debe ser compatible con α. Los términos comparten posiciones y cada ciclo desplaza todos sus términos, de modo que el error de la tasa lo marca el número de ciclos reservados (1.000: desviación típica de unos 0,007 a α = 0,05 y 0,003 a α = 0,01, que subestima el error real). Como criterio heurístico, no demostrado: una tasa entre α − 0,015 y α + 0,015 (a α = 0,05) y entre 0 y 0,015 (a α = 0,01) es compatible con el nivel nominal; una tasa claramente superior (por encima de 0,07 a α = 0,05 o de 0,02 a α = 0,01) indica un estadístico anticonservador y sus p y q no deben interpretarse; una tasa inferior es conservadora y aceptable. Si el estadístico normalizado falla y el simple no, revisar el umbral `min_nonzero_cycles` y la normalización antes de publicar.
- **Nivel 2, observado.** Número de términos probados y no identificables, genes elegibles y excluidos por término, y el suelo de p (1/(N+1) con N = ciclos de calibración). Un término con q < 0,1 solo es una señal de priorización si la validación anterior es satisfactoria.
- **Informe 14.** Renderiza sin errores, con el panel de falsos positivos y las seis secciones.

