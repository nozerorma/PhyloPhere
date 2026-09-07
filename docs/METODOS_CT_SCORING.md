# PhyloPhere — Selección de contrastes → FOP → discovery → null → desambiguación → postprocesado → scoring

Documento de métodos, cotejado línea a línea con el código de la rama `nongreedy_dunn`
(septiembre 2026). Coordenadas de alineamiento **0-based** salvo donde se indique.
`K` = número de pares de contraste independientes. `mdf` = `min_divergent_fraction`.

> **Nota sobre el test hipergeométrico.** El p-valor hipergeométrico **por sitio** que
> tenía CAAStools clásico está **completamente eliminado** de la cadena CT. `discovery.tab`
> ya no emite columna `pvalue` (cabeceras reales en [`disco.py:211-249`](../subworkflows/CT/local/modules/disco.py)),
> y `7.CT_signification.Rmd` **solo** integra el bootstrap. Las referencias a `pvalue` que
> aún aparecen en el código de desambiguación son ramas defensivas muertas. El único
> hipergeométrico que sobrevive en todo el pipeline está a nivel de *conjunto de genes*, no
> de sitio: la **Parte 1 del test de Lachenbruch** en FCS (Fisher exacto = cola superior
> hipergeométrica sobre la prevalencia de un pathway, [`fcs_enrich.R:199-270`](../subworkflows/ENRICHMENT/local/src/fcs_enrich.R))
> y **DOMINO** (módulos PPI, `scipy.stats.hypergeom`).

---

## 0. Ejemplo de trabajo (hilo conductor)

Todo el documento sigue **una sola posición**:

- **Gen** `OPN1` (una opsina hipotética), con 12 posiciones CAAS candidatas en total.
- **Columna** `210` del alineamiento.
- **Rasgo** `depth_m` = profundidad de hábitat (continuo, metros).
- **`K = 3`** pares de contraste independientes tras la selección.
- El árbol tiene, entre otras, estas especies relevantes:
  `sp_A_deep, sp_A_shallow, sp_A2_deep` (clado A);
  `sp_B_deep, sp_B_shallow` (clado B);
  `sp_C_deep, sp_C_shallow, sp_C2_shallow` (clado C).

En la columna 210, los linajes profundos portan **Trp (W)** y los someros **Tyr (Y)**,
salvo `sp_C2_shallow` que porta **Phe (F)**. Es una sustitución aromática recurrente:
convergencia química, no siempre idéntica a nivel de residuo.

---

## A. Selección de contrastes  (`selection_algorithm.R`, `3.CI-composition.Rmd`)

La entrada es **una columna** de `my_traits.tsv`, nunca una lista de grupos. Foreground y
background emergen de los pares seleccionados.

### A.1. PSS — Phylogenetic Shift Score (núcleo BM/OU compartido)

`pss_core.R` (motor phyloq vendorizado). Se ejecuta **una vez** sobre el rasgo observado:

1. `fit_models(tree, depth_m)` ajusta Brownian Motion (BM) y Ornstein-Uhlenbeck (OU).
2. `select_model` elige por **AIC plano**: OU sii `AIC_OU + 2 < AIC_BM`, si no BM.
   *(Ejemplo: `AIC_BM = 141.2`, `AIC_OU = 137.9` → `Δ = 3.3 > 2` → **OU**.)*
3. `covariances_from_fits` produce `cov_bm` / `cov_ou`; se mantienen **fijas** en todas las
   réplicas de permulación (E), de modo que el null aplica exactamente el modelo que usó la
   selección observada.
4. `calculate_pairwise_scores` da, por par no ordenado de especies: `PatristicDistance`,
   `abs_diff` (diferencia de rasgo) y `FinalScore` = **PSS** — cuánto diverge ese par
   respecto de lo que el modelo evolutivo predice a esa distancia filogenética.

El PSS **siempre ordena** los candidatos. Solo en el caso continuo actúa además como *gate*.

### A.2. Gate de elegibilidad por tipo de rasgo (`3.CI-composition.Rmd:340-385`)

| tipo | gate de elegibilidad (test de no-solapamiento) | rol del PSS |
|------|--------------------------------------------------|-------------|
| **ordinal / binaria (categórica)** | una especie en el **nivel codificado máximo**, la otra en el **mínimo**. Todo par así es elegible. **Sin corte por cuantil de valor** (eliminado en todo el pipeline). | solo ordena |
| **count** (`n_trait`/`c_trait`) | intervalos de credibilidad **Jeffreys 95 %** (Beta(0.5,0.5)) que no se solapan: `lb[hi] > ub[lo]` | solo ordena; desempate por `pair_n` |
| **continua** *(← el ejemplo)* | por defecto: top **`pss_top_pct`** (0.05) de los pares hi>lo por PSS. Es un percentil de rango que escala con el tamaño del árbol, no con la señal. `pss_top_pct` fuera de `(0,1)` **desactiva** el gate → solo ordena, como los demás. | gate **y** orden |

*Ejemplo:* `depth_m` es continuo. Con árbol de ~40 tips y `pss_top_pct = 0.05`, de los
~180 pares hi>lo pasan el gate los ~9 de mayor PSS. Entre ellos: `sp_A_deep~sp_A_shallow`
(PSS 8.1), `sp_B_deep~sp_B_shallow` (PSS 6.4), `sp_C_deep~sp_C_shallow` (PSS 5.9),
`sp_A2_deep~sp_A_shallow` (PSS 7.2), `sp_C_deep~sp_C2_shallow` (PSS 4.8), …

> **Categórica vs continua.** Una codificación binaria H/L fijaría `abs_diff = 1` en todos
> los pares, así que el PSS colapsaría a *aislamiento filogenético puro*: enorme para un par
> hermano estrecho, ínfimo para el resto. La codificación continua preserva la brecha real
> de profundidad, lo que sube a las hipótesis conservadas contrastes que la binaria dejaría
> al fondo del harvest.

`candidate_species.tab` se escribe **antes** del filtro Dunn (unión de miembros de todo par
elegible); es lo que consume FADE, no el conjunto CAAS post-Dunn.

### A.3. Ranking unificado (`rank_candidates`, `lean_contrast_selector.R:78`)

Clave primaria: **PSS descendente** cuando hay PSS finito; si no, distancia patrística
ascendente. Desempates: `|abs_diff|` desc, luego `pair_n` desc (count), luego orden de
entrada (estable). **Compartido bit a bit** entre selección observada y null.

### A.4. Ensamblado greedy con gate de Dunn (`greedy_dunn_select`)

Dunn modificado por cluster (`mod_dunn_lean`): *(distancia mínima a cualquier otro cluster)
/ (diámetro propio del cluster)*. Overall Dunn = mínimo sobre clusters.

- Semilla: el candidato mejor rankeado (cluster 1, Dunn = ∞).
- Iterativo: entre candidatos con ambas especies libres, se calcula el Dunn resultante de
  añadirlos; con `enforce_dunn = TRUE` solo se aceptan los que mantienen **todos** los
  clusters con `mod_dunn ≥ 1`; entre esos, el de Dunn más alto (desempate por rango).
- Se para cuando ningún candidato cualifica o el overall Dunn bajaría de 1.

*Ejemplo:* **H1** = `{ sp_A_deep|sp_A_shallow , sp_B_deep|sp_B_shallow ,
sp_C_deep|sp_C_shallow }`, K = 3 pares mutuamente independientes. `sp_A2_deep~sp_A_shallow`
NO entra en H1: comparte `sp_A_shallow` con el par ya seleccionado.

`CHECK_MIN_CONTRASTS`: si K < 3, **toda la cadena CT se salta en silencio** (exit 0).

---

## B. Expansión FOP — hipótesis paralelas múltiples

`fop_pair_sel.f` (observado), espejeado en el null por `lean_fop_harvest`. Dentro del
"dominio" de cada par canónico suele haber varios pares elegibles alternativos; cada
combinación válida es una **hipótesis** `H2..Hn` que se lleva por toda la cadena.

1. **H1 = baseline canónico** (A.4).
2. **Partición de Voronoi.** Cada especie del árbol se asigna al par canónico `k` cuya
   especie más cercana (patrística) minimiza la distancia → **un dominio por par canónico**.
   *Ejemplo:* dominio 1 = región del clado A (incluye `sp_A2_deep`); dominio 2 = clado B;
   dominio 3 = clado C (incluye `sp_C2_shallow`).
3. **Pools de dominio.** `alt_pools[[k]]` = pares candidatos elegibles (lista rankeada de
   A.3) con **ambas** especies dentro del dominio `k`.
   *Ejemplo:* dominio 1 → `{ A_deep|A_shallow , A2_deep|A_shallow }`; dominio 2 →
   `{ B_deep|B_shallow }`; dominio 3 → `{ C_deep|C_shallow , C_deep|C2_shallow }`.
4. **Harvest combinatorio.** Un par por dominio: exhaustivo con `expand.grid` si
   `total_combos ≤ ITER_CAP = max_fop·20`, si no `ITER_CAP` extracciones aleatorias. Se
   descartan combinaciones con especies repetidas o firmas ya vistas; se conservan solo las
   de `overall_dunn ≥ 1`. Si algún dominio tiene pool vacío, el harvest es insatisfacible
   (solo H1).
5. **Ranking y truncado de calidad.** Hipótesis válidas ordenadas por
   **(min-PSS entre los K pares, luego mean-PSS, luego Dunn)** — una hipótesis vale lo que
   su contraste más débil — y se conservan las **top `max_fop − 1`**.

*Ejemplo:* 2 × 1 × 2 = 4 combinaciones, todas Dunn-válidas:

| hyp | dominio 1 | dominio 2 | dominio 3 | min-PSS |
|-----|-----------|-----------|-----------|---------|
| **H1** | A_deep\|A_shallow | B_deep\|B_shallow | C_deep\|C_shallow | 5.9 |
| **H2** | A2_deep\|A_shallow | B_deep\|B_shallow | C_deep\|C_shallow | 5.9 |
| **H3** | A_deep\|A_shallow | B_deep\|B_shallow | C_deep\|C2_shallow | 4.8 |
| **H4** | A2_deep\|A_shallow | B_deep\|B_shallow | C_deep\|C2_shallow | 4.8 |

Con `max_fop = 100` se conservan las 4.

**Salidas:**

- `contrast_hypotheses_pairs.tsv`: `hypothesis_id`, `pair` (= id de **dominio de Voronoi**),
  `species1`, `species2`, `pss_score`. Es el fichero de pesos PSS que consume el
  domain-pooling del scoring (H.2).
- `traitfile_H1.tab … traitfile_H4.tab`: cada una un fichero de contraste de 3 columnas
  (`species`, `1|0`, `pair`). **Todo paso CT posterior corre una vez por hipótesis**; la
  columna `trait` del `discovery.tab` lleva `traitfile_H<n>.tab`.

Restricciones dobles del FOP: **localidad de dominio** (un par alternativo cabe en un solo
dominio canónico) y **corte `max_fop`**.

---

## C. Discovery CAAS/CAAP con hipótesis múltiples  (`caas_id.py`, `caap_id.py`, `disco.py`)

Una fila de `caastools/discovery.tab` por **(gene, position, esquema, hipótesis)**.
Cabeceras reales: `gene, mode, caap_group, trait, position, caas, [amino_encoded],
pattern, ffgn, fbgn, gfg, gbg, mfg, mbg, ffg, fbg, ms [, is_conserved_meta,
conserved_pair]` — **sin `pvalue`**.

### C.1. Filtro de columna (pre-discovery, `alimport.py`)

`filter_position` descarta una columna **solo** si su fracción de gaps supera el umbral. La
antigua regla de diversidad de aminoácidos (`seconds < min(#fg,#bg)`) fue **eliminada**:
escalaba el número de residuos minoritarios exigidos con el número de contrastes y borraba
convergencia genuina de bajo origen.

### C.2. Regla CAAS (`iscaas`, `caas_id.py:236`)

Con `overlap` = min(residuos fg también en bg, residuos bg también en fg), y
`non_overlapping_fg/_bg` = residuos exclusivos de cada lado:

- **estricto** (`max_conserved = 0`): `overlap == 0` **y** (`non_overlapping_fg ≥ 2` o
  `non_overlapping_bg ≥ 2`).
- **con tolerancia** (`max_conserved = floor(K·(1 − mdf))`; p. ej. `floor(3·0.5) = 1`):
  `overlap ≤ max_conserved` **y** (`non_overlapping_* ≥ 2`). Los pares que comparten el
  residuo ancestral se registran en `conserved_pair` (`"overlap:par1,par2"`).

### C.3. Cinco esquemas de codificación

`US` (residuo literal) y `GS1..GS4` (agrupaciones bioquímicas progresivamente más gruesas).
**Cuántos esquemas disparan es una propiedad determinista de qué aminoácidos intervienen**
(distancia bioquímica discretizada), no una medida de fuerza de evidencia. Por eso el
scoring (H.4) agrega esquemas con una **media**, sin pesos.

*Ejemplo — columna 210, patrón fg/bg por hipótesis:*

| hyp | pares cambiados (fg/bg) | patrón US | ¿CAAS US? | ¿CAAS GS3 (aromáticos juntos)? |
|-----|------------------------|-----------|-----------|-------------------------------|
| H1 | W,W,W / Y,Y,Y | `WWW/YYY` | sí | sí (`aaa/bbb`) |
| H2 | W,W,W / Y,Y,Y | `WWW/YYY` | sí | sí |
| H3 | W,W,**F** / Y,Y,Y | `WWF/YYY` | sí (fg = {W,F}, bg = {Y}, overlap 0, non_overlap_fg ≥ 2) | sí |
| H4 | W,W,**F** / Y,Y,Y | `WWF/YYY` | sí | sí |

Así, para la columna 210 el discovery emite filas para `caap_group ∈ {US, GS4, GS3, GS2,
GS1}` × `trait ∈ {traitfile_H1..H4}` (las que superen la regla en cada esquema).

---

## D. Significación  (`7.CT_signification.Rmd`) — **solo bootstrap**

`signification/meta_caas/<scheme>_meta_caas.tsv`. La única cantidad de significación por
sitio es:

**`recovery_boot` = `occurrences / total`** — fracción de **labelings fenotípicos
permulados** (E.1) en los que el CAAS se vuelve a llamar en esa (posición, esquema). Se
comporta como una p empírica: *baja* = raro bajo re-etiquetado del fenotipo = **distintivo
del foreground real**. Es la entrada del `phen_score` en el scoring (H.3).

`gate_all` / `gate_sig` (mencionados en comentarios obsoletos de `scoring_compute.R`) ya
**no existen**: no hay ningún gate de significación por sitio en el scoring; la
priorización por sitio la dan `phen_score`, `pos_perm_p` y (a nivel gen) `gene_caas_pperm`.

---

## E. Bootstrap con FOP + null de permulación

**Dos cosas distintas, ambas construidas sobre el mismo pool de permulaciones de fenotipo**
generado por `permulations.R` (`RESAMPLE` en `ct_resample.nf` — no hay remuestreo de
columnas ni de especies al azar; el modo `random` de `init_bootstrap.py` es legado y no
está cableado).

### E.1. Bootstrap sobre labelings permulados → `recovery_boot`  (`boot.py`, `boot_vec.py`)

El kernel vectorizado (`VectorizedBootstrap`, BLAS, bit-idéntico al bucle escalar
`caasboot`) fija la columna del alineamiento y **reevalúa la regla CAAS bajo cada labeling
permulado** (cada "ciclo" es un vector fg/bg del pool). `F` (ciclos × especies) es la
máscara de foreground por ciclo. `recovery_boot` = (nº de ciclos que llaman CAAS en esa
posición/esquema) / (nº de ciclos).

**Con FOP** (`params.caas_perms_fop && params.multi_hypothesis`, `ct_bootstrap.nf:54`): el
pool de remuestreo es `fop_labelings.tab` (etiquetas `<base>~H<m>`) y los aciertos se
**colapsan a unidades de ciclo base** (`collapse_fop_hits_by_base`): un ciclo base cuenta
sii **cualquiera** de sus `~H<m>` (hipótesis alternativas de ese ciclo nulo) llama un CAAS
ahí. Es un OR a nivel discovery — sin ASR, sin pooling.

*Ejemplo:* con 40 ciclos base, la columna 210 de `OPN1` recibe un CAAS `W/Y` en 3 de ellos
(alguna hipótesis `~H<m>` lo llama) → **`recovery_boot(OPN1,210,US) = 3/40 = 0.075`**.

### E.2. Permulación de fenotipo → *foreground-specificity null*  (`permulations.R`)

Responde: *"¿el rasgo observado apunta a estos linajes más de lo que lo haría un rasgo
aleatorio filogenéticamente emparejado?"*

1. Ajuste único BM/OU por AIC; covarianzas fijas.
2. **Simulación + rank-match** (`simpermvec`): `sim.char` bajo el modelo (árbol
   OU-reescalado si OU) da el **orden**; los **valores** observados se reasignan según ese
   orden → distribución marginal permulada idéntica a la real por construcción.
3. **`evaluate_lean_contrast_selection`**: mismo ordenador + greedy Dunn que el observado,
   pero corre hasta **exactamente `target_pairs = K`** y gradúa la independencia:
   **Tier 1** (K pares con `mod_dunn ≥ 1`), **Tier 2** (exactamente uno por debajo),
   **Tier 0** = rechazo (≥ 2 por debajo, o no se forman K pares).
4. **Cosecha del pool** con escalada de presupuesto *target-aware* (extiende si la
   aceptación Tier-1 es sana pero el pool no se llena; aborta si es < 5 % tras 3 escaladas).
   `resample_000.tab` = `b_0` (etiqueta real); `resample_*.tab` = el pool.
5. **FOP mirror.** Para cada ciclo aceptado, `lean_fop_harvest` corre el **mismo** harvest
   de hipótesis alternativas, tomando `canon_pairs` = el contraste ya aceptado por el null
   (H1 verbatim, Voronoi sembrado desde él). Salidas: `fop_labelings.tab`
   (`<cycle>~H<m>` \t fg \t bg) y `fop_pairs.tsv` (`cycle, hypothesis_id, pair`=dominio,
   `species1/2, pss_score`).
6. **Replay ASR** (`disambiguation_perms_main.py` → `process_all_genes_perms`). La ASR es
   invariante al fenotipo (`load_precomputed_asr` = función pura de alineamiento+árbol), así
   que se **carga una vez por gen** y se re-scorea cada labeling vía
   `compute_asr_path_score` **verbatim** → el null queda calibrado por construcción. El eje
   fenotípico del null es un **rango por ciclo**: `1 − percent_rank(n_detected)` dentro del
   pool que ese ciclo descubrió (mismo `1 − percent_rank` que el observado). Con
   `--postproc-filter` se aplican al null los mismos filtros cluster+gen de CT_POSTPROC.
7. **Salidas del null:**
   - `perm_pos_pval.tsv`: **`pos_perm_p`** por `(Gene, Position, caap_group)` — null
     calibrado a nivel de posición. *Ejemplo:* `pos_perm_p(OPN1, 210, US) = 0.006`.
   - `caas_perms.rds` → **`caas_corStat_byrank`**: matriz densa **genes × N ciclos** en
     tres rankings (`global`, `top`, `bottom`); el null de *exceso* genome-wide que usan
     el `gene_caas_pperm` (H.6) y los tres tests de FCS (I.2). Sello `gene_stat =
     "size_adj_max"` (si no coincide con el observado, el consumidor deja `p.perm = NA` y
     pide reconstruir).

**Diferencia clave `recovery_boot` vs `pos_perm_p`:** mismo pool de permulación; el
primero es un OR a nivel discovery (barato, sin ASR), el segundo es el replay ASR completo
y domain-pooled por hipótesis.

---

## F. Desambiguación — integración de las hipótesis múltiples  (`disambiguate_single.py`, `path_scores.py`)

Una fila de salida por **(gene, position, esquema, hipótesis)**.

### F.1. ASR

`codeml` (modelo empírico `LG`), **una vez por gen**. Por par de contraste se reconstruye
el MRCA del foreground y el camino de él a la raíz, con posteriors por nodo.
`convergence_mode = "mrca"`.

### F.2. Cada fila de metadata se desambigua SOLO contra los pares de SU hipótesis  (`disambiguate_single.py:630-730`)

`parse_trait_pairs` devuelve `{contraste → pares}`: un contraste por hipótesis FOP
(`traitfile_H<n>.tab → n`). `_resolve_contrast(entry)` lee el campo `trait` de la fila
CAAS; si contiene `H<n>`, la fila **pertenece a la hipótesis n** y se desambigua
**exclusivamente contra los pares de esa hipótesis**.

> Esto es deliberado. La unión de los pares de todas las hipótesis (`_flattened_fallback`,
> red de seguridad con warning) **contaminaba** todos los ejes multi-par del ASR path score
> (puntos de fusión LCA, `independence`, `mrca_diversity`) con contrastes de hipótesis no
> relacionadas y destruía la independencia Dunn por-hipótesis que el harvest FOP impone.

*Ejemplo:* la fila `(OPN1, 210, US, trait=traitfile_H3.tab)` se desambigua contra los pares
de H3 = `{A_deep|A_shallow, B_deep|B_shallow, C_deep|C2_shallow}`, **no** contra los de H1.

Invariantes hoisted por gen (una sola vez, no por posición ni por ciclo): lookup del
alineamiento, `node_index`, MRCA de cada par (memoizado), y `gene_walk_cache` que memoiza
los recorridos MRCA→raíz invariantes al labeling.

### F.3. `compute_asr_path_score` — álgebra (por posición, esquema, hipótesis)

Recorrido acotado a dos regiones, nunca al camino MRCA→raíz completo:

- **segmento privado** del par: del nodo justo encima de su MRCA hasta (excluyéndolo) el
  **LCA** más cercano (donde su linaje se fusiona con el de otro par cambiado).
- **nodos LCA**: puntos de fusión de los MRCAs de los pares cambiados (≤ n−1), lo que hace
  que un *producto* sobre LCA sea independiente de la profundidad.

Cinco factores multiplicados:

| factor | qué mide | cómo |
|--------|----------|------|
| **`independence`** | ¿los ancestros compartidos ya llevaban el estado derivado? | `∏_LCA (1 − P(cualquiera del pool derivado en el LCA))`, cota *worst-case* (masa exacta si el residuo está registrado; si no, el remanente no registrado, sumado una vez) |
| **`core`** (replicación) | `P(≥2 cambios independientes)` | inclusión-exclusión sobre los scores de aislamiento privado por par, **por lado fenotípico**: `core_top`, `core_bottom`, y `core = 1 − (1−core_top)(1−core_bottom)`. Un par que cambia solo en top y otro solo en bottom → `core = 0` |
| **`mrca_diversity`** (paralelismo) | ¿los segmentos privados pasaron por fondos ancestrales distintos, o el estado del MRCA de A aparece en el recorrido de B? | ambas direcciones por par de pares, unión "encontrado en algún sitio"; `diversity = 1 − media_pairwise(prob. de fondo compartido)` |
| **`derived_agreement`** (convergencia) | dentro de cada lado con ≥ 2 pares cambiados, ¿qué fracción cae en el residuo derivado de pluralidad? | conteos de pares puros; media sobre lados cualificados; 1.0 si ningún lado tiene ≥ 2 |
| **`conservation_gate`** | ¿los pares conservados sostienen el contraste? | `0.5 + 0.5·media(conservación-a-raíz de pares conservados)`; **1.0 si no hay pares conservados** (posición novel). Solo confirma/socava, nunca aumenta |

```
replication    = independence · core
strength       = (0.75 + 0.25 · diversity) · derived_agreement
asr_path_score = replication · strength · conservation_gate            ∈ [0,1]
```

*Ejemplo — `(OPN1, 210, US)` por hipótesis:*

| hyp | pares cambiados (lado top) | independence | core_top | diversity | derived_agreement | gate | **asr_path_score** |
|-----|---------------------------|-------------|----------|-----------|-------------------|------|--------------------|
| H1 | W, W, W (3 orígenes limpios) | 0.92 | 0.88 | 0.86 | **1.00** | 1.0 | **≈ 0.72** |
| H2 | W, W, W | 0.92 | 0.88 | 0.84 | 1.00 | 1.0 | ≈ 0.71 |
| H3 | W, W, **F** | 0.90 | 0.85 | 0.82 | **0.67** (2 de 3 en W) | 1.0 | **≈ 0.47** |
| H4 | W, W, **F** | 0.90 | 0.84 | 0.80 | 0.67 | 1.0 | ≈ 0.46 |

Bajo **GS3** (W y F en el mismo grupo aromático) `derived_agreement = 1.0` también en
H3/H4, así que `asr_path_score ≈ 0.71` en las cuatro.

### F.4. Bloques columna emitidos (para el pooling FOP)  (`gene_wrapper.py::convert_convergence_result_to_dict`)

Aplanados por par `i` (= dominio de Voronoi por construcción del harvest):

- `mrca_<i>_node` / `_state` / `_posterior` — MRCA reconstruido.
- `mrca_<i>_path_score` — score del par (media top/bottom).
- `mrca_<i>_top_path_score` / `_bot_path_score` — scores **sin promediar** por lado
  (para el core direccional del pooling).
- `mrca_<i>_anc_aa` / `_top_aa` / `_bot_aa` — residuos crudos ancestral y derivado por lado.
- `conserved_<j>_node` / `_cons` — bloque paralelo de pares conservados, `j` por `pair_id`.
- escalares: `independence`, `mrca_diversity`, `derived_agreement`, `conservation_gate`,
  `core`, `convergence_type`, `change_top`, `change_bottom`, `change_side`.

Salida: `ct_disambiguation/caas_convergence_master.csv` (una fila por gene × position ×
esquema × hipótesis).

---

## G. CT_POSTPROC — con hipótesis múltiples

`prepare_postproc_input.py` → `filter_caas_clusters-param.py` → `filter_caas_genes.py`.

### G.1. Precluster (por fila, `prepare_postproc_input.py:37-56`)

Se elimina toda fila con **algún** `mrca_<i>_posterior < mrca_threshold`. Es un filtro
**por fila**, así que actúa por (esquema × hipótesis): si la reconstrucción de H3 tiene un
MRCA de baja confianza pero la de H1 no, se caen solo las filas de H3 de esa posición.

### G.2. Descriptores de residuo a nivel posición (`residue_descriptors.py`)

Calculados **aquí**, aguas arriba de `filtered_discovery.tsv`, agregando sobre **todas** las
filas de la `(Gene, Position)` (todo esquema, toda hipótesis):

- `derived_residues` = `"<top>/<bottom>"`, el lado sancionado por `change_side` muestra sus
  residuos **derivados**, el otro el **ancestral**.
- `top_residue_support` / `bottom_residue_support` = `"W:3,F:1"`: por residuo, el número de
  **pares de contraste CAAS distintos** (bloques `mrca_<i>`) que lo portan. **Acotado por
  el nº de pares del CAAS, NO inflado por el nº de hipótesis descubridoras** (un par cuenta
  una vez aunque aparezca en varias hipótesis).
- `top_residue_support_detail` = igual, pero contando **nodos ancestrales reconstruidos
  distintos** (un par resuelve a nodos distintos bajo hipótesis distintas → mezcla nº de
  pares con multiplicidad de hipótesis; vista secundaria, no un conteo de evidencia).
- `n_conserved_pairs` = nº de `conserved_<j>_node` distintos.
- (`add_species_tally`) `top_species_residues` / `n_top_species` = tally real de especies
  del contraste que portan cada residuo en la columna del alineamiento (chequeo "X de N
  especies", independiente de hipótesis).

*Ejemplo — `(OPN1, 210)`:* pares 1 y 2 portan W en todas las hipótesis; el par 3 porta W
(en H1/H2) o F (en H3/H4). `top_residue_support = "W:3,F:1"` (W: pares {1,2,3}; F: par {3}).
`change_side = "top"` (los cambios están en el linaje profundo). `derived_residues =
"WF/Y"`.

### G.3. Cluster filter (por posición, `filter_caas_clusters-param.py::ctrain`)

**Agnóstico de esquema y de hipótesis.** Toma el **conjunto de posiciones únicas** del gen
(`np.unique`), busca "trenes" (intervalos `[l,r]` con `span ≥ minlen` y
`densidad = count/span ≥ maxcaas`, por defecto `minlen 3`, `maxcaas 0.7`) y marca
`Discarded` **todas** las posiciones dentro de un tren. Una posición marcada elimina **todas
sus filas** (todo esquema, toda hipótesis).

### G.4. Gene filter

`filter_caas_genes.py` (`gene_filter_mode`: `dubious` / `extreme` / `both` / `none`) marca
genes con densidad de CAAS anómala respecto de su longitud. También por gen, no por fila.

**Salida:** `postproc/gene_filtering/filtered_discovery.tsv` — **único input del scoring**,
todavía con **una fila por (gene, position, esquema, hipótesis)** superviviente.

---

## H. Scoring — integración de señales por posición_gen  (`scoring_compute.R`, `fop_pool.R`)

### H.1. Alcance de esquemas y tag de hipótesis (§2a)

Cinco esquemas de scoring (`US, GS4, GS3, GS2, GS1`). `scheme_priority` (US=5…GS1=1) se usa
**solo** para elegir un esquema representativo en columnas de display (`change_side`,
`caap_group`, …), **nunca** en una cantidad scoreada. `hyp_id` se deriva de `trait`
(`H<n>` o `NA`).

### H.2. Domain-pooling FOP (§2b, `fop_pool.R::apply_fop_pooling`) — **cómo se integran las hipótesis**

**Objetivo:** colapsar las filas `H1..Hn` de un `(Gene, Position, caap_group)` a **una sola
fila** con un `asr_path_score` poolado, más los descriptores de posición.

**Por qué no una media simple sobre hipótesis:** `H1..Hn` son diseños de K pares
**solapados** sobre los mismos dominios de Voronoi, no réplicas independientes. Una media
uniforme (a) diluiría el canónico fuerte y (b) dejaría que una posición que cosechó muchas
hipótesis distorsione todo rango genome-wide en el que entre.

El pooling reconstruye **exactamente el álgebra de `path_scores.py`**, pero sobre valores
agregados por dominio en vez de por par:

**Paso 1 — agrupar.** `group_by(Gene, Position, caap_group)`. Si el grupo tiene ≤ 1
hipótesis distinta → *passthrough* (la fila pasa intacta). Con > 1:

**Paso 2 — dos trabajos de peso PSS distintos** (los pesos salen de
`contrast_hypotheses_pairs.tsv` vía `read_hypothesis_pairs`; equipeso si falta el fichero):

- **Job A — pool por dominio `c_i`** (`.pool_domain_col`). Para el dominio `i`: se juntan
  los `mrca_<i>_path_score` de todas las filas del grupo, se **deduplican por
  `mrca_<i>_node`** (= un par físico; un dominio puede aportar 2 pares distintos a lo largo
  del harvest), y se hace media ponderada por **el PSS propio de ese par en el dominio `i`**
  (`pair_pss(hyp, i)`, máximo sobre filas que comparten nodo). → un escalar `c_i` por
  dominio.
- **Job B — pool de ejes** (`independence`, `mrca_diversity`, `core` de fila): media
  ponderada por hipótesis, peso = **media del PSS de los K pares de esa hipótesis** (su
  credibilidad global, no su eslabón más débil).

**Paso 3 — core direccional** (`has_side_path`). En vez de poolear el
`mrca_<i>_path_score` (que ya promedia top/bottom y dejaría que un cambio top en un dominio
y un cambio bottom en otro cuenten como "2 dominios de acuerdo"), se poolean **por
separado** `c_top_i` y `c_bot_i` (de `mrca_<i>_top_path_score` / `_bot_path_score`, misma
dedup + pesos Job A):

```
core_top    = P(≥2 de {c_top_1 … c_top_K})       # inclusión-exclusión
core_bottom = P(≥2 de {c_bot_1 … c_bot_K})
core        = 1 − (1 − core_top)(1 − core_bottom)
```

**Paso 4 — `derived_agreement` harvest-wide y por esquema** (`rebuild_derived_agreement`).
Dos hipótesis pueden ser unánimes cada una pero aterrizar en residuos **distintos**;
poolear sus `da` nunca lo vería. Con el bloque `*_aa` presente, `da` se recalcula sobre el
conjunto de **pares cambiados distintos del grupo** (`.collect_changed_pairs`: una moda por
bloque `mrca_<i>`, indexado por par, no por nodo), por lado, con la **lógica exacta de
pluralidad de `path_scores.py`**, sobre residuos codificados bajo **el esquema de la fila**
(`encode_aa_r`).

**Paso 5 — `conservation_gate`** (`have_cons_cols`): se reconstruye desde los **pares
conservados distintos** del grupo (dedup por `conserved_<j>_node`):
`cg = 0.5 + 0.5·wmean(cons, w)`, con `w` = PSS del par en el dominio cuyo `mrca_<i>_node`
coincide. Sin pares conservados → `cg = 1.0`.

**Paso 6 — recombinar** (misma fórmula que F.3):

```
diversity_mult = 0.75 + 0.25 · diversity_pooled
replication    = independence_pooled · core_pooled
strength       = diversity_mult · da_pooled
asr_pooled     = clamp01( replication · strength · cg_pooled )
```

**Descriptores de posición añadidos aquí** (fuera del pooling por `caap_group`, unidos por
`(Gene, Position)`):

- `convergence_schemes` (`.position_descriptors`): test de **acuerdo entre linajes** sobre
  los residuos derivados reconstruidos. `""` = < 2 pares cambiados en un lado o desacuerdo
  genuino; `"US"` = residuo **idéntico** en todo par cambiado; `"GS1,GS2,…"` = ≥ 2 residuos
  distintos que aún comparten clase fisicoquímica bajo esos esquemas (el caso informativo —
  US ausente por construcción). NO es paralelo a `scheme_set` (que dice qué esquemas
  *descubrieron* la posición).
- `derived_residues`, `{top,bottom}_residue_support`, etc.: se **traen** de G.2 intactos.

**Input no-FOP** (contraste único): la única fila pasa sin cambios (`asr_path_score` intacto).

*Ejemplo — grupo `(OPN1, 210, US)`, filas H1–H4:*

| dominio `i` | pares distintos (dedup por nodo) | `c_i` (pooled, PSS-weighted) |
|-------------|--------------------------------|------------------------------|
| 1 | `A_deep\|A_shallow` (H1, H3) y `A2_deep\|A_shallow` (H2, H4) → **2 pares** | media ponderada de sus `path_score` ≈ 0.85 |
| 2 | `B_deep\|B_shallow` (todas) → **1 par** | 0.88 |
| 3 | `C_deep\|C_shallow` (H1,H2) y `C_deep\|C2_shallow` (H3,H4) → **2 pares** | ≈ 0.80 |

- `core_top = P(≥2 de {0.85, 0.88, 0.80}) ≈ 0.94`; `core_bottom = 0` (no hay cambios en el
  linaje somero) → `core = 0.94`.
- `independence_pooled ≈ 0.91`, `diversity_pooled ≈ 0.84` → `diversity_mult ≈ 0.96`.
- `da` harvest-wide **bajo US**: pares cambiados top = {1:W, 2:W, 3: moda(W,W,F,F)}. El par
  3 empata → la moda coge uno; con residuos {W,W,X} la concentración es 2/3 →
  **`da_US ≈ 0.67`**. **Bajo GS3**: W y F → aromático → {a,a,a} → **`da_GS3 = 1.0`**.
- `cg_pooled = 1.0` (sin pares conservados).
- **`asr_pooled(US)` ≈ 0.91 · 0.94 · 0.96 · 0.67 · 1.0 ≈ 0.55**.
  **`asr_pooled(GS3)` ≈ 0.91 · 0.94 · 0.96 · 1.0 · 1.0 ≈ 0.72**.
- `convergence_schemes = "GS4,GS3,GS2,GS1"` (2 residuos distintos W/F, todos aromáticos;
  US falla). `derived_residues = "WF/Y"`, `top_residue_support = "W:3,F:1"`.

Tras §2b, la columna 210 de `OPN1` tiene **una fila por esquema**: US, GS4, GS3, GS2, GS1.

### H.3. Score de dos ejes por fila (§2f) y null de posición (§2f-bis)

```
phen_score = 1 − percent_rank(recovery_boot)      # percentil genome-wide, uniforme en [0,1]
asr_score  = asr_path_score  (= asr_pooled)        # limpieza en el árbol
caas_row   = phen_score · asr_score               # dos ejes ortogonales
```

`percent_rank` es genome-wide, por eso §2b tuvo que colapsar las hipótesis primero (si no,
una posición con 4 hipótesis contaría 4 veces en el rango). El `recovery_boot` de la fila
es el de la hipótesis representativa (mayor asr) tras el colapso.

`pos_perm_p` (§2f-bis) se une por `(Gene, Position, caap_group)` **antes** del colapso de
esquemas, para que ride en la misma `mean()` de §2g que los demás ejes.

*Ejemplo:* `recovery_boot(OPN1,210) = 0.075`; su percentil genome-wide es bajo (raro bajo
permulación) → `percent_rank ≈ 0.08` → **`phen_score ≈ 0.92`**.

| esquema | `asr_score` | `phen_score` | `caas_row` |
|---------|-------------|--------------|------------|
| US  | 0.55 | 0.92 | **0.51** |
| GS4 | 0.63 | 0.92 | 0.58 |
| GS3 | 0.72 | 0.92 | 0.66 |
| GS2 | 0.70 | 0.92 | 0.64 |
| GS1 | 0.64 | 0.92 | 0.59 |

### H.4. Agregación a Gene×Position (§2g)

`group_by(Gene, Position)` sobre los esquemas que dispararon:

- **`CAAS_score = mean(caas_row)`** *(no un máximo, no un sum: el nº de esquemas es una
  propiedad bioquímica del cambio, no evidencia)*.
- `asr_score`, `mrca_diversity`, `derived_agreement`, `conservation_gate`, `core`,
  `phen_score`, **`pos_perm_p`** = **media** sobre esquemas.
- `n_schemes`, `scheme_set = "GS1+GS2+GS3+GS4+US"`, `n_hypotheses = 4`,
  `supporting_hypotheses = "H1,H2,H3,H4"` — **descriptores**; la recurrencia **nunca**
  multiplica `CAAS_score`.
- `change_side ∈ {top, bottom, both, none}` desde `change_top` / `change_bottom`.
- `pos_perm_p_adj` (§2h): BH sobre exactamente las posiciones con match en el null.

*Ejemplo — `(OPN1, 210)`:*
`CAAS_score = mean(0.51, 0.58, 0.66, 0.64, 0.59) ≈ 0.60`;
`asr_score ≈ 0.65`; `pos_perm_p = mean ≈ 0.006`; `change_side = "top"`.

### H.5. Score de gen (§4a) — `size_adj_max`

```
gene_caas_score = size_adj_max(CAAS_score, pool) = F( max(CAAS_score) ) ^ n
```

`F` = ECDF del pool genome-wide de posiciones; `n` = nº de posiciones del gen. `max(x)`,
como cualquier cuantil alto, crece con `n` por artefacto de estadística de orden
(Spearman(n_positions, raw max) ≈ +0.47). `F(max)^n = P(las n extracciones ≤ max)` mide
**cuán improbable es que un gen aleatorio de tamaño n no superase ese valor** → neutraliza
`n` (Spearman ≈ −0.09) y **preserva el orden dentro de una clase de tamaño**. Pools
*direction-matched*: `.pool_all`, `.pool_top` (`change_side ∈ {top,both}`), `.pool_bottom`.
Un guard numérico recomputa 25 genes y aborta si difieren > 1e-6 del null.

*Ejemplo:* la columna 210 (`CAAS_score = 0.60`) es el máximo de `OPN1`. `F(0.60) = 0.94`
(percentil 94 de todas las posiciones). Con `OPN1` de 12 posiciones →
**`gene_caas_score = 0.94^12 ≈ 0.48`**. Si `OPN1` tuviera solo 3 posiciones → `0.94^3 ≈
0.83`: misma evidencia por posición, pero un gen de 3 posiciones que llega a 0.94 es más
notable que uno de 12 (que casi siempre tiene alguna posición por encima del percentil 94).

### H.6. p.perm de gen (§4f, Tier 1A)

`gene_caas_score{,_top_all,_bottom_all}` contra la fila del gen en
`caas_corStat_byrank[[global|top|bottom]]`: p nominal **de cola derecha**
`(Σ null ≥ obs + 1)/(N + 1)`, BH-ajustada por dirección (`gene_caas_pperm_adj*`). Gen fuera
del universo del null → `NA`, nunca `1/(N+1)` implícito. En `caas_full_perms` de producción
se satura hacia `1/(N+1)`: es **priorización con contexto FDR**, no significación
genome-wide por sí sola. Variante opcional `_pooled` = null poolado dentro del decil de
`n_positions` del gen (válido solo dentro de un estrato-n, porque `size_adj_max` es
n-dependiente).

*Ejemplo:* `gene_caas_score(OPN1) = 0.48`; en 1000 ciclos, 12 dan null ≥ 0.48 →
`gene_caas_pperm = 13/1001 ≈ 0.013`.

### H.7. Ejes independientes (§4b–§4e) — `full_join`, columnas separadas

`gene_caas` se une por **`full_join`** (no `left_join`: cada módulo tiene su propio universo
de genes, más ancho que el CAAS) con:

- **Accumulation**: por dirección (`all`/`top`/`bottom`), Cauchy Combination Test
  (CCT/ACAT) sobre las p empíricas de los 5 esquemas → `accum_cct_p`, BH → `accum_fdr`,
  `accum_significant` (< 0.05) sobre genes con ≥ 1 CAAS.
- **RER**: `rer_min_pval` (usa `p.perm` si existe), `rer_significant` (≤ 0.05), `rer_rho`,
  `rer_acceleration`.
- **FADE**: `fade_max_bf_{top,bottom}` (máx Bayes Factor por dirección), `fade_significant`
  (BF ≥ 100).

La tabla se ordena por `gene_caas_score`; **los ejes NO se funden en un compuesto**.

*Ejemplo — fila `OPN1` de `gene_scores.tsv`:* `gene_caas_score 0.48`,
`gene_caas_score_top 0.55`, `gene_caas_pperm 0.013`, `fade_max_bf_top 300`
(`fade_significant_top = TRUE`), `rer_min_pval 0.02` `rer_rho +0.3`
(`rer_acceleration "accelerated"`), `accum_cct_p 0.04`.

---

## I. Uso del null: caracterización y enriquecimientos

### I.1. En la caracterización (scoring)

- Nivel posición: `pos_perm_p` / `pos_perm_p_adj` en `position_scores.tsv`.
- Nivel gen: `gene_caas_pperm{,_top,_bottom}` + `_adj` en `gene_scores.tsv`.

Ambos del **replay ASR verbatim** de la permulación de fenotipo (E.2) → misma escala que el
score observado.

### I.2. En los enriquecimientos FCS — `p.perm` en los **3 tests**  (`fcs_enrich.R::fcs_run_all`)

Rankings de gen (zero-floored sobre el universo del background limpio): para CAAS
`global` / `top` / `bottom` (cada uno con su matriz `caas_corStat_byrank`); para RER
`global` (con signo, `two.sided`) / `accelerating` / `decelerating` (derivados vía
`corRho`). **Los 3 tests comparten la misma matriz de null** — el mismo modelo nulo de
permulación de fenotipo que el `p.perm` de RER.

| test | estadístico observado | `p.perm` |
|------|----------------------|----------|
| **1. Wilcoxon** (`sig_wilcoxon`, `fastwilcoxGMTall`) | rank-sum / AUC del set vs resto | `fcs_permpvalenrich_vectorized`: `(#{null ≥ obs} + 1)/(N_valid + 1)` (`greater`, rankings de magnitud) o `(#{|null| ≥ |obs|}+1)/(N+1)` (`two.sided`, RER con signo). Universo = genes anotados en el GMT. Gate: `p.adj < fdr_wilcoxon` **y** (`p.perm` NA o `< p_perm_thr`) **y** `stat > 0` |
| **2. Lachenbruch two-part** (`sig_lachenbruch`, `fcs_run_lachenbruch`) | Parte 1: **Fisher exacto** en la 2×2 de prevalencia (`score > 0` vs `= 0`) = cola superior **hipergeométrica**, → χ²(1). Parte 2: Wilcoxon `greater` sobre los genes con score > 0, → χ²(1). Suma → χ²(2) → `lach_pval` | `fcs_compute_lach_p_perm`: vectoriza la Parte 1 con **una llamada `phyper()`** por GMT sobre `set × columna-de-null` (universo = dominio completo del ranking) + Parte 2 vectorizada; combina y compara contra el `lach_chi_total` observado. Gate: `lach_p.adj < fdr_lachenbruch` **y** (`lach_p.perm` NA o `< p_perm_thr`) |
| **3. Path-sum permulation** (`sig_permulation`, `fcs_run_permulation`) | suma de scores de gen del pathway; `NES = (obs − null_mean)/null_sd` | `(rowSums(null ≥ obs) + 1)/(N + 1)`; usa la matriz de null compartida cuando existe, si no un label-shuffle privado. Gate: `perm_p.adj < fdr_permsum` **y** `perm_nes > 0` |

`evidence_count` = suma de los 3 gates → `evidence_label` ("Hard evidence" = 3,
"Supported" = 2, "Exploratory" = 1). Lachenbruch y path-sum se saltan para rankings con
signo (RER global).

*Ejemplo:* `OPN1` (con `score_global = 0.48`) es miembro del pathway
`GO:phototransduction`. En el ranking `global`, ese set obtiene Wilcoxon `stat > 0`,
`p.adj = 0.03`, `p.perm = 0.008` (8 de 1000 columnas de `caas_corStat_byrank[["global"]]`
dan un rank-sum ≥ el observado) → `sig_wilcoxon = TRUE`. Lachenbruch Parte 1 (¿el pathway
está enriquecido en genes con `score > 0`?) `lach_p.adj = 0.04`, `lach_p.perm = 0.02` →
`sig_lachenbruch = TRUE`. Path-sum `perm_p.adj = 0.06` → `sig_permulation = FALSE`.
`evidence_count = 2` → **"Supported"**.

### I.3. Enriquecimiento posicional (`posenrich_enrich.py`)

Test independiente con **permulación propia** (de **posiciones**, no de fenotipo):
*path-sum* por término, `rng.choice(N, K, replace=False)` por permutación,
`M · P` disperso. `perm_nes = (obs − null_mu)/null_sd`;
`p_value = (Σ null ≥ obs + 1)/(n_perms + 1)`. Direcciones `top`/`bottom`/`global`
filtradas por `change_side`, ranking sobre posiciones con `CAAS_score > 0`.

---

## Resumen del ejemplo (columna 210 de `OPN1`)

| etapa | valor |
|-------|-------|
| Selección | K = 3 pares (H1: A\|A, B\|B, C\|C) |
| FOP | 4 hipótesis (H1–H4; alternativas A2 y C2) |
| Discovery | filas para US/GS4/GS3/GS2/GS1 × H1–H4; patrón `WWW/YYY` (H1,H2), `WWF/YYY` (H3,H4) |
| `recovery_boot` (D/E.1) | 0.075 (3/40 ciclos base) |
| ASR por hipótesis (F.3) | asr ≈ 0.72 (H1), 0.47 (H3) bajo US; ≈ 0.71 todas bajo GS3 |
| FOP pooling (H.2) | `asr_pooled` ≈ 0.55 (US), 0.72 (GS3); `da_US ≈ 0.67`, `da_GS3 = 1.0`; `convergence_schemes = "GS4,GS3,GS2,GS1"` |
| `caas_row` (H.3) | 0.51 (US) … 0.66 (GS3) |
| `pos_perm_p` (E.2) | 0.006 |
| `CAAS_score` posición (H.4) | ≈ 0.60 (media sobre esquemas), `change_side = "top"` |
| `gene_caas_score` (H.5) | 0.94^12 ≈ 0.48 |
| `gene_caas_pperm` (H.6) | ≈ 0.013 |
| Ejes independientes (H.7) | FADE BF 300, RER p.perm 0.02 (acc), accum CCT p 0.04 — columnas separadas |
| FCS (I.2) | `GO:phototransduction` "Supported" (Wilcoxon + Lachenbruch, no path-sum) |

---

## Advertencias

- Un `max_fop` ancho recupera más sitios reales pero admite contrastes marginales en el
  harvest (una hipótesis puede meter en el foreground una especie que porta el residuo
  "derivado" por causas no adaptativas), lo que ensucia `top_residue_support` al poolear.
- `pos_perm_p` y `gene_caas_pperm` a `caas_full_perms` de producción se saturan en
  `1/(N+1)`: priorizan, no declaran significación genome-wide por sí solos.
- La maquinaria a nivel gen / FCS (`caas_perms.rds`) está construida para input
  genome-wide; en un fixture de 1–3 genes no lleva señal.
- El gate PSS continuo por defecto es un percentil de rango que escala con el tamaño del
  árbol, no con la señal (cuestión de diseño abierta).
- Comentarios obsoletos en el código (`ctpp_signification.nf:2` "hypergeometric and
  permutation", `scoring_compute.R` "gate_all / gate_sig", docstring de
  `caap_id.py::fetch_caap` con `pvalue`) NO reflejan el comportamiento actual: no hay test
  hipergeométrico ni gate de significación por sitio.
