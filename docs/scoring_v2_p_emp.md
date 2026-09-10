# scoring_v2 · Diseño de `p.emp` — p empírico de posición «detecta Y supera»

Rama `scoring_v2`. Entregable de diseño (sin código). `p.emp` opera sobre el
`CAAS_score` de posición ya agregado (§2g de
[`scoring_compute.R`](../subworkflows/SCORING/local/src/scoring_compute.R)) y no
toca la construcción del `core`.

**Integrado en la secuencia core v3.** El rediseño del `core` de dominio Voronoi
(`docs/scoring_v3_core.md`, Apéndice D) va por commits `V3-0..V3-6`. Estado:
secuencia `V3-0..V3-6` **COMPLETA** en `scoring_v2`. El plumbing del null de
`p.emp` entró en `V3-4a`; el consumidor R (§2f-ter) en `V3-4`. El flip de headline
(§7.3) y el borrado de `null_pvalue_boot` (§7.4) están **HECHOS** (commits
`p.emp §7.3:` / `p.emp §7.4:`), adelantados por decisión del usuario antes de la
corrida PEPC — esa corrida pasa de gate a verificación: si la distribución de
`p.emp` bajo fenotipo nulo-por-construcción se apila en 0, hay que **revertir** el
§7.3 (volver `pos_perm_p_adj` al headline).
Core v3 §4 fija como invariante justo lo que `p.emp` asume: `CAAS_score` = media
§2g de `core_s` sobre los 5 esquemas, `caas_row = asr_path_score = core_s`, shard
de 8 columnas, esquema de `perm_pos_pval.tsv`. Core v3 cambia los **valores** de
`CAAS_score` (Apéndice C: suben), no la estructura.

Cotejado contra
[`gene_wrapper.py`](../subworkflows/CT_DISAMBIGUATION/local/src/utils/gene_wrapper.py)
(`_perms_worker`, `_finalize_perm_scores`, `_finalize_perm_pos_pval`,
`process_all_genes_perms`),
[`reaggregate_perm_scores.py`](../subworkflows/CT_DISAMBIGUATION/local/reaggregate_perm_scores.py),
[`scoring_compute.R`](../subworkflows/SCORING/local/src/scoring_compute.R)
§2f-bis / §2g / §2h / §4a / §4b / §4d / §4f,
[`scoring_caas_perms.R`](../subworkflows/SCORING/local/src/scoring_caas_perms.R),
[`randomize.py`](../subworkflows/CT_ACCUMULATION/local/src/randomization/randomize.py)
(`run_permulation_null`) y
[`fcs_enrich.R`](../subworkflows/ENRICHMENT/local/src/fcs_enrich.R). Estado del
árbol en `665a0b7` (V3-3).

Coordenadas 0-based. `s ∈ {top, bottom}` es la clave de dirección (T4b: `side`
única). `N` = nº de ciclos base del null (tras colapso FOP `~H*` si aplica).

---

## 0. Resumen de decisiones

| decisión | valor |
|---|---|
| `p.emp` | `(k_emp + 1)/(N + 1)`, `k_emp` = ciclos que **re-detectan Y superan** el `CAAS_score` observado, por `(Gene, Position)` — **pooled a posición** (ver amendment V3-4a) |
| regla de score de la superación | **max sobre lados** de la media §2g por esquema — el eje "all"/`.pos_undirected` (amendment V3-4a; era per-lado) |
| dónde se computa la media por ciclo | **Opción C**: el null emite `(caas_sum, n_schemes)` por `(Gene, Position, side, cycle)`; R divide y cuenta (§5) |
| secuenciación | plumbing null: V3-4a · consumidor R: V3-4 · flip de headline: commit `p.emp §7.3:` (HECHO, verificación PEPC pendiente) · borrado `null_pvalue_boot`: commit `p.emp §7.4:` (HECHO) |
| headline de posición | `p.emp` / `p.emp_adj` — **sustituye** a `pos_perm_p_adj`; `pos_perm_p` fuera de `position_scores.tsv`, solo en `perm_pos_pval.tsv` |
| `pos_perm_p` | pooled a `(Gene, Position, caap_group)`; baja a diagnóstico (`perm_pos_pval.tsv`), fuera del headline y de `position_scores.tsv` |
| `null_pvalue_boot` | **borrado** (commit `p.emp §7.4:`) — inerte bajo el `percent_rank` de aguas abajo, sin consumidor externo |
| `p.emp_score_only`, `gene_perm_p_detect` | **descartados** (ver §3) |
| p de gen | sin cambios: `gene_caas_pperm` (§4f, magnitud) + `accum_cct_p` (§4b, conteo) — dos ejes, nunca combinados (§6e) |
| corrección LOO en `p.emp` | **no** — add-one puro (§6d) |

---

## 1. Objetivo y definición

Recuperar un p-valor de posición que combine **detección** y **superación de
score** sobre el mismo eje de remuestreo (relabelado de fenotipo) que ya usan los
p-valores del null existentes:

> Se cuenta un ciclo `i` del null a favor de la hipótesis nula **sii** en ese
> ciclo la posición `(Gene, Position, s)`
> **(a)** vuelve a detectarse como CAAS en el lado `s`, **y**
> **(b)** su `CAAS_score` en ese ciclo —media sobre los esquemas detectados,
> idéntica regla que §2g— iguala o supera el `CAAS_score` observado.

Con `k_emp` = nº de ciclos que cumplen (a)∧(b):

```
p.emp = (k_emp + 1) / (N + 1)              add-one (Davison & Hinkley), right-tailed
```

Acumulación exactamente la de `pos_perm_p` (add-one, no leave-one-out; §6d);
comparación de score exactamente la agregación de §2g (media sobre esquemas, por
lado; §4).

---

## 2. Los dos p-valores del null que ya existen

### 2.1 Nivel posición, solo-detección — `_perms_worker` (bloque de detección ≈1299-1366 tras V3-3)

`process_all_genes_perms` reproduce, gen a gen, `N` relabelados del fenotipo sobre
la ASR cacheada y vuelve a correr la disambiguación (`analyze_gene_disambiguation`
sobre las entradas de `perm_discovery_file`, es decir el *bootstrap discovery*
por-ciclo de `caas_permulation.nf` → `ct` bootstrap full-pool con
`--export_perm_discovery`). Para cada `(Gene, Position, caap_group)` cuenta en
cuántos ciclos re-apareció:

- `k` = `len(cycles_set)` — nº de ciclos que re-detectaron la candidata.
  **Hoy side-agnóstico** (≈1332-1345, intacto tras V3-3): se cuenta por
  `(pos, caap_group)` y el `side` solo entra en `sides_seen` para difundir el
  **mismo** `k` a una fila de `perm_pos_pval` por lado. Las filas que recorre
  vienen de `_expand_pooled` (V3-3, sustituto de `_expand_sides`), que ya emite
  ≤2 `PositionAxes` por `(cyc, pos, grp)`, cada una con su `side` autoritativo y
  el colapso `side="none"`. Ver §3: pasa a side-consciente por presencia de fila.
- ~~`null_pvalue_boot = (k − 1) / max(N − 1, 1)`~~ — **borrado** en el commit
  `p.emp §7.4:`. Era la fracción de replicación
  leave-one-out, análogo null del `recovery_boot` observado. **Numéricamente
  inerte** bajo el `percent_rank` de aguas abajo —`(k−1)/(N−1)` y `k/N` ordenan
  igual—; solo sobrevivía porque `perm_pos_pval.tsv` se lee como tabla de
  contraste. Candidato a borrado (§3).
- `pos_perm_p = (k + 1) / (N + 1)` — p permulacional calibrado, add-one, **no**
  LOO. Por `(Gene, Position, caap_group)`.

**Ninguno mira el score.** `asr_path_score` viaja en el shard
`perm_pos_detail/<Gene>.tsv.gz` (8 columnas: `Gene, cycle, Position, caap_group,
asr_path_score, n_detected, clust, side`) pero hoy solo lo consumen los ejes de
gen (§2.2).

**Ruta en R.** §2f-bis (≈325-392) hace `left_join` de `pos_perm_p` sobre `df`
—granularidad `(Gene, Position, caap_group, side)`, con *fallback* de 3 claves
`pos_perm_p_3key` cuando la etiqueta `side` difiere entre observado y null—
**antes** de §2g, de modo que `pos_perm_p` cabalga la misma `mean()` sobre
esquemas que cualquier otro eje per-esquema (línea 441). §2h (≈473-485) hace
`p.adjust(..., "BH")` sobre exactamente las posiciones con match →
`pos_perm_p_adj`. La clave del join **ya incluye `side`** aunque el valor no
varíe por lado hoy: side-consciente (§3) le da contenido.

### 2.2 Nivel gen — dos ejes distintos

**`gene_caas_pperm` (§4f, magnitud).** `_finalize_perm_scores` agrega el shard a
`gene_cycle_scores.tsv` (`_size_adj_max_null` = `(bisect_right(pool_ciclo, max) /
|pool|) ** n`, espejo exacto de `size_adj_max` observada, §4a);
`scoring_caas_perms.R` pivota a `caas_corStat_byrank` (genes×N, denso: ceros
estructurales conservados) y estampa `gene_stat = "size_adj_max"`. §4f: por
fila-gen, right-tailed contra su propia fila null, `(Σ[nr ≥ obs] + 1)/(N + 1)`,
BH dentro de dirección → `gene_caas_pperm_adj*`. Gen fuera del universo del null
→ `NA`, nunca `1/(N+1)` (break-point #4a). §4d es el guard numérico que mantiene
`size_adj_max` (R) y `_size_adj_max_null` (Py) bit-compatibles (`F(max)^n`
amplifica deriva de ~1e-16). `fcs_enrich.R` consume la misma matriz para el
`p.perm` de vías.

**`accum_cct_p` (§4b, conteo).** `randomize.py` cuenta, por gen y esquema, las
posiciones detectadas y saca `PValueEmpirical = (#{null_count ≥ obs_count} +
1)/(N + 1)`. §4b combina los 5 esquemas con Cauchy CCT/ACAT → `accum_cct_p`, BH
→ `accum_fdr`. El null depende de `accumulation_randomization_type`
(`conf/ct_accumulation.config:10`, **default `cons_decile`**):

- `cons_decile` (default): null de **ocupación**, redibuja el conteo observado de
  un pool decil-de-conservación. Controla longitud/cobertura + conservación;
  **no** controla el árbol-fenotipo.
- `permulation` (opt-in): lee `run_permulation_null` del **mismo**
  `perm_pos_detail/` que §4f. Controla el árbol-fenotipo; **no** es
  decil-de-conservación.

### 2.3 Dónde encaja `p.emp`

`p.emp` es, por construcción, una cantidad **post-§2g**: la comparación (b) es
contra el `CAAS_score` observado, que solo existe tras promediar los 5 esquemas
por lado. No puede cabalgar la `mean()` de §2g como hace `pos_perm_p`. Su lugar
es **junto a §2h**, unido a `pos_scores` por `(Gene, Position, side)`
—cardinalidad ya 1:1— y BH-ajustado en el mismo bloque (§5, §6a-b).

---

## 3. El conjunto de p-valores (recortado)

La rejilla conceptual es 2×2: {posición, gen} × {solo-detección/conteo,
detección+score}.

| | solo detección / conteo | detección + score |
|---|---|---|
| **posición** | `pos_perm_p` (side-consciente, diagnóstico) | **`p.emp` / `p.emp_adj`** (headline) |
| **gen** | `accum_cct_p` / `accum_fdr` (§4b) | `gene_caas_pperm(_adj)` (§4f) |

Tres de las cuatro celdas **ya están pobladas**. El diseño puebla una nueva
(`p.emp`), hace side-consciente otra (`pos_perm_p`) y no añade nada al nivel gen.

**Decisiones:**

- **(a) `p.emp` sustituye a `pos_perm_p_adj` como headline de posición.**
  Justificación mecanística, no heurística: `k_emp ≤ k_det` por definición ⇒
  `p.emp ≤ pos_perm_p` **siempre**. `p.emp` domina en potencia: no declara
  significativa una posición que un relabelado tanto re-detecta como iguala en
  score. Headline = columna ordenadora en `gene_lists/*_pos*`,
  `position_scores.tsv`, Rmd de SCORING/POSENRICH.
  *Límite de la inferencia:* la mayor potencia solo es mejora si la distribución
  null del `CAAS_score` está bien calibrada. Si los relabelados producen
  `asr_path_score` sistemáticamente bajos por una razón estructural (topología:
  fg permutados que rara vez recrean la profundidad de MRCA real), `p.emp` será
  anticonservador. Contraste Tier 1 PEPC pendiente: bajo fenotipo
  nulo-por-construcción, la distribución de `p.emp` debería salir ~uniforme; si
  se apila en 0, quedarse con `pos_perm_p` como headline.

- **`pos_perm_p` pasa a side-consciente y baja a diagnóstico.** La detección en
  el null ya se emite per-lado (`_expand_pooled` tras V3-3; el shard lleva `side`
  nativo); lo único side-agnóstico es el conteo. Cambio: llavear `n_detected` por
  `(pos, group, side)` en el bloque de detección de `_perms_worker` (≈1332-1345)
  y emitir un `k` por lado;
  `_finalize_perm_pos_pval` / `reaggregate_perm_scores.py` recalculan
  `k_side` = nº de filas de shard que casan `(Gene, Position, caap_group, side)`
  (una fila por ciclo/lado/esquema/posición). Efecto: `k_side ≤ k` agrupado ⇒
  `pos_perm_p` per-lado sube (el valor actual poolea lados, es optimista). Queda
  **solo en `perm_pos_pval.tsv`**, fuera de `position_scores.tsv` y de los
  gene_lists. Sirve como descomposición apples-to-apples de `p.emp` en el mismo
  `(Gene, Position, side)`:
  - `p.emp` sig. **y** `pos_perm_p` sig. → ambos canales;
  - `p.emp` sig. **y** `pos_perm_p` no → significancia solo por score (aviso de
    calibración).
  Si por lo que sea NO se hace side-consciente, entonces **quitarlo**: un
  companion de detección en distinta unidad que `p.emp` no es chequeo válido.

- **`null_pvalue_boot`**: inerte por diseño. Candidato a borrado en una pasada
  de limpieza aparte; no se hace side-consciente.

- **Descartados del plan:**
  - `p.emp_score_only` (p de posición solo-score, sin gate de detección) — era
    puro diagnóstico de descomposición. Recuperable en cualquier momento desde
    `perm_pos_cycle_caas.tsv.gz`. No es columna estándar.
  - `gene_perm_p_detect` (companion solo-detección a nivel gen) — **redundante
    con `accum_cct_p`**, que ya ocupa la celda "gen / conteo". No se añade.

---

## 4. La regla de score de (b): idéntica a §2g

`CAAS_score` observado de `(Gene, Position, s)` (§2g, línea 409):

```
CAAS_score_obs(G,P,s) = mean_{k ∈ schemes_obs(G,P,s)}  asr_score_k
```

media sobre las filas per-esquema presentes (`caap_group`), `na.rm=TRUE`.

`CAAS_score` del null en el ciclo `i`:

```
CAAS_score_null(G,P,s,i) = mean_{k ∈ schemes_null(G,P,s,i)}  asr_path_score_k^{(i)}
```

misma operación sobre las filas del shard con `cycle == i`, `Position == P`,
`side == s`, agrupando por `caap_group`. Coincide término a término con §2g
porque:

- `caap_group` en el shard **es** el esquema (US/GS4/GS3/GS2/GS1);
- el shard solo trae los esquemas que ese ciclo detectó, igual que `df` en §2g
  solo trae los que el observado detectó — ambas medias son sobre "esquemas
  presentes", no sobre 5 fijos;
- ya es la media que `_build_cycle_score_pools._drain_side` calcula de forma
  transitoria (`entry[0]/entry[1]` por `(cyc,pos,side)`) para el pool de §4f — no
  es fórmula nueva.

**Detección (a) por lado.** Ciclo `i` detecta `(G,P,s)` sii existe ≥1 fila de
shard `(i, P, *, s)`. La presencia de fila ya codifica (a); no hace falta columna
nueva.

---

## 5. (b) ¿Dónde ocurre la comparación observado-vs-null por ciclo?

La media por ciclo (§4) hay que formarla en algún sitio. Tres opciones.

### 5.1 Opción A — en `_perms_worker`, antes de escribir el shard

- **Contra, decisivo:** el worker **no** conoce `CAAS_score_obs` (se calcula en R,
  §2g). Solo podría emitir la media null por ciclo, no el conteo `k_emp`. Y
  duplicaría la regla de media de §2g en Python → nuevo break-point tipo #11, con
  guard.

### 5.2 Opción B — pase finalizador que ingiere `position_scores.tsv`

- **Contra, decisivo:** invierte el DAG. `CAAS_PERMS_DISAMBIGUATE` es upstream /
  paralelo a `SCORING`; hacer que el null dependa de un artefacto del observado
  obliga a un segundo proceso post-SCORING. Complejidad de orquestación real por
  una columna. También re-implementa la media de §2g → guard.

### 5.3 Opción C — emitir `(caas_sum, n_schemes)` por ciclo; **R divide y cuenta** (recomendada)

`_finalize_perm_scores` gana un output compacto `perm_pos_cycle_caas.tsv.gz` con,
por `(Gene, Position, side, cycle)`:

| col | significado |
|---|---|
| `caas_sum` | `Σ_k asr_path_score_k` sobre esquemas de ese ciclo/lado |
| `n_schemes` | `#k` (>0 ⇒ ciclo detectó ese lado) |

Sale **gratis**: es el `entry = [asr_sum, n, caas_sum]` que `_flush` ya acumula
por `(cyc, pos, side)` en pass B2 (líneas 1795-1800); solo hay que volcarlo en
vez de descartarlo tras el pool / q90.

Luego, en R, nueva **§2f-ter** (tras §2g, dentro o antes de §2h):

```
para cada (Gene, Position, s) de pos_scores con CAAS_score_obs no-NA:
    filas_null ← perm_pos_cycle_caas[Gene, Position, s]        # ≤ N filas
    k_emp ← #{ filas_null : n_schemes > 0
                            ∧ (caas_sum / n_schemes) ≥ CAAS_score_obs }
    p.emp ← (k_emp + 1) / (N + 1)
```

- **Media única.** `caas_sum / n_schemes` lo hace R, con la misma aritmética IEEE
  que `mean(caas_row)` de §2g (mismo orden de suma si el emisor recorre
  `caap_group` en el orden de prioridad de §2g). **No hay segunda implementación
  de la media ⇒ no hace falta guard nuevo tipo §4d.** Argumento decisivo frente a
  A/B.
- **DAG intacto.** Emisor = null; consumidor = R; `p.emp` vive donde
  `pos_perm_p_adj`.
- **Tamaño.** `N` × (posiciones detectadas × lados). Del orden del shard actual
  ÷ 5 (sin la dimensión esquema), gz. Medir en el primer `caas_full_perms` real
  (§8).
- **`reaggregate_perm_scores.py`** gana la misma llamada al finalizador (ya reusa
  `_finalize_perm_scores`), así que un rebuild sin ASR reproduce
  `perm_pos_cycle_caas.tsv.gz` desde los shards.

**Veredicto:** Opción C.

---

## 6. Puntos de interacción

### 6a. Cardinalidad del join

`pos_perm_p` se une en §2f-bis a `df` (per-esquema) y depende de la `mean()` de
§2g para colapsar. `p.emp` **no**: se une a `pos_scores` (post-§2g, ya 1:1 en
`(Gene, Position, side)`).

- Join cardinalidad-neutro por construcción; no necesita el `stopifnot(nrow(df)
  == .n_rows_before)` de §2f-bis.
- Reusa el mismo *fallback* de clave: intento por `(Gene, Position, side)`,
  `coalesce` con `(Gene, Position)` 3-key para filas cuyo `side` difiera entre
  observado y null (post-T3 debería ser 0). Con `pos_perm_p` ya side-consciente,
  ambos p de posición usan exactamente la misma clave y el mismo fallback.
- Mismo guard de tasa de match < 50 % de §2f-bis (coordenadas
  `filtered_discovery.tsv` vs shard): si el join de `pos_perm_p` está sano,
  `p.emp` hereda la salud.
- Posición "both": dos filas en `pos_scores` (`side` top/bottom), cada una con su
  `CAAS_score` ⇒ `p.emp_top` ≠ `p.emp_bottom` en general. Con `pos_perm_p`
  side-consciente los dos p de posición son per-lado y consistentes.

### 6b. `p.adjust` en §2h

`p.emp_adj = p.adjust(p.emp[tested], "BH")`, `tested = !is.na(p.emp)`, mismo
patrón que `pos_perm_p_adj`. Cada `(Gene, Position, side)` es un test. `p.emp` y
`pos_perm_p` se BH-ajustan **por separado** (dos familias de hipótesis; el
usuario elige headline, no se penaliza por reportar ambas). Sin
`p.emp_score_only` (§3).

### 6c. Guard `size_adj_max` (§4d) y `_size_adj_max_null`

`p.emp` **no interactúa** con `size_adj_max` ni `_size_adj_max_null`: es
nivel-posición, compara `CAAS_score` crudo, no `F(max)^n`. El guard §4d sigue
cubriendo solo el eje de gen.

El requisito de paridad que `p.emp` introduce (Opción C) es más débil:
`caas_sum/n_schemes` (en R) debe igualar `mean(caas_row)` de §2g. Ambas
divisiones en R, mismos sumandos ⇒ exacto salvo orden de suma. Mitigación: el
emisor C acumula `caas_sum` recorriendo `caap_group` en el orden de prioridad de
§2g (`US > GS4 > GS3 > GS2 > GS1`). Con eso, sin guard. Cinturón-y-tirantes
opcional: check ligero «re-deriva `CAAS_score_null` de 10 tuplas desde el shard
de 8 col y compara».

### 6d. ¿La corrección LOO (`k − 1`) se generaliza a «detecta Y supera»?

**No, y no debe.** El LOO de `null_pvalue_boot` corrige auto-inclusión: al puntuar
el ciclo `i` como muestra de replicación, `i` no debe contar en su propia
evidencia. `pos_perm_p` ya renuncia al LOO (add-one Davison-Hinkley) porque es un
p calibrado, no una fracción de replicación.

`p.emp` es un p de permutación estándar: el estadístico observado
(`CAAS_score_obs`) **no** es uno de los `N` ciclos del null (viene del análisis
real). No hay auto-inclusión que corregir. Add-one puro `(k_emp + 1)/(N + 1)`,
que además lo alinea con `pos_perm_p` y con el p de gen §4f. Aplicar
`(k_emp − 1)/(N − 1)` sería un error categórico (trataría el observado como un
ciclo null).

### 6e. Los dos p de gen — qué debe quedar clarísimo

> **`gene_caas_pperm` (§4f) y `accum_cct_p` (§4b) son dos ejes. Nunca se
> combinan.**
>
> - **Pregunta distinta.** §4f: magnitud n-ajustada — `F(max)^n` de la mejor
>   posición del gen, calibrada por cuántas posiciones tuvo oportunidad de
>   sortear. §4b: conteo — cuántas posiciones detectadas, vs lo que predice el
>   null, CCT sobre los 5 esquemas.
> - **Null distinto en el default.** `accumulation_randomization_type =
>   "cons_decile"` ⇒ §4b es null de ocupación decil-de-conservación: controla
>   longitud/cobertura de gen + conservación, **no** el árbol-fenotipo. §4f
>   controla el árbol-fenotipo, **no** la conservación. En este modo son casi
>   ortogonales y aportan evidencia parcialmente independiente. Un gen con muchos
>   CAAS débiles puntúa en §4b y no en §4f; uno con un CAAS muy fuerte al revés.
>   Ambos son modos de descubrimiento reales.
> - **Solo comparten null en modo `permulation`** (opt-in,
>   `run_permulation_null` lee el mismo `perm_pos_detail/` que §4f). Ahí quedan
>   correlacionados (conteo vs máximo; `Spearman(n_positions, raw max) = +0.47`
>   según el comentario de §4a, ~50 % compartido como mucho) y hay que leerlos
>   como «dos estadísticos del mismo null», no dos experimentos.
> - **Nunca un p combinado de los dos.** Se reportan lado a lado y se correlacionan
>   en §5 de `scoring_compute.R`, igual que RER/FADE (líneas 1134-1135:
>   "represented by their native significance"). Combinar dos p correlacionados
>   del mismo null exigiría CCT/ACAT y seguiría sin ser evidencia independiente.

**`p.emp` (posición) y `gene_caas_pperm` (gen) son el par coherente por nivel de
decisión**, misma familia empírica (add-one, right-tailed, mismos `N` relabelados,
BH-en-tested), difieren solo en la agregación —que *debe* diferir—:

| | `p.emp` (posición) | `gene_caas_pperm` (§4f) |
|---|---|---|
| unidad | `(Gene, Position, side)` | `(Gene, direction)` |
| estadístico obs | `CAAS_score` (media esquemas) | `size_adj_max` = `F(max)^n` sobre `CAAS_score` de las posiciones del gen |
| comparación | `CAAS_score_null(i) ≥ obs`, por ciclo | `size_adj_max_null(i) ≥ obs`, por ciclo |
| gate de detección | **sí** — el ciclo debe re-detectar el lado | **no** — un gen «tiene» siempre sus posiciones (la no-detección entra por los ceros estructurales de `gene_cycle_scores.tsv`) |
| calibración de pool | ninguna (comparación directa) | ECDF por-ciclo (`_build_cycle_score_pools`) dentro de `F` |

**No son evidencia independiente para un mismo locus.** `gene_caas_pperm` es
función determinista de los mismos scores null de posición que alimentan `p.emp`
(vía `_build_cycle_score_pools` → `_size_adj_max_null`). Están anidados: un gen
significativo + una posición significativa dentro de él son *una* línea de
evidencia a dos resoluciones, no dos confirmaciones.

`p.emp` **no** es «el §4f llevado a posición»: el §4f-a-posición sería un `p.emp`
solo-score, descartado (§3).

---

## 7. Cambios concretos y secuenciación

> **Amendment V3-4a (implementado).** Decisión del usuario en el checkpoint de
> V3-4: `p.emp` y `pos_perm_p` **NO** van per-lado. Ambos se poolean a
> `(Gene, Position)` sobre el eje **max-sobre-lados** (el mismo que
> `scoring_compute.R` `.pos_undirected` / `_build_cycle_score_pools` `pc["all"]`):
>
> - **detección (a)**: un ciclo re-detecta la posición si re-detecta **cualquier
>   lado**.
> - **score del ciclo**: `max_lado( Σ_esquemas asr / n_esquemas )` — la media §2g
>   por lado, y luego el máximo sobre los lados que ese ciclo detectó.
> - **score observado**: `max_lado( CAAS_score )` = `.pos_undirected`.
> - `p.emp` / `p.emp_adj` se **broadcastean idénticos** a las dos filas de lado de
>   una posición "both" en `position_scores.tsv`.
> - `perm_pos_pval.tsv` **pierde la columna `side`** (una fila por
>   `(Gene, Position, caap_group)`); `pos_perm_p` sigue per-`(pos, group)` (ya lo
>   era en V3-3, solo se retira el broadcast por lado).
> - **`pos_perm_p_adj` se QUEDA en `position_scores.tsv`** (headline de posición)
>   hasta el flip post-V3-6 — resuelve la contradicción §7.2 vs §7.3 a favor de
>   §7.3.
>
> El plumbing del null (§7.1) se implementó en el commit `V3-4a`: nuevo
> `perm_pos_cycle_caas.tsv.gz` (columnas `Gene, Position, side, cycle, caas_sum,
> n_schemes`; `side` se mantiene **en este fichero** para que R tome el max),
> `perm_pos_pval.tsv` sin `side`, y el cableado Nextflow
> (`--caas_pos_cycle_caas`). El consumidor R (§7.2) va en `V3-4`.

`V3-3` ya está en el árbol (`665a0b7`) y **no** llevó nada de `p.emp` (reescribió
la expansión per-lado `_expand_sides` → `_expand_pooled`; el bloque de detección
y los finalizadores quedaron intactos). El plumbing de `p.emp` es aditivo sobre
ese estado. Reparto en tres tramos:

### 7.1 Commit `p.emp-null` — ahora, sobre V3-3 (independiente de V3-4)

**`gene_wrapper.py`:**

- `_perms_worker`, bloque de detección (≈1322-1345): llavear `n_detected` por
  `(pos, caap_group, side)` en vez de `(pos, caap_group)`; emitir un `k` (y
  `pos_perm_p`, `null_pvalue_boot`) **por lado**. Retirar la difusión `sides_seen`
  del mismo `k` a los dos lados. Las filas ya vienen per-lado de `_expand_pooled`,
  así que "detecta el lado `s`" = existe registro con ese `side` (fork 3 §8
  resuelto por V3-3).
- `_finalize_perm_pos_pval` (≈1830-1880) y `reaggregate_perm_scores.py`: contar
  `k_side` = nº de filas de shard que casan `(Gene, Position, caap_group, side)`
  (una fila por ciclo/lado/esquema/posición).
- `_finalize_perm_scores` (`_flush`, ≈1680-1735): volcar `(caas_sum, n_schemes)`
  por `(cyc, pos, side)` a `perm_pos_cycle_caas.tsv.gz`
  (`Gene, Position, side, cycle, caas_sum, n_schemes`); acumular `caas_sum` por
  `caap_group` en el orden de prioridad de §2g.
- `process_all_genes_perms` / `reaggregate_perm_scores.py`: declarar el output
  nuevo (el segundo sale con la misma llamada al finalizador).
- Docstrings: `perm_pos_cycle_caas` es la entrada de `p.emp`; la media por-ciclo
  replica §2g término a término.

**Nextflow:** publicar `perm_pos_cycle_caas.tsv.gz` desde `caas_permulation.nf` y
cablearlo como `--caas_pos_cycle_caas` a SCORING, junto a `--caas_pos_pval`.

**Tests:** extender `test_null_domain_pool_wiring.py` (o uno nuevo) para fijar
`k_side` per-lado y el volcado `caas_sum/n_schemes` contra los escenarios
Apéndice B.

Este commit es válido en aislamiento: añade un output, hace `pos_perm_p`
per-lado, y no toca R. `perm_pos_cycle_caas.tsv.gz` ya sale con valores core v3
(V3-1/-2/-3 ya migraron `compute_domain_scores` / `pool_domains`).

### 7.2 Dentro de V3-4 — el consumidor R

`V3-4` ya reescribe `scoring_compute.R` (schema downstream, `position_scores.tsv`).
`p.emp` entra ahí, no en un pase aparte, para no migrar el schema dos veces:

- Nueva §2f-ter (tras §2g): cargar `perm_pos_cycle_caas.tsv.gz`, join a
  `pos_scores` por `(Gene, Position, side)` con fallback 3-key, `k`-loop →
  `p.emp`. Guard de tasa de match análogo a §2f-bis.
- §2h: `p.emp_adj` vía `p.adjust("BH")` sobre tested; BH separado de `pos_perm_p`.
- Escritura de `position_scores.tsv`: añadir `p.emp`, `p.emp_adj`; **quitar**
  `pos_perm_p` / `pos_perm_p_adj` de esa tabla (pasan a vivir solo en
  `perm_pos_pval.tsv`).
- `scoring_caas_perms.R`, `randomize.py`: sin cambios.

### 7.3 Flip de headline — HECHO (commit `p.emp §7.3:`)

Adelantado antes de la corrida PEPC por decisión del usuario ("cerrarlo todo ya").
Cambios:

- `scoring_compute.R`: §2f-bis reducido a leer **solo `n_cycles`** de
  `perm_pos_pval.tsv` (el N del add-one de `p.emp`); el join de `pos_perm_p` a
  `df` desaparece. §2g deja de agregar `pos_perm_p`. §2h deja de calcular
  `pos_perm_p_adj`. `pos_out` (§6) escribe `position_scores.tsv` **sin**
  `pos_perm_p` / `pos_perm_p_adj` — solo `p.emp` / `p.emp_adj` (decisión del
  usuario: "quitar de position_scores.tsv", resuelve la contradicción §7.2 vs
  amendment V3-4a a favor de §7.2).
- `11.Scoring_report.Rmd`: sección headline reescrita sobre `p.emp` / `p.emp_adj`
  (leídos de `position_scores.tsv`), con nota de que bajo fenotipo nulo debería
  salir ~uniforme. `pos_perm_p` baja a sub-sección diagnóstica que lee
  `perm_pos_pval.tsv` directo y BH-ajusta ahí mismo (ya no viene ajustado de
  aguas arriba). Nuevo param `scoring_p_emp_thr` (default 0.1); `scoring_pos_perm_p_thr`
  se conserva para la sub-sección diagnóstica.
- `conf/scoring.config` + `scoring_report.nf`: `scoring_p_emp_thr` cableado.
- `position_scores.tsv` y `position_lists/` ya ordenaban por `CAAS_score`, no por
  `pos_perm_p_adj` — no hay reordenación de ficheros que hacer, solo la prosa de
  reporte y el esquema de columnas.

**Verificación pendiente (NO gate):** la corrida Tier 1 PEPC de v3. Si `p.emp`
bajo fenotipo nulo-por-construcción se apila en 0 → anticonservador → **revertir**
este commit (`pos_perm_p_adj` vuelve al headline, `pos_perm_p` vuelve a
`position_scores.tsv`).

### 7.4 Borrado de `null_pvalue_boot` — HECHO (commit `p.emp §7.4:`)

Confirmado por el usuario que ningún consumidor externo lo lee. Borrado de
`_perms_worker` (cálculo `(k-1)/loo_denom`), `_finalize_perm_pos_pval`, esquema de
`perm_pos_pval.tsv` (ahora `Gene, Position, caap_group, n_detected, n_cycles,
pos_perm_p`), y de la prosa/plots de `11.Scoring_report.Rmd` (el marcador de
esquema "ranked scale" pasa a mirar `pos_perm_p`; el panel RAW usa
`n_detected / n_cycles`).

---

## 8. Forks abiertos

1. ~~**Flip de headline (§7.3).**~~ — **HECHO** (commit `p.emp §7.3:`), adelantado
   antes de la corrida PEPC. Pasa de decisión de merge a verificación: si `p.emp`
   se apila en 0 bajo fenotipo nulo en la corrida Tier 1 PEPC, se **revierte**.
2. **Tamaño de `perm_pos_cycle_caas.tsv.gz`** en `caas_full_perms` — medir en el
   primer run real; si molesta, emitir solo `(Gene, Position, side)` detectados
   alguna vez (no requiere el observado). **Único fork que queda; se resuelve con
   la corrida PEPC, no antes.**
3. ~~`pos_perm_p` side-consciente vs quitarlo~~ — **resuelto**. Amendment V3-4a lo
   pooló a `(Gene, Position, caap_group)`; el §7.3 flip lo sacó de
   `position_scores.tsv` (queda solo en `perm_pos_pval.tsv`, diagnóstico).
4. ~~**Borrado de `null_pvalue_boot`**~~ — **HECHO** (commit `p.emp §7.4:`).
   Usuario confirmó que ningún consumidor externo lo lee.
