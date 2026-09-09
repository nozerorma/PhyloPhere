# scoring_v2 · T3-doc — Diseño del `core` pareado por lado

Rama `scoring_v2`. Entregable de diseño del tramo T3 (sin código). Fija la fórmula
**per-lado** del nuevo `asr_path_score` con la precisión necesaria para recalcular
a mano los goldens de T0. Cotejado contra
[`path_scores.py`](../subworkflows/CT_DISAMBIGUATION/local/src/convergence/path_scores.py)
tras T1 (`b1e8f87`) y contra `path_scores_golden.json` (12 escenarios).

Coordenadas 0-based. `enc(x)` = residuo codificado en el esquema activo (US → residuo
crudo; GS → etiqueta de grupo). `s ∈ {top, bottom}`.

---

## 1. Objetivo

Sustituir el producto de tres factores globales de T1

```
asr_path_score = independence · core · derived_agreement          (T1, por posición·esquema·hipótesis)
core           = 1 − (1 − core_top)(1 − core_bottom)
```

por **un `core` por `(Gene, Position, scheme, side)`**, construido como agregación
**pareada** que:

1. absorbe `independence` (entra en cada pareja, sobre su nodo de fusión concreto y
   su residuo compartido, no como factor global);
2. absorbe el acuerdo (entra por pareja: `agree(c,d) ∈ {0,1}`, no como media de
   posición);
3. contabiliza los pares conservados **en el denominador** (un diseño de `k` pares
   donde solo `j` participan en una convergencia puntúa ~`j/k` de lo que puntuaría
   con `k` participantes) — sin factor de descuento ni constante nueva;
4. no recombina lados: una posición "both" son dos filas independientes.

---

## 2. Qué cambia respecto a T1

| # (roadmap) | Problema T1 | Corrección T3 |
|---|---|---|
| 1 | `independence` global multiplica todo el score con el pool derivado completo en cada nodo LCA | `independence` por pareja: `1 − P_wc(enc(der_c) @ LCA(c,d))` — solo el residuo compartido de esa pareja, solo su nodo de fusión |
| 2 | `conservation_gate` (T1: ya latente) vivía fuera del `core` y solo penalizaba | par conservado = miembro con `score_c = 0` de `D_s`; entra en el denominador de la media |
| 3 | dos heurísticas paralelas: `derived_agreement` continuo (score) vs etiquetas categóricas de `convergence.py` (`convergence_type`) | el acuerdo entra por pareja en el score; `convergence_type` se deriva de `(agree_num_s, agree_den_s)` (§7); T4a borra `convergence.py::classify_change_and_parallelism` |
| 4 | `1 − (1−core_top)(1−core_bottom)` infla las posiciones "both" | eliminado: `side` es clave de primera clase, `"both"` → dos filas `(Gene, Position, top)` y `(Gene, Position, bottom)` |
| 5 | `models.py` arrastra payloads muertos; sin red de tests numérica por eje | fuera de alcance de T3-doc (T4a / T0) |

**Desaparecen del cálculo:** `derived_agreement` como factor global,
`_p_at_least_2` / la descomposición `p0/p1/p≥2`, `conservation_gate`,
`mrca_diversity`, `diversity_floor`, el factor `independence` de bloque, y
`1 − (1−t)(1−b)`.

**Se conserva emitido** (diagnóstico / consumo aguas abajo): `concentration_s`
(= `agree_num_s / agree_den_s`, la reencarnación por lado de `derived_agreement`),
`conserved_pair_scores` / `conserved_pair_nodes`, `pair_ancestral` /
`pair_derived_{top,bot}`.

---

## 3. Primitivas reutilizadas (sin primitiva nueva salvo noisy-OR)

- `side_path_score(..., is_changed=True, stop_at_id=…)` → `s_c^s`, aislamiento del
  segmento privado del par `c` en el lado `s`. **Sin cambios.**
- `worst_case_group_probability(dist, enc, scheme)` → `P_wc(enc @ nodo)`: masa
  exacta si el residuo está registrado en el posterior; si no, el remanente no
  registrado (`1 − Σ registrado`) como cota superior.
- `find_lca`, `path_to_root_ids`, `build_node_index`, `encode_aa`.
- **noisy-OR** `1 − ∏(1 − x_i)`: función nueva, trivial (§5.5).

`_p_at_least_2` deja de usarse en el observado (sigue en el null hasta T3c).

---

## 4. Clave de salida

`compute_asr_path_score` devuelve `{"top": {…}, "bottom": {…}}` (T3a). Cada
sub-dict es una fila `(Gene, Position, scheme, side)`. Un lado sin pares que
cambiaron devuelve `asr_path_score = 0.0`, `n = |D_s|`, `n_participating = 0`. La
fontanería de T3b decide si la fila se emite (según `change_side`).

---

## 5. Definición del `core` por lado

Para una posición, un esquema y un lado `s`:

### 5.1. `D_s` — el conjunto de pares del lado

```
P_s = pares que cambiaron en el lado s
    = { par c : enc(tip_c^s) definido y ≠ enc(focal_state_c) }

C   = pares conservados en metadata (fg y bg AMBOS con el residuo ancestral),
      definidos en esta posición                                    # conserved_pair

D_s = P_s ∪ C
n   = |D_s|
```

**`D_s` NO incluye** los pares que cambiaron **solo en el otro lado**: esos son
evidencia *a favor* de la convergencia del otro lado, no *en contra* de la de `s`.
Un par que cambió en top y no en bottom no está en `D_bottom` (salvo que sea
conservado en metadata, que es otra cosa: fg y bg iguales).

`★ Insight ─────────────────────────────────────`
- Distinción clave: un **conservado de metadata** (fg = bg = ancestral) SÍ diluye
  `core_s` — es un linaje con el fenotipo que no hizo la sustitución, evidencia de
  que la sustitución no es universalmente necesaria. Un **par que cambió solo en
  el otro lado** NO diluye `core_s` — su señal pertenece al otro lado.
- Los lados nunca se mezclan. `D_top` y `D_bottom` se construyen y puntúan
  independientes (decisión A del roadmap).
`─────────────────────────────────────────────────`

### 5.2. `s_c^s` — aislamiento del segmento privado (sin cambios respecto a T1)

```
L_s        = { LCA(mrca_a, mrca_b) : a, b ∈ P_s }        # puntos de fusión de los pares del lado
stop_c     = primer nodo de path_to_root(mrca_c) que está en L_s
s_c^s      = ∏_{k ∈ path_to_root(mrca_c), antes de stop_c} (1 − P_wc(enc(der_c^s) @ k))
```

- Segmento privado vacío (fusión de hermanos) → `s_c^s = 1.0`.
- `path_to_root` vacío (MRCA en la raíz) → `s_c^s = EMPTY_PATH_SCORE = 0.5`.
- Único cambio vs T1: `L_s` es **por lado** (LCAs entre pares del lado `s`), no el
  conjunto global sobre ambos lados. Coincide con el global salvo cuando hay un
  par de un lado cuyo MRCA induce una fusión más cercana en el camino de un par
  del otro lado; en los 12 goldens de T0 coinciden.

### 5.3. `contrib(c, d)` — fuerza de una pareja convergente

Para cada pareja **no ordenada** `{c, d}` de pares que **cambiaron en `s`**
(`c, d ∈ P_s`):

```
agree(c,d)   = 1  si  enc(der_c^s) = enc(der_d^s) ,  0  si no
contrib(c,d) = s_c^s · s_d^s · agree(c,d) · ( 1 − P_wc(enc(der_c^s) @ LCA(mrca_c, mrca_d)) )
```

Lectura: probabilidad de que `c` y `d` formen una pareja convergente **genuina**
(ambos cambiaron de verdad: `s_c·s_d`), **al mismo residuo** (`agree`), con el
**ancestro compartido limpio** (`1 − P` del residuo en su nodo de fusión).

`contrib(c,d) = 0` si `agree(c,d) = 0` (residuos distintos) — no se calcula el
resto.

### 5.4. `agree` duro, y de dónde sale la suavidad bioquímica

`agree(c,d) ∈ {0, 1}` sobre residuos **codificados**. Dos residuos distintos bajo
US dan `0`; dos residuos que co-codifican bajo un esquema GS dan `1`
automáticamente.

La "media tinta" para residuos distintos-pero-parecidos **no entra aquí**, entra
en `scoring_compute.R` §2g, que promedia `caas_row` sobre los **5 esquemas**
(US, GS4, GS3, GS2, GS1):

| par `{V, L}` | US (20 grupos) | GS4 (12) | GS3 (6) | GS2 (7) | GS1 (6) |
|--------------|------|------|------|------|------|
| `agree` | 0 | 0 | 0 ó 1 | 1 | 1 |
| `core_s` | 0 | 0 | … | > 0 | > 0 |

`CAAS_score = mean(0, 0, …, >0, >0)` → un valor intermedio = "en cuántas de las 5
granularidades bioquímicas V y L son la misma solución". Es una **distancia
bioquímica discretizada** (nota de `METODOS_CT_SCORING.md` §F.3), no una constante.
Si V y L difieren en los 5 esquemas → `mean = 0`: soluciones genuinamente
distintas, cero convergencia.

`★ Insight ─────────────────────────────────────`
- Dentro de un esquema, `agree` es un hecho crudo: ¿misma clase o no?
  "Convergencia" = "misma solución"; dos residuos distintos = dos soluciones = `0`
  natural, no forzado. Añadir un `0.5` sería **añadir una regla** y premiaría la
  labilidad de sitio (el confundente principal de la hipótesis CAAS: un sitio
  rápido donde cada linaje hace lo suyo acumularía muchos `0.5`).
- El gradiente bioquímico existe, cuantizado en 5 niveles, y se aplica al
  promediar esquemas — no dentro de uno. **Requiere que los 5 esquemas se corran
  siempre** (hoy `scoring_schemes` en §2a los corre).
`─────────────────────────────────────────────────`

### 5.5. `score_c` — la mejor evidencia de que `c` converge (noisy-OR)

```
partners(c) = { d ∈ P_s : d ≠ c , agree(c,d) = 1 }

score_c = 1 − ∏_{d ∈ partners(c)} ( 1 − contrib(c,d) )        # 0 si partners(c) = ∅
```

**noisy-OR** = probabilidad de que `c` forme una pareja convergente genuina **con
al menos uno** de sus compañeros del mismo residuo. `∏(1 − contrib)` es "ningún
compañero convergió con `c`"; `1 −` eso, "al menos uno".

| `partners(c)` | `score_c` |
|---------------|-----------|
| `∅` (conservado, o residuo único en el lado) | `0` |
| `{d}` | `contrib(c,d)` — recupera el valor pareado |
| `{d, e}`, `contrib` = 0.8, 0.8 | `1 − 0.2·0.2 = 0.96` |
| `{d, e}`, `contrib` = 0.9, 0.1 | `1 − 0.1·0.9 = 0.91` (el fuerte manda) |

Acotado en `[0,1]`, monótono (añadir compañero nunca baja `score_c`).

**Alternativa (fork, §14):** `score_c = max_{d} contrib(c,d)`. Más simple; no
premia testigos extra dentro de un grupo de residuo (`3→V` empataría con `2→V`).

### 5.6. `core_s` — media sobre el diseño

```
core_s = ( Σ_{c ∈ D_s} score_c ) / n            # n = |D_s| ;  core_s = 0 si n = 0
```

Media de `score_c` sobre **todos** los pares relevantes del lado (participantes +
conservados). Los conservados y los que cambiaron a un residuo huérfano aportan
`score_c = 0` y engordan `n`: la media se reparte entre el diseño completo.

### 5.7. Salida

```
asr_path_score[Gene, Position, scheme, side = s] = clamp01(core_s)
```

Sin `independence` de bloque, sin `derived_agreement` global, sin drag, sin
recombinación de lados.

---

## 6. Racional de cada decisión

### 6.1. Por qué media sobre `D_s` y no inclusión-exclusión / unión / `max` global

- **Comparabilidad entre hipótesis (petición del usuario).** Una hipótesis FOP con
  `k` pares "tiene esos `k` pares y hay que contabilizarlos". Si diseñaste 4
  contrastes esperando convergencia y solo 2 la produjeron, tu `core` debe ser
  ~½ del de un 4-de-4. El denominador `n = |D_s|` lo implementa. Unión / `max`
  global ignoran los pares que no participan.
- **Penalización lineal, no combinatoria.** Promediar `score_c` **por par** (donde
  `score_c` ya es el noisy-OR de *sus* compañeros) mete **un** cero por par muerto.
  Promediar sobre las `C(n,2)` parejas metería `n−1` ceros por par muerto →
  castigo cuadrático (un conservado entre 2 cambiados → `1/6` del score). Lineal:
  un par muerto entre `k` que convergen → factor `≈ k/(k+1)`.
- **Multi-modal.** `AAVV` (dos convergencias, A×2 y V×2): con `agree` por pareja,
  `{P1,P2}` y `{P3,P4}` tienen `contrib > 0`, las cruzadas `0`. Cada par encuentra
  su compañero; `core_s` promedia las dos convergencias. Con pluralidad se perdía
  una entera.

### 6.2. Por qué noisy-OR dentro de `score_c`

`c` puede tener varios compañeros del mismo residuo. `mean` de sus `contrib`
castiga a `c` por tener un compañero débil aunque tenga uno fuerte; `sum` no está
acotado; `max` ignora que varias parejas mediocres juntas convencen más. noisy-OR
= "≥1 pareja genuina, y más compañeros suben la confianza", que es la semántica de
replicación. Supuesto de independencia entre `contrib(c,d)` y `contrib(c,e)`:
comparten `s_c` y, si `d` y `e` se fusionan con `c` en el mismo nodo, comparten
el factor de independencia → doble conteo suave (pesimista) en topología en
estrella. Aceptable.

### 6.3. Por qué los conservados en el denominador y no un factor de descuento

Un factor `∏(1 − λ_c · cons_k)` necesitaba una constante nueva `λ_c` (justo lo que
la decisión B del roadmap prohíbe) y su magnitud dependía de combinatoria no
controlada (un conservado se comía el 70% del score en el ejemplo de trabajo).
`score_c = 0` en el denominador es "el par conservado sencillamente como un 0"
(decisión C, lectura de "media ponderada"), sin parámetros, con castigo lineal y
proporcional al peso del conservado en el diseño.

### 6.4. Por qué `independence` por pareja concreta y no global

En el nodo 1 (fusión de `c` y `d`, **ambos → V**) la pregunta correcta es "¿estaba
V ya en el nodo 1?" = `P(V @ 1)`, no `P({todo el pool} @ 1)`. Que un tercer par
`e → L` se fusione más arriba (nodo 0) es irrelevante para si `c` y `d`
convergieron independientemente. El factor global de T1 aplicaba el pool completo
en cada nodo LCA y arrastraba nodos ajenos a cada relación. Por pareja: `indep`
aparece **una vez** por `{c,d}`, keyed a su nodo y su residuo → sin doble conteo,
y la lectura de probabilidad conjunta se mantiene (segmentos privados y nodos de
fusión son regiones disjuntas del árbol).

### 6.5. Por qué `agree` duro

Ver §5.4: "misma solución" es binario dentro de un esquema; el gradiente
bioquímico lo da el promedio de los 5 esquemas en §2g; un `0.5` ad-hoc premiaría
labilidad de sitio.

---

## 7. `convergence_type` desde `(agree_num_s, agree_den_s)` (spec para T4a)

```
agree_den_s = |P_s|                                             # nº de líneas que cambiaron en s
agree_num_s = max_r |{ c ∈ P_s : enc(der_c^s) = r }|            # tamaño del mayor grupo de residuo
concentration_s = agree_num_s / agree_den_s                     # emitido como `derived_agreement` (diagnóstico)
```

| condición | `convergence_type` (lado `s`) |
|-----------|-------------------------------|
| `agree_den_s ≥ 2` y `agree_num_s ≥ 2` | `convergent` |
| `agree_den_s ≥ 2` y `agree_num_s < 2` (todas a residuos distintos) | `divergent` |
| `agree_den_s = 1` | `single` |
| `agree_den_s = 0` | `no_change` |

`AAVV` → `den = 4`, `num = 2` (max de {A:2, V:2}) → `convergent` (hay una
convergencia de ≥2 vías). La etiqueta es el **conteo físico** de la decisión B
(sin umbral de probabilidad, sin parámetro), y **no entra en el score continuo**.

---

## 8. Comportamiento con `n` y multi-modal

- **`< 2` pares con residuo compartido en el lado** → todos los `score_c = 0` →
  `core_s = 0`. El gate "≥2" del roadmap emerge, no se codifica aparte.
- **`k` pares, todos al mismo residuo, limpios** (`contrib ≈ c`): `score_c =
  1 − (1−c)^{k−1}` (noisy-OR de `k−1` compañeros) → sube con `k`. `core_s ≈` ese
  valor. Un 4-de-4 supera a un 2-de-2. Con `max` empatarían (fork §14).
- **`k` participan de un diseño de `n > k`** (resto conservado o huérfano):
  `core_s ≈ [convergencia] · k/n`. Comparabilidad lineal.
- **multi-modal `A×j / V×l`**: cada grupo de residuo se auto-empareja; `core_s`
  promedia `score_c` de los `j + l` participantes (los de cada grupo con su
  noisy-OR interno) sobre `n`. Ninguna convergencia se pierde.

---

## 9. Casos borde

| caso | resultado |
|------|-----------|
| `P_s = ∅` y `C = ∅` | `n = 0` → `core_s = 0` |
| `P_s = {c}`, sin compañero de residuo | `score_c = 0` → `core_s = 0` |
| `P_s = {c, d}`, `enc(der_c) ≠ enc(der_d)` | ambos `score = 0` → `core_s = 0` (acantilado deliberado: dos soluciones distintas ≠ convergencia; gradiente vía §2g) |
| MRCA en la raíz | `s_c^s = 0.5` |
| fusión de hermanos | `s_c^s = 1.0` |
| `LCA(c,d)` = raíz o nodo sin posterior | `P_wc = 0` → factor `1.0` (o `EMPTY` si procede) |
| par conservado de metadata | en `D_s`, `score = 0`, cuenta en `n` |
| par que cambió solo en el otro lado | **no** en `D_s` |

---

## 10. Recálculo de los 12 goldens de T0

`s_c^s` e `indep` se leen de los valores actuales del golden (los walks no
cambian). Todos los escenarios tienen los pares cambiados en `top` salvo donde se
indica.

| escenario | `D_top` (`n`) | `contrib` / `score_c` | **`core_top` nuevo** | `asr` T1 | Δ |
|-----------|---------------|-----------------------|----------------------|----------|---|
| `single_changed_pair` | {P1} (1) | `score_P1 = 0` (sin compañero) | **0** | 0 | = |
| `two_same_side_converge` | {P1,P2} (2) | `contrib(P1,P2)=0.9025·0.9025·1·0.95=0.7738`; `score` ambos `0.7738` | **0.7738** | 0.7738 | = |
| `opposite_sides_only` | top {P1} (1), bot {P2} (1) | `score = 0` en cada lado | **top 0 / bot 0** | 0 | = (2 filas) |
| `pair_changes_both_sides` | top {P1,P2} (2), bot {P1} (1) | top: `contrib=0.7738`; bot: `score_P1=0` | **top 0.7738 / bot 0** | 0.7738 | schema; `.pool_bottom` deja de recibir 0.7738 (H5) |
| `with_conserved_pair` | {P1,P2,**K3**} (**3**) | `score_P1=score_P2=0.7738`, `score_K3=0` | **`1.5476/3 = 0.5159`** | 0.7738 | **↓** conservado en el denominador (2 de 3) |
| `contaminated_hop1` | {P1,P2} (2) | `contrib=0.285·0.9025·1·0.95=0.2444` | **0.2444** | 0.2444 | = |
| `n_gt_2_mixed_residues` | {P1,P2,P3(→V), **P4(→T)**} (**4**) | `contrib(P1,P2)=0.857`, `(P1,P3)=0.857`, `(P2,P3)=0.95`; `score_{P1,P2,P3}≈{0.980, 0.993, 0.993}`, `score_P4=0` | **`2.965/4 = 0.741`** | 0.6075 | **↑** convergencia V fuerte; el par→T solo mete un 0 (no arrastra una media a 0.75) |
| `gs3_coencoded_agreement` | {P1,P2} (2) | V,I co-encode bajo GS3 → `agree=1` → `contrib=0.7738` | **0.7738** | 0.7738 | = |
| `mrca_at_root` | {P1,P2} (2) | `s_c=0.5` cada; `contrib=0.5·0.5·1·0.95=0.2375` | **0.2375** | 0.2375 | = |
| `sibling_merge` | {P1,P2} (2) | `s_c=1.0`; `contrib=1·1·1·0.95=0.95` | **0.95** | 0.95 | = |
| `no_changed_pairs` | ∅ (0) | — | **0** | 0 | = |
| `soft_posteriors_midrange` | {P1,P2} (2) | `contrib=0.49·0.49·1·0.70=0.1681` | **0.1681** | 0.1681 | = |

**10 de 12 idénticos.** Se mueven `with_conserved_pair` (↓: el conservado ahora
cuenta) y `n_gt_2_mixed_residues` (↑: la convergencia mayoritaria deja de ser
castigada por una `derived_agreement` de posición). **Cero constantes nuevas.**

Razón de que los `n = 2` no se muevan: `contrib(P1,P2) = s_A·s_B·1·(1 − P(V@LCA))`
es exactamente `core_T1 · independence_T1 · da_T1` cuando hay un solo par de
pares, un solo nodo LCA y residuo unánime — que es la estructura de esos goldens.

### Golden nuevo para T3a — `both_sides_two_rows`

Diseño: `pair 1: top→V, bottom→L`; `pair 2: top→V`; `pair 3: bottom→L`. Árbol con
3 MRCAs en subclados separados.

- `D_top = {P1, P2}` (`n=2`), ambos → V → `core_top = contrib(P1,P2) > 0`.
- `D_bottom = {P1, P3}` (`n=2`), ambos → L → `core_bottom = contrib(P1,P3) > 0`.

Dos filas con scores **positivos e independientes**. Congela el fin de
`1 − (1−t)(1−b)` (T1 habría dado un único `core` combinado). `gen_golden.py`
fija los números exactos en T3a.

---

## 11. Deltas esperados en `position_scores.tsv` / `gene_scores.tsv`

1. **Posiciones "both" → 2 filas.** `n_distinct(Gene, Position)` no cambia;
   `nrow(position_scores)` sube. Downstream ya agrupa por `(Gene, Position, side)`
   desde T2a.
2. **Posiciones con par conservado bajan.** Ahora cuenta en el denominador. Baja
   proporcional al peso del conservado en el diseño (lineal).
3. **Posiciones con un residuo minoritario suben un poco.** El minoritario mete un
   `0` en vez de arrastrar una media de posición (`derived_agreement`) a la baja.
4. **`.pool_bottom` / `.pool_top` cambian de composición (hazard H5).** Ver §12.

La validación Tier 1 (PEPC en Marvin2) al cerrar T3d comprueba que estos 4 deltas
aparecen **en la dirección prevista** y ningún otro.

---

## 12. Pooling — nivel lado (roadmap bullet 7, H5)

`pos_scores` tiene una fila por `(Gene, Position, side)` con `CAAS_score` = media
sobre esquemas de `core_s`.

- `.pool_top` ← filas `side == "top"`; `.pool_bottom` ← `side == "bottom"`;
  `.pool_all` ← todas las filas-lado (una posición "both" aporta **2 entradas**).
- `size_adj_max` = `F(max)^n`: `n` = nº de filas-lado del gen. Premia genes con
  muchos eventos de convergencia direccional.
- Modo acumulación "top/bottom/all": filtro directo `side %in% c("top")` /
  `c("bottom")` / todas — sin `c(dir, "both")` (simplificación prevista T3d §4a).

**Hazard H5.** En T1 una posición "both" fuerte-en-conjunto pero débil en un lado
entraba en `.pool_<ese lado>` con su score combinado. Con nivel lado entra con su
`core_<lado>` real. Los rankings genome-wide se mueven: es la corrección buscada
(evidencia de un solo lado o de un lado débil deja de hacerse pasar por
bidireccional). Changelog, no regresión.

---

## 13. Columnas emitidas por `(Gene, Position, scheme, side)` (T3a/b)

| columna | contenido |
|---------|-----------|
| `asr_path_score` | `core_s` |
| `core` | alias de `asr_path_score` (compat) |
| `n_pairs_side` | `n = |D_s|` |
| `n_participating` | `|P_s|` |
| `n_conserved` | `|C|` |
| `derived_agreement` | `concentration_s = agree_num_s / agree_den_s` (diagnóstico; num/den para T4a) |
| `agree_num`, `agree_den` | enteros (T4a `convergence_type`) |
| `convergence_type` | derivado per §7 (T4a) |
| `side` | `top` / `bottom` (clave, ya en T2a) |
| `pair_scores` / `pair_top_scores` / `pair_bottom_scores` | `{pid: s_c^s}` del lado |
| `pair_partner_scores` | `{pid: score_c}` (diagnóstico) |
| `conserved_pair_scores` / `conserved_pair_nodes` | sin cambios (latente; ya no multiplican, solo informan) |
| `pair_ancestral` / `pair_derived_{top,bot}` | sin cambios (pooler FOP) |

Se eliminan del dict de retorno: `replication`, `independence` (de bloque),
`mrca_diversity`, `conservation_gate`, el `core` escalar combinado de lados.

---

## 14. Forks abiertos — a cerrar antes de T3a

| # | fork | recomendación |
|---|------|---------------|
| 1 | `score_c` = **noisy-OR** (§5.5) vs `max` | noisy-OR (premia testigos extra dentro de un grupo de residuo) |
| 2 | `both_sides_two_rows` golden — diseño exacto del escenario (§10) | confirmar árbol y residuos antes de regenerar `path_scores_golden.json` |
| 3 | Pooling nivel lado (§12) — confirmar antes de T3d | nivel lado |

Confirmados en la discusión de diseño (no relitigar): `independence` por pareja
(§6.4); `agree` duro con gradiente vía §2g (§5.4); `D_s` = participantes ∪
conservados-metadata, sin los que cambiaron solo en el otro lado (§5.1); lados
siempre separados; sin factor de drag ni constante nueva (§6.3).

Cerrados esos 3 puntos, T3a reescribe la sección de agregación de
`compute_asr_path_score` (retorno `{"top": …, "bottom": …}`) y `gen_golden.py`
regenera los 12 + 1 goldens; el diff debe coincidir con la columna "nuevo" de §10.

---

## Apéndice — walkthrough completo

Un ejemplo que ejerce cada pieza, calculado nodo a nodo.

### A.1. Árbol y entradas

```
                          root (0)
                    ┌────────┴────────┐
                  (1)                 (2)
             ┌─────┴─────┐       ┌─────┴─────┐
           (10)        (11)    (20)        (21)
           ┌─┴─┐       ┌─┴─┐   ┌─┴─┐       ┌─┴─┐
          3   4       5   6   7   8       9   12
        MRCA P1     MRCA P2  MRCA P3     MRCA P4
```

`path_to_root`: P1(3)=`[10,1,0]`, P2(5)=`[11,1,0]`, P3(7)=`[20,2,0]`, P4(9)=`[21,2,0]`.

`focal_state = A` para los 4. Esquema **US**.

| par | `top_tip` | `bottom_tip` | qué es |
|-----|-----------|--------------|--------|
| P1 | V | A | cambió solo en top, a V |
| P2 | V | V | cambió en top y bottom, a V |
| P3 | L | A | cambió solo en top, a **L** (discrepa) |
| P4 | A | A | **conservado** (metadata) |

Posteriores del sitio:

| nodo | posterior |
|------|-----------|
| 0 | `{A:.90, V:.05, L:.05}` |
| 1 | `{A:.80, V:.15, L:.05}` |
| 2 | `{A:.90, V:.05, L:.05}` |
| 10 | `{A:.70, V:.30}` |
| 11 | `{A:.95, V:.05}` |
| 20 | `{A:.90, L:.10}` |
| 21 | `{A:.98, V:.02}` |

### A.2. Lado TOP

**`D_top`** = `P_top ∪ C` = `{P1, P2, P3}` ∪ `{P4}` = **`{P1, P2, P3, P4}`**, `n = 4`.

**`L_top`** = LCAs entre `{3, 5, 7}`: `LCA(3,5)=1`, `LCA(3,7)=0`, `LCA(5,7)=0` →
**`{1, 0}`**.

**Walks privados `s_c^top`:**

| par | camino | `stop` (∈ `L_top`) | nodos caminados | factores | `s_c^top` |
|-----|--------|--------------------|-----------------|----------|-----------|
| P1 | `[10,1,0]` | 1 | `[10]` | `1−P_wc(V@10)=1−.30=.70` | **0.70** |
| P2 | `[11,1,0]` | 1 | `[11]` | `1−.05=.95` | **0.95** |
| P3 | `[20,2,0]` | 0 | `[20,2]` | `(1−P_wc(L@20)=.90)·(1−P_wc(L@2)=.95)` | **0.855** |

**Acuerdo:** P1–P2 (V,V)→1 · P1–P3 (V,L)→0 · P2–P3 (V,L)→0.

**`contrib`:**

```
contrib(P1,P2) = 0.70 · 0.95 · 1 · (1 − P_wc(V @ LCA(3,5)=1))
              = 0.70 · 0.95 · 1 · (1 − 0.15)
              = 0.665 · 0.85 = 0.565
contrib(P1,P3) = 0        (agree 0)
contrib(P2,P3) = 0
```

**`score_c`:**

```
score_P1 = 1 − (1 − contrib(P1,P2)) = 0.565            partners(P1) = {P2}
score_P2 = 1 − (1 − contrib(P1,P2)) = 0.565            partners(P2) = {P1}
score_P3 = 0                                           partners(P3) = ∅  (residuo L huérfano)
score_P4 = 0                                           partners(P4) = ∅  (conservado)
```

**`core_top`:**

```
core_top = (0.565 + 0.565 + 0 + 0) / 4 = 1.130 / 4 = 0.2825
```

### A.3. Lado BOTTOM

**`D_bottom`** = `{P2}` ∪ `{P4}` = `{P2, P4}`, `n = 2`.
`partners(P2) = ∅` (ningún otro par cambió en bottom) → `score_P2 = 0`.
`score_P4 = 0`.

```
core_bottom = (0 + 0) / 2 = 0
```

### A.4. Resultado y comparación con T1

| | T1 (una fila) | T3 |
|---|---------------|-----|
| filas | 1 | 2 |
| score | **0.451** | top **0.2825**, bottom **0** |

Descomposición de la diferencia en top (`0.565` de convergencia P1–P2, diluido):

```
convergencia P1–P2 pura:                      0.565
· (2 de 4 pares del diseño participan):        · 0.5    →  0.2825
   ├─ P3 cambió, pero a L (huérfano)  → score 0, en el denominador
   └─ P4 conservado                   → score 0, en el denominador
```

T1 daba `0.451` porque: (a) ignoraba P4 por completo; (b) penalizaba P3 con una
`derived_agreement = 2/3` de posición en vez de un `0` en el denominador; (c) su
`independence` global `0.72` metía el nodo 0 y el pool `{V,L}` en la relación
P1–P2, que aquí solo debería ver `P(V @ nodo 1) = 0.15`.

### A.5. Variantes (para intuición)

- **P3 también → V** (`P_top = {P1,P2,P3}` todos V, P4 conservado, `n = 4`):
  `contrib(P1,P2)=0.565`, `contrib(P1,P3)` y `contrib(P2,P3)` con sus LCAs;
  `score_P1 = 1−(1−contrib(P1,P2))(1−contrib(P1,P3))` (noisy-OR de 2 compañeros) →
  sube. `core_top = (score_P1+score_P2+score_P3+0)/4`. Más testigos → cada
  `score_c` sube; el conservado sigue diluyendo por `1/4`.
- **`AAVV`** (P1,P2→A, P3,P4→V, sin conservados, `n = 4`): `{P1,P2}` y `{P3,P4}`
  con `contrib > 0`, cruzadas `0`. `score` de los 4 = su `contrib` intra-grupo.
  `core_top = (c_AA + c_AA + c_VV + c_VV)/4` — las **dos** convergencias contadas.
- **Comparabilidad**: hipótesis X = `{P1,P2→V}` (`n=2`) da `core = contrib`;
  hipótesis Y = `{P1,P2→V, K3,K4 conservados}` (`n=4`) da `core = contrib/2`.
  Mitad de contrastes efectivos → mitad de score.
