# scoring_v2 · T3-doc — Diseño del `core` pareado por lado

Rama `scoring_v2`. Entregable de diseño del tramo T3 (sin código). Fija la fórmula
**per-lado** del nuevo `asr_path_score` con la precisión necesaria para recalcular
a mano los goldens de T0. Cotejado contra
[`path_scores.py`](../subworkflows/CT_DISAMBIGUATION/local/src/convergence/path_scores.py)
tras T1 (`b1e8f87`) y contra `path_scores_golden.json` (12 escenarios).

Coordenadas de alineamiento 0-based. `m_s` = nº de pares que cambiaron en el lado
`s`. `enc(·)` = residuo codificado en el esquema activo (US → residuo crudo; GS →
etiqueta de grupo).

---

## 1. Qué cambia y por qué

El `core` actual (T1) se construye así, por `(Gene, Position, scheme, hypothesis)`:

```
core_top    = P(≥2)  sobre  { s_c : c cambió en top }          # _p_at_least_2, incl-excl
core_bottom = P(≥2)  sobre  { s_c : c cambió en bottom }
core        = 1 − (1 − core_top)(1 − core_bottom)              # unión de lados
asr         = independence · core · derived_agreement          # independence y da GLOBALES
```

con `independence = ∏_{L ∈ LCA} (1 − P_wc(pool_derivado @ L))` un único factor global
y `derived_agreement` la media, a nivel posición, de la concentración de residuo por
lado.

Problemas (numeración del roadmap):

| # | Problema | Corrección en T3 |
|---|----------|------------------|
| 1 | Doble ponderación del paralelismo (`mrca_diversity` ya se fue en T1; queda la redundancia conceptual entre `independence` global y limpieza de camino privado) | `independence` pasa a **por lado**, sobre los LCA de los pares de ese lado |
| 2 | `conservation_gate` solo penalizaba y vivía fuera del `core` | par conservado entra como **drag multiplicativo** dentro de la construcción de `core_s` (§5.6) |
| 3 | Dos heurísticas paralelas (`derived_agreement` continuo vs etiquetas categóricas de `convergence.py`) | el acuerdo entra **por par** dentro de `g_c` (§5.2); `convergence_type` se deriva del numerador/denominador de acuerdo (§5.3), T4a borra la clasificación categórica |
| 4 | `1 − (1−core_top)(1−core_bottom)` infla las posiciones que cambian en ambos lados | **desaparece**: `side` es clave de primera clase, una posición "both" son dos filas independientes `(Gene, Position, top)` y `(Gene, Position, bottom)` |
| 5 | (models.py bloat, tests) | fuera de alcance de T3-doc |

Resultado: **un `core_s` por `(Gene, Position, scheme, side)`**, que absorbe
`independence` y el acuerdo, con los pares conservados empujando a la baja.
`asr_path_score[·, ·, ·, s] = core_s`. No hay recombinación de lados.

---

## 2. Primitivas reutilizadas (sin primitiva nueva)

Todo lo que sigue se calcula con funciones que ya existen en `path_scores.py`:

- `side_path_score(..., is_changed=True, stop_at_id=LCA)` → `s_c^s`, el score de
  aislamiento del segmento privado del par `c` en el lado `s` (∏ de `1 − P_wc(derived)`
  sobre los nodos privados). **Sin cambios.**
- `_conserved_side_score(...)` → `cons_k`, conservación-a-raíz del par conservado `k`
  (media de `P(ancestral)` sobre el camino MRCA→raíz). Ya se emite en
  `conserved_pair_scores`. **Sin cambios.**
- `worst_case_any_group_probability(dist, pool, scheme)` → `P_wc(pool @ nodo)`.
- `find_lca`, `path_to_root_ids`, `build_node_index`.
- `_p_at_least_2(lista)` → inclusión-exclusión exacta `P(≥2 éxitos)` sobre
  Bernoullis independientes. **Sin cambios.**
- `encode_aa` / `modal_encoded`.

La reescritura de T3a es de **agregación**, no de recorridos.

---

## 3. Clave de salida

`compute_asr_path_score` deja de devolver un dict plano y devuelve
`{"top": {...}, "bottom": {...}}` (T3a). Cada sub-dict corresponde a una fila
`(Gene, Position, scheme, side)`. Un lado sin pares cambiados devuelve el sub-dict
con `asr_path_score = 0.0` y `n_participating = 0` (no se omite: la fontanería de
dos filas de T3b decide si se emite o no según `change_side`).

---

## 4. Definición: pares participantes y conservados

Para una posición, esquema y lado `s ∈ {top, bottom}`:

- **participantes** `P_s` = pares `c` con `enc(tip_c^s)` definido y `≠ enc(focal_state_c)`
  (el tip de ese lado divergió del modal en el MRCA del par). Es la clasificación
  que ya hace la Fase 1 de `compute_asr_path_score`, restringida al lado.
- **conservados** `C` = pares en `conserved_pair` (metadata), comunes a ambos lados.
  Un par conservado nunca está en `P_s`.

`m_s = |P_s|`.

---

## 5. Fórmula per-lado

### 5.1. Éxito por par participante `g_c`

Para cada `c ∈ P_s`:

```
plur_s   = argmax_r |{ d ∈ P_s : enc(tip_d^s) = r }|          # residuo codificado de pluralidad
                                                              # empate → el de mayor s_d (determinista)
agree_c  = 1  si  enc(tip_c^s) = plur_s   ,   0  en otro caso
g_c      = s_c^s · agree_c                                     # ∈ [0, 1]
```

`s_c^s` es el score de aislamiento privado ya calculado (`side_path_score`,
`stop_at_id` = LCA más cercano del par en ese lado).

`agree_c` es **duro (0/1) sobre residuos codificados**, no una similitud bioquímica
continua: el esquema GS ya absorbe la bioquímica en `enc(·)`, así que dos residuos
que co-codifican (p.ej. V/I bajo GS3) dan `agree = 1` automáticamente y bajo US
dan `agree = 0`. Esto es coherente con la nota del docstring actual ("agreement is
tested on the scheme-encoded residues").

Consecuencia deliberada: **`m_s = 2` con residuos codificados distintos → `g` = [s₁, 0]
→ `core_s = 0`**. Dos cambios a cosas distintas, uno cada uno, no son convergencia.
El acantilado solo ocurre cuando genuinamente no hay dos pares que compartan residuo.

### 5.2. Numerador / denominador de acuerdo (spec para T4a)

```
agree_den_s = m_s
agree_num_s = |{ c ∈ P_s : enc(tip_c^s) = plur_s }|
concentration_s = agree_num_s / agree_den_s        si m_s ≥ 1 ,  indefinido si m_s = 0
```

`concentration_s` es la reencarnación **por lado** del antiguo `derived_agreement`
(que era la media de `concentration` sobre lados cualificados). Se emite como
columna diagnóstica. T4a deriva de `(agree_num_s, agree_den_s)`:

| condición | `convergence_type` (lado `s`) |
|-----------|-------------------------------|
| `agree_den_s ≥ 2` y `agree_num_s ≥ 2` | `convergent` |
| `agree_den_s ≥ 2` y `agree_num_s < 2` | `divergent` |
| `agree_den_s = 1` | `single` |
| `agree_den_s = 0` | `no_change` |

(La etiqueta de "lado cambiado" del roadmap decisión B — conteo físico de pares con
tip ≠ ancestral ≥ 2 — es exactamente `agree_den_s ≥ 2`.)

### 5.3. Probabilidad de replicación `p_ge2`

Inclusión-exclusión exacta sobre los `g_c` (Bernoullis independientes), **solo
participantes**, conservados excluidos:

```
p0    = ∏_{c ∈ P_s} (1 − g_c)                                  # P(0 cambios genuinos y de acuerdo)
p1    = Σ_{c ∈ P_s} [ g_c · ∏_{d ≠ c} (1 − g_d) ]              # P(exactamente 1)
p_ge2 = max(0, 1 − p0 − p1)                                    # P(≥ 2)
```

- `m_s < 2` → `p_ge2 = 0` (no puede haber ≥2). Cubre `single_changed_pair`,
  `opposite_sides_only`, y cualquier lado con un único par.
- Es literalmente `_p_at_least_2([g_c for c in P_s])`. Sin primitiva nueva.
- Un par conservado NO entra aquí (entraría como `g = 0`, que en `_p_at_least_2`
  aporta factor `(1−0) = 1` a cada término: **efecto nulo**, no drag). Los
  conservados se tratan en §5.4.

### 5.4. Independencia por lado `independence_s`

```
L_s            = { LCA(a, b) : a, b ∈ MRCAs de P_s }           # puntos de fusión pareados, cada nodo una vez
pool_derivado_s = { enc(tip_c^s) : c ∈ P_s }
independence_s = ∏_{L ∈ L_s} (1 − P_wc(pool_derivado_s @ L))
```

Es el `independence` actual **restringido a los pares del lado** (hoy el pool y los
LCA mezclan ambos lados). Se mantiene como **un único factor por lado**, no como
producto por par.

> **Desviación consciente de la decisión C.** La decisión C dice "`independence`
> entra como producto sobre el/los LCA **de cada pareja concreta**". Tomado al pie
> de la letra, cada `g_c` llevaría su propio `indep_c`, y en el término `p2` del
> par `{A, B}` el nodo `LCA(A, B)` aparecería en `g_A` y en `g_B` → **se elevaría
> al cuadrado**. El nodo de fusión de `{A, B}` debe contar una vez por ese par, no
> una por cada miembro. Mantener `independence_s` como factor único sobre `L_s`
> (cada nodo una vez) evita el doble conteo y conserva la lectura de probabilidad
> conjunta: `p_ge2` vive en los segmentos privados (nodos disjuntos de `L_s`) e
> `independence_s` en los de fusión, condicionalmente independientes → el producto
> es `P(≥2 privado-limpio ∩ fusión-limpia)`. **Requiere visto bueno del usuario
> para cerrar la decisión C con esta forma** (§12).

### 5.5. Drag de pares conservados `conserved_drag_s`

```
conserved_drag_s = ∏_{k ∈ C} (1 − λ_c · cons_k)
```

`cons_k ∈ [0, 1]` = conservación-a-raíz del par `k` (`conserved_pair_scores`, ya
calculada). `λ_c` = constante nueva en `[0, 1]`.

Racional: un par conservado es un contraste fg/bg que **no** adquirió el residuo
esperado pese al fenotipo. Es evidencia en contra de que la convergencia en los
pares que sí cambiaron esté guiada por el fenotipo, luego descuenta.

Un par conservado no puede entrar como Bernoulli en `_p_at_least_2` (§5.3): p=0
tiene efecto nulo, no drag. Por eso entra multiplicativo, análogo al antiguo
`conservation_gate` pero (i) por lado, (ii) construido con la primitiva pareada
(`cons_k`), (iii) dentro de `core_s`, no bolted-on fuera.

**`λ_c` — decisión abierta (§12).** Candidatos:

| `λ_c` | efecto de 1 par conservado con `cons_k = 0.9`, `m_s = 2` | comentario |
|-------|--------------------------------------------------------|------------|
| `1.0` | `drag = 0.10` | agresivo; un conservado casi anula el lado |
| `0.5` | `drag = 0.55` | fijo, simple |
| `m_s / (m_s + 1)` | `m_s=2 → 2/3 → drag = 0.40` ; `m_s=5 → 5/6 → drag = 0.25` | escala con el nº de testigos que sí cambiaron; recomendado |

Recomendación: **`λ_c = m_s / (m_s + 1)`** (un conservado entre muchos cambiados
pesa poco; entre 2, pesa ~⅔). Es la lectura de "como 0 en una media ponderada" de
la decisión C trasladada a la forma multiplicativa.

### 5.6. `core_s` y `asr_path_score`

```
core_s = p_ge2 · independence_s · conserved_drag_s
asr_path_score[Gene, Position, scheme, side=s] = clamp01(core_s)
```

Sin `derived_agreement` global (absorbido en `g_c`), sin recombinación de lados,
sin `mrca_diversity`, sin `conservation_gate` externo.

---

## 6. Comportamiento con `n`

`p_ge2` sobre `_p_at_least_2` es **monótona creciente en el nº de participantes de
buena calidad**: cada par adicional con `g_c` alto sube `P(≥2)`. Esto es
"replicación como existencia de ≥2 testigos", no "calidad media" — es lo correcto
para un score de replicación y coincide con la semántica actual.

No infla con `n` de baja calidad: `_p_at_least_2` es lineal en cada `g_c`, no
cuadrática en `C(n,2)`. Un par que discrepa de la pluralidad entra con `g_c = 0` y
no cuenta.

Se descartaron:

- **unión sobre contribuciones pareadas** `1 − ∏_{A<B}(1 − contrib(A,B))`: el
  exponente crece `C(n,2)`, satura demasiado rápido, sobre-pondera `n`.
- **media sobre contribuciones pareadas**: pierde la monotonía en testigos (3
  pares al 90% dan 0.81 en vez de subir), y no es una probabilidad de conteo, así
  que la descomposición `p0/p1/p≥2` (necesaria para T4a) no sale de forma natural.

---

## 7. Casos borde

| caso | resultado |
|------|-----------|
| `m_s = 0` | `asr[·, s] = 0`, `n_participating = 0` |
| `m_s = 1` | `p_ge2 = 0` → `asr[·, s] = 0` |
| `m_s ≥ 2`, residuos codificados todos distintos | `plur_s` gana con 1; todos menos uno `agree = 0` → como mucho un `g_c > 0` → `p_ge2 = 0` |
| MRCA en la raíz (camino privado vacío) | `s_c^s = EMPTY_PATH_SCORE = 0.5` (sin cambios) |
| fusión de hermanos (segmento privado vacío) | `s_c^s = 1.0` (sin cambios) |
| `L_s = ∅` (p.ej. `m_s = 1`, ya cae por `p_ge2 = 0`) | `independence_s = 1` (producto vacío) |
| `C = ∅` | `conserved_drag_s = 1` |

---

## 8. Recálculo de los 12 goldens de T0

`independence_s` y `s_c^s` se leen de los valores actuales del golden (los
recorridos no cambian). `da` global se elimina; entra `agree_c` por par.
Escenarios todos con pares en `top` salvo donde se indica.

| escenario | `P_top` `g_c` | `p_ge2` | `ind_top` | `drag` | **`core_top` (nuevo)** | `asr` T1 | Δ |
|-----------|---------------|---------|-----------|--------|------------------------|----------|---|
| `single_changed_pair` | [0.857·1] | 0 (m<2) | — | 1 | **0** | 0 | = |
| `two_same_side_converge` | [0.902, 0.902] | 0.8145 | 0.95 | 1 | **0.7738** | 0.7738 | = |
| `opposite_sides_only` | top [0.902], bottom [0.902] | 0 / 0 | — | 1 | **top 0, bottom 0** | 0 | = (ahora 2 filas) |
| `pair_changes_both_sides` | top [0.902, 0.902], bottom [0.902] | 0.8145 / 0 | 0.95 / — | 1 | **top 0.7738, bottom 0** | 0.7738 | schema: 2 filas; `.pool_bottom` deja de recibir 0.7738 (H5) |
| `with_conserved_pair` | [0.902, 0.902], `C={3}` `cons₃≈?` | 0.8145 | 0.95 | `1 − λ_c·cons₃` | **0.7738 · drag** | 0.7738 | ↓ por el conservado (magnitud según `λ_c`, `cons₃`) |
| `contaminated_hop1` | [0.285, 0.902] | 0.2572 | 0.95 | 1 | **0.2444** | 0.2444 | = |
| `n_gt_2_mixed_residues` | V,V,V,T → `g` = [0.902, 1.0, 1.0, 0] | 1.0 | 0.81 | 1 | **0.81** | 0.6075 | ↑ el par discrepante (→T) se excluye en vez de arrastrar una media a 0.75 |
| `gs3_coencoded_agreement` | V,I co-codifican bajo GS3 → `agree` = [1, 1] → `g` = [0.902, 0.902] | 0.8145 | 0.95 | 1 | **0.7738** | 0.7738 | = |
| `mrca_at_root` | [0.5, 0.5] | 0.25 | 0.95 | 1 | **0.2375** | 0.2375 | = |
| `sibling_merge` | [1.0, 1.0] | 1.0 | 0.95 | 1 | **0.95** | 0.95 | = |
| `no_changed_pairs` | ∅ | 0 | — | 1 | **0** | 0 | = |
| `soft_posteriors_midrange` | [0.49, 0.49] | 0.2401 | 0.70 | 1 | **0.1681** | 0.1681 | = |

**Se mueven solo 2 valores** (`n_gt_2_mixed_residues` sube, `with_conserved_pair`
baja) más el cambio de cardinalidad de las posiciones "both". Todos los demás
goldens son **idénticos** porque, para `n = 2` con residuo unánime y camino limpio,
`asr_T1 = (s_A·s_B) · independence · 1 = p_ge2 · independence_s · 1 = core_s`.

Nota: `_p_at_least_2([p, q]) = p·q` exacto para `n = 2` (≥2 de 2 = ambos), así que
el `core` de T1 para `n = 2` ya era `s_A·s_B`.

### Golden nuevo que T3a debe añadir

`gen_golden.py` gana un escenario `both_sides_two_rows`: una posición con pares
`{1: top→V, bottom→V}` y `{2: top→V}` y `{3: bottom→V}` → fila `top` con
`P_top = {1, 2}` (`core_top > 0`) y fila `bottom` con `P_bottom = {1, 3}`
(`core_bottom > 0`), dos filas con scores independientes. Congela el fin de
`1 − (1−t)(1−b)`.

---

## 9. Deltas esperados a nivel `position_scores.tsv` / `gene_scores.tsv`

1. **Posiciones "both" → 2 filas.** `n_distinct(Gene, Position)` no cambia;
   `nrow(position_scores)` sube por cada posición "both". Downstream que agrupa por
   `(Gene, Position, side)` (ya preparado en T2a) lo absorbe.
2. **`n_gt_2` y análogos suben.** Un residuo minoritario deja de castigar una
   convergencia mayoritaria real. Efecto: cola alta de `CAAS_score` un poco más
   poblada.
3. **Posiciones con par conservado bajan.** Nuevo drag. Magnitud a fijar con `λ_c`.
4. **`.pool_bottom` / `.pool_top` cambian de composición (hazard H5).** Ver §10.

El punto de la validación Tier 1 (PEPC en Marvin2) al cerrar T3d es que estos
cuatro deltas aparezcan **en la dirección prevista** y ningún otro.

---

## 10. Decisión de pooling — nivel lado vs nivel posición (roadmap bullet 7, H5)

Tras el split, `pos_scores` tiene una fila por `(Gene, Position, side)` con su
`CAAS_score` (media sobre esquemas de `core_s`). Para `gene_caas_score` /
`.pool_all` / `.pool_top` / `.pool_bottom`:

- **Nivel lado (recomendado).**
  - `.pool_top` = `CAAS_score` de las filas `side == "top"`.
  - `.pool_bottom` = `CAAS_score` de las filas `side == "bottom"`.
  - `.pool_all` = todas las filas-lado. Una posición "both" aporta **2 entradas**
    (una por dirección), porque tiene dos piezas independientes de evidencia
    direccional.
  - `size_adj_max` (`F(max)^n`): `n` = nº de filas-lado del gen. Premia a un gen
    con muchos eventos de convergencia direccional, que es la intención.
- **Nivel posición (dedup).** Colapsar las 2 filas-lado a 1 por `(Gene, Position)`
  (max o media) antes de poolear. Reintroduce el colapso que T3 elimina; solo si
  la validación muestra que el nivel lado infla `gene_caas_score`.

**Hazard H5 explícito.** Hoy una posición "both" fuerte-en-conjunto pero débil en
`bottom` entra en `.pool_bottom` con su score combinado (fuerte). Con nivel lado
entra con su `core_bottom` (débil). Los rankings genome-wide se mueven; es la
corrección buscada (evidencia unidireccional o de un lado débil deja de hacerse
pasar por bidireccional). Va al changelog, no es regresión.

`.pool_all` "top/bottom/all" del modo de acumulación (`change_side != 'none'` etc.)
se simplifica a `side %in% c("top")` / `c("bottom")` / todas — filtro directo sin
`c(dir, "both")` (previsto en T3d §4a).

---

## 11. Columnas emitidas por `(Gene, Position, scheme, side)` (T3a/b)

| columna | contenido |
|---------|-----------|
| `asr_path_score` | `core_s` |
| `core` | alias de `asr_path_score` (compat; `core_top`/`core_bottom` como columnas separadas se van con las filas) |
| `p_ge2` | `P(≥2)` sobre `g_c` (diagnóstico) |
| `independence` | `independence_s` (diagnóstico) |
| `conserved_drag` | `conserved_drag_s` (diagnóstico) |
| `derived_agreement` | `concentration_s` (diagnóstico; num/den del acuerdo) |
| `agree_num`, `agree_den` | numerador/denominador (T4a `convergence_type`) |
| `convergence_type` | derivado de `(agree_num, agree_den)` per §5.3 (T4a) |
| `side` | `top` / `bottom` (clave, ya en T2a) |
| `n_participating` | `m_s` |
| `pair_scores` / `pair_top_scores` / `pair_bottom_scores` | `{pid: s_c}` del lado |
| `conserved_pair_scores` / `conserved_pair_nodes` | sin cambios (latente → ahora consumido por `conserved_drag_s`) |
| `pair_ancestral` / `pair_derived_{top,bot}` | sin cambios (pooler FOP) |

Se eliminan del dict de retorno: `replication` (= `p_ge2 · independence_s`, se puede
recomponer), el `core` escalar único combinado de lados.

---

## 12. Forks abiertos — a cerrar antes de T3a

1. **`independence` por lado como factor único sobre `L_s`** (§5.4), en vez del
   "producto por pareja concreta" literal de la decisión C. Recomiendo cerrar la
   decisión C con la forma de factor único (evita el doble conteo del nodo de
   fusión). **Necesita OK.**
2. **`λ_c`** del drag de conservados (§5.5). Recomiendo `λ_c = m_s / (m_s + 1)`.
   Alternativas: `1.0`, `0.5`. **Necesita elección.**
3. **`agree_c` duro 0/1 sobre residuo codificado** (§5.1), con el acantilado
   `m_s = 2` residuos distintos → `core_s = 0`. Recomiendo mantenerlo (coherente
   con "GS encoding handles biochemistry"). Alternativa: `agree_c = concentration_s`
   (suave, mismo valor para todos los pares del lado) — reintroduciría una media.
   **Confirmar.**
4. **Pooling nivel lado** (§10). Recomiendo nivel lado. **Confirmar antes de T3d.**
5. **Golden `both_sides_two_rows`** (§8): confirmar el diseño del escenario antes
   de regenerar `path_scores_golden.json` en T3a.

Cerrados esos cinco puntos, T3a reescribe la sección de agregación de
`compute_asr_path_score` (retorno `{"top": ..., "bottom": ...}`) y `gen_golden.py`
regenera los 12 + 1 goldens; el diff debe coincidir con la columna "nuevo" de §8.
