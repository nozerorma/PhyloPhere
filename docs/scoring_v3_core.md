# scoring_v2 · core v3 — el score de convergencia CAAS sobre el dominio de Voronoi

Rama `scoring_v2`, desde `9de3aef` (post-T4b). Sucesor de
[`scoring_v2_T3_core_pareado.md`](scoring_v2_T3_core_pareado.md) (marcado *superseded*).
Especificación **hand-computable**: fija la fórmula v3 con la precisión necesaria para
recalcular a mano cada golden del apéndice. Verificación local únicamente; el usuario corre
Tier 1 PEPC en Marvin2 antes del merge.

Coordenadas 0-based. `enc(x)` = residuo codificado en el esquema activo (US → residuo crudo;
GS → etiqueta de grupo). `s ∈ {top, bottom}`. Sin em dashes (estilo de la casa).

**Estado: DRAFT (V3-0).** Se marca `FINALIZED` en V3-6.

---

## 0. Qué cambia respecto al `core` pareado (T3a–T4b)

El `core` pareado operaba sobre **pares fg/bg individuales**: cada par recorría un segmento
privado MRCA→raíz (`s_c^s`), los pares se combinaban vía `contrib(c,d)` + noisy-OR, y la
cosecha FOP multi-hipótesis se pooleaba deduplicando la unión de participantes **por nodo
MRCA**. Ese dedup por nodo era un bug latente (T4c): la unidad real es el **dominio de
Voronoi**, no el nodo, y el proxy de nodo a la vez parte dominios dispersos y fusiona
dominios distintos.

core v3 mueve todo el score al **dominio de Voronoi** como unidad. La cosecha FOP ya extrae
exactamente **un par candidato por dominio** (K dominios fijados por la selección Dunn
canónica de H1), los pesos PSS ya son por `(hipótesis, dominio)`, y la capa por-par era
precisión aparente contra la que el pooler llevaba peleando. Colapsar a dominios:

- **disuelve T4c entero** — el índice de dominio *es* la clave, siempre K de ellos;
- **borra** el walk del segmento privado, `L_s`, `walk_cache`, `iso_override`, y deja el
  pooler FOP **sin árbol** (promedia escalares ya calculados por hipótesis);
- **unifica** observado (M=1) y null-FOP (M≥1) en un solo code path;
- el factor de aislamiento `s_c^s · s_d^s` **desaparece** de `contrib` (era el principal
  motivo por el que los pares `n=2` no se movían en T3; ahora `contrib` es solo
  `agree · (1 − P_wc_any)`).

**Fuera de alcance:** el p-valor empírico `p.emp` (score-exceedance) — conversación aparte.
RER/FADE. Nada en `main`.

---

## 1. Definiciones

Por hipótesis `h`, lado `s`, esquema, posición. **K dominios de Voronoi fijos** (selección
Dunn canónica de H1). `M` = nº de hipótesis en la cosecha (`M = 1` en el observado sin FOP).

### 1.1. Dominio `d`

Un dominio es un par `(sp1_fg, sp2_bg)` con:

- `nodo_d^h = MRCA(sp1, sp2)`;
- `anc_d = enc(modal ASR @ nodo_d^h)` — en el código, `focal_state` (el estado ancestral
  modal en el MRCA del par) codificado en el esquema activo;
- `r_d^{h,s} = enc(tip fg observado en el lado s)`, si el tip está definido.

**Dominio cambiado en el lado `s`** ⟺ `r_d^{h,s}` definido **y** `r_d^{h,s} ≠ anc_d`.

No hay concepto de *conservado de metadata*: un dominio que no converge es simplemente
`score_d = 0`. Un dominio **sin reconstrucción válida** (sin ASR modal en `nodo_d`, o tip
no mapeable, o `focal_state` ausente → `anc_d` indefinido) no cuenta como cambiado y aporta
`score_d = 0`, **pero permanece en el denominador** (§3, decisión del usuario: el dato
faltante penaliza `core_s`).

### 1.2. Contraste `{d_i, d_j}` (ambos cambiados en `s`)

```
set          = { enc(r_i^{h,s}), enc(r_j^{h,s}) }          # ≤ 2 elementos; 1 si co-codifican
agree(i,j)   = 1  si  d_i y d_j comparten clase bajo el esquema  (⟺ |set| == 1)
             = 0  si no
contrib_h(d_i, d_j) = agree(i,j) · ( 1 − P_wc_any( set @ LCA(nodo_i^h, nodo_j^h) ) )
```

`P_wc_any` = `path_scores.worst_case_any_group_probability(dist, set, scheme)`: masa exacta
registrada para los targets presentes, más el remanente no registrado **una vez** si algún
target falta. Cuando `agree = 1` el `set` es un singleton `{enc(r)}` y `P_wc_any` coincide
con `worst_case_group_probability`.

**Sin factor de aislamiento. Sin PSS aquí.** El PSS entra solo en el pooling (§3).

`agree` es duro 0/1 sobre residuos codificados. El gradiente bioquímico lo da la media de
los 5 esquemas en `scoring_compute.R` §2g, no una regla suave aquí (racional idéntico al de
`scoring_v2_T3_core_pareado.md` §5.4 / §6.5).

### 1.3. `score_d^{h,s}` — noisy-OR sobre compañeros del mismo residuo

```
partners(d) = { d' ≠ d : d' cambió en s  ∧  agree(d, d') = 1 }
score_d^{h,s} = noisy_or_{d' ∈ partners(d)} contrib_h(d, d')          # 0 si partners(d) = ∅
```

`noisy_or(x) = 1 − ∏(1 − x_i)`, `path_scores.noisy_or`, cada `x_i` clampado a `[0,1]`,
entrada vacía → `0.0`.

| `partners(d)` | `score_d` |
|---|---|
| `∅` (residuo huérfano en el lado, o único dominio cambiado) | `0` |
| `{d'}` | `contrib_h(d, d')` |
| `{d', d''}`, contribs `0.8, 0.8` | `1 − 0.2·0.2 = 0.96` |

---

## 2. `compute_domain_scores` — capa por hipótesis (V3-1)

`compute_domain_scores(pair_details, per_node_dist, node_index, scheme)` devuelve, para
**una** `(Gene, Position, scheme, hipótesis)`:

```
{
  "top":    { "domain_scores": {d: score_d^{h,top}},        # solo dominios cambiados en top
              "domain_der":     {d: raw_tip_top},
              "domain_der_enc": {d: enc(raw_tip_top)},        # interno: pooling sin esquema
              "domain_anc":     {d: raw_focal_state},
              "agree_num": int, "agree_den": int, "n_changed": int,
              "convergence_type": str },                      # _convergence_type(agree_num, agree_den)
  "bottom": { ... idéntico para el lado bottom ... },
  "domain_meta": {d: {"mrca_id": int,
                      "state":    raw_focal_state | null,     # null si focal_state ausente
                      "posterior": group_probability(dist@mrca_id, anc_enc, scheme)}}
}
```

- `d` = `pair_id`; `mrca_id` = `node_id` del par.
- `agree_den == n_changed` = nº de dominios cambiados en el lado.
- `agree_num` = tamaño del mayor grupo de un mismo `domain_der_enc` entre los dominios
  cambiados (`0` si ninguno).
- `_convergence_type`: `convergent` si `den≥2 ∧ num≥2`; `divergent` si `den≥2 ∧ num<2`;
  `single` si `den==1`; `no_change` si `den==0`.
- `domain_meta` lleva **los K dominios** (todo `pair_detail`), incluidos los sin
  reconstrucción (`state: null`, `posterior: 0.0`).

**Función `score_domains_side(domains, node_index, per_node_dist, scheme)`** — `domains` =
`list[{d, mrca_id, anc_enc, der_enc}]` de los dominios cambiados en ese lado. Empareja sobre
`der_enc` igual; `p_shared = worst_case_any_group_probability(node_dist(per_node_dist,
find_lca(node_index, a.mrca_id, b.mrca_id)), {a.der_enc, b.der_enc}, scheme)`;
`contrib = 1 − p_shared`; `score_d = noisy_or`. `per_node_dist` se consulta **solo en los
nodos LCA** — no hay walk.

`compute_domain_scores` **no** recibe args de conservados. Borra de `path_scores.py`:
`side_path_score`, `_changed_side_walk`, `_apply_changed_stop`, `_conserved_side_score`,
`EMPTY_PATH_SCORE`, `parse_conserved_ids`, `_p_at_least_2`, `group_probability` (queda
`encoded_distribution` + `worst_case_*`), `modal_encoded`, `aggregate_core_side`,
`compute_asr_path_score`. **Conserva** (revividas, ahora en vivo):
`worst_case_group_probability`, `worst_case_any_group_probability`, `noisy_or`, `find_lca`,
`path_to_root_ids`, `build_node_index`, `encode_aa`, `encoded_distribution`, `node_dist`.

---

## 3. `pool_domains` — pooling = media de medias; el PSS entra solo aquí (V3-2)

`pool_domains(hyp_records, pss_by_hyp_domain=None) -> {"top": agg, "bottom": agg,
"n_hypotheses": M}`, **sin árbol** (promedia escalares).

- `hyp_records` = `list[{"hyp": str, "sides": <retorno de compute_domain_scores>}]`, `M ≥ 1`.
- `pss_by_hyp_domain` = `{(hyp, d): pss}`; `None` o clave ausente → peso `1.0`.

### 3.1. Universo de dominios (decisión V3-0, refina el sketch de V3-2)

```
universo = ⋃_h keys(hyp_records[h].sides["domain_meta"])   ∪   keys(pss_by_hyp_domain)
```

`domain_meta` lleva siempre los K dominios fijos (todo `pair_detail`), así que el universo
**siempre es K**, exista o no fichero PSS. (El sketch de V3-2 decía "union of
`domain_scores` keys ∪ `pss` keys"; eso deja fuera un dominio que no cambió en ninguna
hipótesis y no tiene PSS — se añade `domain_meta` para cerrar ese hueco y honrar la regla
"denominador = todos los K".)

### 3.2. Por lado `s`

```
M          = len(hyp_records con clave "hyp")
para cada d en universo:
    s̄_d   = ( Σ_h  sides_h["s"]["domain_scores"].get(d, 0.0) ) / M
    w̄_d   = ( Σ_h  pss(h, d, default 1.0) ) / M
core_s     = clamp01( Σ_d  w̄_d · s̄_d  /  Σ_d  w̄_d )          # 0 si Σ_d w̄_d == 0
```

Una hipótesis donde `d` no cambió (o `d` sin reconstrucción) aporta `0.0` al numerador de
`s̄_d` y su peso `pss(h,d)` (o `1.0`) al de `w̄_d`.

**`M = 1` degenera exactamente** a `core_s = Σ_d PSS_d · score_d / Σ_d PSS_d` sobre los K
dominios fijos. Sin PSS → `w̄_d = 1.0` → media aritmética sobre K.

### 3.3. `agg`

```
agg = { "asr_path_score": core_s, "core": core_s,        # alias, idénticos
        "domain_scores":  {d: s̄_d}   ∀ d en universo,
        "domain_weights": {d: w̄_d}   ∀ d en universo,
        "domain_der":     {d: _modal_str(der crudo de d sobre h)}   # solo d cambiados en ≥1 h
        "domain_anc":     {d: _modal_str(anc crudo de d sobre h)}   # solo d cambiados en ≥1 h
        "agree_num": int, "agree_den": int, "n_participating": int,
        "convergence_type": str }
```

- `agree_den == n_participating` = nº de dominios cambiados en ≥1 hipótesis.
- `agree_num` = mayor grupo de `_modal_str(domain_der_enc de d sobre h)` entre esos dominios.
- `convergence_type = _convergence_type(agree_num, agree_den)` (harvest-wide).
- El pooling **no recibe `scheme`**: la codificación ya está horneada en cada
  `domain_scores` / `domain_der_enc` aguas arriba. Un dominio que "se parte" V/I/L entre
  hipótesis bajo US puntúa bajo; bajo GS3 puntúa alto; `pool_domains` solo promedia.

`fop_pool.py` conserva `_modal_str`, `base_cycle`, `_cands`, `_dget`; borra `pool_hypotheses`,
`pool_hypotheses_pairwise`, `_domain_side_dists`, `_da_frac_from_dists`,
`_rebuild_derived_agreement`, `p_at_least_2`, `_encode_aa`, y (V3-5, ya sin caller)
`_num` / `_wmean`.

### 3.4. Colapso `side = "none"`

Si **ningún** dominio cambió en ningún lado en ninguna hipótesis, `_emit_pooled_side_rows`
(en `disambiguate_single.py`) emite una única fila `side = "none"` con `core = 0` en vez de
dos filas top/bottom a 0. `pool_domains` en sí devuelve ambos `agg` con `core_s = 0`; el
colapso vive en el caller.

---

## 4. Invariantes que **no** cambian

- Lados siempre separados; una posición `"both"` son **2 filas**, `side` es la única clave
  de dirección.
- `CAAS_score` = media (§2g de `scoring_compute.R`) de `core_s` sobre los 5 esquemas
  (US, GS4, GS3, GS2, GS1), todos corridos siempre.
- `caas_row = asr_score = asr_path_score = core_s`.
- `gene_stat == "size_adj_max"`; el guard §4d de consistencia interna (`max|delta| = 0`).
- Cardinalidad de filas por posición (≤ 2); qué posiciones se detectan (discovery +
  CT_FILTER aguas arriba); hipergeométrico `pvalue` / `gate_all` / `gate_sig`;
  `null_pvalue_boot` / `pos_perm_p` (detección-only); esquema del detail shard (8 columnas);
  esquema de `perm_pos_pval.tsv`.

---

## Apéndice A — goldens `compute_domain_scores`, nodo a nodo

Fixture de 8 nodos, dos clados (`golden/gen_golden.py::_tree_two_clades`), extendida a 3
clados donde se indica:

```
        0                       0                (extensión e3: + nodo 30 → {7,8})
       / \                     /|\
      1   2                   1 2 30
      |   |
     10   20                 path_to_root: 3→[10,1,0]  5→[20,2,0]  7→[30,0]
    / \   / \                LCA(3,5)=0  LCA(3,4)=10  LCA(5,6)=20  LCA(3,7)=0
   3   4 5   6
```

Posteriores base: `A = {A:.90, V:.05, T:.05}`, `V_here = {V:.85, A:.10, T:.05}`,
`T_here = {T:.85, A:.10, V:.05}`, `L_here = {L:.85, A:.10, V:.05}`. Nodos internos
(0,1,2,10,20,30) → `A` salvo override.

`domain_meta[d].posterior = group_probability(dist@mrca_id, enc(focal_state), scheme)`.
Con `focal_state = A` y `dist@MRCA` un `*_here` (que lleva `A:.10`): `posterior = 0.10`
bajo US. Bajo GS3, `enc(A)='n'` y el `*_here` aporta `A:.10 + (T:.05 → 'n')` cuando aplica.

Bajo **US** `worst_case_group_probability(A_dist, "V", US) = 0.05` (V registrada);
`("L") = 1 − (.90+.05+.05) = 0.0` (no registrada, remanente 0).

---

### A.1. `one_domain_changed` — sin compañero → 0

`pair_details`: `P1 {node 3, focal A, top V, bottom A}`. Esquema US. `K = 1`.

- `domain_meta`: `{1: {mrca_id: 3, state: "A", posterior: 0.10}}`.
- **top**: cambiados `{1}` (V ≠ A), `der_enc = "V"`. `partners(1) = ∅` → `score_1 = 0`.
  `domain_scores = {1: 0.0}`, `domain_der = {1: "V"}`, `domain_der_enc = {1: "V"}`,
  `domain_anc = {1: "A"}`. `agree_den = n_changed = 1`, `agree_num = 1` →
  `convergence_type = "single"`.
- **bottom**: cambiados `∅`. `domain_scores = {}`. `agree_den = agree_num = n_changed = 0`
  → `"no_change"`.

### A.2. `two_domains_same_residue_clean` — contrib limpio en la raíz

`P1 {node 3, A, top V}`, `P2 {node 5, A, top V}`; bottom ambos `A`. US. `K = 2`.

- `domain_meta`: `1 → {3, "A", 0.10}`, `2 → {5, "A", 0.10}`.
- **top**: cambiados `{1, 2}`, ambos `der_enc "V"`. `partners(1) = {2}`.
  `LCA(3,5) = 0`; `dist@0 = A`; `P_wc_any({"V"} @ 0) = 0.05`;
  `contrib_h(1,2) = 1 − 0.05 = 0.95`.
  `score_1 = noisy_or([0.95]) = 0.95`; `score_2 = 0.95`.
  `domain_scores = {1: 0.95, 2: 0.95}`. `agree_den = 2`, `agree_num = 2` →
  `"convergent"`.
- **bottom**: cambiados `∅` → `{}`, `"no_change"`.

### A.3. `two_domains_diff_residue_US` — `agree = 0`, acantilado

`P1 {node 3, A, top V}`, `P2 {node 6, A, top T}`; bottom `A`. `post[6] = T_here`. US.
`K = 2`.

- `domain_meta`: `1 → {3, "A", 0.10}`, `2 → {6, "A", 0.10}`.
- **top**: cambiados `{1 der V, 2 der T}`. `enc("V") ≠ enc("T")` → `partners(1) = partners(2)
  = ∅` → `score_1 = score_2 = 0`. `domain_scores = {1: 0.0, 2: 0.0}`,
  `domain_der = {1: "V", 2: "T"}`, `domain_der_enc = {1: "V", 2: "T"}`.
  `agree_den = 2`, `agree_num = 1` (grupos `{V:1, T:1}`) → `"divergent"`.
- **bottom**: `{}`, `"no_change"`.

### A.4. `two_domains_coencoded_GS3` — `agree = 1` vía grupo

`P1 {node 3, A, top V}`, `P2 {node 6, A, top I}`; bottom `A`.
`post[6] = {I:.85, A:.10, V:.05}`. Esquema **GS3** (`V, I, L → 'l'`; `A, T → 'n'`).
`K = 2`.

- `anc_enc` ambos = `enc("A", GS3) = 'n'`.
- `domain_meta.posterior = group_probability(dist, 'n', GS3)`:
  - `d1` @ node 3 `V_here = {V:.85, A:.10, T:.05}` → `'n'` = `A:.10 + T:.05` = **0.15**.
  - `d2` @ node 6 `{I:.85, A:.10, V:.05}` → `'n'` = `A:.10` = **0.10**.
  - `domain_meta`: `1 → {3, "A", 0.15}`, `2 → {6, "A", 0.10}`.
- **top**: `d1` tip `V → 'l' ≠ 'n'` cambiado, `der_enc 'l'`, raw `"V"`.
  `d2` tip `I → 'l' ≠ 'n'` cambiado, `der_enc 'l'`, raw `"I"`. `partners(1) = {2}`.
  `LCA(3,6) = 0`; `dist@0 = A`; enc GS3 `{n: .95, l: .05}`; `P_wc_any({'l'} @ 0) = 0.05`;
  `contrib = 0.95`. `score_1 = score_2 = 0.95`.
  `domain_scores = {1: 0.95, 2: 0.95}`, `domain_der = {1: "V", 2: "I"}`,
  `domain_der_enc = {1: "l", 2: "l"}`, `domain_anc = {1: "A", 2: "A"}`.
  `agree_den = 2`, `agree_num = 2` → `"convergent"`.
- **bottom**: `{}`, `"no_change"`.

### A.5. `two_domains_opposite_sides` — cada lado bajo umbral

`P1 {node 3, A, top V, bottom A}`, `P2 {node 5, A, top A, bottom V}`. US. `K = 2`.

- `domain_meta`: `1 → {3, "A", 0.10}`, `2 → {5, "A", 0.10}`.
- **top**: cambiados `{1}` (der V). `partners(1) = ∅` → `score_1 = 0`.
  `domain_scores = {1: 0.0}`. `agree_den = 1` → `"single"`.
- **bottom**: cambiados `{2}` (der V). `score_2 = 0`. `domain_scores = {2: 0.0}`.
  `agree_den = 1` → `"single"`.

### A.6. `three_domains_VVL_majority` — mayoría `V`, el `L` mete un 0

Extensión e3. `P1 {node 3, A, top V}`, `P2 {node 5, A, top V}`, `P3 {node 7, A, top L}`;
bottom `A`. `post[3] = post[5] = V_here`, `post[7] = L_here`, resto `A`. US. `K = 3`.

- `domain_meta`: `1 → {3, "A", 0.10}`, `2 → {5, "A", 0.10}`, `3 → {7, "A", 0.10}`.
- **top**: cambiados `{1 V, 2 V, 3 L}`.
  - `partners(1) = {2}`, `partners(2) = {1}`, `partners(3) = ∅`.
  - `LCA(3,5) = 0`; `contrib(1,2) = 1 − P_wc_any({"V"}@0) = 1 − 0.05 = 0.95`.
  - `score_1 = score_2 = noisy_or([0.95]) = 0.95`; `score_3 = 0`.
  - `domain_scores = {1: 0.95, 2: 0.95, 3: 0.0}`,
    `domain_der = {1: "V", 2: "V", 3: "L"}`, `domain_der_enc` idem, `domain_anc` todo `"A"`.
  - `agree_den = 3`, `agree_num = 2` (`{V:2, L:1}`) → `"convergent"`.
- **bottom**: `{}`, `"no_change"`.

### A.7. `shared_lca_contaminated` — el LCA ya lleva el derivado

`P1 {node 3, A, top V}`, `P2 {node 4, A, top V}`; bottom `A`. `post[3] = post[4] = V_here`,
`post[10] = {V:.70, A:.30}` (override). US. `K = 2`.

- `domain_meta`: `1 → {3, "A", 0.10}`, `2 → {4, "A", 0.10}`.
- **top**: cambiados `{1, 2}` der V. `partners(1) = {2}`. `LCA(3,4) = 10`;
  `dist@10 = {V:.70, A:.30}`; `P_wc_any({"V"} @ 10) = 0.70`;
  `contrib(1,2) = 1 − 0.70 = 0.30`. `score_1 = score_2 = noisy_or([0.30]) = 0.30`.
  `domain_scores = {1: 0.30, 2: 0.30}`. `agree_den = 2`, `agree_num = 2` → `"convergent"`.
- **bottom**: `{}`, `"no_change"`.

### A.8. `no_domain_changed` — `no_change` en ambos lados

`P1 {node 3, A, top A, bottom A}`, `P2 {node 5, A, top A, bottom A}`. US. `K = 2`.

- `domain_meta`: `1 → {3, "A", 0.10}`, `2 → {5, "A", 0.10}`.
- **top** y **bottom**: cambiados `∅` → `domain_scores = {}`,
  `agree_den = agree_num = n_changed = 0` → `"no_change"`.

### A.9. `domain_without_reconstruction` — dominio sin `anc`

Extensión e3. `P1 {node 3, A, top V}`, `P2 {node 5, A, top V}`,
`P3 {node 7, focal_state "", top V, bottom A}`; bottom `A`. `post[3] = post[5] = V_here`,
resto `A`. US. `K = 3`.

- `d3`: `focal_state ""` → `anc_enc = encode_aa("") = None` → no clasificable como cambiado.
  `domain_meta`: `3 → {mrca_id: 7, state: null, posterior: 0.0}`.
- `domain_meta`: `1 → {3, "A", 0.10}`, `2 → {5, "A", 0.10}`, `3 → {7, null, 0.0}`.
- **top**: cambiados `{1, 2}` (d3 excluido). `partners(1) = {2}`. `LCA(3,5) = 0`;
  `contrib(1,2) = 0.95`. `score_1 = score_2 = 0.95`.
  `domain_scores = {1: 0.95, 2: 0.95}` (**sin `3`**), `domain_der = {1: "V", 2: "V"}`,
  `domain_anc = {1: "A", 2: "A"}`. `agree_den = 2`, `agree_num = 2` → `"convergent"`.
- **bottom**: `{}`, `"no_change"`.

El efecto del slot vacío se ve en `pool_domains` (§ Apéndice B.5), no aquí.

---

## Apéndice B — goldens `pool_domains`, aritmética completa

Cada `hyp_record` = `{"hyp": "H<n>", "sides": <retorno de compute_domain_scores>}`.
`pss` como lista `[[hyp, domain, weight], ...]` → `{(hyp, domain): weight}`; ausente →
peso `1.0`. Resultado: `{"top": agg, "bottom": agg, "n_hypotheses": M}`.

### B.0. `m1_passthrough` — `M = 1`, sin PSS, degenera a media

`hyp_records = [{"hyp": "H1", "sides": <A.2>}]`. Universo (`domain_meta` de A.2) = `{1, 2}`.

- **top**: `s̄_1 = 0.95/1 = 0.95`, `s̄_2 = 0.95`. `w̄_1 = w̄_2 = 1.0/1 = 1.0`.
  `core_top = (1.0·0.95 + 1.0·0.95) / (1.0 + 1.0) = 1.90/2 = 0.95`.
  `domain_scores = {1: 0.95, 2: 0.95}`, `domain_weights = {1: 1.0, 2: 1.0}`,
  `domain_der = {1: "V", 2: "V"}`, `domain_anc = {1: "A", 2: "A"}`.
  `agree_den = n_participating = 2`, `agree_num = 2` → `"convergent"`.
  `asr_path_score = core = 0.95`.
- **bottom**: `s̄_1 = s̄_2 = 0.0`, `w̄ = 1.0`. `core_bottom = 0 / 2 = 0.0`.
  `domain_scores = {1: 0.0, 2: 0.0}`, `domain_weights = {1: 1.0, 2: 1.0}`,
  `domain_der = {}`, `domain_anc = {}`. `agree_den = n_participating = 0` → `"no_change"`.
- `n_hypotheses = 1`.

### B.1. `m1_pss_shifts_to_d1` — `M = 1`, PSS reordena el diseño

`hyp_records = [{"hyp": "H1", "sides": <A.6>}]` (VVL: `score` top `{1: 0.95, 2: 0.95,
3: 0.0}`). `pss = [["H1",1,0.9], ["H1",2,0.9], ["H1",3,0.1]]`. Universo = `{1, 2, 3}`.

- **top**: `s̄ = {1: 0.95, 2: 0.95, 3: 0.0}`. `w̄ = {1: 0.9, 2: 0.9, 3: 0.1}`.
  `Σ w̄·s̄ = 0.9·0.95 + 0.9·0.95 + 0.1·0 = 0.855 + 0.855 = 1.71`.
  `Σ w̄ = 1.9`. `core_top = 1.71 / 1.9 = 0.90`.
  (Equipeso daría `(0.95 + 0.95 + 0) / 3 = 0.6333`; el PSS empuja hacia los dominios V.)
  `domain_scores = {1: 0.95, 2: 0.95, 3: 0.0}`, `domain_weights = {1: 0.9, 2: 0.9,
  3: 0.1}`, `domain_der = {1: "V", 2: "V", 3: "L"}`, `domain_anc` todo `"A"`.
  `agree_den = 3`, `agree_num = 2` → `"convergent"`. `n_participating = 3`.
- **bottom**: `s̄ = 0` ∀ d. `Σ w̄·s̄ = 0`, `Σ w̄ = 1.9` → `core_bottom = 0.0`.
  `domain_der = {}`, `agree_den = 0` → `"no_change"`.
- `n_hypotheses = 1`.

### B.2. `m2_domain_changed_only_in_H1` — `s̄_d = score / M`

`hyp_records = [{"hyp": "H1", "sides": <A.2>}, {"hyp": "H2", "sides": <A.8>}]`.
Sin PSS. Universo = `{1, 2}`. `M = 2`.

- **top**: `Σ_h score_1 = 0.95 (H1) + 0.0 (H2) = 0.95` → `s̄_1 = 0.95/2 = 0.475`.
  `s̄_2 = 0.475`. `w̄_1 = w̄_2 = (1.0 + 1.0)/2 = 1.0`.
  `core_top = (1.0·0.475 + 1.0·0.475) / 2 = 0.475`.
  `domain_scores = {1: 0.475, 2: 0.475}`, `domain_weights = {1: 1.0, 2: 1.0}`.
  `domain_der`: `d1` sobre `h` = `["V"]` (H2 no lo cambió) → `_modal_str = "V"`; `d2 = "V"`.
  `agree_den = 2` (cambiados en ≥1 h), `agree_num = 2` → `"convergent"`.
  `n_participating = 2`.
- **bottom**: todo 0 → `core_bottom = 0.0`, `"no_change"`.
- `n_hypotheses = 2`.

### B.3. `m2_asymmetric_pss_same_domain` — `w̄_d` = media de PSS por hipótesis

`hyp_records = [{"hyp": "H1", "sides": <A.6>}, {"hyp": "H2", "sides": <A.6>}]`.
`pss = [["H1",1,0.1], ["H2",1,0.9], ["H1",2,0.9], ["H2",2,0.9], ["H1",3,0.5],
["H2",3,0.5]]`. Universo = `{1, 2, 3}`. `M = 2`.

- **top**: `s̄_1 = (0.95 + 0.95)/2 = 0.95`, `s̄_2 = 0.95`, `s̄_3 = (0 + 0)/2 = 0`.
  `w̄_1 = (0.1 + 0.9)/2 = 0.5`, `w̄_2 = (0.9 + 0.9)/2 = 0.9`, `w̄_3 = (0.5 + 0.5)/2 = 0.5`.
  `Σ w̄·s̄ = 0.5·0.95 + 0.9·0.95 + 0.5·0 = 0.475 + 0.855 = 1.33`. `Σ w̄ = 1.9`.
  `core_top = 1.33 / 1.9 = 0.70`.
  `domain_scores = {1: 0.95, 2: 0.95, 3: 0.0}`, `domain_weights = {1: 0.5, 2: 0.9,
  3: 0.5}`, `domain_der = {1: "V", 2: "V", 3: "L"}`. `agree_den = 3`, `agree_num = 2`
  → `"convergent"`. `n_participating = 3`.
- **bottom**: `core_bottom = 0.0`, `"no_change"`.
- `n_hypotheses = 2`.

### B.4. `side_none_collapse` — nada cambió; ambos `agg` a 0

`hyp_records = [{"hyp": "H1", "sides": <A.8>}]`. Sin PSS. Universo = `{1, 2}`. `M = 1`.

- **top** y **bottom**: `s̄_1 = s̄_2 = 0.0`, `w̄_1 = w̄_2 = 1.0`.
  `Σ w̄·s̄ = 0`, `Σ w̄ = 2.0` → `core_s = 0 / 2 = 0.0`.
  `domain_scores = {1: 0.0, 2: 0.0}`, `domain_weights = {1: 1.0, 2: 1.0}`,
  `domain_der = {}`, `domain_anc = {}`. `agree_den = agree_num = n_participating = 0`
  → `"no_change"`.
- `n_hypotheses = 1`. El caller (`_emit_pooled_side_rows`) colapsa esto a una fila
  `side = "none"`.

### B.5. `no_reconstruction_in_denominator` — el slot vacío tira `core` abajo

`hyp_records = [{"hyp": "H1", "sides": <A.9>}]`. Sin PSS.
Universo = `domain_meta` de A.9 = `{1, 2, 3}`. `M = 1`.

- **top**: `s̄_1 = 0.95`, `s̄_2 = 0.95`, `s̄_3 = get(3, 0.0) / 1 = 0.0` (d3 nunca en
  `domain_scores`). `w̄_1 = w̄_2 = w̄_3 = 1.0` (PSS por defecto para el slot vacío también).
  `Σ w̄·s̄ = 0.95 + 0.95 + 0 = 1.90`. `Σ w̄ = 3.0`. `core_top = 1.90 / 3 = 0.633333…`.
  (Si d3 se excluyera del denominador: `1.90 / 2 = 0.95`. La regla "denominador = todos los
  K" penaliza el dato faltante.)
  `domain_scores = {1: 0.95, 2: 0.95, 3: 0.0}`, `domain_weights = {1: 1.0, 2: 1.0,
  3: 1.0}`, `domain_der = {1: "V", 2: "V"}` (d3 nunca cambió), `domain_anc = {1: "A",
  2: "A"}`. `agree_den = 2`, `agree_num = 2` → `"convergent"`. `n_participating = 2`.
- **bottom**: `core_bottom = 0.0`, `"no_change"`.
- `n_hypotheses = 1`.

### B.6. `m3_split_scheme_resolved_US` vs `_GS3` — el pooling es agnóstico de esquema

3 hipótesis, dominios `d1` (node 3) y `d2` (node 5). El tip fg de `d1` varía por hipótesis
(H1 → V, H2 → I, H3 → L); `d2 → V` en las tres. `focal A`, `post[3] = post[5] = V_here`
(y `{I}` / `{L}` según haga falta para el tip; el tip no lee del posterior, se pasa en
`pair_details`).

**`compute_domain_scores` por hipótesis, esquema US:**

| h | `d1` der | `d2` der | `domain_scores` top | nota |
|---|---|---|---|---|
| H1 | V | V | `{1: 0.95, 2: 0.95}` | `partners` mutuos, `contrib(1,2)@0 = 0.95` |
| H2 | I | V | `{1: 0.0, 2: 0.0}` | `enc(I) ≠ enc(V)` → sin `partners`; `agree_den 2, agree_num 1` → `divergent` |
| H3 | L | V | `{1: 0.0, 2: 0.0}` | idem |

**`pool_domains` US** (`M = 3`, sin PSS, universo `{1, 2}`):

- `s̄_1 = (0.95 + 0 + 0)/3 = 0.316667`. `s̄_2 = (0.95 + 0 + 0)/3 = 0.316667`.
  `w̄_1 = w̄_2 = 1.0`. `core_top = (0.316667 + 0.316667)/2 = 0.316667`.
- `domain_der`: `d1` sobre h = `["V", "I", "L"]` → `_modal_str` (empate, first-seen) = `"V"`;
  `d2 = "V"`. `domain_der_enc` modal: `d1 → "V"`, `d2 → "V"`.
  `agree_den = 2`, `agree_num = 2` (`{V: 2}`) → `"convergent"` (etiqueta = conteo físico;
  el score es bajo — son cosas distintas). `n_participating = 2`.
- `domain_scores = {1: 0.316667, 2: 0.316667}`, `domain_weights = {1: 1.0, 2: 1.0}`.
- **bottom**: `core_bottom = 0.0`. `n_hypotheses = 3`.

**`compute_domain_scores` por hipótesis, esquema GS3** (`V, I, L → 'l'`, `A → 'n'`;
`dist@0` enc GS3 `{n: .95, l: .05}`):

| h | `d1` der_enc | `d2` der_enc | `domain_scores` top |
|---|---|---|---|
| H1 | l | l | `{1: 0.95, 2: 0.95}` |
| H2 | l | l | `{1: 0.95, 2: 0.95}` |
| H3 | l | l | `{1: 0.95, 2: 0.95}` |

**`pool_domains` GS3** (`M = 3`, universo `{1, 2}`):

- `s̄_1 = (0.95·3)/3 = 0.95`. `s̄_2 = 0.95`. `w̄ = 1.0`.
  `core_top = (0.95 + 0.95)/2 = 0.95`.
- `domain_der` (crudo): `d1` → `_modal_str(["V","I","L"]) = "V"`, `d2 = "V"`.
  `domain_der_enc` modal: `d1 → "l"`, `d2 → "l"`. `agree_den = 2`, `agree_num = 2` (`{l: 2}`)
  → `"convergent"`. `n_participating = 2`.
- `domain_scores = {1: 0.95, 2: 0.95}`, `domain_weights = {1: 1.0, 2: 1.0}`.
- **bottom**: `core_bottom = 0.0`. `n_hypotheses = 3`.

Mismo harvest, mismo `pool_domains` (sin `scheme`), resultado `0.317` vs `0.95` según el
esquema bajo el que se corrió `compute_domain_scores`. La resolución bioquímica vive aguas
arriba; el pooling solo promedia.

---

## Apéndice C — deltas Tier-1 (PEPC) esperados

**Cambian (changelog, no regresión):**

- `CAAS_score` por posición **sube** de forma amplia — el factor de aislamiento
  `s_c^s · s_d^s` (cada uno `∈ [0,1]`, normalmente `< 1`) desaparece de `contrib`.
- Semántica del denominador — v2 era `|P_s| + n_conserved`; v3 es `Σ_{d=1..K} w̄_d` sobre
  los K dominios fijos. Posiciones cuyos dominios cambiaron **solo en el otro lado**
  puntúan más bajo en el lado `s`.
- Posiciones FOP multi-hipótesis re-escoran — media de medias reemplaza el
  `noisy_or` / `p_at_least_2` de pooling por dominio; el PSS pesa dominios, no hipótesis.
- Etiquetas `convergence_type` cambian para algunas posiciones (recomputadas harvest-wide).
- `gene_caas_score`, rankings de `position_scores.tsv` / `gene_scores.tsv` se reordenan.
- `position_scores.tsv` pierde `convergence_schemes`, `n_conserved_pairs`;
  `caas_convergence_master.csv` renombra `mrca_<i>_*` → `domain_<d>_*`, deja caer
  `conserved_<j>_*` y `independence`.
- Conjuntos de posiciones VEP PrimateAI / COSMIC crecen — el skip-gate `convergence_schemes`
  se retira.

**NO cambian:** cardinalidad de filas por posición (≤ 2); asignación de `side`; qué
posiciones se detectan; estructura de la media de 5 esquemas §2g; hipergeométrico `pvalue` /
`gate_all` / `gate_sig`; `null_pvalue_boot` / `pos_perm_p`; esquema del detail shard
(8 columnas); `perm_pos_pval.tsv`; `gene_stat == "size_adj_max"` y el guard §4d; RER / FADE.

---

## Apéndice D — secuencia de commits

`V3-0` (este documento + red golden xfail) · `V3-1` (`compute_domain_scores` en
`path_scores.py`; regen de `core_v3_golden.json` con diff vacío) · `V3-2` (`pool_domains` +
`disambiguate_single` + `models.py` + writers) · **`V3-3` SHIPPED** (null path:
`gene_wrapper.py::_perms_worker` reescrito sobre `pool_domains`; `_expand_pooled` reemplaza
`_expand_sides`; se retiran `pool_hypotheses_pairwise` / `_nss_node_index` /
`_nss_per_node_dist` / `build_node_index` del worker — el pooler null es sin árbol.
FOP: una llamada `pool_domains` por `(base cycle, pos, scheme)`; no-FOP: `M = 1` por
record. Detail shard 8-col y `perm_pos_pval.tsv` sin cambios. Test nuevo
`test_null_domain_pool_wiring.py` fija la identidad observado == null sobre los goldens
Apéndice B) · `V3-4` (`scoring_compute.R` + `residue_descriptors.py` +
VEP + schema downstream) · `V3-5` (borrar `fop_pool.R` + poda de huérfanos) · `V3-6`
(borrar `dunn_modified.R` + unificar parsers TSV + prosa final; marcar este doc `FINALIZED`).

`git branch --show-current` antes de cada commit (hazard de worktree concurrente).
