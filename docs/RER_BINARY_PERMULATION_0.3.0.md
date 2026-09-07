# Permulación binaria de RERconverge 0.3.0 vs. el null de contrastes de PhyloPhere

Documento de investigación, cotejado línea a línea con el código fuente de RERconverge en el
commit **`2bd328f7530b4aca9b48c0b3997875c9b77a7026`** (2026-04-17, `DESCRIPTION: Version 0.3.0`)
y con la rama `nongreedy_dunn` de PhyloPhere (septiembre 2026).

`K` = número de pares de contraste independientes = `N_pairs_obs`. `Q` = matriz de tasas de
transición Mk. `masterTree` = árbol maestro del `treesObj` de `readTrees()`.

---

## 0. Aclaración de versiones — cuál es "0.3.0"

Hay tres cosas distintas llamadas "0.3.0" y no coinciden:

| Referencia | Qué es | Motor CC binario | `categoricalPermulations` |
|---|---|---|---|
| `environment/phylophere.yml:631` → `rerconverge=0.3.0=r44h7b50bb2_3` (bioconda) | tag `v0.3.0` de git = commit `f6ea2ac`, **2021-03-02** | `generatePermulatedBinPhen` → `simBinPhenoCC` (BM + enraizado + rejection) | **NO existe** |
| `environment/install_env.sh:79` → `install_github("nclark-lab/RERconverge@2bd328f7…")`, **2026-04-17** | master reciente, `DESCRIPTION` sigue diciendo `0.1.0`/`0.3.0` | `getPermsBinary(permmode="cc")` → **`categoricalPermulations`** (Mk + stochastic mapping) | **SÍ** |
| entorno instalado ahora en `correfoc` (`~/.conda/envs/phylophere`) | build intermedio; `packageVersion` reporta `0.1.0` | legacy `simBinPhenoCC` | **NO** (`exists("categoricalPermulations") == FALSE`) |

**`install_env.sh` es la fuente de verdad del proyecto** y fija el commit `2bd328f7`, que
tiene `categoricalPermulations`. La entrada `rerconverge=0.3.0` del `.yml` es *dead weight*:
`install_env.sh` la sobrescribe con `install_github`. El entorno de la cuádrícula está
**desactualizado** respecto a ese pin y hay que reinstalarlo antes de que las permulaciones
binarias del RER puedan correr; hasta entonces el guardián de versión en
[`binary_rer.R`](../subworkflows/RERCONVERGE/local/binary_rer.R) hace `stop()` de forma
correcta.

El resto de este documento se refiere por "0.3.0" al commit **`2bd328f7`**, salvo mención
explícita del motor *legacy*.

---

## 1. `categoricalPermulations` 0.3.0 — algoritmo exacto

Ruta del código: `getPermsBinary(permmode="cc")`
([`PermulationFuncs.R:250-304`](https://github.com/nclark-lab/RERconverge/blob/2bd328f7530b4aca9b48c0b3997875c9b77a7026/R/PermulationFuncs.R#L250))
→ `categoricalPermulations`
([`:2688-2730`](https://github.com/nclark-lab/RERconverge/blob/2bd328f7530b4aca9b48c0b3997875c9b77a7026/R/PermulationFuncs.R#L2688))
→ `getNullTips` `:2424` / `getNullTrees` `:2485` / `shuffleInternalNodes` `:2461` / `improveTree` `:2512`
→ `getAncLiks` ([`RERfuncs.R:2064`](https://github.com/nclark-lab/RERconverge/blob/2bd328f7530b4aca9b48c0b3997875c9b77a7026/R/RERfuncs.R#L2064)).

Dependencias externas: `castor::fit_mk`, `castor::simulate_mk_model`, `castor::map_to_state_space`,
`phytools::to.matrix`, `Matrix::expm`.

### 1.1 Modelo de transición

```
tree     <- pruneTree(treesObj$masterTree, species presentes en phenvals)
intlab   <- map_to_state_space(phenvals)          # 0/1 (o n-ario) -> enteros 1..Nstates
Q        <- fit_mk(tree, Nstates, intlab$mapped_states,
                   rate_model = "ER", root_prior = "auto")$transition_matrix
```

- **Mk de tasas iguales (`rate_model="ER"`)**: un único parámetro de tasa, matriz Q
  simétrica, ajustada por **máxima verosimilitud** sobre los estados observados en las
  puntas. Es exactamente el modelo que se usaría para una reconstrucción de estados
  ancestrales de un carácter discreto.
- `root_prior="auto"` deja que `castor` elija el prior de raíz (típicamente el
  estacionario para ER).
- Para binario, `phenotypeVector` es `0` (background) / `1` (foreground) para **todas** las
  puntas del `masterTree` (`getPermsBinary:263-265`), así que no hay poda efectiva.

### 1.2 Simulación de puntas nulas — `getNullTips`

```
true_counts <- table(estados observados en puntas)          # p.ej. {bg: 34, fg: 6}
repeat hasta acumular N=ntrees aceptadas:
    sim <- simulate_mk_model(tree, Q, root_probabilities = "stationary")
    # CTMC hacia adelante desde una raíz muestreada del estacionario de Q;
    # devuelve estados en TODAS las puntas y en TODOS los nodos internos
    sim_counts <- table(sim$tip_states)
    reject si falta algún estado en las puntas
    reject salvo que  |sim_counts - true_counts| <= true_counts * percent_relax   (elemento a elemento)
    # getPermsBinary llama con percent_relax = 0  =>  sim_counts == true_counts EXACTO
    accept: guardar sim$tip_states y sim$node_states
```

> **Consecuencia clave:** con `percent_relax = 0` el null **preserva exactamente el número
> de puntas foreground** (y de background). Lo que se aleatoriza es *qué* puntas son
> foreground, con la distribución inducida por el CTMC Mk (que tiende a agrupar
> filogenéticamente, en la medida en que Q lo permita).

### 1.3 Estados internos de arranque — `getNullTrees` + `shuffleInternalNodes`

```
node_states_obs <- getStatesAtNodes( getAncLiks(tree, estados_obs, Q) )   # MAP por nodo, datos OBSERVADOS
for i in 1..N:
    tips_i   <- null_tips[i, ]
    ancliks_i <- getAncLiks(tree, tips_i, Q)          # verosimilitud marginal por estado y nodo, dado el set nulo
    shuffled <- sample(node_states_obs)               # permuta el MULTISET de estados MAP observados
    for state in shuffled (sin reemplazo):
        node <- sample(nodos_disponibles, prob = ancliks_i[, state])
        asigna 'state' a 'node'; retira 'node' del pool
```

- Se toma el **multiset de estados MAP de los nodos internos en los datos observados**
  (p.ej. "5 nodos foreground, 30 background") y se **baraja** sobre los nodos internos.
- La colocación no es uniforme: cada estado se coloca en un nodo muestreado ∝ a la
  verosimilitud ancestral de ese estado en ese nodo, calculada con Q para *ese* conjunto de
  puntas nulas.
- Resultado: un historia de arranque `(tips_i, internal_states_i)` que respeta (a) el conteo
  de puntas, (b) el conteo de nodos internos por estado, (c) sesgada hacia colocaciones
  verosímiles.

### 1.4 Pulido por *simulated annealing* — `improveTree`

```
improveTree(tree, Q, P = {expm(Q * bl_e) para cada arista e}, nodes, tips,
            T0 = 10, Nk = 10, cycles = 100, alpha = 0.9)
```

```
L(nodes) <- prod_e  P[e][ estado(anc_e), estado(desc_e) ]        # verosimilitud de transición de la historia
Tk <- T0
repeat cycles veces, Nk iteraciones por ciclo:
    elige nodo n1 y un par de estados (s1 actual, s2 alternativo)  ~  ratios de verosimilitud ancestral
    elige nodo n2 que tenga s2 y quiera s1
    swap: nodes[n1] <- s2 ; nodes[n2] <- s1                       # PERMUTACIÓN: el multiset no cambia
    r <- L(nueva) / L(vieja)   (recalculada solo en las aristas afectadas)
    si r >= 1: acepta
    si r < 1 : acepta con prob  u = exp( -(-log(L*r)+log L) / Tk )   (Metropolis)
    al final de cada ciclo:  Tk <- T0 / (1 + alpha * k)              (enfriamiento hiperbólico)
return nodes  (multiset intacto; verosimilitud de transición aumentada)
```

- Es un **recocido simulado sobre permutaciones del multiset de estados internos**, con
  criterio de energía = log-verosimilitud de transición bajo Q.
- **No es aceptación/rechazo del árbol nulo entero**: siempre devuelve un resultado. El
  "polish" solo reorganiza los estados internos de arranque hacia una historia más
  verosímil, sin tocar las puntas ni el multiset.
- `rp="auto"` de `categoricalPermulations` se pasa a `fit_mk` como `root_prior`; NO controla
  `improveTree`. `percent_relax` solo controla `getNullTips` (§1.2).

### 1.5 De historia nula a correlación — `getPermsBinary`

```
para cada historia nula  x = (tips, nodes):
    tr <- masterTree ;  tr$edge.length <- c(x$tips, x$nodes)[ tr$edge[,2] ] - 1
    # arista = 1 (foreground) sii el estado de su nodo hijo es 2 (fg);  0 en otro caso
    fg_null <- getForegroundsFromBinaryTree(tr)                 # puntas con arista terminal = 1
    path_null <- foreground2Paths(fg_null, treesObj, clade = "all")   # LONGITUD == ncol(RERmat)
    cor_i <- correlateWithBinaryPhenotype(path_null, RERmat)    # Kendall (no ponderada) por defecto
devuelve list(corP, corRho, corStat)  con  nrow = nrow(RERmat),  ncol = numperms
```

`correlateWithBinaryPhenotype` se llama **con sus valores por defecto** dentro de
`getPermsBinary`: `clade="all"`, `weighted="auto"` (→ Kendall porque los pesos son 0/1),
**sin winsorización**, `min.sp=10`, `min.pos=2`
([`PermulationFuncs.R:288`](https://github.com/nclark-lab/RERconverge/blob/2bd328f7530b4aca9b48c0b3997875c9b77a7026/R/PermulationFuncs.R#L288)).
El wrapper de PhyloPhere recalcula por eso una correlación observada de referencia
`res_ref` con esos mismos ajustes para que `permpvalcor` compare like-with-like.

### 1.6 Qué preserva y qué aleatoriza

| Preserva exactamente | Aleatoriza | NO controla |
|---|---|---|
| nº de puntas foreground (`percent_relax=0`) | qué puntas son foreground | nº de orígenes independientes / transiciones 0→1 |
| multiset de estados MAP de nodos internos | disposición de estados internos (multiset fijo, pulido SA) | estructura de clados del subárbol foreground |
| matriz de tasas Q (ER, ajuste ML observado) | historia de transición | nada relativo a pares de contraste |

---

## 2. ¿Cualitativo Y cuantitativo? — API de permulación en 0.3.0 (`2bd328f7`)

| Función (exportada) | Tipo de rasgo | Generador del null | Estadístico observado | p.perm |
|---|---|---|---|---|
| `getPermsContinuous(type="simperm")` | continuo | `simpermvec`: BM `ratematrix`+`sim.char` en `mastertree`, **rank-match** de los valores observados sobre el orden simulado | `correlateWithContinuousPhenotype` (Pearson winsorizada) | `permpvalcor` |
| `getPermsBinary(permmode="cc")` | binario (0/1) | `categoricalPermulations(rm="ER")` → Mk + stochastic mapping + pulido SA (§1) | `correlateWithBinaryPhenotype` (Kendall / Pearson ponderada) | `permpvalcor` |
| `getPermsBinary(permmode="ccLegacy")` | binario | `simBinPhenoCC`: BM del `pathvec`, umbral top-K, rejection por `blsum` + orden de profundidad; **requiere `root_sp`** | idem | `permpvalcor` |
| `getPermsBinary(permmode="ssm"/"ssmLegacy")` | binario, por árbol de gen | permulación *species-specific mapping* (una historia nula por topología de gen única) | `calculateCorPermuted` | `permpvalcor` |
| `categoricalPermulations` + `getPermPvalsCategorical(method="kw")` | **categórico n-ario / ordinal** | mismo motor Mk que §1, `percent_relax` por estado configurable | `correlateWithCategoricalPhenotype` (Kruskal–Wallis / ANOVA) + tests por pares de niveles | `getPermPvalsCategorical` |
| `getPermsContinuousExtantOnly` / `getPermsBinaryExtantOnly` | continuo / binario, solo extantes | permuta/simula solo estados de puntas, sin paths internos | `getAllCorExtantOnly` | idem |

Observaciones:

1. **El motor Mk categórico es unificado para binario y multiestado.** El binario es
   simplemente `Nstates = 2`. El ordinal se trata como categórico sin ordenación (KW), no
   hay un motor "ordinal" propio en RERconverge.
2. **El continuo sigue siendo un camino separado** (`getPermsContinuous` / `simpermvec`).
   0.3.0 **no** reencauza el continuo por el motor categórico.
3. **La permulación continua no cambió entre v0.2.0, v0.3.0 y `2bd328f7`**: las firmas y
   los cuerpos de `getPermsContinuous`, `getNullCor`, `simulatevec` y `simpermvec` son
   byte-idénticos en los tres (verificado con `diff`). El único cambio 0.1→0.3 en la
   familia binaria fue mover el CC a `categoricalPermulations` y renombrar el antiguo a
   `ccLegacy`.

---

## 3. Contraste con `permulations.R` de PhyloPhere

Código: [`subworkflows/CT/local/scripts/permulations.R`](../subworkflows/CT/local/scripts/permulations.R),
[`lean_contrast_selector.R`](../subworkflows/CT/local/scripts/lean_contrast_selector.R),
[`pss_core.R`](../subworkflows/CT/local/scripts/pss_core.R).

### 3.1 Qué hace el null de contrastes, paso a paso

```
# --- UNA VEZ, sobre el rasgo observado ---
obs_fits   <- fit_models(pruned.tree, starting.values)      # geiger::fitContinuous BM y OU
sel_model  <- select_model(obs_fits)                        # OU sii AIC_OU + 2 < AIC_BM  (y alpha OU no pegada a la cota)
cov_bm, cov_ou <- covariances_from_fits(...)                # matrices de covarianza, FIJAS para todo el harvest
simulation_tree <- if (sel_model == "OU") rescale(pruned.tree, "OU", alpha)  else  pruned.tree

# --- POR CADA DRAW ---
sim_v  <- simulatevec(starting.values, simulation_tree)     # BM (o BM sobre árbol OU-reescalado) + ratematrix + sim.char
pvec   <- rec_value[ order(order(sim_v)) ]                  # RANK-MATCH: valores observados en el orden simulado
                                                           # (para count data se arrastran también CI_lb, CI_ub, n)
e <- evaluate_lean_contrast_selection(pvec, D, target_pairs = K, cov_bm, cov_ou, sel_model, ...)
```

`evaluate_lean_contrast_selection` (idéntico en política de candidatos al selector observado
`selection_algorithm.R::pair_sel.f`, vía el core compartido `rank_candidates` /
`greedy_dunn_select`):

```
1. PSS: calculate_pairwise_scores(pvec, tr, cov_bm, cov_ou, sel_model)   # sobre covarianzas OBSERVADAS fijas
       s = -expm1(log 2 + pnorm(|Δtrait| / sqrt(var_par), lower=F, log=T))   # prob. de una diferencia >= |Δ| bajo el modelo
       FinalScore = s * (|Δtrait|/max) / (patristic/max)
2. Puerta por tipo de rasgo:
       count    : CI_lb[hi] > CI_ub[lo]        (no solapamiento de Jeffreys)
       ordinal  : trait[hi] == nivel_max  &  trait[lo] == nivel_min
       continuo : todas las parejas con Δtrait > 0
3. (solo continuo) quedarse con el top `pss_top_pct` (def. 1%) de las supervivientes de (2) por PSS
4. rank_candidates: PSS desc -> |Δtrait| desc -> pair_n desc -> orden estable
5. greedy_dunn_select(target = K, enforce_dunn = FALSE):
       semilla = mejor candidato; añade greedy el par con mayor Dunn modificado disponible
       hasta llegar a K pares (NO se detiene en Dunn < 1, a diferencia del observado)
6. graduación:
       n_below = nº de pares con mod_dunn < 1
       n_below == 0  -> Tier 1   (acepta al pool)
       n_below == 1  -> Tier 2   (solo para rellenar déficit de Tier 1)
       n_below >= 2  -> rechazo (tier 0)
       no se pueden formar K pares no solapados -> rechazo
```

El *pool* se llena con Tier 1; Tier 2 solo entra si Tier 1 se agota. Escalado de presupuesto
según la tasa de aceptación Tier 1 (`MIN_VIABLE_TIER1_RATE = 0.05`, tope duro
`HARVEST_HARD_CAP = 50 * pool_size`).

### 3.2 Tabla comparativa

| Eje | RERconverge 0.3.0 categórico (§1) | RERconverge *legacy* `simBinPhenoCC` | PhyloPhere `permulations.R` |
|---|---|---|---|
| Qué se simula | historia de estados discretos (CTMC Mk) | latente continuo (BM del `pathvec`) | latente continuo (BM, o BM sobre árbol OU-reescalado) |
| Modelo del null | Mk **ER**, raíz = estacionaria; ajuste ML sobre el 0/1 observado | BM, `ratematrix` estimada del vector de paths foreground real | **BM u OU**, seleccionado por AIC sobre el rasgo observado (`OU sii ΔAIC>2`) |
| Reescalado temporal | no (Q absorbe la tasa) | no | **sí**: `geiger::rescale(tree,"OU",alpha)` cuando gana OU |
| De simulación a foreground | emerge del historia Mk (nodo hijo en estado 2) | top-`K` valores del latente | emerge de PSS + Dunn sobre `pvec` |
| Preserva **exactamente** | nº de puntas foreground; multiset de estados de nodos internos | nº de aristas foreground (`blsum`); orden de profundidad del subárbol foreground | **el multiset completo de valores observados del rasgo** (rank-match) |
| Condición de rechazo | conteos de puntas = observados (± `percent_relax`; aquí 0 → exacto) | `blsum == fgnum` **y** `setequal` del orden de profundidad | se pueden formar K pares **y** ≤ 1 par con Dunn < 1 (Tier 1/2) |
| ¿Condiciona sobre la **unidad del estadístico**? | **No** — RER es rama-a-rama, el path pesa todas las ramas foreground | **Parcial** — el orden de profundidad fija el nº de linajes foreground independientes y su anidamiento | **Sí** — K pares de contraste Dunn-independientes = la unidad exacta del estadístico CAAS |
| Distribución marginal del rasgo nulo | binaria por construcción; conteo fijo | binaria; conteo fijo | **idéntica a la observada** por construcción (rank-match) |
| Coste | 1 ajuste `fit_mk` + rejection de `getNullTips` (barato si el conteo fg no es extremo) + SA por árbol | rejection sampling doble (`blsum` + estructura); puede ser caro | rejection sampling condicionado a Tier 1; tasa 5–90% según geometría Dunn del rasgo |

### 3.3 Análisis de los criterios de aceptación

**Los tres condicionan el null sobre una característica de diseño, pero distinta en cada caso.**

- RERconverge categórico condiciona sobre el **conteo marginal de foregrounds**. Es el
  nuisance mínimo: si el nº de ramas foreground variara entre draws, la varianza del
  estadístico de correlación cambiaría y el p.perm quedaría mal calibrado. `getNullTips`
  con `percent_relax = 0` lo fija exactamente.
- `simBinPhenoCC` condiciona además sobre el **nº de orígenes independientes** (vía el orden
  de profundidad del subárbol foreground). Esto es *más* restrictivo: preserva no solo
  cuántas ramas son foreground sino cómo se reparten en linajes. Los mantenedores de
  RERconverge lo relajaron a propósito al pasar al motor categórico (no necesita outgroup,
  admite n estados, es más simple).
- PhyloPhere condiciona sobre la **existencia de K pares de contraste filogenéticamente
  independientes** (Dunn ≥ 1). Ésta es la unidad sobre la que se construye el estadístico
  CAAS: la selección observada produjo K pares independientes, el estadístico de
  convergencia se calcula *sobre esos K pares*, y por tanto el null tiene que replicar esa
  misma condición. Si no lo hiciera, el p-valor compararía "estadístico observado | K pares
  independientes" contra "estadístico nulo | estructura arbitraria", que es una comparación
  mal planteada.

> **Argumento mecanicista (no heurístico):** el condicionamiento Dunn de PhyloPhere y el
> condicionamiento de conteo de `getNullTips` son *el mismo tipo de operación* — restringir
> el soporte del null a las etiquetas que admiten el mismo diseño experimental que produjo
> el estadístico observado. Difieren en *qué* diseño: RER no tiene "pares", su diseño es
> "N ramas foreground sobre el árbol"; CAAS tiene "K pares independientes". Cada null
> condiciona correctamente sobre su propio diseño.

**¿El condicionamiento Tier 1/2 sesga la calibración?**

- El pseudo-conteo `(x·N + 1)/(N + 1)` usa `N = pool_size` (nº de ciclos aceptados), no
  `total_draws`. El p.perm mínimo detectable es `1/(pool_size + 1)`, **independiente de la
  tasa de aceptación**. El rejection sampling encarece el harvest pero no degrada la
  resolución ni la calibración, *siempre que el pool aceptado sea una muestra insesgada de*
  `P(etiqueta | Tier 1)`. Lo es por construcción: cada draw es un `simulatevec`
  independiente, aceptado sii cae en Tier 1.
- **No es importance sampling** (no hay reponderación). Es **rejection sampling de la
  distribución condicional** `P(etiqueta | K pares Dunn-independientes)`. La nota de diseño
  `perms_lambda_conditional_null` (aceptación FGBG 92% > lambda 36% > BM 15%) trata de
  elegir una *distribución propuesta* mejor (una que reproduzca la dispersión filogenética
  λ̂ del split observado) para subir la aceptación **sin cambiar la distribución objetivo**
  (que sigue siendo el filtrado a Tier 1). Propuesta y objetivo son ortogonales: λ elige la
  propuesta, el gate Tier define el objetivo.
- **Sesgo real, menor:** cuando Tier 1 se agota y entran registros Tier 2 (exactamente un
  par con Dunn < 1), el null pasa a ser una mezcla "todos independientes" + "uno casi
  dependiente". Un par semi-dependiente comparte ancestro y puede *inflar* el recuento de
  convergencia por herencia compartida → el estadístico nulo medio sube → **p.perm
  ligeramente conservador** (más difícil rechazar). Solo ocurre con rasgos cuya geometría
  Dunn no llena el pool desde Tier 1. Merece un aviso en el log (ya lo emite:
  `"Tier 1 exhausted … topping up with Tier 2"`).

---

## 4. Recomendación

**El RER binario de PhyloPhere debe usar `categoricalPermulations` (vía
`getPermsBinary(permmode="cc")`) tal cual, SIN añadirle el condicionamiento Dunn/Tier del
null de contrastes.**

Justificación:

1. **Las dos pruebas tienen unidades distintas.** RER correlaciona la tasa evolutiva
   relativa *por rama* contra un indicador foreground *por rama* (`foreground2Paths` +
   `correlateWithBinaryPhenotype`). No hay ninguna estructura de "K pares independientes"
   en el estadístico RER que preservar. CAAS sí: su estadístico se define sobre K pares de
   contraste. El condicionamiento Dunn preserva la unidad del estadístico CAAS; para RER
   sería condicionar sobre una característica (`K` pares de especies) que no interviene en
   el cálculo.
2. **`categoricalPermulations` ya preserva el nuisance relevante para RER** — el número de
   puntas foreground (`percent_relax = 0`) — y genera foregrounds filogenéticamente
   plausibles mediante el modelo Mk ajustado. Es el nuisance que afecta a la varianza de la
   correlación rama-a-rama.
3. **"Coherencia con el resto del pipeline" significa que el null replica el mismo diseño
   que usó el estadístico observado, no que todos los módulos compartan el mismo generador.**
   El diseño observado del RER es `foreground2Paths(clade="all")` + Kendall; su null replica
   exactamente eso. El diseño observado de CAAS es PSS + Dunn; su null replica eso. Ambos
   son coherentes *cada uno con su propio observado*.
4. **Añadir el gate Dunn al RER haría su `p.perm` incomparable con la literatura RERconverge**
   (Kowalczyk et al. 2019; Saputra et al. 2021), que es el marco de referencia con el que se
   contrastarán los resultados.

**Nivel de confianza:**

- *Alto* en "no añadir el condicionamiento Dunn al RER binario": es un argumento mecanicista
  sobre la unidad del estadístico, no una heurística.
- *Moderado* en "`categoricalPermulations` tal cual es suficiente": es la recomendación
  actual de los mantenedores de RERconverge y cubre el nuisance del RER, pero **relaja**
  respecto al `simBinPhenoCC` legacy en un punto (no preserva la estructura de linajes
  independientes del foreground, solo el conteo). Esa relajación no está validada
  empíricamente para el árbol de mamíferos concreto de este proyecto. **Chequeo barato que
  la zanjaría:** correr `getPermsBinary(cc)` sobre un rasgo binario neutral simulado y
  verificar que la distribución de `p.perm` es aprox. uniforme (calibración) antes de
  confiar en los resultados de producción.

**Precondición operativa:** reinstalar el entorno de la cuádrícula desde
`install_env.sh` (commit `2bd328f7`, que trae `categoricalPermulations`). El
`rerconverge=0.3.0` de `phylophere.yml` no basta y debería eliminarse o alinearse con el
pin de `install_env.sh` para evitar esta confusión en el futuro.

**Detalle de implementación ya resuelto:** `getPermsBinary` puntúa sus nulls con
`correlateWithBinaryPhenotype` en valores por defecto (sin winsorización, `clade="all"`,
`min.sp=10`, `min.pos=2`); el wrapper de [`binary_rer.R`](../subworkflows/RERCONVERGE/local/binary_rer.R)
calcula una correlación observada de referencia con esos mismos ajustes para `permpvalcor`,
manteniendo la `res` reportada con el `rer_binary_clade` y `winsorizeRER` del usuario.

---

## Anexo — pseudocódigo compacto de `categoricalPermulations` 0.3.0

```
FUNCTION categoricalPermulations(treesObj, phenvals[0/1], rm="ER", rp="auto", ntrees, percent_relax=0):
    tree  <- prune(treesObj$masterTree to names(phenvals))
    Q     <- fit_mk(tree, states(phenvals), rate_model="ER", root_prior=rp).transition_matrix

    # (a) puntas nulas: CTMC hacia adelante, aceptar solo si el conteo por estado coincide
    null_tips[1..ntrees], sim_node_states[1..ntrees] <- []
    WHILE accepted < ntrees:
        sim <- simulate_mk_model(tree, Q, root="stationary")     # castor
        IF all states present in sim.tips AND |count(sim.tips) - count(obs.tips)| <= count(obs.tips)*percent_relax:
            accept (sim.tip_states, sim.node_states)

    # (b) estados internos de arranque: barajar el multiset MAP observado, colocar por verosimilitud
    obs_node_MAP <- argmax_state( getAncLiks(tree, obs.tips, Q) )
    FOR i in 1..ntrees:
        anc_i <- getAncLiks(tree, null_tips[i], Q)
        start_nodes[i] <- placeByLikelihood( shuffle(obs_node_MAP), anc_i )

    # (c) pulido: recocido simulado sobre PERMUTACIONES de start_nodes[i]
    FOR i in 1..ntrees:
        improved_nodes[i] <- improveTree(tree, Q, P=expm(Q*bl),
                                         nodes=start_nodes[i], tips=null_tips[i],
                                         T0=10, Nk=10, cycles=100, alpha=0.9)
        # swaps s1<->s2 entre nodos, acepta si sube prod_e P[e][anc,desc], si no Metropolis con Tk=T0/(1+alpha*k)

    RETURN trees = [ (tips=null_tips[i], nodes=improved_nodes[i]) for i in 1..ntrees ]

FUNCTION getPermsBinary(..., permmode="cc"):
    phenotypeVector <- 1[species in fg_vec] over ALL masterTree tips
    hist <- categoricalPermulations(trees, phenotypeVector, "ER", "auto", numperms)
    FOR each hist.trees[i]:
        edge_len <- (c(tips, nodes)[child_of_edge]) - 1           # 1 = rama foreground
        fg_null  <- tips whose terminal edge == 1
        path     <- foreground2Paths(fg_null, trees, clade="all")  # length == ncol(RERmat)
        cor_i    <- correlateWithBinaryPhenotype(path, RERmat)     # defaults: Kendall, sin winsor
    RETURN list(corP, corRho, corStat)   # nrow = genes, ncol = numperms
```
