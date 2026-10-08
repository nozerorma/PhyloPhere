# `p.emp`: p empírico de posición «detecta y supera»

`p.emp` es el p-valor de posición del nulo de permulación CAAS. Combina en un
único estadístico la **detección** de una posición y la **magnitud** de su
`CAAS_score`, sobre el mismo remuestreo (relabelado del fenotipo) que usan los
p-valores de gen. Su ajuste por comparaciones múltiples, `p.adj_bh`, se escribe
junto a él en `scoring/position_scores.tsv`.

Código de referencia:
[`scoring_compute.R`](../subworkflows/SCORING/local/src/scoring_compute.R)
§2f-ter (cómputo) y §2h (ajustes);
[`gene_wrapper.py`](../subworkflows/CT_DISAMBIGUATION/local/src/utils/gene_wrapper.py)
`_finalize_perm_scores` (emisión del nulo);
[`reaggregate_perm_scores.py`](../subworkflows/CT_DISAMBIGUATION/local/reaggregate_perm_scores.py)
(reconstrucción del nulo sin ASR).

Coordenadas 0-based. `N` = número de ciclos del nulo.

---

## 1. Definición

Para cada posición `(Gene, Position)` se define el estadístico

```
T = max_lado CAAS_score   si la posición se detecta en algún lado
T = -Inf                  si no se detecta
```

donde `CAAS_score` de un lado es la agregación §2g sobre los esquemas (suma sobre los cinco esquemas, o media sobre los
detectados con `caas_score_aggregation = mean`; el nulo y el observado usan la misma)
(US, GS4, GS3, GS2, GS1). El score de la posición es el máximo sobre los lados
detectados (el eje "all", `.pos_undirected` de §4a).

Un ciclo nulo `i` cuenta a favor de la hipótesis nula si **(a)** re-detecta la
posición en cualquier lado **y** **(b)** su score agrupado, calculado con la
misma regla, iguala o supera el observado:

```
k_emp = #{ i : ciclo i detecta la posición  y  T_null(i) >= T_obs }
p.emp = (k_emp + 1) / (N + 1)          add-one, cola derecha
```

- **Sin corrección leave-one-out.** El estadístico observado procede del
  análisis real, no es uno de los `N` ciclos; no hay auto-inclusión que
  corregir. `(k_emp − 1)/(N − 1)` trataría al observado como un ciclo nulo.
- **Una posición que ningún ciclo re-detecta** tiene `k_emp = 0` y
  `p.emp = 1/(N + 1)`, el mínimo alcanzable.
- **Una posición con `T_obs = 0`** (detectada, sin convergencia: ningún par de dominios comparte
  residuo derivado) tiene `p.emp = 1`. Sin evidencia de convergencia el p no puede ser menor, y
  contado como «detecta y supera» solo mediría la frecuencia con que el nulo detecta la columna.
  Para `T_obs > 0` el resultado no cambia: un ciclo que no detecta la posición (−Inf, o 0 como
  en el nulo de gen) nunca iguala un observado positivo. Entra en la familia de BH con `p = 1`.
- **Las dos filas de lado** de una posición detectada en ambos lados llevan el
  mismo `p.emp` y `p.adj_bh`.

---

### 1.1 Convención de reporte

El mínimo de `p.emp` es `1/(N + 1)`. En tablas y gráficos de los informes un valor en
ese suelo se muestra «< 1/N» («<0.001» con 1000 ciclos), como en Saputra et al. 2021 (con
1000 permulaciones el mínimo reportable es 0.001 y un p calculado como 0 se escribe
«<0.001»). Los ficheros TSV conservan el valor numérico.

---

## 2. Cómputo

### 2.1 Emisión del nulo

`_finalize_perm_scores` escribe `perm_pos_cycle_caas.tsv.gz` con una fila por
`(Gene, Position, side, cycle)` detectado:

| columna | contenido |
|---|---|
| `caas_score` | score de posición del ciclo en ese lado (`core.scores`); vacío si ningún esquema lo puntuó |
| `n_schemes` | número de esquemas que puntuaron la posición |

`caas_score` es la media de los scores por esquema, con suma correctamente
redondeada (`math.fsum`), la misma función que da el `CAAS_score` observado. Se
calcula una sola vez, en Python; R y el enriquecimiento de posiciones lo leen sin
recalcularlo, y un fichero anterior sin esa columna se rechaza (hay que regenerar
el nulo). Con el espejo
FOP, las etiquetas `<ciclo>~H*` se colapsan al ciclo base antes de agregar.
`reaggregate_perm_scores.py` regenera el fichero desde los shards
`perm_pos_detail/` sin volver a la ASR.

El fichero llega a SCORING como `--caas_pos_cycle_caas` (desde
`CAAS_CORE_MERGE` en `main.nf`, o desde `params.caas_pos_cycle_caas_file` en
`workflows/scoring.nf` al reutilizar un nulo previo).

### 2.2 Consumo en R (§2f-ter)

1. Score nulo agrupado por `(Gene, Position, cycle)`: máximo sobre lados de
   `caas_score`.
2. Score observado agrupado por `(Gene, Position)`: máximo sobre lados de
   `CAAS_score`.
3. `k_emp` por posición mediante `left_join` desde lo observado, de modo que
   las posiciones no re-detectadas quedan con `k_emp = 0`.
   La comparación `>=` cuenta como empate lo que difiere menos de
   `TIE_TOL = 1e-12` (la misma constante que `core.scores`): las medias iguales en
   aritmética exacta pueden diferir en 1 ulp según cómo se sumaron o se leyeron.
4. `N` es el roster de ciclos de `caas_perms.rds` (incluye los ciclos que no
   detectan nada), siempre que contenga todos los ciclos presentes en el
   fichero; si no, `N` = ciclos presentes, lo que sólo puede aumentar `p.emp`.

**Guard de coordenadas.** Si hay al menos 10 posiciones observadas y menos de
la mitad aparecen en el nulo, se asume un desajuste de coordenadas entre
`filtered_discovery.tsv` y el nulo: las posiciones sin correspondencia quedan
en NA y se emite un aviso. Con menos de 10 posiciones la tasa no es
informativa y las posiciones sin correspondencia se puntúan con `k_emp = 0`.

### 2.3 `p.emp_fact`: p factorizado de posición

«Detecta y supera» se descompone exactamente en
`P(detecta) × P(score ≥ s | detecta)`. `p.emp_fact` estima el primer factor con los
ciclos de la propia posición y el segundo con las detecciones de las posiciones de su
clase, de modo que no queda limitado por el suelo `1/(N + 1)`:

```
p.emp_fact = (nd + 1) / (N + 1)  ×  (1 + #{detecciones de la clase con score >= s}) / (1 + #{detecciones de la clase})
```

- `nd` = ciclos del nulo que puntúan la posición (en cualquier lado). La **clase** es la
  propensión de la posición, cortada en percentiles de las detecciones del nulo: cada una de
  las `FACT_PROP_CLASSES` = 20 clases reúne en torno a 1/20 de todas las detecciones (pares
  posición-ciclo con score), y las posiciones con el mismo `nd` no se parten, de modo que una
  clase puede tener más. La posición que ningún ciclo puntúa va a la primera clase. Los cortes
  salen del propio nulo, sin valores fijados a mano. La posición observada cuenta como una
  detección más de sí misma, de modo que se clasifica por `nd + 1`, el número de detecciones
  que tiene una posición de su clase en la reserva. Agrupar posiciones de propensión distinta no es válido: las
  que el nulo puntúa a menudo también alcanzan scores altos más a menudo por azar.
- Un score observado 0 da `p.emp_fact = 1`. Donde el guard de coordenadas deja `p.emp` en
  NA, `p.emp_fact` también.
- `p.adj_bh_fact` es BH sobre la misma familia que `p.adj_bh` (§3.1), con las posiciones
  sólo-nulo en `p = 1`.
- **Supuestos y comprobación.** Los ciclos nulos son intercambiables con el observado, y
  las clases se definen con el propio nulo. Tomando cada ciclo nulo como observado frente
  a los demás, la tasa de pares (posición, ciclo) con `p.emp_fact ≤ α` no supera α en
  ninguna clase de propensión. La prueba está en `PhyloPhere_validation`
  (`test_factorized_p_calibration.py`).

---

## 3. Comparaciones múltiples (§2h)

### 3.1 `p.adj_bh`

BH con un test por `(Gene, Position)` sobre la familia

```
familia = posiciones detectadas por >= 1 ciclo nulo  ∪  posiciones observadas con p.emp
```

Las posiciones que el nulo detecta y el observado no entran con `p = 1`: el
estadístico es `-Inf` para ellas, y restringir BH a lo observado filtraría
sobre el propio estadístico. Las posiciones que no detectan ni el nulo ni el
observado quedan fuera.

**Limitación.** m cuenta sólo las columnas que la muestra finita del nulo
alcanzó. Una columna que el nulo nunca detecta pero el observado sí entra con
`p = 1/(N + 1)` en una familia que no incluye las demás columnas de baja tasa
de detección, con lo que m queda infraestimado. El control negativo nc13 de
Tier 1 (§5) es un ejemplo: una única posición observada, nunca re-detectada,
en una familia de 17 posiciones, da `p.adj_bh = 0.017`.

### 3.2 Sin FDR por permutación sobre el score

No se calcula un FDR por permutación sobre el score agrupado (la idea de SAM,
Tusher et al. 2001). Estima la fracción de falsas entre las llamadas por encima
de un corte de score con una distribución nula agrupada sobre todas las
posiciones, de modo que es un FDR de conjunto: las llamadas falsas se
concentran en las posiciones que el nulo detecta muchas veces, y su valor
depende de la elección entre media y mediana del número de falsos y de cómo se
fija el corte. Las afirmaciones por posición descansan en `p.emp`.

### 3.3 Uso de los p-valores

`p.emp`, `p.emp_fact` y sus ajustes `p.adj_bh` y `p.adj_bh_fact` son informativos: se reportan tal cual,
cada uno a su nivel, y no se deriva de ellos ningún umbral ni columna de significatividad. Los informes
muestran los cuatro. Ninguno filtra `position_scores.tsv`.

---

## 4. Nivel de gen

| | detección / conteo | detección + score |
|---|---|---|
| **posición** | | `p.emp` / `p.adj_bh`, `p.emp_fact` / `p.adj_bh_fact` |
| **gen** | `accum_cct_p` / `accum_fdr` (§4b) | (sin p de gen) |

- **Sin p de gen del CAAS.** El CAAS se contrasta a nivel de posición. El gen lleva
  `gene_caas_score` (`size_adj_max = F(max)^n` sobre los `CAAS_score` de sus posiciones, calculado una
  vez por `core.scores` y reproducido por `_size_adj_max_null` en el nulo), que es el ranking que entra al
  enriquecimiento (FCS) con su nulo por permulación en la misma escala (`caas_perms.rds`,
  `scoring_caas_perms.R`).
- **`accum_cct_p` (§4b).** Conteo de posiciones detectadas por gen y esquema,
  combinado con Cauchy (CCT/ACAT). Con `accumulation_randomization_type =
  "cons_decile"` (por defecto) su nulo es de ocupación por deciles de
  conservación y no controla la relación árbol-fenotipo; con `"permulation"`
  lee `perm_pos_detail/`.

---

## 5. Comportamiento observado (Tier 1 PEPC)

Resultados en `validation/tier1/reports/pepc_genotypic_vs_phenotypic.md` §6.

- **Sub-uniformidad en el conjunto candidato.** Para una posición nula que el
  observado detecta, `p.emp ≈ d × U`, con `d` la tasa de detección del nulo;
  `p.emp` no puede superar `d`. Filtrar por detección observada es filtrar
  sobre el estadístico (Bourgon et al. 2010): BH o Storey aplicados sólo a los
  candidatos son anticonservadores, y Storey sobre ese conjunto estima
  π0 ≈ 0.1 y declara significativos todos los candidatos. Por eso la familia de
  `p.adj_bh` incluye las posiciones sólo-nulo con `p = 1`.
- **Controles negativos.** Veinte rasgos generados con el propio generador del
  nulo y pasados por el camino observado: el número de detecciones es el de un
  ciclo nulo más (mediana P = 0.40); las ejecuciones con alguna llamada son
  1/20 a `p.adj_bh < 0.05` y 4/20 a `p.adj_bh < 0.1`.
  Sobre las familias agrupadas, P(p ≤ 0.01 / 0.05 / 0.10) = 0.028 / 0.069 /
  0.093: exceso en la cola extrema, repartido entre controles, sin explicar
  por las posiciones de baja detección ni por el modelo (OU/BM) del nulo.

---

## 6. Limitaciones abiertas

- **Familia de BH (§3.1).** m depende de qué columnas alcanzó la muestra finita
  del nulo. En los controles negativos el ritmo de falsos positivos por
  ejecución es nominal, pero el mecanismo existe.
- **Exceso en la cola de `p.emp`** en los controles negativos (§5), de causa no
  identificada. `p.emp_fact` lo hereda: en los controles de PEPC el número de posiciones con
  p ≤ 0.01 es de 11 a 13 frente a 5.2 esperadas como máximo (mismo orden que `p.emp`: 14).
- **Resolución.** El mínimo de `p.emp` es `1/(N + 1)`; con familias de decenas
  de posiciones, una llamada BH a 0.05 exige estar en ese suelo o cerca, y los
  rangos entre las posiciones más fuertes dependen de unos pocos ciclos.
- **Tamaño de `perm_pos_cycle_caas.tsv.gz`** a escala genómica: una fila por
  `(Gene, Position, side, cycle)` detectado; sin medir fuera de Tier 1.

---

## Referencias

- Saputra E, Kowalczyk A, Cusick L, Clark N, Chikina M. 2021. Phylogenetic permulations: a statistically rigorous approach to measure confidence in associations in a phylogenetic context. Mol Biol Evol 38(7):3004–3021. doi:10.1093/molbev/msab068.
- Bourgon R, Gentleman R, Huber W. 2010. Independent filtering increases detection power for high-throughput experiments. Proc Natl Acad Sci USA 107(21):9546–9551. doi:10.1073/pnas.0914005107.
- Davison AC, Hinkley DV. 1997. Bootstrap Methods and Their Application. Cambridge University Press.
