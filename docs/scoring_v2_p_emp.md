# `p.emp`: p empírico de posición «detecta y supera»

`p.emp` es el p-valor de posición del nulo de permulación CAAS. Combina en un
único estadístico la **detección** de una posición y la **magnitud** de su
`CAAS_score`, sobre el mismo remuestreo (relabelado del fenotipo) que usan los
p-valores de gen. Sus dos ajustes por comparaciones múltiples, `p.adj_bh` y
`p.adj_sam`, se escriben junto a él en `scoring/position_scores.tsv`.

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

donde `CAAS_score` de un lado es la media §2g sobre los esquemas detectados
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
- **Las dos filas de lado** de una posición detectada en ambos lados llevan el
  mismo `p.emp`, `p.adj_bh` y `p.adj_sam`.

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

### 3.2 `p.adj_sam`

FDR por permutación (Tusher et al. 2001) sobre el score agrupado. Para un
umbral `t`:

```
FDR(t) = [ #pares (ciclo, posición) nulos con score >= t ] / N
         / #{ posiciones observadas con score >= t }
```

con π0 = 1. El valor de una posición es el mínimo de `FDR(t)` sobre los
umbrales `t` iguales o inferiores a su score. El conteo nulo esperado recorre
todas las posiciones que detecta cada ciclo, incluidas las que el observado no
detecta, así que no necesita definir una familia.

### 3.3 Diferencias entre ambos

`p.adj_bh` ordena por el nulo propio de cada posición (`p.emp`); `p.adj_sam`
ordena por el score observado frente a la distribución nula de todo el run.
Una posición con score alto cuyo nulo lo alcanza a menudo sale mejor en
`p.adj_sam`; una posición con score bajo que su nulo casi nunca alcanza, mejor
en `p.adj_bh`. En PEPC (§5) ambos coinciden en el núcleo de llamadas y
difieren en dos o tres posiciones por run.

### 3.4 Umbral

`scoring_p_emp_thr` (por defecto 0.05, `conf/scoring.config`) se aplica a los
dos ajustes en los informes. `position_significant` en los informes 15 y 16
sigue a `p.adj_bh`; `position_significant_sam` acompaña a SAM. El mismo
parámetro controla `gene_caas_pperm_adj`. Ninguno filtra `position_scores.tsv`.

---

## 4. Relación con los p-valores de gen

| | detección / conteo | detección + score |
|---|---|---|
| **posición** | | `p.emp` / `p.adj_bh` / `p.adj_sam` |
| **gen** | `accum_cct_p` / `accum_fdr` (§4b) | `gene_caas_pperm(_adj)` (§4f) |

- **`gene_caas_pperm` (§4f).** Estadístico `size_adj_max = F(max)^n` sobre los
  `CAAS_score` de las posiciones del gen, contra su fila en
  `caas_corStat_byrank` (`scoring_caas_perms.R`); add-one, cola derecha, BH
  dentro de cada dirección. No tiene compuerta de detección: la no detección
  entra como ceros estructurales de `gene_cycle_scores.tsv`. El guard numérico
  de §4 mantiene `size_adj_max` (R) y `_size_adj_max_null` (Python) idénticos.
- **`accum_cct_p` (§4b).** Conteo de posiciones detectadas por gen y esquema,
  combinado con Cauchy (CCT/ACAT). Con `accumulation_randomization_type =
  "cons_decile"` (por defecto) su nulo es de ocupación por deciles de
  conservación y no controla la relación árbol-fenotipo; con `"permulation"`
  lee el mismo `perm_pos_detail/` que §4f.
- **§4f y §4b son dos ejes y no se combinan.** Responden a preguntas distintas
  (magnitud de la mejor posición frente a número de posiciones) y, en el modo
  por defecto, con nulos distintos.
- **`p.emp` y `gene_caas_pperm` no son evidencia independiente para un mismo
  locus.** `gene_caas_pperm` es función de los mismos scores nulos de posición
  que alimentan `p.emp`; un gen significativo y una posición significativa
  dentro de él son una sola línea de evidencia a dos resoluciones.

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
- **Calibración de `p.adj_sam`.** Tomando cada ciclo nulo como observado frente
  a los 999 restantes, P(≥ 1 llamada) a q < 0.1 es 0.104 (genotípico) y 0.097
  (fenotípico). BH sobre la familia a 0.1 da 0.011 y 0.067: con `N = 1000` el
  suelo `1/(N + 1)` limita la potencia de BH.
- **Controles negativos.** Veinte rasgos generados con el propio generador del
  nulo y pasados por el camino observado: el número de detecciones es el de un
  ciclo nulo más (mediana P = 0.40); las ejecuciones con alguna llamada son
  1/20 a `p.adj_bh < 0.05`, 4/20 a `p.adj_bh < 0.1` y 1/20 a `p.adj_sam < 0.1`.
  Sobre las familias agrupadas, P(p ≤ 0.01 / 0.05 / 0.10) = 0.028 / 0.069 /
  0.093: exceso en la cola extrema, repartido entre controles, sin explicar
  por las posiciones de baja detección ni por el modelo (OU/BM) del nulo.

---

## 6. Limitaciones abiertas

- **Familia de BH (§3.1).** m depende de qué columnas alcanzó la muestra finita
  del nulo. En los controles negativos el ritmo de falsos positivos por
  ejecución es nominal, pero el mecanismo existe.
- **Exceso en la cola de `p.emp`** en los controles negativos (§5), de causa no
  identificada.
- **Resolución.** El mínimo de `p.emp` es `1/(N + 1)`; con familias de decenas
  de posiciones, una llamada BH a 0.05 exige estar en ese suelo o cerca, y los
  rangos entre las posiciones más fuertes dependen de unos pocos ciclos.
- **Tamaño de `perm_pos_cycle_caas.tsv.gz`** a escala genómica: una fila por
  `(Gene, Position, side, cycle)` detectado; sin medir fuera de Tier 1.

---

## Referencias

- Bourgon R, Gentleman R, Huber W. 2010. Independent filtering increases detection power for high-throughput experiments. Proc Natl Acad Sci USA 107(21):9546–9551. doi:10.1073/pnas.0914005107.
- Davison AC, Hinkley DV. 1997. Bootstrap Methods and Their Application. Cambridge University Press.
- Tusher VG, Tibshirani R, Chu G. 2001. Significance analysis of microarrays applied to the ionizing radiation response. Proc Natl Acad Sci USA 98(9):5116–5121. doi:10.1073/pnas.091062498.
