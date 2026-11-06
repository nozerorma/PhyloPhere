# Plan: nulo génico sobre todas las posiciones puntuables, U solo para posiciones

Estado: implementado en el repo local (reductor, tablas, merge, report, parámetro y GUI), sin commit y sin ejecutar en el cluster. Quedan los pasos de validación 3 a 5 de la sección 6, que exigen regenerar tiles.

## 1. Problema

El nulo de permulación solo puntúa las posiciones del universo posicional U (positions que `b_0` probó y que pasan el test de CAAS bajo alguna orientación de los pares observados). El observado cae dentro de U por construcción; un ciclo nulo, con otras especies, solo conserva las detecciones que coinciden con U.

Consecuencia medida en Carn_complete (1000 ciclos):

| Medida | `b_0` | Nulo (mediana) |
|---|---|---|
| Detecciones crudas por unidad (`n_pos`) | 1411 | 1391 |
| Posiciones puntuadas | 996 | 95 |
| Genes con puntuación (global) | 862 | 91 |

Los rankings génicos del nulo no son extracciones del mismo proceso que el ranking observado. Afecta a:

- el FCS de CAAS (`gene_stat = max`): 477 de 2964 sets con q.perm < 0,15, 20 % con p.perm < 0,05;
- el p génico (`gene_p`, `scoring_evidence.R`): el máximo nulo se toma sobre menos posiciones que el observado.

No afecta al p.emp por posición, a la familia del BH ni a POSENRICH: ahí condicionar en U es correcto.

## 2. Objetivo y criterio de éxito

Que el estadístico génico (máximo de las puntuaciones de posición del gen, por etiquetado) se calcule igual en `b_0` y en cada ciclo nulo, sobre todas las posiciones puntuadas del etiquetado, de modo que:

1. `p_gen = (1 + #{nulo >= obs}) / (n_valid + 1)` sea un p de permulación calibrado.
2. El FCS de CAAS compare rankings intercambiables.

"p_gen adecuado" se verifica con el control de la sección 6, no se da por supuesto.

## 3. Diseño

Se separan dos usos de la misma unidad (gen x ciclo).

| Pieza | Contenido | Universo | Lo usan |
|---|---|---|---|
| `units`, columnas actuales (`n_rows`, `n_top/bottom/all`, `max_*`) | sin cambios | U | cortes de `scores`, muestra y cuantiles del nulo |
| `units`, columnas nuevas (`gn_rows`, `gn_top/bottom/all`, `gmax_top/bottom/all`) | conteos y máximos sobre todas las filas tras quitar trenes | todas | matrices génicas (FCS, p génico), conteos totales por ciclo |
| `scores` | posición y puntuación exacta (float64) | U | p.emp, BH, POSENRICH |
| `scores_diag` (nueva) | gen, ciclo, bin de 0,01 de `top` y de `bottom` (uint8), sin posición | fuera de U | plots diagnósticos del report de scoring |

- El máximo génico se guarda exacto: el FCS compara rangos con empates en cero y un bin los alteraría.
- `scores_diag` guarda solo las filas fuera de U de los primeros `caas_null_diag_cycles` ciclos (por defecto 1000, o todos si hay menos). Los ciclos son intercambiables, así que el subconjunto vale para diagnóstico.
- `b_0` pasa por el mismo reductor: sus valores observados no cambian (todas sus posiciones están en U).
- No cambian: `n_pos`, las unidades que se eliminan por gen, el universo, `background.output`, el contrato de p.emp.

## 4. Pasos de implementación

Todo en el repo local; la sincronización con el cluster la hace el usuario.

1. `CT_DISAMBIGUATION/local/src/utils/null_tables.py`
   - `CycleReducer`: calcular, además de lo actual, los conteos y máximos sobre todas las filas que quedan tras quitar trenes, y los bins de las filas fuera de U si el ciclo está dentro de `caas_null_diag_cycles`.
   - `UNITS_SCHEMA`: añadir las columnas `g*`. Nuevo esquema y fichero `scores_diag.parquet`. Lo mismo para el pooling estricto (`universe_strict`).
2. `CT_DISAMBIGUATION/local/src/utils/gene_wrapper.py` y `disambiguation_perms_main.py`: pasar `caas_null_diag_cycles`; incluirlo en los metadatos del tile.
3. `CT_DISAMBIGUATION/local/src/utils/null_merge.py`
   - `gene_statistics`: leer `gn_*` y `gmax_*`.
   - Aplicar la eliminación de genes a `scores_diag` igual que `kept_scores`.
   - Escribir `cycle_summary_total.tsv` y `cycle_hist_total.tsv` sin alterar los ficheros actuales.
   - `sample_and_quantiles` y `cycle_summaries` siguen con las columnas de U.
4. `CT/caas_permulation.nf` y `conf/ct.config`: parámetro `caas_null_diag_cycles`; clave de caché de `CAAS_CORE_BATCHED`.
5. `SCORING/local/10.Scoring_report.Rmd` y `scoring_report.nf`: leer los ficheros totales, superponer el nulo total y el condicional a U, y rotular cada uno.
6. GUI (`gui/models`, `gui/templates`): exponer `caas_null_diag_cycles`.
7. Tests en `PhyloPhere_validation` (no en el repo): ver sección 6.

## 5. Qué hay que regenerar

Los tiles del nulo y de `b_0` (cambian las columnas de `units`). Es la etapa cara (`CAAS_CORE_BATCHED`). No se repite el ASR por gen, solo la reducción por unidad, pero el proceso entero se vuelve a lanzar por la clave de caché. Se recomienda ejecutarlo primero en la toy de Carn.

## 6. Validación (puntos de parada)

Cada uno decide si se continúa.

1. **Prueba unitaria del reductor.** Con `universe = None`: `gn_* == n_*` y `gmax_* == max_*`. Con U: `gn_* >= n_*`, `gmax_* >= max_*`. En `b_0`: columnas `g*` iguales a las actuales.
2. **Prueba unitaria del merge.** Que el corte de `scores` por unidad siga cuadrando y que los ficheros actuales salgan idénticos.
3. **Predicción falsable en la toy de Carn** (por `sbatch`, no en el login node). La hipótesis es que el nulo pierde detecciones por restringirse a U. Se espera que el número de genes con puntuación por ciclo nulo suba hacia el nivel del observado (en Carn completo, de ~91 hacia el orden de las 1000 posiciones puntuables). Si sigue en ~100, la hipótesis es falsa y se para el plan.
4. **Calibración de `p_gen`** (control con ciclo nulo como observado, ahora factible porque habrá ciclos de nivel comparable). Para varios ciclos nulos, usar cada uno como observado contra los demás y comprobar:
   - la distribución del p génico entre los genes puntuados es compatible con la uniforme;
   - la fracción con p <= 0,05 es próxima a 5 %;
   - no hay exceso sistemático de genes en el suelo (1/(N+1)).
5. **Observado.** Comparar el número de genes con p <= 0,01 y 0,05 frente al esperado, y el FCS de CAAS frente al control del paso 4.

Nota: con N = 1000 el BH génico no puede dar aciertos salvo que haya muchos genes en el suelo del p (del orden de 240 para q < 0,05). Un p_gen adecuado significa calibrado, no necesariamente con aciertos.

## 7. Efectos esperados

- Observado: sin cambios.
- Nulo génico: el máximo sube; el p génico se vuelve más conservador respecto al actual.
- FCS de CAAS: deberían desaparecer la mayor parte de los significativos que vienen de la diferencia de nivel; la magnitud no está medida.
- p.emp por posición, BH y POSENRICH: sin cambios.
- Plots del report: los condicionales a U se mantienen; se añaden los totales.
- Tamaño: `units` gana columnas; `scores` no cambia; `scores_diag` suma unas pocas decenas de MB en Carn y del orden de GB en cáncer con todos los ciclos (por eso se acota con `caas_null_diag_cycles`).

## 8. Riesgos y supuestos sin comprobar

- No se ha medido cuánto sube el nulo génico; la predicción del paso 3 es una hipótesis.
- Quedan otras asimetrías entre `b_0` y el nulo que este plan no corrige: distancia filogenética de los pares (H1 observada 16,5 frente a mediana nula 10,3) y el filtro de genes por longitud, que hoy elimina menos en el nulo. El p génico puede seguir con exceso si alguna pesa.
- No se ha leído el resto de `10.Scoring_report.Rmd` más allá de la parte que consume el resumen del nulo.
- Si en el futuro el estadístico génico pasa a `F(max)` o `size_adj_max`, harán falta puntuaciones exactas (no bins) para el ECDF de cada ciclo.
