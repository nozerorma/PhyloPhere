# WP9. Estadístico base del p de gen

Estado: A1 hecho (2026-10-21); control en el FCS y A2 pendientes. Puede ejecutarse en paralelo con las demás fases; A2 depende de WP5.
Conversación independiente. Antes de empezar, leer `docs/fases/README.md`.

## Objetivo

Decidir con qué estadístico se construye el p de permulación del gen que sirve de ranking al FCS de CAAS: máximo crudo de posición, `F(max)` o `F(max)^n`, y medir el peso específico de `^n`. Ya está decidido que el ranking es −log10 de ese p; no se confirma `F(max)^n` frente al score crudo.

## Diseño del experimento

- Variantes por (gen, ciclo, dirección): máximo crudo; `F(max)^(n^β)` con β ∈ {0; 0,25; 0,5; 0,75; 1} (β = 0 es `F(max)`, β = 1 es `F(max)^n`). `F` es la ECDF de los scores de posición del mismo ciclo y dirección (con `TIE_TOL` = 1e-12) y n el número de posiciones puntuadas del gen.
- Métricas: uniformidad del score dentro del ciclo (KS), Spearman con n y con la propensión del gen, ICC del efecto de gen; p de gen contra los otros ciclos (genes con p ≤ 0,01 y 0,05 por ciclo nulo y percentil del observado entre ciclos nulos de nivel parecido, 0,8 a 1,5 veces los genes puntuados del observado); sensibilidad al nivel; fracción de genes cuyo p cambia más de una unidad en −log10 p al pasar de β = 0 a β = 1; control "ciclo nulo como observado" del FCS con la regla nueva.
- Criterio de elección: `F(max)^n` solo si su calibración no es peor que `F(max)` y la dependencia de n que quita es relevante; si no, `F(max)`; máximo crudo solo si está mejor calibrado. Las tres variantes las calcula el merge de WP5, de modo que la elección es un parámetro (`gene_stat`), y la estampa de `fcs_enrich.R` y `scoring_compute.R` (hoy `size_adj_max`) debe seguirla.

## A1: resultado (neoplasia estricto, 1000 ciclos nulos, 208 de nivel parecido)

Datos y código: `/home/miguel/IBE-UPF/PhD/PhyloPhere_validation/golden_581b999/a1_gene_statistic/` (`a1_variants.R`, `a1_variants.out`, `a1_variants_out.rds`); entradas locales `~/.cc/calib_ctrl/pos.rds` (posiciones con score de los 1000 ciclos nulos) y `b0_pos.tsv.gz` (observado). 11.344 genes; el observado puntúa 2480 (73 % con n = 1).

| Variante | KS | ρ con n | ρ con propensión | ICC | p ≤ 0,01: observado / nulo / % nulos ≥ obs | p ≤ 0,05: observado / nulo / % nulos ≥ obs | ρ nivel |
|---|---|---|---|---|---|---|---|
| máximo crudo | | 0,229 | 0,406 | 0,080 | 309 / 194,2 / 19 % | 1129 / 844,3 / 15 % | 0,853 |
| β = 0 (`F(max)`) | 0,060 | 0,229 | 0,377 | 0,098 | 234 / 224,0 / 36 % | 922 / 890,0 / 35 % | 0,989 |
| β = 0,25 | 0,049 | 0,171 | 0,340 | 0,097 | 229 / 221,4 / 37 % | 918 / 881,4 / 34 % | 0,989 |
| β = 0,5 | 0,038 | 0,108 | 0,295 | 0,096 | 222 / 218,9 / 38 % | 912 / 872,7 / 33 % | 0,989 |
| β = 0,75 | 0,026 | 0,041 | 0,244 | 0,096 | 218 / 216,9 / 39 % | 900 / 864,7 / 34 % | 0,989 |
| β = 1 (`F(max)^n`) | 0,017 | −0,026 | 0,189 | 0,097 | 216 / 215,0 / 38 % | 894 / 857,1 / 34 % | 0,989 |

Lectura (medida, una sola corrida):

- **El exceso de magnitud del observado lo quita la normalización por el pool del ciclo (`F`), no el exponente.** Con el máximo crudo el observado queda en el 15 al 19 % superior de los nulos de nivel parecido; con cualquier variante de `F`, incluida β = 0, queda en el 33 al 39 %.
- **Peso de `^n`:** mejora la uniformidad (KS 0,060 → 0,017) y quita la dependencia de n (ρ 0,229 → −0,026) y parte del acople a la propensión (0,377 → 0,189), sin cambiar el ICC. Sobre el p del gen observado casi no se nota: ningún gen de los 2480 puntuados cambia más de una unidad en −log10 p entre `F(max)` y `F(max)^n` (tampoco entre los de n > 1; solo 4 genes tienen n > 8); la correlación de Spearman de −log10 p es 0,971 entre `F(max)` y `F(max)^n`, 0,992 entre crudo y `F(max)` y 0,959 entre crudo y `F(max)^n`.
- **Nivel:** en los ciclos nulos, el número de genes con p ≤ 0,05 sigue casi por completo al número de genes puntuados (ρ 0,989 con `F`, 0,853 con el crudo): cualquier comparación con el nulo exige emparejar por nivel.
- **Límites:** una corrida (neoplasia estricto); falta el control a nivel del FCS (métrica 5) y la repetición a 10.000 ciclos (A2), que es la que cambia la resolución del p (suelo 1/1001).

## Pendiente

1. **Control en el FCS** ("ciclo nulo como observado", con la regla de WP1) para el candidato elegido y para el máximo crudo; hace falta construir el ranking con `F(max)` (hoy hay scripts para el crudo, `adhoc_scripts/fcs_rawmax_ranking_control.sbatch`, y para `F(max)^n`, `fcs_calibrated_ranking_control.sbatch`; el de `F(max)` no existe). Repetir en las otras tres corridas (neoplasia completo, maligno estricto y completo) si el resultado en neoplasia estricto no es concluyente.
2. **A2 (tras WP5):** repetir a 10.000 ciclos con las tablas de tesela.
3. Fijar el parámetro `gene_stat` y la estampa; documentar la decisión.

## Registro de la fase

```
- Conversación y fechas: A1 en la conversación del 2026-10-21 (plan y mediciones)
- Commit de partida / commit de cierre: 581b999 / pendiente
- Qué se hizo (ficheros y cambios): análisis local sin cambios de código (script en /home/miguel/IBE-UPF/PhD/PhyloPhere_validation/golden_581b999/a1_gene_statistic/)
- Comprobaciones ejecutadas (comando, conjunto de datos, tolerancia, resultado): ver tabla de A1
- Mediciones (recurso antes y después): no aplica
- Desviaciones del plan y decisiones tomadas (con su justificación): ninguna todavía; la elección entre variantes queda abierta
- Pendiente que pasa a otra fase: control en el FCS y A2
```
