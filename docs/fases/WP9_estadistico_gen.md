# WP9. Estadístico base del p de gen

Estado vigente (2026-10-23): A1 histórico y control FCS dirigido local realizados; A2 a 10.000 y elección definitiva pendientes. Se conserva provisionalmente el default size_adj_max, sin declararlo decisión final. A2 depende de nulos reales de escala; WP5–WP7 ya están implementados.
Conversación independiente. Antes de empezar, leer `docs/fases/README.md`.

## Objetivo

Decidir con qué estadístico se construye el p de permulación del gen que sirve de ranking al FCS de CAAS: máximo crudo de posición, `F(max)` o `F(max)^n`, y medir el peso específico de `^n`. Ya está decidido que el ranking es −log10 de ese p; no se confirma `F(max)^n` frente al score crudo.

## Diseño del experimento

- Variantes por (gen, ciclo, dirección): máximo crudo; `F(max)^(n^β)` con β ∈ {0; 0,25; 0,5; 0,75; 1} (β = 0 es `F(max)`, β = 1 es `F(max)^n`). `F` es la ECDF de los scores de posición del mismo ciclo y dirección (con `TIE_TOL` = 1e-12) y n el número de posiciones puntuadas del gen.
- Métricas: uniformidad del score dentro del ciclo (KS), Spearman con n y con la propensión del gen, ICC del efecto de gen; p de gen contra los otros ciclos (genes con p ≤ 0,01 y 0,05 por ciclo nulo y percentil del observado entre ciclos nulos de nivel parecido, 0,8 a 1,5 veces los genes puntuados del observado); sensibilidad al nivel; fracción de genes cuyo p cambia más de una unidad en −log10 p al pasar de β = 0 a β = 1; control "ciclo nulo como observado" del FCS con la regla nueva.
- Criterio de elección: `F(max)^n` solo si su calibración no es peor que `F(max)` y la dependencia de n que quita es relevante; si no, `F(max)`; máximo crudo solo si está mejor calibrado. Las tres variantes las calcula el merge de WP5, de modo que la elección es un parámetro (`gene_stat`), y la estampa de `fcs_enrich.R` y `scoring_compute.R` (hoy `size_adj_max`) debe seguirla.

## A1 histórico: resultado (neoplasia estricto, 1000 ciclos nulos, 208 de nivel parecido)

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

## Pendiente vigente

1. Ampliar los controles a familias/direcciones completas y a neoplasia/maligno completo/strict. El control local nuevo incluye fmax; ya no falta ese camino.
2. **A2:** repetir con nulos reales de 10.000 ciclos, sin reciclar columnas. El merge ya conserva max/fmax/n; no se necesitan las tablas de posición globales retiradas para comparar estas variantes.
3. Revisar calibración condicionada al nivel y decidir gene_stat con el usuario tras A2. La estampa y las parejas observado/nulo ya están parametrizadas; no cambiar defaults sin una decisión científica explícita.

## Registro histórico de A1

```
- Conversación y fechas: A1 en la conversación del 2026-10-21 (plan y mediciones)
- Commit de partida / commit de cierre: 581b999 / pendiente
- Qué se hizo (ficheros y cambios): análisis local sin cambios de código (script en /home/miguel/IBE-UPF/PhD/PhyloPhere_validation/golden_581b999/a1_gene_statistic/)
- Comprobaciones ejecutadas (comando, conjunto de datos, tolerancia, resultado): ver tabla de A1
- Mediciones (recurso antes y después): no aplica
- Desviaciones del plan y decisiones tomadas (con su justificación): ninguna todavía; la elección entre variantes queda abierta
- Pendiente que pasa a otra fase: control en el FCS y A2
```


## Control FCS dirigido y preparación A2 (2026-10-23, vigente)

Continuación de WP8 en la misma conversación a petición del usuario. Partida `d6edd5b2`, rama remove_caap; cambios previos de WP4/validation/prompts preservados. Implementación del control y oráculos: `e8c0d483`; documentación en el commit siguiente. Este avance no cierra WP9 ni la aceptación integrada.

Resultados y límites completos: `WP9_resultados_dirigidos_2026-10-23.md`; resumen numérico versionado `WP9_holdout_2026-10-23.tsv`. Evidencia y snapshots estables: `PhyloPhere_validation/.tmp/wp9/`. Entorno phylophere, TMPDIR en IBE. Ningún trabajo ni edición de cluster; ningún nulo generado ni corrida completa a 10.000.

### Qué se hizo

- Nuevo control de validación `adhoc_scripts/wp9_gene_stat_control.R`: máximo, fmax, β=0,25/0,5/0,75 y size_adj_max, a partir de parts canónicas o posiciones históricas contrastadas con RDS de referencia. Ceros y roster íntegros, universo fijo del consumidor, validación de índices/metadatos/b0, límites explícitos de memoria, fuentes/MD5 y sessionInfo. No altera el pipeline ni sus parámetros científicos.
- Usa gene_p.R productivo (tolerancia 1e-12) y AUC centrado por rutas (tolerancia exactamente cero), BH por familia y Redlich independiente. Observado y pseudo-observados usan el mismo estadístico base por variante.
- Se evalúan todos los ciclos como pseudo-observados; además, dos mitades disjuntas de referencia/evaluación, intercambiadas, semilla 1998. El resumen incluye conteos válidos por base, dirección y nivel. La tasa marginal LOO cercana a 5 % está casi impuesta por los rangos: no es un certificado de calibración. El holdout tampoco convierte rutas/folds dependientes en réplicas independientes ni estima FDR.
- Se amplía `fcs_calibrated_ranking_control.sbatch` con WP9_CONTROL=true: un trabajo por combinación/ranking, snapshot nuevo de HEAD, todas las variantes/ciclos, EXPECTED_CYCLES y MAX_CELLS explícitos. Se conserva su modo antiguo. Ambos launchers de control calibrado/crudo pasan a max_g=0 por defecto y guardan q.perm/sig_perm/redlich_level, eliminando la selección de columnas retiradas. El FCS productivo permanece intacto.

### Verificación y medidas

- `pytest tests/test_wp9.py`: 2 PASS. El oráculo R calcula las variantes elemento a elemento, AUC con la función nativa, recuentos LOO/BH exhaustivos y holdout con calibración manual/AUC nativo. Cubre selección de universo, pools históricos antes de esa selección, roster con ciclos vacíos, b0, metadatos, límites e integridad. El test del launcher usa un Git mínimo y un ejecutor simulado: verifica snapshot, argumentos y rechazo de reutilización/arrays, sin enviar jobs.
- Marino completo y strict: 120 genes, 100 ciclos, tres rankings, 17 GMT leídos y 399 rutas por variante/modo. Los 399 resultados completos de size_adj_max coinciden con WP7 a 1e-12 en stat, p/q, claves y sig_perm (`compare_wp7.log`).
- Neoplasia estricto histórico: 12.473 genes del universo depurado, 1000 ciclos, global, Hallmark+KEGG fijados a priori (364 rutas). F(max)^n reconstruido sin redondeo adicional reproduce el RDS con delta máximo 2,22e-16. Son datos históricos, no teselas regeneradas con el checkout actual.
- El universo usado por A1 era el de todos los genes detectados del almacén; el FCS usa el universo depurado. La selección actual conserva los pools originales y los p individuales se calculan con la regla vigente. El grupo emparejado cambia de 208 a 137 ciclos; no se atribuye ese cambio al exponente. Los recuentos observados a p≤0,05 y la correlación F(max)/F(max)^n a tres decimales reproducen A1.
- Con holdout: marino completo 43,22 s y 616.564 KiB RSS; strict 43,38 s y 605.176 KiB; neoplasia global/two-GMT 93,47 s y 3.492.284 KiB (3,33 GiB). Son controles deliberadamente de memoria completa, no memoria productiva ni medidas de escala 10k.
- bash -n de ambos launchers y parse R del control y cuerpos heredoc: PASS. Ayuda CLI y smoke final de una variante/base comprobados. No se repite Nextflow: no cambian sus procesos, entradas, formato ni cálculo productivo.

### Lectura científica y recomendación provisional

En neoplasia, ρ(score,n) pasa de 0,225 a −0,020 y KS de 0,064 a 0,023 al pasar de F(max) a F(max)^n. La correlación entre p observados transformados es 0,9711; ningún gen cambia >1 unidad en −log10 p y solo cuatro observados tienen n>8. La tasa holdout a p≤0,05 en ciclos de nivel observado es 8,25 % frente a 8,26 % respectivamente; Hallmark ronda 10,3 % y KEGG 7,9 %. Es asociación con nivel en el subconjunto evaluado, no FDR ni prueba de desajuste biológico. Ninguna variante elimina esa asociación; no elegir por números de significativos.

Recomendación: conservar provisionalmente size_adj_max, que elimina la dependencia de n sin diferencias apreciables en este FCS, y mantener A2/decisión definitiva abiertos. No se cambia gene_stat ni se declara calibración no inferior demostrada por un margen preespecificado. No hay soporte para seleccionar un β intermedio ni máximo crudo.

### Lanzamiento A2 preparado (solo el usuario, después de sincronizar)

Partir de las corridas reales preparadas en `WP8_aceptacion_integrada.md`; STORE y UNIVERSE deben corresponder exactamente a la misma corrida/modo. El ejemplo reserva recursos para un control de memoria completa, no propone cambiar reservas del pipeline. Ajustar el presupuesto con las medidas disponibles.

```bash
export STORE="/ruta/a/la/corrida_10000/strict/caas_permulation"
export UNIVERSE="/ruta/a/la/corrida_10000/universes/cleaned_background_main_neoplasia_prevalence.txt"
export EXPECTED_CYCLES=10000 RANKINGS=global MAX_CELLS=200000000
sbatch --mem=64G --export=ALL,WP9_CONTROL=true,TRAIT=neoplasia,MODE=strict   adhoc_scripts/fcs_calibrated_ranking_control.sbatch
```

Repetir por neoplasia/maligno, completo/strict y rankings global/top/bottom, cada trabajo con salida propia. No usar un array del modo antiguo para WP9_CONTROL. Conservar defaults num_g=10, max_g=0, todas las familias y variantes. EXPECTED_CYCLES evita etiquetar como A2 un almacén de 1000; el control rechaza un roster/matriz incompleto y guarda los identificadores, no fabrica ciclos. Guardar Slurm ID, salida de time -v, metadata/sessionInfo y resultados por base/nivel. Revisar entonces la elección definitiva, separada de la aceptación integrada y de la calibración biológica del modelo nulo.
