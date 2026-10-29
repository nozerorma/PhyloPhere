# WP9: evaluación dirigida del estadístico de gen — 23/10/2026

La evidencia disponible favorece **conservar provisionalmente `size_adj_max = F(max)^n`**. Reduce la dependencia del score con el número de posiciones sin mostrar una diferencia apreciable del FCS frente a F(max) en los controles examinados. No hay evidencia aquí para sustituirlo por máximo crudo ni para elegir un β intermedio. **No se ha cambiado el default ni cerrado la decisión definitiva: A2 a 10.000 ciclos sigue pendiente.**

## Qué se comparó

Se añadieron máximo, F(max), β=0,25/0,5/0,75 y F(max)^n a un control común de validación, sin modificar el cálculo productivo. El ranking sigue siendo −log10(p del gen); gene_p.R conserva empates inclusivos y tolerancia 1e-12. Rutas usan AUC−0,5 unilateral, comparación exacta (sin esa tolerancia) y pseudo-cuenta. BH por base/ranking y Redlich permanecen separados.

- Marino completo y strict: bloques canónicos, 120 genes del universo FCS, 100 ciclos y tres rankings. Se cargan los 17 GMT; resultan 133 rutas por ranking (399 por variante y modo). ADAM3A queda fuera del universo depurado, igual que en la transposición productiva.
- Neoplasia estricto histórico: 12.473 genes del universo depurado, 1.000 ciclos, ranking global y dos bases fijadas antes de mirar resultados: Hallmark y KEGG (50+314 rutas). No representa todas las bases, direcciones ni las cuatro corridas.
- Las piezas históricas se reconstruyen sin redondeo adicional; F(max)^n reproduce el RDS guardado con diferencia máxima 2,22e-16. Los pools se forman antes de seleccionar el universo del FCS, preservando la definición del score. El roster completo procede del RDS, incluidos ciclos sin detección.

Se calcula el diagnóstico leave-one-out y un holdout de dos mitades disjuntas, intercambiadas, semilla 1998. Cada gen de evaluación se calibra contra la otra mitad; el ranking de los ciclos de referencia usa leave-one-out interno. Las rutas de evaluación se comparan solo con AUC de referencia. La resolución es 1/51 en marino y 1/501 en neoplasia, no la de una referencia de 10.000.

La tasa marginal leave-one-out próxima al 5 % está casi impuesta por los rangos; no certifica calibración del modelo. El holdout elimina la inclusión del ciclo evaluado en su referencia, pero comparte ciclos entre folds y conjuntos entre rutas: **no son réplicas independientes ni una estimación de FDR**.

## Neoplasia estricto, global, Hallmark+KEGG

| Variante | KS | ρ con n | p≤0,05 en holdout (todos) | p≤0,05 en ciclos de nivel observado |
|---|---:|---:|---:|---:|
| Máximo crudo | — | 0.225 | 4.965 % | 8.248 % |
| F(max) | 0.064 | 0.225 | 4.961 % | 8.250 % |
| β=0,25 | 0.053 | 0.169 | 4.960 % | 8.234 % |
| β=0,5 | 0.043 | 0.108 | 4.966 % | 8.256 % |
| β=0,75 | 0.032 | 0.045 | 4.965 % | 8.258 % |
| F(max)^n | 0.023 | -0.020 | 4.965 % | 8.258 % |

Hay 364.000 evaluaciones válidas por variante; 49.868 corresponden a los 137 ciclos cuyo número de genes puntuados está entre 0,8 y 1,5 veces el observado. Esa condición se aplica al resumen, nunca modifica el nulo usado para calcular el p. Hallmark da aproximadamente 10,3 % y KEGG 7,9 % en ese subconjunto. La asociación con el nivel persiste con todas las variantes; no la corrige cambiar el exponente y no demuestra por sí sola un fallo de control de FDR.

F(max) frente a F(max)^n: Spearman de −log10(p observado)=0,9711; ninguno de los 2.480 genes observados cambia más de una unidad. Solo cuatro tienen n>8. Los recuentos de genes con p≤0,05 son 922 y 894; BH por gen sobre la familia completa da cero a q<0,15 en todas las variantes. Los recuentos de rutas significativas se guardan como diagnóstico y no se usan para elegir el estadístico.

Los resultados no sustituyen mecánicamente los de A1: aquel análisis medía los genes detectados de todo el almacén; aquí se usa el universo fijo del FCS. Por eso cambian KS/correlaciones y el grupo emparejado (208→137 ciclos). Los recuentos observados a p≤0,05 y el efecto pequeño del exponente coinciden con A1.

## Marino

En global, el holdout da p≤0,05 en 3,47/3,46/3,35 % para máximo/F(max)/F(max)^n en completo, y 3,44/3,34/3,62 % en strict. La resolución de 50 referencias y los empates limitan estas tasas. En completo, ρ(score,n) pasa de 0,627 a 0,035; en strict, de 0,620 a 0,011, al pasar de F(max) a F(max)^n. No se exige que un score dentro de ciclo sea uniforme como criterio suficiente de calibración del p por gen o ruta.

El modo F(max)^n reproduce los 399 resultados completos de WP7: claves, stat, pval, p.adj, p.perm, q.perm y sig_perm a 1e-12.

## Validación, recursos y pendientes

Oráculos versionados: fórmulas calculadas elemento a elemento, AUC nativo por ciclo, conteos exhaustivos de p/BH, holdout con calibración manual y AUC nativo, universo fijo, ciclos vacíos, metadatos/partes y rechazo de almacenes incompletos. Launcher existente ampliado, comprobado con repositorio Git mínimo y ejecutor simulado; bash -n y parse R aprobados. Se conserva el resultado completo de pruebas en PhyloPhere_validation/.tmp/wp9/.

Controles en memoria completa con holdout, medidos con /usr/bin/time -v: marino completo 43,22 s / 616.564 KiB RSS; strict 43,38 s / 605.176 KiB; neoplasia global sobre dos bases 93,47 s / 3.492.284 KiB (3,33 GiB). Son controles de validación, no el camino por bloques productivo ni medidas a 10.000. No se ha ejecutado una corrida completa ni trabajo de cluster.

A2 requiere nulos reales de 10.000, las cuatro combinaciones neoplasia/maligno×completo/strict, direcciones y familias completas, y la misma definición de universo. El launcher `fcs_calibrated_ranking_control.sbatch` admite `WP9_CONTROL=true`, variantes, referencias canónicas, snapshot de HEAD y `EXPECTED_CYCLES=10000`. Rechaza reutilizar directorios y archivos incompletos; no sustituye el nulo por columnas recicladas. La memoria de estos controles es global y se presupuesta aparte del pipeline.

La elección final de gene_stat permanece abierta hasta revisar A2. La recomendación provisional conserva el valor vigente; no se basa en maximizar significativos ni afirma que ^n resuelva el efecto de nivel.
