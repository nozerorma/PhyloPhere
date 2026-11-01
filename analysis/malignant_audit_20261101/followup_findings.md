# ASR: efecto sobre scores y causa de la discrepancia b0/master

Comprobación del 2 de noviembre de 2026. No se ha ejecutado PAML ni modificado código de producción o archivos remotos. Se han leído reconstrucciones ya existentes y recalculado scores exclusivamente en local.

## Conclusiones

1. En el primer lote, 121 genes y 1.095 registros gen/posición/lado, la mayoría de los scores observados son estables al cambiar el ASR antiguo por el nuevo. Hay 39 cambios mayores de 1e-9, distribuidos entre 13 genes; hay efectos importantes en posiciones concretas y una reclasificación.
2. La discrepancia de A3GALT2 entre el replay b0 y el master se explica exactamente por dos reconstrucciones PAML diferentes publicadas bajo la misma clave de caché. No es una discrepancia entre las fórmulas de las dos rutas.
3. No se ha estimado la variación de p-valores/FDR, ni la estabilidad de todos los genes o de todos los ciclos nulos. Los scores y rankings descritos son del master sin filtrado posterior.

## Comparación descriptiva de 121 genes

Se compara el master antiguo publicado en malignant_def2 con los 121 shards del primer lote b0 nuevo. Los IDs se resuelven contra los genes del master antiguo: se reconocen sufijos de especie productora y, cuando corresponde, versiones de IDs como AC005551.1. Sin esa resolución, 25 registros parecían nuevos por una diferencia de nombre. Tras resolverla, hay exactamente 1.095 registros compartidos, sin pérdidas ni incorporaciones de posiciones/lados.

- 39/1.095 scores cambian más de 1e-9 (3,56%).
- 13/121 genes tienen algún cambio de score; 7/121 cambian su máximo más de 1e-9.
- Ningún registro cambia de cero a positivo ni de positivo a cero.
- Una clasificación cambia: A1BG, posición cero 229, top, convergent → divergent, score 0,0511930 → 0,0356203.
- El máximo por gen tiene correlación de rangos 0,99994743.
- Los top 20 y top 50 contienen los mismos genes; el top 10 conserva 9: sale A1BG y entra ACSS1.
- A1BG pasa de rango 8 a 11 y su máximo baja de 0,6607798 a 0,6137109.
- La mediana de |delta| es 0; el percentil 99 es 0,00159205.

| Gen / posición cero / lado | Score antiguo | Score nuevo | Delta | Efecto sobre máximo por gen |
|---|---:|---:|---:|---|
| ACAD10 / 495 / top | 0,335996 | 0,191946 | -0,144050 | No cambia el máximo |
| AASDH / 169 / top | 0,169990 | 0,065353 | -0,104637 | No cambia el máximo |
| A1BG / 186 / top | 0,660780 | 0,613711 | -0,047069 | Baja el máximo y sale del top 10 |
| ABCA6 / 292 / top | 0,107272 | 0,088738 | -0,018534 | El máximo cambia mucho menos: -0,000374 |

El conjunto procede de un lote ordenado por nombre (genes A...), no de un muestreo aleatorio de todo el genoma. Estos porcentajes no se deben extrapolar como tasas poblacionales. Los rankings son internos a esos 121 genes.

## Prueba causal: 24 genes y 278 registros

Se incluyeron todos los 13 genes con diferencias del conjunto de 121, A3GALT2 y controles repartidos por el lote. Se dispuso de dos cachés para 24 genes. La selección está enriquecida deliberadamente en discrepancias; su porcentaje de cambios no estima la frecuencia global.

Para cada gen:

- Código actual + ASR antiguo, con discovery, diseño, pesos y umbral fijos.
- Código actual + ASR nuevo, manteniendo las demás entradas.
- Código b112c1fb + exactamente el mismo contexto ASR nuevo serializado.

Resultado: los 39 deltas de score y el cambio de clasificación se reproducen sustituyendo solo el ASR. Los masters antiguos y nuevos publicados se reproducen dentro de tolerancia 1e-9 en todos los genes. Al cambiar solo la implementación hay cero diferencias de score en los 278 registros y cero cambios de clasificación. No hay cambios de claves de registro.

En 23 de los 24 pares, el ASR nuevo incorpora taxid 1268362 (Trachypithecus melamera). La prueba no separa los efectos de incorporación de especie, normalización de missingness, serialización del árbol y ajuste numérico de PAML; confirma el efecto del contexto ASR completo.

Resultados completos: controlled_subset24.json; tablas de 121 genes: master_comparison_121.tsv y gene_comparison_121.tsv; métricas: master_summary_121.json.

## A3GALT2: discrepancia resuelta

Clave compartida de caché:

`validated_v1/19a32fb54b7260c18e34553304be0089ee27bb9abc2c6dab622d5e522b23f420/`

Se conservan dos entradas:

- `.invalid-5fb933a837c34daeae21b61100fa1417/`: rst de las 13:12:09 del 1 de noviembre; usado en el replay b0 de la tarea.
- `asr_A3GALT2.Homo_sapiens/`: rst de las 13:12:13; usado por el master observado posterior.

Los provenance.json declaran entradas idénticas. alignment_paml.phy, tree_paml.nwk, codeml.ctl y lg.dat coinciden byte a byte. Los rst/rst1 difieren. Los dos ajustes publicados muestran:

| Reconstrucción | alpha | lnL |
|---|---:|---:|
| Primera, sustituida | 0,13502 | -4940,736444 |
| Segunda, conservada | 0,06740 | -4882,728025 |

Esto demuestra diferente resultado del ajuste de PAML pese a entradas idénticas. No se ha aislado por qué el optimizador llegó a dos resultados (no se ha relanzado PAML); no se atribuye a un cambio de fórmula, de especies ni de parámetros declarados entre estas dos ejecuciones.

Recalculando con el primer rst se reproducen exactamente los cinco scores positivos guardados del replay b0 de A3GALT2, incluidos:

| Posición cero / lado | Primer ASR = replay b0 | Segundo ASR = master |
|---|---:|---:|
| 222 / bottom | 0,18481481481340437 | 0,1851851851837719 |
| 294 / top | 0,2351730098184759 | 0,23900865523855294 |

Con un mismo contexto ASR, la ruta full y la ruta axes-only dan igualdad exacta, tanto usando el discovery fijo como repitiendo el kernel CAAS de b0. Los labelings y pesos b0 reconstruidos para el experimento coinciden con los originales de la tarea.

El chequeo completo del tile encuentra cuatro discrepancias gen/posición/lado entre su replay b0 y master: las dos de A3GALT2 y una en A4GNT y ABRA respectivamente. Solo las dos de A3GALT2 se han atribuido causalmente a una entrada de caché sustituida en este análisis; no se presupone el mismo origen para las otras dos. Valores en tile_internal_mismatches.json.

## Implicación para recuperación

Los scores son mayoritariamente estables, pero algunos deltas locales y el cambio de clasificación son reales. No es defendible declarar que ambos conjuntos son idénticos; tampoco hay evidencia de un cambio generalizado de fórmula del score observado.

Antes de recuperar el conjunto nuevo, fijar un resultado ASR por clave y evitar reemplazar una entrada válida cuando termina otro cálculo concurrente. La corrección de recuperación tras FileExistsError no elimina por sí sola la rama de publicación que mueve un target existente a .invalid y publica otro. La publicación debería volver a validar y conservar al ganador válido, con exclusión mutua/atomicidad, y rechazar mezclas en los tiles mediante identidad de los artefactos ASR utilizados. Esto es una recomendación: no se ha modificado aquí la implementación.

Recalcular los tiles que usaron otra reconstrucción es suficiente en principio para esta causa, manteniendo la elegida y sin volver a ejecutar PAML. No se han inventariado todos esos tiles ni verificado la identidad rst por gen en toda la ejecución: la procedencia de tiles actual no registra directamente el checksum del rst usado, por lo que detectar todo el alcance requiere logs/cachés históricas o un replay con la caché fijada.
