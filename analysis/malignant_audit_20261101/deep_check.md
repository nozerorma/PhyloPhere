# Comprobación controlada del score y de las entradas ASR

Fecha: 2 de noviembre de 2026. Cómputo exclusivamente local sobre copias de archivos; no se ha ejecutado PAML, modificado cachés remotas ni lanzado jobs.

## Diseño

Se han probado A1BG, A1CF, A3GALT2 y A4GALT. Para cada gen se cruzaron dos implementaciones de desambiguación y dos contextos ASR: código del commit b112c1fb (26 octubre) y código local actual (7c10bfda, equivalente científicamente al del run nuevo en estas funciones), con reconstrucción antigua asr_<gen> y reconstrucción nueva validated_v1. Se mantuvieron fijos el discovery nuevo, el diseño de 72 hipótesis, los pesos PSS y el umbral 0,1.

El contexto incluye posteriors y árbol de cada reconstrucción. Se desactivó explícitamente la posibilidad de ejecutar PAML en el arnés. Los contextos se construyen con los lectores actuales y se serializan para compartir exactamente los mismos objetos entre los dos scorers. Por tanto, esta prueba aísla la implementación posterior al contexto ASR; no compara por separado los antiguos/nuevos lectores de posteriors. Esos lectores no presentan cambios en el diff consultado.

El código antiguo se ha obtenido de git, no de un snapshot acreditado del task original. Como comprobación independiente, los resultados antiguos publicados se reproducen en todos los registros de los cuatro genes con diferencias máximas inferiores a 1e-14. Los nuevos publicados se reproducen exactamente.

## Resultados

| Gen | Registros gen/posición/lado | Cambios de score entre implementaciones, ASR fijo | Cambios al sustituir solo ASR |
|---|---:|---:|---:|
| A1BG | 25 | 0 con ambos ASR | 8 |
| A1CF | 2 | 0 con ambos ASR | 0 |
| A3GALT2 | 9 | 0 con ambos ASR | 0 |
| A4GALT | 7 | 0 con ambos ASR | 0 |

Son 43 registros con cada contexto, 86 comparaciones entre implementaciones: igualdad exacta del score, de las claves y del número de registros. No constituye una demostración global para todo el dataset ni para strict, residuos ambiguos en pares focales u otros casos límite.

En A1BG, cambiar solo ASR reproduce exactamente los ocho deltas de los masters publicados, incluido 0,6607797606785337 → 0,6137109154590084 en posición cero 186, lado top. El máximo delta es 0,04706884521953125.

## Qué cambia en los ASR y por qué

| Gen | Especies ASR antiguo | Especies ASR nuevo | Longitud en ambos | Sustituciones en secuencias compartidas |
|---|---:|---:|---:|---|
| A1BG | 225 | 226 | 466 | 30 X → - |
| A1CF | 225 | 226 | 579 | ninguna |
| A3GALT2 | 208 | 209 | 316 | 12 X → - |
| A4GALT | 197 | 198 | 350 | 28 X → - |

En los cuatro casos la especie adicional es taxid 1268362, Trachypithecus melamera. No hay cambios de aminoácidos estándar ni desplazamientos de longitud en las secuencias compartidas de estos cuatro genes.

Se ejecutaron ambos preparadores (load_alignment_and_mappings y load_and_match_tree) con los mismos archivos actuales de A1BG: ambos producen 226 secuencias, 226 tips e incluyen 1268362. Por tanto, la diferencia de especies no está causada por una nueva regla de estos preparadores. La reconstrucción antigua guardada se generó con un conjunto diferente; reutilizarla por nombre de gen conservaba esa diferencia. Su alignment_paml.phy y rst tienen fecha remota de diciembre de 2025; la antigüedad por sí sola no establece incompatibilidad, pero la ausencia comprobada de la especie sí demuestra que las entradas no son idénticas.

El cambio X → - es explícito en el commit 3defc510 (31 octubre, 03:55), asr_single.py: todas las secuencias se pasan a mayúsculas y cualquier símbolo fuera de STANDARD_AAS se sustituye por '-'. Antes no se aplicaba esa normalización en esta entrada ASR. La operación representa residuos no resueltos como missingness; no se ha realizado un experimento PAML que separe el efecto de X → - del de añadir una especie. No se debe atribuir causalidad numérica a cada cambio por separado.

En A1BG, tras podar la especie añadida de la copia del árbol nuevo, la topología de clados coincide con la antigua; las longitudes de ramas no raíz coinciden dentro de 1e-4. Hay también un cambio de precisión de serialización del árbol en el código nuevo. La rama raíz se representa de forma distinta; no se ha cuantificado su efecto en PAML.

## Comprobación adicional: ruta observada frente a replay b0

El master antiguo y su nulo b0 coinciden numéricamente en las 27.080 posiciones compartidas al colapsar ambos lados por máximo (tolerancia 1e-9).

En los cuatro genes del tile nuevo, el master y null_tables/scores.parquet b0 coinciden en A1BG, A1CF y A4GALT. En A3GALT2 hay dos discrepancias en las posiciones compartidas:

- Posición 222: master 0,1851851851837719; replay b0 0,18481481481340437.
- Posición 294: master 0,23900865523855294; replay b0 0,2351730098184759.

Esto no es un delta entre las dos versiones de la fórmula observada: el experimento de masters de A3GALT2 da igualdad exacta para código/ASR cruzados. El origen de esta discrepancia entre rutas no está aislado en esta comprobación. Es una razón adicional para no certificar todavía el conjunto nuevo observado/nulo como intercambiable por simple cambio de metadatos. Los números de comparaciones anteriores basados en Parquet b0 mezclaban esta cuestión con los cambios ASR del master.

## Implicaciones

La evidencia controlada no respalda afirmar que la fórmula central de asr_path_score haya cambiado entre estas dos versiones para los casos probados. Sí respalda que reutilizar ASR anteriores con entradas distintas cambia los scores de A1BG. No es correcto convertir automáticamente esas entradas en compatibles añadiendo provenance.json.

Para decidir qué reutilizar, verificar especie/seq/árbol/modelo por gen, conservar los ASR que coincidan con las entradas elegidas y fijar un conjunto coherente para observado y nulo. La ampliación del universo posicional y los cambios de filtros siguen siendo cuestiones independientes del score observado. La discrepancia adicional en el replay b0 requiere una comprobación específica antes de dar por válida una recuperación global.

Arnés: controlled_score_check.py. Resumen numérico: controlled_summary.json. Los contextos y copias de entrada están en .tmp/asr_deep_audit/ (datos temporales, no forman parte de la implementación).
