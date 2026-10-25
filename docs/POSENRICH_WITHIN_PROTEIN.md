# POSENRICH: pregunta principal dentro de proteína

Decisión del usuario, 25/10/2026: comparar posiciones **dentro de cada proteína**. Una proteína puede acumular señal CAAS por razones generales; queremos saber si sus sitios/regiones anotados concentran señal adicional respecto a otras posiciones de esa misma proteína.

La asociación global entre todas las proteínas es complementaria. No sustituye esta pregunta. PP2A es un piloto del método general.

## Contraste que debe preservarse

Para cada proteína, comparar posiciones pertenecientes al término con posiciones elegibles de la misma proteína que no pertenecen a él. Una proteína sin ambos grupos informativos no aporta ese contraste; informar su exclusión, sin imputar efecto ni p. Resumir después los efectos por proteína, declarando sus pesos y mostrando si el resultado depende de pocos genes.

La AUC dentro de proteína ya implementada es descriptiva. Su ranking CAAS actual se limita a detecciones con puntuación canónica positiva: compara intensidad entre esas detecciones, y no reconstruye el universo completo de posiciones probadas. No presentar esta pregunta condicionada como una prueba de frecuencia de detección en todas las posiciones.

Para contrastar frecuencia, recuperar y verificar la elegibilidad desde el productor CAAS. La ausencia de puntuación no se convierte en score=0 ni en negativo sin demostrar que esa posición fue efectivamente evaluada. Para PP2A se debe respetar además la cobertura experimental: el fondo medido disponible de Sapkota no es automáticamente el fondo de Kauko/Stewart ni todas las posiciones S/T/Y.

## Proximidad

Contrastar si las detecciones están más próximas a los sitios anotados de su misma proteína de lo esperado bajo un control que respete proteína, elegibilidad, cobertura y densidad de sitios. Mantener radios predefinidos y una sensibilidad excluyendo las coincidencias exactas. No mezclar distancia proteica en residuos con distancia genómica ni presentarlo como cercanía 3D.

## Orden de trabajo

Actualización operativa, 25/10/2026: implementación descriptiva operativa. Se
desplegó el contrato Pfam/UCR corregido en los checkouts remotos del productor y
consumidor, con verificación de hashes y respaldo de las fuentes. Las tablas
históricas no se sobrescribieron: para reutilizarlas se requiere reexportación
nativa. La primera optimización de caracterización está implementada, probada
y desplegada: acumulación de términos en una pasada, carga de columnas usadas y
liberación progresiva de anotaciones. La validación inferencial siguiente es una
línea científica separada; no bloquea producir las salidas descriptivas actuales.
La proximidad genómica actual se restringe al mismo gen y cromosoma.

1. Contrato productor→consumidor para UCR/Pfam y pruebas de correspondencia con el MSA/MAP; PP2A ya comprobado hasta péptido y residuo.
2. Definir y verificar elegibilidad/fondo por colección; conservar por separado frecuencia e intensidad.
3. Validar controles y multiplicidad para el contraste dentro de proteína y la proximidad. p/q siguen retenidos mientras tanto.
4. Optimizar el módulo de caracterización al final, como pidió el usuario.
5. Regenerar registro, módulos e informes con las coordenadas y contratos validados, conservando el cálculo CAAS/nulos reutilizables.
