# WP10 propuesto: POSENRICH por posición

Estado: PREPARADO, sin implementación ni aceptación. El usuario eligió preparar la opción A el 2026-10-23 en WP8. Posteriormente pidió investigar en un fork si un contraste de recuentos como acumulación sería preferible. Esa comparación sigue abierta; preparar A no autoriza activarla ni resolver por cuenta propia las elecciones científicas restantes.

## Entrada y pregunta científica

Leer `POSENRICH_PORTE_EVALUACION.md`, los contratos vigentes de WP5/WP6/WP7 y `CT_ACCUMULATION/local/accumulation_test.py`. Distinguir enriquecimiento de detecciones (recuentos) y enriquecimiento del ranking (AUC). Los contadores marginales `nd` permiten reconstruir el universo de posiciones positivas del nulo, pero no el score ni la coocurrencia de posiciones de un mismo ciclo. Un test hipergeométrico analítico necesita tamaños y solapamiento observado una vez definido el universo; un p de permulación requiere evaluar el estadístico de cada ciclo o conservar sus recuentos suficientes por conjunto. No son intercambiables sin cambiar la hipótesis nula.

## Propuesta concreta para la opción A

1. Fijar antes de implementar: universo probado por dirección, unidad `gen:posición`, duplicados, posiciones no detectadas (cero) frente a no evaluables, ranking exacto por posición y lateralidad. Propuesta a contrastar: AUC centrado sobre score de posición en cada dirección, usando el mismo tratamiento de empates entre observado y nulo; no adoptar automáticamente el ranking calibrado por gen para posiciones.
2. Corregir la difusión de `charset` a todos los lotes (`.collect().map { it[0] }`), y convertir el inventario de fuentes/GMT en un contrato. Ninguna fuente desaparece porque un canal se agote; cada ausencia exige motivo (sin posiciones, sin anotación, fuera de tamaño).
3. Prototipo: consumir las teselas/bloques de WP5 para calcular el AUC de cada conjunto por ciclo y reducir contadores aditivos. Evaluar si debe almacenarse un formato disperso de posiciones positivas por bloque o calcular los estadísticos durante el merge y descartarlo; no reintroducir por defecto la tabla global de posiciones × N. Medir por separado costes del universo, del bloque y de los conjuntos.
4. Mantener observado interno RDS, identificadores estables y contadores exactos de ciclos válidos. BH por familia explícita y p inclusivo con pseudo-cuenta si se aprueba el mismo contrato de WP6. No usar los p históricos Lachenbruch como oráculo de equivalencia numérica de AUC.
5. Conservar POSENRICH omitido hasta que el nuevo módulo y sus consumidores estén probados. No interpretar una salida omitida como ausencia de enriquecimiento.

## Comparación con recuentos antes de activar

Comparar AUC con hipergeométrico sobre exactamente el mismo universo, anotaciones y datos: qué efecto detecta cada uno, dependencia de detectabilidad por gen, sensibilidad al nivel observado/nulo, redundancia con acumulación, coste y resolución. Separar p analítico de p calibrado por ciclos. No elegir por cantidad de significativos. El fork «Acumulación y alternativas para POSENRICH» investiga esta puerta; incorporar su conclusión tras revisión del usuario.

## Validación y criterio de entrega

- Oráculos independientes pequeños: conjuntos de posiciones con ceros, empates, NA, vacíos, conjuntos de todo el fondo y posiciones repetidas; AUC manual y test nativo concordantes.
- Particiones {1,7,N}: mismo observado y recuentos enteros de excedencias/ciclos válidos, p/q a 1e-12; cobertura de todas las fuentes y lotes (regresión del fallo 5/12).
- Nulo como observado con referencia separada; tasas descriptivas por fuente y nivel, sin interpretar réplicas dependientes como independientes. No certificar FDR solo por recuentos de BH.
- Toy marino/cáncer con posiciones anotadas, completo/strict; no usar PEPC de un solo gen como equivalencia positiva de rutas.
- Medir RSS, tiempo y E/S en pequeño; las corridas 1000/4000/10000 se preparan para el cluster. Comparar métodos científicos antes de habilitar un default.
- Decisión final del usuario sobre método, ranking y lateralidad; implementación, GUI, ayuda, informe 14 y aceptación se registran por separado.
