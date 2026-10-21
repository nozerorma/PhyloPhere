# WP8. Cierre: limpieza, documentación y sincronía de GUI y plantillas

Estado: pendiente. Depende de: WP5 (y de lo que hayan cerrado WP6 y WP7 si ya están).
Conversación independiente. Antes de empezar, leer `docs/fases/README.md`.

## Objetivo

Cerrar las decisiones que quedaron condicionadas a los cambios del nulo y dejar el repositorio coherente.

## Decisiones y puertas

- **Destino de `remove_intermediates`:** una vez el nulo no tiene exports ni detalle de texto que podar, retirar la opción, o redefinirla como "no publicar el almacén masivo por posición" si aún se publica. Los `rm` residuales de las plantillas (higiene de FADE, selección, `__pycache__`, enlaces) son incondicionales y no los gobierna la opción.
- **Validación pendiente de `remove_work`** en una corrida real (comparación de huellas y ahorro de disco), si no se hizo en WP2.
- **Reducción de redundancia:** calibración del umbral y la similitud por defecto con las tablas existentes (recomendaciones de EnrichmentMap: solapamiento 0,5 o Jaccard 0,25, de su documentación y no de un estudio independiente).
- **Porte de POSENRICH:** decidir con `docs/POSENRICH_PORTE_EVALUACION.md` (opciones A, B, C) y cerrar antes la anomalía de 5 de 12 fuentes de posiciones probadas en la corrida del 19/10.

## Alcance

README, `workflows/help.nf`, `docs/FCS_REVISION_RERCONVERGE.md` (añadir los controles del ranking calibrado, el máximo crudo, las tres pruebas de posición y gen y el BH por gen: hoy solo están en `docs/PLAN_NULO_10K_PERMULACIONES.md`, secciones 2 y 5), `docs/PLAN_NULO_10K_PERMULACIONES.md` (mediciones, decisiones y fases), `docs/fases/README.md`, sincronía de `gui/templates/*.json` con `validation/tier1/templates` y `validation/tier2/templates` (la plantilla tier 1 y `gui/templates/tier1_pepc_c4.json` deben quedar idénticas), campos de la GUI de todos los parámetros nuevos o retirados, `gui/i18n.py`.

## Validación (criterios de paso)

- Plantillas y GUI coherentes: abrir y guardar cada plantilla; scripts generados sin variables retiradas.
- Documentos sin referencias a parámetros o ficheros retirados.
- Una corrida de extremo a extremo a 10.000 ciclos de neoplasia completo con memoria por tarea de pocos GB y los mismos resultados que la versión de referencia donde ambas son comparables.

## Registro de la fase

```
- Conversación y fechas:
- Commit de partida / commit de cierre:
- Qué se hizo (ficheros y cambios):
- Comprobaciones ejecutadas (comando, conjunto de datos, tolerancia, resultado):
- Mediciones (recurso antes y después):
- Desviaciones del plan y decisiones tomadas (con su justificación):
- Pendiente que pasa a otra fase:
```
