# Sustitución completa de percentiles por evidencia empírica

Plan solicitado el 23 de octubre de 2026. Inventario sobre el checkout actual; este documento no implementa la migración ni autoriza lanzar trabajos. El árbol contiene cambios concurrentes de anotaciones de posiciones: coordinar el contrato de entrada sin sobrescribirlos.

## Contrato de resultados

- Tres familias de AMI: CAAS, RER y FADE. Conservar las direcciones pertinentes; «una familia» no significa mezclar aceleración y desaceleración ni eliminar las direcciones CAAS.
- CAAS: genes activos con p.gen < 0,05; resaltar p.gen.adj < 0,05. RER: activos con p.perm < 0,05; resaltar p.perm.adj < 0,05. FADE: BF >= 100, conservando el umbral configurado.
- Una ejecución por método/dirección definida. No generar otra ejecución solo por el ajuste BH. Las tablas deben permitir identificar/filtrar ambos subconjuntos. Filtrar una tabla por BH no equivale a recalcular módulos con semillas BH: la segunda ejecución queda fuera de este contrato.
- Cada nodo distingue: semilla nominal; semilla que también cumple BH; nodo incorporado por la red. Un nodo incorporado no adquiere significación individual por pertenecer al módulo. Conservar la significación del módulo como resultado distinto.
- Posiciones: sustituir TODOS los cortes por p.emp < 0,05 y mostrar p.adj_bh < 0,05 como evidencia ajustada. Por decisión del usuario (25/10/2026), se retiran p.emp_fact y p.adj_bh_fact de productores, tablas e informes activos. Las auditorías históricas se conservan.
- FCS CAAS sigue usando el ranking continuo -log10(p.gen), con todos los genes del universo y sus nulos calibrados. No filtrar el ranking ni sustituirlo por BH. Mantener AUC, reglas de colas, empates, pseudocuentas, LOO y significación de pathways.
- Los valores continuos CAAS_score y gene_stat siguen disponibles como estadísticos de origen. Se eliminan las selecciones/anotaciones por percentiles, no el cálculo del estadístico.

## 1. Calibración de genes como productor compartido

Archivos actuales: `subworkflows/ENRICHMENT/fcs_blocks.nf`, `fcs_inputs.nf`, `local/src/fcs_blocks.R`, `caas_blocks.R`, `gene_p.R`; integración en `workflows/scoring.nf`, `workflows/enrichment.nf` y `main.nf`.

Hoy FCS_GENE_P calcula los p-valores de genes dentro de FCS y exporta internamente -log10(p) por batch. Extraer esa etapa a un productor compartido independiente de GMT, FCS y AMI. Debe estar disponible con FCS apagado y funcionar con nulos canónicos vivos o reutilizados, complete y strict. Reutilizar el cálculo existente y sus bloques, evitando repetir transposición y calibración para cada consumidor.

Publicar una tabla canónica por gen/ranking con IDs, estadístico observado, gene_stat, p.gen, p.gen.adj, flags nominal/ajustado, N y familia de ajuste. Ajustar BH después de reunir el universo completo del ranking, nunca por batch, sobre hits solamente o tras el mapeo a STRING. Documentar global/top/bottom como familias explícitas: una llamada dentro de una familia no implica control conjunto sobre las tres. Ausencias sin señal del universo conservan su tratamiento canónico; no eliminar genes para favorecer BH.

Conservar los bloques -log10(p) del nulo como entrada de FCS. Resolver la dependencia sin ciclos: estadísticos observados -> calibración -> tabla de evidencia -> informes/AMI/FCS. SCORING_COMPUTE no debe esperar a un informe que a su vez requiera la calibración. Las rutas históricas que no permiten reconstruir p.gen deben fallar de forma informativa o indicar evidencia no disponible; no usar percentiles como fallback.

## 2. Productor de posiciones y familia BH

Archivos: `subworkflows/SCORING/local/src/scoring_compute.R`, `null_position_p.R`, `null_position_stats.py`, `subworkflows/SCORING/scoring_null_positions.nf`, `scoring_compute.nf`, `workflows/scoring.nf`.

Reutilizar p.emp y p.adj_bh existentes desde contadores, roster y universo. No hacen falta tablas completas posición x ciclo para seleccionar posiciones por estos valores. Publicar flags nominal/ajustado con los mismos criterios usados por los informes.

El p.emp actual usa el máximo entre lados por (Gene, Position): sus filas top/bottom comparten el mismo p. Deduplicar la posición para BH, conteos globales y cruces; separar por lado solo como anotación/selección de las filas observadas. No presentarlo como dos pruebas direccionales independientes. Mantener la familia posicional existente, incluidos los miembros requeridos por su definición aunque no sean hits observados.

Si falta el nulo, los p quedan no disponibles. Producir listas vacías/no disponibles con un estado explícito, sin seleccionar por score ni recuperar percentiles. Mantener nd, null_unseen, N y diagnósticos del nulo.

## 3. Retirada de listas, flags y análisis por cortes

`subworkflows/SCORING/local/src/scoring_compute.R`:

- Retirar las 12 listas de genes y 12 de posiciones de 25/10/5/1%.
- Sustituirlas por listas nominales por dirección, con columnas p y BH y flags ajustados; las listas son proyecciones de las tablas canónicas.
- Retirar la tabla común 100/50/25/10/5/1 y sus barridos de enriquecimiento por genes y posiciones. Los análisis de solapamiento conservados usan conjuntos nominales/ajustados identificados y sus universos de comparación explícitos.
- Eliminar los campos pct_* y las extracciones por tamaño fijo de candidatos en toda salida de genes/posiciones. No mantener rutas que reconstruyan percentiles desde scores.

`subworkflows/ENRICHMENT/local/src/percentile_flags.R`: retirar los flags y sus fallbacks. Actualizar importaciones, staging y hashes en `fcs_enrich.R`, `fcs_compute.R`, `fcs_blocks.R`, `fcs_blocks.nf` y demás consumidores encontrados por búsqueda completa. Sustituir por lectura de evidencia canónica.

## 4. AMI: entradas, orquestación y visualización

Archivos: `subworkflows/RERCONVERGE/rer_gene_lists.nf`, `subworkflows/FADE/fade_gene_lists.nf`, `subworkflows/ENRICHMENT/domino.nf`, `scoring_enrichment.nf`, `local/src/run_domino_modules.py`, `workflows/enrichment.nf`.

CAAS deja de consumir slice_*.tsv y recibe listas de actividad nominal. RER mantiene una familia nominal con flags BH. FADE utiliza su BF configurado y una definición coherente con SCORING. Revisar el staging de listas .txt para que no incorpore por accidente ficheros de background u otras listas auxiliares.

La red de cada método se construye sobre su universo evaluado, restringido por el mapeo/conexiones de STRING; el ajuste BH se hace antes de esa restricción. Las listas nominal y ajustada se calculan sobre la misma familia; no contar RER nominal y BH como evidencias independientes. Respetar ramas sin hits, sin aristas o sin módulos.

## 5. Informes afectados

| Informe | Sustitución prevista |
| --- | --- |
| 5, RER | Mostrar p.perm y p.perm.adj, llamada nominal y ajustada; sincronizar definiciones con SCORING y AMI. Corregir el alcance/nombre del selector para no sugerir que una lista ajustada inicia otro AMI. |
| 11, SCORING | Retirar umbrales/colores/tablas por porcentajes en genes y posiciones, incluidos los seis cortes de composición de convergencia. Mostrar p/BH y conteos de conjuntos nominales/ajustados. |
| 12, FCS general | Mantener ranking y estadística FCS. Cambiar anotaciones pct_* de leading edges por p.gen/BH y flags; sincronizar composición/exportaciones. |
| 13, AMI | Tres familias de método, direcciones existentes; sustituir listas por cortes, tablas top10/top5/top1 y descripciones de módulos por evidencia nominal/BH. Mostrar ambos estados en cada red, tabla, leyenda y exportación. |
| 14, posiciones | Sustituir listas y anotaciones por p.emp/BH, incluidas tablas de solapamiento con FADE/VEP y resúmenes por gen. Coordinar con el porte de POSENRICH; este cambio no activa por sí mismo el enriquecimiento por permulaciones de conjuntos de posiciones. |
| 15, comparación | Retirar columnas de percentiles, .pct_tier_score y sus contribuciones a la puntuación compuesta de genes y posiciones. Mostrar p nominal y ajustado por método, dirección y llamadas. No inventar otra ponderación numérica ni sumar nominal y BH como pruebas independientes. Mantener las anotaciones de FCS/AMI separadas de las llamadas individuales. |
| 16, significación CAAS | Reutilizar p.emp/BH de posición como ya hace; verificar consistencia de flags, familia, joins y conteos con 11/14/15. |

Fuentes: `subworkflows/RERCONVERGE/local/5.RERconverge_report.Rmd`, `subworkflows/SCORING/local/11.Scoring_report.Rmd`, `subworkflows/ENRICHMENT/local/{12.FCS_general_report,13.AMI_analysis,14.Position_enrichment_report,15.Comparison_report}.Rmd`, `subworkflows/CT_META_CAAS/local/16.CAAS_significance_report.Rmd` y sus wrappers Nextflow, incluido `fcs.nf`, `scoring_enrichment.nf`, `posenrich.nf`, `scoring_report.nf` y `ctpp_meta_caas.nf`.

## 6. Parámetros, GUI, migración y generación

Retirar scoring_gene_top_pct y scoring_position_top_pct y los argumentos derivados top25/top5/top1. Revisar `conf/scoring.config`, `gui/models/modules.py`, `gui/models/serialization.py`, `gui/widgets/tabs/scoring_tab.py`, `gui/generation/context.py`, ambos templates shell y `gui/generation/report_registry.py`, traducciones, ayuda y todos los proyectos/fixtures versionados.

Proponer dos umbrales compartidos de actividad nominal (genes/posiciones, default 0,05) y umbral BH de resaltado (default 0,05), con nombres y alcance inequívocos. Evitar selectores que permitan que AMI, listas y reportes diverjan. Cargar proyectos antiguos eliminando parámetros retirados con aviso; nunca convertir un percentil en un p-valor numéricamente igual. Los parámetros de RER se migran al contrato común conservando preferencias explícitas compatibles.

Conservar PSS_TOP_PCT/pss_top_pct del diseño de permulaciones y controles de contrastes, y DOMINO slice_threshold del algoritmo: son parámetros distintos, no percentiles de selección del score observado. Los límites de filas/etiquetas de visualización tampoco son pruebas de significación, pero deben quedar identificados como presentación.

## 7. Orden de implementación y aceptación

1. Congelar el contrato de tablas/listas y casos de referencia completos/strict; coordinar entradas de posiciones con el trabajo concurrente.
2. Extraer calibración compartida de genes y publicar tabla canónica. Comprobar identidad numérica del ranking/FCS anterior para entradas idénticas.
3. Construir listas nominales y flags ajustados de genes y posiciones desde las tablas canónicas.
4. Actualizar AMI y staging para tres familias; luego informes 5/11/12/13/14/15/16 y exportaciones.
5. Retirar configuración/GUI y migrar fixtures. Actualizar regeneración independiente de informes.
6. Retirar funciones/importaciones antiguas; búsqueda final de consumidores activos de pct_*, slice_<dir><pct>, top100/top50/top25/top10/top5/top1 y parámetros scoring_*_top_pct. Las coincidencias históricas deben distinguirse de código ejecutable; no borrar documentación científica de controles pasados.
7. Aceptación pequeña integrada con datos reales reutilizados, antes del run 10k. Verificar FCS apagado + AMI encendido, solo informes, complete/strict, nulo ausente, cero hits y red vacía. No lanzar trabajos grandes como parte de esta migración sin la decisión de ejecución correspondiente.

## 8. Verificación necesaria

- Oráculo de p.gen observado y LOO conservado; BH idéntico frente a un cálculo sobre el universo completo, independiente del tamaño de batch.
- Identidad de p.emp/BH posicional frente al cálculo actual; posiciones duplicadas por lado no duplican tests ni conteos.
- Identidad de estadísticas/p-valores FCS para las mismas entradas; las selecciones AMI no filtran FCS.
- Listas activas iguales a los IDs nominales de su tabla; flags BH como subconjunto; anotaciones, reportes y exportaciones consistentes.
- Tres familias de AMI, sin ejecuciones extras por BH; las direcciones y sus nombres son coherentes. Los nodos incorporados por la red no se etiquetan como hits.
- Casos con ningún ajustado pero sí nominales muestran ambos recuentos y leyenda correcta. No variar umbrales automáticamente para obtener hits.
- Ausencia de nulo no produce actividad inventada. Reutilización correcta después de borrar work, con artefactos públicos suficientes.
- Regenerar scripts/proyectos e informes representativos; comprobar sintaxis, joins, vacíos y que no queden columnas o narrativas de percentiles.

El cambio requiere recomputar calibraciones/listas/reportes/AMI, no volver a generar los alineamientos ni los ciclos canónicos existentes si sus partes, roster, contadores y metadatos son compatibles. La decisión científica de POSENRICH sobre conjuntos de posiciones sigue siendo un trabajo separado.
