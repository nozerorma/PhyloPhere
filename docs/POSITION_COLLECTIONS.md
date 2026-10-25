# Colecciones posicionales y preparación de PP2A

Decisión del usuario, 23/10/2026: conservar Pfam dominios/clanes, UCR core/flancos globales, FUBAR positiva/purificadora globales, FADE top/bottom globales, ventanas genómicas de 1 Mb y conjuntos personalizados. Retirar colecciones basadas en ortogrupos, también sus intersecciones con UCR, FUBAR, COSMIC y PrimateAI. Los eggNOG de **FCS por gen** se conservan. COSMIC, PrimateAI y Ensembl VEP siguen como anotaciones descriptivas del informe.

`POSITION_ANNOTATIONS` construye estas colecciones y `position_characterization.tsv` en `posenrich/gmts/`, sin nulo, contraste ni informe 14. El informe 15 recibe esta tabla aunque POSENRICH esté omitido. Las tablas Pfam/UCR indicadas se reutilizan. Con la casilla de la capa activada y la ruta vacía, se genera automáticamente sobre `--alignment`; desactivar la casilla omite la capa. Pfam usa hmmscan y descarga/cachea Pfam-A si falta; UCR calcula Valdar y las ventanas globales sliding. Las salidas y manifiestos se publican en `core_inputs/`. La selección de prueba y nulo posicional sigue pendiente de WP10. No se reconstruye el TSV antiguo del nulo ni se altera acumulación/FCS.

## Coordenadas y contratos

El identificador final es `Gene:Position`, con `Position` **cero-based en el alineamiento recortado de CAAS**. El MAP se lee por encabezados: `ori_codon_col`, `status`, `prot_ali_col`, `hg38_aa_pos`, `hg38_nt_coord`. Solo filas selected; `prot_ali_col - 1` es Position. Un directorio con dos MAP del mismo gen se rechaza, sin escoger una referencia arbitraria.

| Fuente | Entrada y proyección |
|---|---|
| Pfam generado por `bin/compute_domain_variability.py` | El nuevo productor añade `alignment_positions` (lista separada por comas) y `coordinate_space=CAAS_0based`. Conserva las columnas exactas del alineamiento suministrado al productor, excluyendo gaps y X del query hmmscan. **Ese alineamiento debe ser el mismo alineamiento recortado de CAAS.** |
| Pfam externo | También admite `coordinate_space=reference_protein_1based`: ali_start/ali_end deben corresponder a la misma proteína de referencia del MAP. Se recuperan los residuos interiores conservados aunque se hayan eliminado los extremos. Las tablas sin coordinate_space se rechazan. En el PPGA verificado de 20260615_200027, los extremos hmmscan corresponden a la consulta humana recortada sin gaps; no equivalen a hg38_aa_pos original. Se convierten primero mediante el FASTA productivo a alignment_positions explícitas. |
| UCR | Requiere coordinate_space explícito. CAAS_0based se conserva; trimmed_alignment_1based se convierte únicamente restando uno; raw_alignment_1based solo se proyecta con ori_codon_col si su procedencia cruda está verificada. En el PPGA 20260615_200027, position indexa el MSA ya recortado, uno-based: se convierte a position − 1 y se declara CAAS_0based. Solo método sliding. La generación automática ya emite CAAS_0based. |
| FUBAR | `hg38_aa_pos` se busca en el número explícito de residuo del MAP; no se cuenta la sucesión de filas selected. |
| FADE | CSV `gene,position,max_bf,target_aa`, producido por FADE_TABLES. position ya es CAAS. Se respeta la selección upstream, sin otro umbral de BF. |
| Ventanas | Coordenada genómica del MAP; bins de 1.000.000 nt, separados por cromosoma, con extremos uno-based inclusivos (1–1.000.000, 1.000.001–2.000.000, etc.). Es proximidad genómica, distinta de proximidad en proteína. |
| Personalizados | `--posenrich_custom_file`: TSV sin encabezado, gene, Position, term y descripción opcional. Comentarios #. Position ya debe estar validada en coordenadas CAAS; no acepta posiciones negativas ni fuera del MAP disponible. |

La anotación Pfam conserva todos los dominios/clanes solapantes, separados por `;`. Ausencia de una anotación no demuestra ausencia biológica. El universo de genes limita todas las fuentes. Sin MAP no se pueden proyectar UCR/FUBAR ni producir ventanas; sí se pueden conservar anotaciones que ya llevan coordenadas CAAS verificadas.

No se borran resultados históricos publicados. Una carpeta antigua puede contener GMT retirados: no debe usarse su inventario como salida nueva. El flujo nuevo solo emite sus archivos declarados, y permanece separado del motor estadístico histórico.

## Recuperación histórica de las referencias PPGA (25/10/2026)

La solución productiva es que el generador de alineamientos exporte `alignment_positions`/`coordinate_space` para Pfam y `trimmed_alignment_1based` para UCR; PhyloPhere consume esas tablas directamente. La utilidad siguiente es una recuperación histórica, no un paso requerido del flujo productivo.

`python bin/build_reference_tables.py project-ppga --base <results_ppga/20260615_200027> --output <nueva_version>` selecciona explícitamente el contrato del productor PPGA verificado. No es un detector genérico del sistema de coordenadas. Conserva los originales y rechaza sobrescribir una versión existente.

Valida que cada columna del MSA tenga MAP y entropía, convierte UCR sin segunda proyección, y contrasta cada extremo/query length Pfam con el domtblout preservado y la secuencia humana del MSA. En ese productor hmmscan eliminaba solo `-`: X consume posición. La máscara Pfam conserva residuos de consulta, excluye gaps de referencia y recalcula sus resúmenes de variabilidad. FUBAR conserva su sitio productivo y se contrasta con las coordenadas originales del MAP. Las tablas `.tsv.gz` y sus manifiestos vinculan fuentes, hashes e inventarios completos de alineamientos y MAP; solo `reference_bundle_manifest.json` con status complete acredita que terminó el paquete.

La tabla `ucr_positions.tsv.gz` exporta el método sliding, el mismo que consume POSENRICH. En una fuente externa de otro origen se debe declarar y verificar su contrato o regenerar sobre el alineamiento productivo; nunca asumir el contrato PPGA por el nombre del archivo.

Decisión del usuario, 25/10/2026: la pregunta principal de localización es si la señal se concentra en los sitios anotados **dentro de cada proteína**, más allá de su carga general de señal. Los controles inferenciales siguen pendientes. Orden: corregir coordenadas y reconstruir resultados dependientes; validar fondos y controles; optimizar eficiencia al final. Esta reparación no activa inferencia ni recalcula nulos CAAS.

## FADE y VEP en informes

El constructor ahora lee los CSV seleccionados FADE y produce las capas globales y sus flags por posición. La regeneración GUI de informes 14/15 busca esos CSV, en lugar del TSV de BF con otro delimitador. La regeneración de 15 incorpora también position_characterization.tsv.

La entrada Ensembl VEP del informe 15 ya estaba conectada. Se conserva el join por gen, posición y dirección, y la elección de transcrito MANE Select, después canonical y después consecuencia más severa. Se prefiere su archivo canónico `vep/ensembl_vep_mapped.tsv` al descubrirlo en GUI. Las anotaciones VEP no necesitan reactivar POSENRICH.

## Tres estudios PP2A

`bin/prepare_pp2a.py --source-dir <carpeta Excel> --output-dir <salida>` requiere openpyxl. No modifica los originales. Produce:

- `source_records.tsv.gz`: fila de origen completa, incluyendo encabezados repetidos y abundancias, con ID de registro.
- `site_candidates.tsv`: candidatos por proteína/residuo, localización reportada, estudio, hoja, fila, condición y papel de la tabla (hits/fondo medido).
- `unresolved.tsv`: registros sin sitio numérico o con localización incompleta, conservando el texto original.
- `manifest.json`: huellas de los Excel, recuentos y límites de interpretación.

| Estudio | Evidencia conservada |
|---|---|
| Kauko | Dos hojas distintas: disminución tras activación por depleción de inhibidores (CIP2A/PME-1/SET, respetando las celdas combinadas) e incremento tras siPP2A-A. No se mezclan sus perturbaciones. La columna Modulation tiene muchas celdas vacías sin combinar: permanecen unspecified, sin propagar automáticamente la etiqueta anterior. Esas asignaciones por bloque requieren revisar el original. |
| Stewart | Incremento asociado a PPP2R1A R183P; abundancias de WT/mutante y probabilidades de localización permanecen en la fila original. No se inventa un p-valor a partir de ellas. |
| Sapkota | Hoja Hyper-phosphorylated y Entire Data separadas. La segunda aporta cobertura de péptidos medidos; no equivale a todas las serinas/treoninas del proteoma. |

No se imponen nuevos cortes de p, FDR o probabilidad de localización. No se equipara fosforilación asociada a una perturbación con sustrato directo de PP2A. Las filas multipéptido/multiposición o con mappings alternativos no son observaciones independientes; `ambiguous_mapping` impide tratarlas como sitios inequívocos. Incluso `reported_site` solo significa que la hoja da una localización: falta validar secuencia y referencia.

**Puerta pendiente para probar solapamiento/cercanía con CAAS:** validar accession/isoforma y secuencia del péptido contra la referencia efectivamente usada en cada gen, conservar o resolver mappings alternativos y proyectar por su MAP real. La carpeta local suministrada contiene los estudios, pero no los MAP de la referencia productiva: la plantilla neoplasia apunta a `/homes/users/mramon/scratch/4.Generate_alignments_from_codons/results_ppga/20260615_200027/MAP`. No se sustituye ese MAP por otro toy ni se iguala residuo UniProt con Position. Los archivos preparados todavía **no son GMT listos para CAAS**. El test inferencial, el fondo por estudio y el criterio de cercanía se deben fijar antes de interpretar significación.

## Validación dirigida de esta implementación

Evidencia local: `PhyloPhere_validation/.tmp/posenrich_redesign_20261023/`.

- `tests/test_position_collections.py` y carga/guardado/render de todas las plantillas: **7 PASS**. Incluye coordenadas con columnas eliminadas, gaps, primer residuo distinto de 1, frontera genómica de 1 Mb, dominios solapantes, Pfam con X/gaps, ausencia de opcionales, FADE, retirada de colecciones y migración del selector eggNOG.
- Una ejecución real de `POSITION_ANNOTATIONS` en Nextflow, con FADE y el resto de capas ausentes: **COMPLETED**, tabla publicada y sin llamar al motor inferencial. El primer harness puso outdir después del include y produjo la tabla bajo null/; corregido pasando outdir antes de evaluar el módulo. No era un fallo del proceso.
- `Rscript tests/test_position_report.R <repo> <scratch>`: render completo del informe 15 y aserciones sobre preferencia MANE, VEP distinto por dirección, dominios solapantes y FADE top sin contagio a bottom.
- Compilación de main.nf mediante ayuda: termina con código 1 deliberado de HELP, sin error de compilación; no es una ejecución integrada de enriquecimiento.
- Normalización de los tres Excel reales; SHA-256 de originales preservados. Las tablas derivadas y sus recuentos se entregan bajo outputs/pp2a_prepared de esta conversación.

No se han ejecutado nuevos nulos, una corrida completa, trabajos en cluster ni una prueba biológica de enriquecimiento PP2A. Las modificaciones ajenas del checkout se conservan.

## Registro común y componentes modulares (WP10, `--posenrich_modular`)

`subworkflows/ENRICHMENT/posenrich_modular.nf` orquesta cinco procesos sobre el registro
común (`posenrich/registry/`), que es la tabla normalizada única; los GMT de
`POSITION_ANNOTATIONS` son una exportación compatible (el registro exporta los mismos GMT
byte a byte, verificado en `tests/test_posenrich_registry.py`).

| Componente | Script | Publica | Estado |
|---|---|---|---|
| Registro | `local/src/posenrich_registry.py` | `posenrich/registry/`: positions.tsv (estados de evaluación + coordenadas), position_evidence.tsv (score direccional + p.emp compartido), collections.tsv, members.tsv, coverage.tsv, pending_mappings.tsv, manifest.json (SHA-256 de fuentes) | implementado y validado dirigido |
| Caracterización | `local/src/posenrich_characterize.py` | `posenrich/characterization/`: position_characterization.tsv (propietario único para informes 14/15), collection_summary.tsv, ranking_summary.tsv, universe_summary.tsv | implementado y validado dirigido |
| Localización | `local/src/posenrich_localization.py` | `posenrich/localization/`: localization_results.tsv (AUC por gen, agregación mean principal + pair_weighted, control condicionado dentro de proteína, p Monte Carlo unilateral (k+1)/(B+1)), localization_by_gene.tsv | implementado y validado dirigido (oráculos AUC, calibración sintética del control) |
| Asociación global | `local/src/posenrich_global.py` | `posenrich/global/`: global_results.tsv (AUC vs fondo global, composición por gen, LOO; **descriptivo, sin p**: la calibración por ciclo no es reproducible desde contadores marginales) | implementado y validado dirigido |
| Proximidad | `local/src/posenrich_proximity.py` | `posenrich/proximity/`: distancias proteicas (residuos de referencia) y genómicas (nt), separadas; descriptivo, control de distancias no validado | implementado y validado dirigido |
| Informe | `local/posenrich_modular_report.Rmd` | `posenrich/POSENRICH_modular_report_<trait>.html` + `html_reports/` | render verificado dirigido |

`build_position_gmt.py` quedó refactorizado en builders reutilizables
(`collect_*_members` con membresía normalizada y mapeos pendientes opcionales); su salida
no cambió (`tests/test_position_collections.py` sigue pasando). La inferencia histórica
POSENRICH (Lachenbruch) sigue omitida: `--posenrich` permanece saltado con aviso y el
informe 14 inferencial no se produce.

Pruebas dirigidas: `tests/test_posenrich_registry.py`, `test_posenrich_characterize.py`,
`test_posenrich_localization.py` (incluye calibración sintética del control bajo nulo),
`test_posenrich_global.py`, `test_posenrich_proximity.py`,
`tests/test_posenrich_modular_report.R` (render del informe) y
`tests/test_posenrich_modular_smoke.py` (smoke Nextflow del subworkflow completo).

Pendiente: validación del piloto PP2A (el MAP productivo de la plantilla neoplasia no es
accesible desde este entorno), aceptación integrada 1000/4000/10000 en cluster, y la
familia BH de conjuntos (analysis × scope × ranking × collection) cuando se active la
inferencia a escala.


## Revisión del contrato modular

Cuando `posenrich_modular=true`, el informe 15 recibe la tabla de
`posenrich/characterization/position_characterization.tsv`; el constructor antiguo
solo se ejecuta con el módulo apagado. La regeneración GUI prefiere esa tabla.
El registro conserva discovery sin score en `position_evidence.tsv`, con
`raw_asr_path_score`, `score_status` y CAAS_score=NA si el productor canónico no lo
puntúa. Nunca convierte MAP sin evaluación, valores no positivos o faltantes en
score cero ni en p=1. AUC y proximidad son descriptivos sobre scores disponibles.
FUBAR restringe el fondo de AUC a `fubar_tested`; la cobertura experimental de otras
fuentes se declara desconocida. La inferencia de localización permanece retenida.

## Tablas Pfam antiguas sin sistema de coordenadas

El lector vuelve a aceptar el formato original de PPGA sin modificar la tabla ni
exigir una columna nueva. Se recupera la proyección histórica del intervalo de
residuos de referencia mediante MAP, manteniendo la lectura por nombres de columna
y la numeración explícita de residuos (no se cuentan las filas después del trimming).
Se conservan los dominios solapantes y los dominios sin clan.

Los formatos explícitos nuevos siguen admitidos. `bin/migrate_pfam_coordinates.py`
sigue disponible como conversión opcional para queries verificables; no es un
requisito para leer las tablas originales.

Con `posenrich_domains=true` / `posenrich_ucr=true`, un archivo vacío activa la
generación; un archivo indicado se reutiliza y se valida su encabezado. FUBAR
permanece precomputado y no se vuelve a ajustar. La generación UCR automática usa
solo las ventanas globales sliding consumidas por las colecciones, sin estadísticas
por clado. Las tablas externas deben corresponder a su referencia/MAP: aceptar un
formato antiguo no establece equivalencia entre los alineamientos PPGA y marine.

La especificación de la pregunta principal está en [POSENRICH_WITHIN_PROTEIN.md](POSENRICH_WITHIN_PROTEIN.md).
