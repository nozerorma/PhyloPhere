# Colecciones posicionales y preparación de PP2A

Referencia de las fuentes de posiciones que lee POSENRICH y de la preparación de los estudios PP2A. El diseño del módulo está en `POSENRICH_DISENO.md`.

Colecciones: Pfam dominios y clanes, UCR core y flancos, FUBAR positiva y purificadora (más la cobertura `fubar_tested`), FADE top y bottom, ventanas genómicas y listas personalizadas. COSMIC, PrimateAI y Ensembl VEP son anotaciones descriptivas del informe 15, no colecciones.

El registro (`posenrich_registry.py`, `--posenrich` con `--scoring`) lee las fuentes con `subworkflows/POSENRICH/local/src/position_sources.py`, las proyecta a coordenadas CAAS y escribe las colecciones en `posenrich/registry/`; `posenrich/characterization/position_characterization.tsv` lo produce la caracterización. El informe 15 de ENRICHMENT recibe esa tabla como entrada opcional; sin ella omite las columnas de caracterización. Las tablas Pfam y UCR indicadas se reutilizan. Con la casilla de la capa activada y la ruta vacía, se genera sobre `--alignment`; desactivar la casilla omite la capa. Pfam usa hmmscan y descarga y cachea Pfam-A si falta; UCR calcula Valdar y las ventanas globales sliding. Las salidas y manifiestos se publican en `core_inputs/`.

## Coordenadas y contratos

El identificador final es `Gene:Position`, con `Position` **cero-based en el alineamiento recortado de CAAS**. El MAP se lee por encabezados: `ori_codon_col`, `status`, `prot_ali_col`, `hg38_aa_pos`, `hg38_nt_coord`. Solo filas `selected`; `prot_ali_col - 1` es Position. Un directorio con dos MAP del mismo gen se rechaza, sin escoger una referencia arbitraria.

| Fuente | Entrada y proyección |
|---|---|
| Pfam generado por `bin/compute_domain_variability.py` | Añade `alignment_positions` (lista separada por comas) y `coordinate_space=CAAS_0based`. Conserva las columnas exactas del alineamiento suministrado al productor, excluyendo gaps y X del query hmmscan. **Ese alineamiento debe ser el mismo alineamiento recortado de CAAS.** |
| Pfam externo | Admite `coordinate_space=reference_protein_1based`: ali_start y ali_end deben corresponder a la misma proteína de referencia del MAP. Se recuperan los residuos interiores conservados aunque se hayan eliminado los extremos. Las tablas sin `coordinate_space` se rechazan; las originales de PPGA se convierten antes con `bin/migrate_pfam_coordinates.py` o `bin/build_reference_tables.py` (ver abajo). |
| UCR | Requiere `coordinate_space` explícito. `CAAS_0based` se conserva; `trimmed_alignment_1based` se convierte restando uno; `raw_alignment_1based` solo se proyecta con `ori_codon_col` si su procedencia cruda está verificada. Solo método sliding. La generación automática emite `CAAS_0based`. |
| FUBAR | `hg38_aa_pos` se busca en el número explícito de residuo del MAP; no se cuenta la sucesión de filas selected. |
| FADE | CSV `gene,position,max_bf,target_aa`, producido por FADE_TABLES. `position` ya es CAAS. Se respeta la selección upstream, sin otro umbral de BF. |
| Ventanas | Coordenada genómica del MAP; ventanas de `--window_size_bp` (valor de `scoring_window_size_bp`), separadas por cromosoma, con extremos uno-based inclusivos. Es proximidad genómica, distinta de la proteica. |
| Personalizados | `--posenrich_custom_file`: TSV sin encabezado: gene, Position, term y descripción opcional. Comentarios `#`. Position ya debe estar validada en coordenadas CAAS; no acepta posiciones negativas ni fuera del MAP disponible. |

La anotación Pfam conserva todos los dominios y clanes solapantes. Un dominio sin clan solo genera la región del dominio. La ausencia de una anotación no demuestra ausencia biológica. El universo de genes limita todas las fuentes. Sin MAP no se pueden proyectar UCR ni FUBAR ni producir ventanas; sí se pueden conservar anotaciones que ya llevan coordenadas CAAS verificadas. Los registros que no se proyectan quedan en `pending_mappings.tsv`, nunca como negativos.

Las tablas externas deben corresponder a su referencia y su MAP: aceptar un formato antiguo no establece equivalencia entre alineamientos de distinto origen. FUBAR es una entrada precomputada y no se vuelve a ajustar.

## Referencias PPGA existentes

La solución productiva es que el generador de alineamientos exporte `alignment_positions` y `coordinate_space` para Pfam y `trimmed_alignment_1based` para UCR, de modo que PhyloPhere consuma esas tablas directamente. Para tablas PPGA ya generadas:

`python bin/build_reference_tables.py project-ppga --base <results_ppga/20260615_200027> --output <nueva_version>` selecciona el contrato del productor PPGA verificado. No es un detector genérico del sistema de coordenadas. Conserva los originales y rechaza sobrescribir una versión existente.

Valida que cada columna del MSA tenga MAP y entropía, convierte UCR sin segunda proyección y contrasta cada extremo y longitud de consulta Pfam con el domtblout preservado y la secuencia humana del MSA. En ese productor hmmscan eliminaba solo `-`: X consume posición. La máscara Pfam conserva residuos de consulta, excluye gaps de referencia y recalcula sus resúmenes de variabilidad. FUBAR conserva su sitio productivo y se contrasta con las coordenadas originales del MAP. Las tablas `.tsv.gz` y sus manifiestos vinculan fuentes, hashes e inventarios completos de alineamientos y MAP; solo `reference_bundle_manifest.json` con status complete acredita que el paquete terminó.

`bin/migrate_pfam_coordinates.py` convierte tablas Pfam con consultas verificables. Las tablas Pfam y UCR sin `coordinate_space` se rechazan en la entrada del workflow (`workflows/posenrich.nf`).

## FADE y VEP en los informes

La regeneración de los informes 14 y 15 desde la GUI busca los CSV de sitios FADE seleccionados (no el TSV de BF). La del informe 15 incorpora además `position_characterization.tsv`. La entrada Ensembl VEP del informe 15 se une por gen, posición y dirección, y elige el transcrito MANE Select, después el canónico y después la consecuencia más severa; se prefiere `vep/ensembl_vep_mapped.tsv`. Las anotaciones VEP no dependen de POSENRICH.

## Tres estudios PP2A

`bin/prepare_pp2a.py --source-dir <carpeta Excel> --output-dir <salida>` requiere openpyxl. No modifica los originales. Produce:

- `source_records.tsv.gz`: fila de origen completa, incluyendo encabezados repetidos y abundancias, con ID de registro.
- `site_candidates.tsv`: candidatos por proteína y residuo, localización reportada, estudio, hoja, fila, condición y papel de la tabla (hits o fondo medido).
- `unresolved.tsv`: registros sin sitio numérico o con localización incompleta, conservando el texto original.
- `manifest.json`: huellas de los Excel, recuentos y límites de interpretación.

| Estudio | Evidencia conservada |
|---|---|
| Kauko | Dos hojas distintas: disminución tras activación por depleción de inhibidores (CIP2A, PME-1, SET; respetando las celdas combinadas) e incremento tras siPP2A-A. No se mezclan sus perturbaciones. La columna Modulation tiene muchas celdas vacías sin combinar: permanecen `unspecified`, sin propagar la etiqueta anterior. Esas asignaciones por bloque requieren revisar el original. |
| Stewart | Incremento asociado a PPP2R1A R183P; abundancias de WT y mutante y probabilidades de localización permanecen en la fila original. No se inventa un p-valor a partir de ellas. |
| Sapkota | Hojas Hyper-phosphorylated y Entire Data separadas. La segunda aporta cobertura de péptidos medidos; no equivale a todas las serinas y treoninas del proteoma. |

No se imponen cortes de p, FDR ni probabilidad de localización. No se equipara fosforilación asociada a una perturbación con sustrato directo de PP2A. Las filas multipéptido, multiposición o con mapeos alternativos no son observaciones independientes; `ambiguous_mapping` impide tratarlas como sitios inequívocos. Incluso `reported_site` solo significa que la hoja da una localización: falta validar secuencia y referencia.

**Puerta pendiente para contrastar solapamiento o cercanía con CAAS:** validar accession, isoforma y secuencia del péptido contra la referencia usada en cada gen, conservar o resolver mapeos alternativos y proyectar con su MAP real (`bin/map_pp2a_to_caas.py`). No se iguala el residuo UniProt con Position. Los archivos preparados no son listas listas para CAAS. El control de proximidad, el fondo por estudio y el criterio de cercanía se fijan antes de interpretar significación (`POSENRICH_DISENO.md`, sección 11).
