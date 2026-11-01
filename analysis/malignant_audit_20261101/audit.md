# Comparación malignant_def / malignant_def2 — 2026-11-01

## Conclusión

Hay diferencias metodológicas y diferencias numéricas observadas en CAAS_CORE_BATCHED. Para obtener resultados con el código actual no basta con reutilizar el CAAS antiguo y regenerar solo downstream: observado y nulo deben proceder del mismo contrato científico. Conservar malignant_def2 como versión histórica completa.

## Fuentes consultadas

Inspección remota de solo lectura en mramon@correfoc-vs.s.upf.edu:

- Nuevo work: `/homes/users/mramon/scratch/3.Work_dirs_malignant_def/malignant_prevalence_complete`.
- Nuevo log y params: directorio hermano `malignant_prevalence_complete_nxf_run`.
- Resultados antiguos: `/homes/users/mramon/scratch/2.Primates/2.Primates_results/CAAS_RESULTS/malignant_def2/malignant_prevalence_complete`.
- Resultados nuevos: misma raíz CAAS_RESULTS, `malignant_def/malignant_prevalence_complete`.
- Repo remoto HEAD: `0df06f0d005321f7236a5ca186a6c255e69e0332`, coincidente con el local. Hay cambios sin commit en driver.py de disambiguation; no se han modificado.
- Historia local: especialmente `5dc078ce` (29 octubre, 18:17) y `3defc510` (31 octubre, 03:55).

Copias consultadas guardadas en `.tmp/malignant_audit/` del workspace. No se han lanzado análisis ni modificado archivos remotos.

## Estado de las ejecuciones

malignant_def2 contiene discovery, background, nulo b0 y 20 bloques, scoring, FCS, POSENRICH, comparación e informes. Su metadata declara 15.818 genes y 10.000 ciclos. counts.json declara score_rule `caas_pos_species_pss`, 92.052.791 unidades y 1.455.593 unidades eliminadas. Los bloques CAAS tienen fecha del 26 de octubre; algunos productos globales heredados son del 19 de octubre. El log conservado corresponde a un relanzamiento downstream del 29 de octubre: finaliza con 1 proceso completado, 354 cacheados y 0 fallidos. Las rutas internas todavía apuntan a malignant_def; por tanto, el nombre actual malignant_def2 no identifica una ejecución independiente ni demuestra cuándo se generó cada producto.

El nuevo run arranca el 1 de noviembre a las 12:57, ignora resume por ausencia de historial y vuelve a generar CAAS. La copia del trace inspeccionada registra 282 lotes CAAS_CORE_BATCHED y 62 CAAS_UNIVERSE completados, sin merge final de CAAS ni scoring posterior. El log termina a las 16:57 sin cierre de workflow. Esto acredita ejecución parcial, no su estado actual en Slurm: squeue/sacct no están disponibles en el PATH de la sesión SSH.

## Diferencias relevantes

1. Los nuevos tiles declaran `caas_pos_species_pss_v4` en null_tables/meta.json. El nulo anterior declara `caas_pos_species_pss`. `subworkflows/SCORING/local/src/null_position_stats.py` rechaza la regla antigua y exige regenerar observado y nulo juntos.
2. El 31 de octubre se introduce CAAS_UNIVERSE: posiciones testeadas por b0 capaces de pasar el test bajo alguna orientación de los pares de alguna hipótesis observada. Ese universo restringe el score nulo antes del máximo por gen, tiene una versión strict y se incorpora a la procedencia de los tiles. No es solo un cambio de formato.
3. Se separa el solapamiento de residuos entre foreground/background del número e identidad de pares realmente conservados. Se corrige la asociación por ID de par ante gaps/missingness asimétricos. No se ha cuantificado cuánto afecta esa corrección específica a este dataset.
4. Desde el 29 de octubre cambia el filtro de genes a umbrales estratificados por longitud. El run nuevo solicita stratified, 10 estratos, k_far=3 y k_near=1,5. La eliminación de genes se resuelve en merge, pero el filtrado de clusters también opera en los tiles; el nuevo tile registra minlen=3, min_positions=3, maxcaas=0,7. El antiguo nulo publicado precede a este cambio.
5. La generación/validación de permulaciones cambia: reutilización del modelo observado, comprobaciones del diseño y rechazo de relajaciones no autorizadas. Los manifests antiguos/nuevos ya difieren en b1 y b2; esa diferencia por sí sola no mide un efecto biológico. El nuevo ct_preparation.json registra 10.000 ciclos, 72 hipótesis, OU reutilizado del observado y ninguna relajación. match_pss está solicitado pero no aplicado, coherente con el tratamiento de traits de conteo del código.
6. También cambian resolución de nombres de genes/longitudes/MAP, validaciones de procedencia y tratamiento de fallos de scoring. No todas estas modificaciones implican cambios de resultado para todos los genes.

## Comprobación empírica del primer lote observado

Tile nuevo: `5e/7dd4474f818c1afed35128b4b08230`, caas_core_batch_00001, b0, 200 genes.

Comparación local de su null_tables/scores.parquet con el antiguo caas_permulation/b0/positions.parquet. Se normaliza únicamente el sufijo `.Homo_sapiens` de los IDs nuevos y se calcula el score nuevo como máximo de top/bottom finitos. Este corresponde al colapso de direcciones de la función collapse_sides actual. Se comparan solo claves gen/posición presentes en ambos productos, ya que el producto antiguo está después del merge y el nuevo es un tile previo al filtrado global de genes.

- 406 posiciones compartidas.
- 43 scores cambian en más de 1e-9; 363 coinciden con esa tolerancia.
- A1BG, posición 186: 0,6607797606785337 → 0,6137109154590084.
- A1BG, posición 229: 0,05119302195235501 → 0,03562031947034417.
- Las 700 filas de discovery de A1BG coinciden exactamente en las columnas compartidas tras normalizar el sufijo del gen. El nuevo discovery añade shared_residue_count. Por tanto, igualdad de detecciones iniciales no garantiza igualdad del score posterior.

Todas las posiciones de esta comparación son índices cero del alineamiento recortado. Tabla reproducible de valores: `sample_score_comparison.tsv`. No es una comparación global de rankings, p-valores o conclusiones biológicas, y no permite atribuir cada delta a un cambio de código concreto.

## Recomendación

Conservar malignant_def2 como referencia histórica. Para trabajar con la metodología actual, completar o reanudar de forma coherente la generación nueva de CAAS (universo, tiles observado/nulo y merges), y después regenerar sus pasos posteriores. No mezclar scores observados v4 con el nulo antiguo. Esta auditoría no ha lanzado ni cancelado ningún job.
