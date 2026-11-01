# Reutilización para completar CAAS v4

Revisión de solo lectura: 1 de noviembre de 2026. No se ha preparado ni ejecutado un relanzamiento remoto.

## Material reutilizable comprobado

Raíz del work nuevo: `/homes/users/mramon/scratch/3.Work_dirs_malignant_def/malignant_prevalence_complete`.

- Diseño observado completo y pesos: `6a/2d5221770e3724145e48b0f3860cc7/traitfiles_ok_dir/`. Contiene traitfile_H*.tab, contrast_hypotheses_pairs.tsv y observed_selection.rds.
- Harvest nuevo: `13/baa9ba680b8e374d24e7384da5dfaa/nw_tree.resampled.output/`. Contiene fop_labelings.tab, fop_pairs.tsv, manifests y permulation_provenance.json. Usar el directorio completo: caastools/resample.tab por sí solo no conserva todas las hipótesis y pesos requeridos.
- Bloques ya preparados: `c0/f5544c*/blocks/`, además de ct_preparation.json. Diseño/harvest han pasado ya la preparación v4: 10.000 ciclos, 72 hipótesis, sin relajación.
- CAAS_UNIVERSE y CAAS_CORE_BATCHED completados: registrados en execution_trace.txt y en la base Nextflow del run.
- ASR: `/homes/users/mramon/scratch/2.Primates/1.Primates_data/7.ASR_primates/validated_v1/`. Unas 9.000 carpetas de claves al inspeccionarlo. Cada entrada debe pasar las comprobaciones de procedencia y checksum que ya hace el código; no se ha auditado individualmente todo el directorio. El log del primer lote confirma uso de resultados de ese namespace.
- FADE de malignant_def2: el run nuevo ya está configurado para consumirlo como precomputado.
- Trait, árboles y matriz RER nuevos: ya generados; preservar su caché. No se ha certificado aquí equivalencia científica del resultado RER antiguo con la implementación actual.

## Lo que no permite evitar el recálculo

Los scores/matrices antiguos son reducciones con el contrato anterior. No conservan toda la información necesaria para reconstruir los scores v4 por una simple conversión, filtrado de posiciones o cambio de metadatos. Hay deltas de score aun con discovery idéntico. El cambio de umbrales de eliminación de genes se resuelve en merge si se conservan tiles compatibles; cambios de score, universo y clusters requieren revisar/recalcular los tiles afectados.

## Orden recomendado

1. Mantener coherentes diseño, harvest y scoring del run nuevo; evitar importar el nulo final antiguo.
2. Recuperar los tiles v4 ya terminados y recalcular solo los faltantes o realmente incompatibles. Para ello hay que resolver primero la invalidación de caché descrita abajo. Esto exige un ajuste revisado de caché o una importación validada; no se ha implementado ni probado aquí.
3. Mantener la caché ASR validada para evitar repetir reconstrucciones compatibles, incluso si se repiten tiles.
4. Ejecutar merges v4, postprocesado y downstream a partir de ese conjunto coherente. Conservar FADE precomputado y cachés de módulos independientes cuando sean compatibles.

Si se decide iniciar un workflow de recuperación separado, la interfaz actual admite `caas_config` con el diseño guardado, `resample_from` con el directorio completo del harvest y `ct_tool=discovery` sin `resample`, con `contrast_selection=false`. Hay que mantener árbol, alineamientos, MAP, taxid y parámetros científicos del run nuevo. Cambiar a esa interfaz puede cambiar las claves de los tiles: no sustituye por sí sola la recuperación validada de los ya completados. No se ha probado un lanzamiento de esa configuración.

## Por qué un resume ciego puede desperdiciar trabajo

CAAS_CORE_BATCHED incorpora coreScoringSourceKey(), que recorre todos los .py y .dat de CT_DISAMBIGUATION/local, incluyendo tests y driver.py.

SHA256 de la parte de fuentes observadas/replay del primer tile:

`d05c0e84499783eeb8071fb0f59fdf75062eb9ae9a8ca7269295ae67118fa79e`

La misma parte calculada con el workspace actual:

`f56cb76f2eaa11628edac05e9cdf48f617baee592786dc2003ec4a7e77c7945a2b`

Al excluir el nuevo test_asr_cache_race.py y usar la versión HEAD de driver.py, la huella vuelve exactamente a la del tile. El diff local de driver.py introduce recuperación de una entrada ASR validada publicada por otro worker tras una carrera; no modifica la fórmula de scoring. Esto demuestra el origen de esa diferencia de huella local, pero no basta por sí solo para certificar todos los tiles remotos o cambiar sus claves sin validación.

La clave ASR incluye los archivos asr_single.py, reconstruct.py, posterior.py, tree_parser.py y cache.py, no driver.py ni ese test. Por tanto, esa corrección concreta no invalida por sí misma las entradas ASR.

No se recomienda editar los stamps antiguos ni rebautizar la regla del nulo para eludir controles.
