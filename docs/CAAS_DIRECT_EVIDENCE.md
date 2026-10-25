# Evidencia CAAS directa — decisión del 25/10/2026

El usuario ha retirado el factorizado del flujo activo. La evidencia principal es
p.emp por posición y p.gen por gen, con nominal <0,05 y BH <0,05. FCS conserva su
corte propio; no se cambia el estadístico génico size_adj_max por esta retirada.

SCORING deja de producir p.emp_fact/p.adj_bh_fact. El contador por bloque calcula
solo k_emp; el reductor añade nd. Los nombres de staging null_classes.parquet y
null_fit.tsv se conservan para compatibilidad de Nextflow: contienen una sola
clase inerte, sin pooling ni inferencia condicional. El argumento legacy fit del
helper R se ignora. Los contadores históricos con c1..cK siguen siendo legibles.
Los informes SCORING y CAAS eliminan las columnas retiradas al importar archivos
antiguos. POSENRICH no las propaga en su registro. La ausencia de p sigue siendo
NA; nunca se rellena con evidencia factorizada.

Las herramientas bin/audit_caas_evidence.py y bin/audit_caas_gene_evidence.R
conservan explícitamente el método histórico para reproducir la auditoría T2.
Los planes/documentos históricos que lo describen no definen el contrato activo.
No se incorpora ningún p.gen_fact productivo.

La comparación bin/compare_caas_gene_statistics.R evalúa máximo crudo, F(max) y
F(max)^n contra los mismos 10.000 nulos propios de cada gen, con BH por ranking
y universo limpio. No cambia la selección productiva. La evaluación en ciclos
reservados es un control de calibración bajo el modelo nulo guardado; no mide
potencia frente a alternativas biológicas ni demuestra calibración del FCS.
