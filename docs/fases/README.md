# Ejecución por fases: arranque, cierre y registro

Cada paquete de trabajo (WP) se ejecuta en una conversación independiente. Este directorio es lo que debe bastar para arrancar una conversación nueva sin contexto previo y para dejar constancia de lo hecho.
Plan aprobado (copia versionada): `docs/fases/PLAN_APROBADO_2026-10-21.md`. Diagnóstico y mediciones: `docs/PLAN_NULO_10K_PERMULACIONES.md`. Revisión del FCS: `docs/FCS_REVISION_RERCONVERGE.md`. Porte de POSENRICH: `docs/POSENRICH_PORTE_EVALUACION.md`.

## Estado

| WP | Fichero | Fase del plan de 10.000 ciclos | Depende de | Estado |
|---|---|---|---|---|
| WP0 | `WP0_linea_base.md` | | commit del usuario (hecho: `581b999`) | HECHO. PEPC y cancer_b0_toy regenerados con `a92fd2d`; marine toy completada el 2026-10-22 (WP5) |
| WP1 | `WP1_fcs_rer.md` | | WP0 | EN REVISIÓN (RER: NA excluidos, bilateralidad explícita y ajustes observado/nulo compartidos; pruebas y procesos Nextflow pequeños validados; pendientes restantes y decisión a 10.000 en el registro) |
| WP2 | `WP2_exports_posenrich_limpieza.md` | | WP0 | HECHO (PEPC c4 y toy de cancer idénticos a la referencia con su caché de ASR; marine toy validado en WP5) |
| WP3 | `WP3_bloques_de_ciclos.md` | Fase 1 | WP2 | HECHO (PEPC c4 a 1000 ciclos con bloque de 250 y toy de cáncer con bloques de 25 y 40 idénticos a la referencia, delta 0.0; marine toy validado en WP5) |
| WP4 | `WP4_replay_en_flujo.md` | Fase 2 | WP3 | IMPLEMENTADO; VALIDACIÓN DE ESCALA ABIERTA (90c06a4; equivalencia CAAS validada; referencia escalar reparada, 89 pruebas focalizadas pasan; faltan RP1 B=500 y resolver el criterio de memoria) |
| WP5 | `WP5_tablas_y_merge_por_bloque.md` | Fase 3 | WP4 | REVISADO en local: guardas de ciclos corregidas; 68 controles dirigidos pasan; reconstrucciones PEPC/cáncer/marine completo y strict reproducen matrices y contadores; cableado 38 PASS y 7 fallos ya registrados. Pendientes adaptar esos tests y validar escala real en el cluster; scripts renderizados, sin lanzar en esta revisión |
| WP6 | `WP6_fcs_por_bloques.md` | Fase 4 | WP1, WP5 | IMPLEMENTADO y validado dirigido: CAAS calibrado y RER por bloques conectados; RER genera lotes con semillas reproducibles. Equivalencia local con memoria completa y oráculos AUC aprobada; aceptación integrada y recursos reales a escala pendientes |
| WP7 | `WP7_caas_perms_por_bloques.md` | Fase 5 | WP5 | IMPLEMENTADO y validado dirigido: bloques + roster + índice canónicos; RDS único opcional (false por defecto). Seis reconstrucciones y Nextflow con/sin exportación aprobados; SCORING y 399 resultados FCS conservados. Aceptación integrada y recursos reales a escala pendientes |
| WP8 | `WP8_cierre.md` | | WP5–WP7 | IMPLEMENTADO y validado dirigido: limpieza inerte retirada, remove_work mínimo real, GUI/plantillas y DOMINO sin aristas; aceptación integrada y escala pendientes |
| WP9 | `WP9_estadistico_gen.md` | | A1 ya; A2 tras WP5 | A1 hecho, A2 pendiente |
| WP10 | `WP10_posenrich_auc.md` | | WP8, evaluación recuentos/AUC | PREPARADO a petición del usuario; sin implementación ni activación |

Actualiza esta tabla (estado y, en el fichero de la fase, el registro) al cerrar cada fase.
Paralelismo posible: WP1 y WP2 son independientes (ficheros disjuntos); WP9-A1 ya está hecho; WP4 y WP5 son seriales porque comparten el formato de tabla de la tesela.

## Validación y commits por WP (decisión del usuario, 2026-10-22)

Cada WP se guarda en un commit temático con sus pruebas dirigidas, referencias guardadas y pendientes explícitos. No se exige una ejecución completa del pipeline por cada WP. La regeneración del nulo, la corrida completa y la aceptación a escala se planifican conjuntamente al terminar la integración de los WP, usando las referencias congeladas para contrastar los resultados. Las mediciones específicas de recursos pueden reutilizar esa corrida integrada o ejecutarse como benchmarks focalizados.

Un commit de implementación o de avance parcial no implica haber superado la aceptación integrada ni haber medido 10.000 ciclos. Los estados y registros deben distinguir esas situaciones. Esta decisión sustituye los requisitos históricos de volver a ejecutar una corrida completa como condición de cada commit de WP. El usuario autoriza los commits temáticos y el push de esta revisión; no se reescribe el historial.

## Reglas de trabajo comunes

- **Nada cambia un número sin comparación con lo actual.** Referencias congeladas en `/home/miguel/IBE-UPF/PhD/PhyloPhere_validation/golden_581b999/` (local) y, para el cluster, en las salidas de `~/scratch/2.Primates/2.Primates_results/CAAS_RESULTS/{neoplasia,malignant}_def/*_prevalence_complete` con la huella de `~/scratch/0.Phylophere/.tmp/golden/*_manifest.sha256`. Código de referencia del FCS: `PhyloPhere_validation/unification/legacy_fcs/` (con `SHA256SUMS`).
- **Commits:** el usuario autoriza commits temáticos y push de WP8; no incluir cambios ajenos. La prohibición histórica queda sustituida por esa autorización. El commit de partida de la línea base es `581b999` (rama `remove_caap`); el checkout del cluster está en `d246e82`, y entre ambos la cadena del nulo es idéntica (RESAMPLE solo añade `pair_center_gate`, sin efecto con `none`; la acumulación pasa a la prueba exacta).
- **Cluster:** nada de cómputo en el nodo de entrada (`sbatch`, `salloc` o `srun --pty`); solo lecturas baratas. Sin lecturas pesadas (`wc -l` de una tabla de 1 GB colgó el nodo). Nunca editar código en el cluster: el usuario lo sincroniza por git. Directorio temporal: `~/scratch/0.Phylophere/.tmp`. Tras una caída de la VPN, `ssh correfoc` se cuelga con un socket obsoleto: usar `ssh -o ControlPath=none -o ControlMaster=no correfoc`.
- **Local:** nunca `/tmp`. Temporales grandes en la partición IBE (`/home/miguel/IBE-UPF/PhD/PhyloPhere_validation/.tmp`, también como `TMPDIR` de R y pytest); `/home` está casi lleno. Entorno con pandas, R y RERconverge: `~/micromamba/envs/phylophere/bin/{python,Rscript}`.
- **Código:** los comentarios describen solo el comportamiento actual (sin historia, fechas ni "antes era"). Todo parámetro nuevo lleva campo en la GUI y se sincroniza en `gui/templates/*.json` y `validation/tier*/templates`. El modo toy no fija parámetros. Cada desviación de RERconverge se justifica (referencia: commit 2bd328f y el walkthrough de Saputra).
- **Medidas:** distinguir medido, leído en el código y supuesto. No citar extrapolaciones como medidas. No usar la cifra de 61 GB del merge (anterior a su corrección).
- **Memoria:** al cerrar la fase, actualizar `project_fcs_10k_implementation_plan_2026-10-21.md` de la memoria del proyecto con el estado.

## Cómo arrancar una fase

1. Leer `docs/fases/README.md` (este fichero), el fichero de la fase y `docs/fases/PLAN_APROBADO_2026-10-21.md`.
2. Comprobar `git log --oneline -5` y `git status`: la fase parte del commit anotado en el registro de la fase anterior (preservar los cambios ajenos y aislar los del WP; no exigir un árbol limpio).
3. Confirmar que las dependencias están cerradas en la tabla de estado.
4. Si hay duda sobre el estado del cluster o de la VPN, comprobarlo antes de planificar.

## Cómo cerrar una fase

1. Ejecutar las comprobaciones de la sección "Validación" del fichero de la fase y anotar comando, conjunto de datos, tolerancia y resultado.
2. Rellenar la sección "Registro de la fase" (plantilla abajo).
3. Marcar el estado en la tabla de arriba y en la memoria del proyecto.
4. Dejar lista la entrada de la fase siguiente: interfaz producida y qué debe comprobar la próxima conversación primero.
5. Con autorización vigente, hacer commits temáticos y push; anotar sus hashes en el registro.

## Plantilla del registro de fase

```
- Conversación y fechas:
- Commit de partida / commit de cierre:
- Qué se hizo (ficheros y cambios):
- Comprobaciones ejecutadas (comando, conjunto de datos, tolerancia, resultado):
- Mediciones (recurso antes y después):
- Desviaciones del plan y decisiones tomadas (con su justificación):
- Pendiente que pasa a otra fase:
```

## Resumen de lo hecho antes de la primera fase (2026-10-21)

- Plan aprobado y copiado aquí; commit `581b999` hecho por el usuario.
- `PhyloPhere_validation/unification/legacy_fcs/` con `fcs_enrich.R`, `fcs_compute.R` y `percentile_flags.R` idénticos a los del commit (verificado por SHA-256).
- Trabajo 7568379 en el cluster: huellas SHA-256 de las salidas de neoplasia y maligno completo en `~/scratch/0.Phylophere/.tmp/golden/` y revisión del checkout (`d246e82`).
- WP9-A1 ejecutado (resultados en `WP9_estadistico_gen.md`; ficheros en `/home/miguel/IBE-UPF/PhD/PhyloPhere_validation/golden_581b999/a1_gene_statistic/`).
- `docs/POSENRICH_PORTE_EVALUACION.md` redactado (entregable de WP2).
- Medidas de recursos a 1000 a 10.000 ciclos (FCS, `scoring_caas_perms.R`, RESAMPLE, parseo de etiquetas): `docs/PLAN_NULO_10K_PERMULACIONES.md`, secciones 3, 5 y 8.

## Hallazgos pendientes de asignar a una fase

- **DOMINO sin aristas: resuelto en WP8.** Salidas vacías, slicer omitido, CLI/proceso de módulos toleran vacío y el informe 13 produce un grafo vacío. Prueba Nextflow mínima con dos procesos completados; ver `WP8_cierre.md`.

- **Referencia del marine toy (2026-10-21/22).** La corrida murió con el reinicio del equipo (21:07); se completó el 2026-10-22 al montar `/media/miguel/external` (`golden_581b999/marine_toy`, 17 min). Validación de WP3 y WP4 en marine: registro de WP4. Referencia local con el código de WP4: `PhyloPhere_validation/golden_marine_wp4/`.
