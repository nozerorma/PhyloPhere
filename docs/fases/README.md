# Ejecución por fases: arranque, cierre y registro

Cada paquete de trabajo (WP) se ejecuta en una conversación independiente. Este directorio es lo que debe bastar para arrancar una conversación nueva sin contexto previo y para dejar constancia de lo hecho.
Plan aprobado (copia versionada): `docs/fases/PLAN_APROBADO_2026-10-21.md`. Diagnóstico y mediciones: `docs/PLAN_NULO_10K_PERMULACIONES.md`. Revisión del FCS: `docs/FCS_REVISION_RERCONVERGE.md`. Porte de POSENRICH: `docs/POSENRICH_PORTE_EVALUACION.md`.

## Estado

| WP | Fichero | Fase del plan de 10.000 ciclos | Depende de | Estado |
|---|---|---|---|---|
| WP0 | `WP0_linea_base.md` | | commit del usuario (hecho: `581b999`) | HECHO con un pendiente explícito: marine toy (falta el disco externo con los alineamientos toy). PEPC y cancer_b0_toy regenerados con `a92fd2d` |
| WP1 | `WP1_fcs_rer.md` | | WP0 | HECHO (p.perm de CAAS idéntico a la referencia; RER global firmado; pendiente medir RER_CONT a 10.000 y los recuentos 1, 0, 6, 6) |
| WP2 | `WP2_exports_posenrich_limpieza.md` | | WP0 | HECHO (PEPC c4 y toy de cancer idénticos a la referencia con su caché de ASR; marine toy pendiente) |
| WP3 | `WP3_bloques_de_ciclos.md` | Fase 1 | WP2 | pendiente |
| WP4 | `WP4_replay_en_flujo.md` | Fase 2 | WP3 | pendiente |
| WP5 | `WP5_tablas_y_merge_por_bloque.md` | Fase 3 | WP4 | pendiente |
| WP6 | `WP6_fcs_por_bloques.md` | Fase 4 | WP1, WP5 | pendiente |
| WP7 | `WP7_caas_perms_por_bloques.md` | Fase 5 | WP5 | pendiente |
| WP8 | `WP8_cierre.md` | | WP5 | pendiente |
| WP9 | `WP9_estadistico_gen.md` | | A1 ya; A2 tras WP5 | A1 hecho, A2 pendiente |

Actualiza esta tabla (estado y, en el fichero de la fase, el registro) al cerrar cada fase.
Paralelismo posible: WP1 y WP2 son independientes (ficheros disjuntos); WP9-A1 ya está hecho; WP4 y WP5 son seriales porque comparten el formato de tabla de la tesela.

## Reglas de trabajo comunes

- **Nada cambia un número sin comparación con lo actual.** Referencias congeladas en `/home/miguel/IBE-UPF/PhD/PhyloPhere_validation/golden_581b999/` (local) y, para el cluster, en las salidas de `~/scratch/2.Primates/2.Primates_results/CAAS_RESULTS/{neoplasia,malignant}_def/*_prevalence_complete` con la huella de `~/scratch/0.Phylophere/.tmp/golden/*_manifest.sha256`. Código de referencia del FCS: `PhyloPhere_validation/unification/legacy_fcs/` (con `SHA256SUMS`).
- **Commits:** los hace el usuario. No hacer commit ni push. El commit de partida de la línea base es `581b999` (rama `remove_caap`); el checkout del cluster está en `d246e82`, y entre ambos la cadena del nulo es idéntica (RESAMPLE solo añade `pair_center_gate`, sin efecto con `none`; la acumulación pasa a la prueba exacta).
- **Cluster:** nada de cómputo en el nodo de entrada (`sbatch`, `salloc` o `srun --pty`); solo lecturas baratas. Sin lecturas pesadas (`wc -l` de una tabla de 1 GB colgó el nodo). Nunca editar código en el cluster: el usuario lo sincroniza por git. Directorio temporal: `~/scratch/0.Phylophere/.tmp`. Tras una caída de la VPN, `ssh correfoc` se cuelga con un socket obsoleto: usar `ssh -o ControlPath=none -o ControlMaster=no correfoc`.
- **Local:** nunca `/tmp`. Temporales grandes en la partición IBE (`/home/miguel/IBE-UPF/PhD/PhyloPhere_validation/.tmp`, también como `TMPDIR` de R y pytest); `/home` está casi lleno. Entorno con pandas, R y RERconverge: `~/micromamba/envs/phylophere/bin/{python,Rscript}`.
- **Código:** los comentarios describen solo el comportamiento actual (sin historia, fechas ni "antes era"). Todo parámetro nuevo lleva campo en la GUI y se sincroniza en `gui/templates/*.json` y `validation/tier*/templates`. El modo toy no fija parámetros. Cada desviación de RERconverge se justifica (referencia: commit 2bd328f y el walkthrough de Saputra).
- **Medidas:** distinguir medido, leído en el código y supuesto. No citar extrapolaciones como medidas. No usar la cifra de 61 GB del merge (anterior a su corrección).
- **Memoria:** al cerrar la fase, actualizar `project_fcs_10k_implementation_plan_2026-10-21.md` de la memoria del proyecto con el estado.

## Cómo arrancar una fase

1. Leer `docs/fases/README.md` (este fichero), el fichero de la fase y `docs/fases/PLAN_APROBADO_2026-10-21.md`.
2. Comprobar `git log --oneline -5` y `git status`: la fase parte del commit anotado en el registro de la fase anterior (pedir al usuario que haga commit si el árbol tiene cambios de la fase previa).
3. Confirmar que las dependencias están cerradas en la tabla de estado.
4. Si hay duda sobre el estado del cluster o de la VPN, comprobarlo antes de planificar.

## Cómo cerrar una fase

1. Ejecutar las comprobaciones de la sección "Validación" del fichero de la fase y anotar comando, conjunto de datos, tolerancia y resultado.
2. Rellenar la sección "Registro de la fase" (plantilla abajo).
3. Marcar el estado en la tabla de arriba y en la memoria del proyecto.
4. Dejar lista la entrada de la fase siguiente: interfaz producida y qué debe comprobar la próxima conversación primero.
5. Pedir al usuario el commit de cierre y anotar su hash en el registro.

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

- **DOMINO debe saltarse el paso (skip silencioso) si no hay aristas suficientes** (decisión del usuario, 2026-10-21). Hoy `DOMINO_BUILD_NETWORK` (`subworkflows/ENRICHMENT/local/src/build_domino_network.py:223`) lanza `RuntimeError("no edges survived filtering ...")` y, con la estrategia de error por defecto, aborta toda la corrida. Visto en el toy de cancer (45 genes, 31 de fondo, 0 aristas con puntuación >= 700 sobre 13.715.404 enlaces de STRING). Comportamiento deseado: sin aristas, escribir las salidas vacías (como ya hace la rama de fondo ausente, que crea `network.sif`, `network_edge_scores.tsv` y `slices.txt`), avisar en el log y dejar que los pasos siguientes tomen el centinela `NO_*`. Comprobar además que `run_domino_modules.py` y `13.AMI_analysis.Rmd` toleran redes vacías. Sin código cambiado todavía; candidato a WP2 (limpieza) o WP8.
