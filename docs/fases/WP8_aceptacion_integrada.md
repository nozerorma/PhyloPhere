# Aceptación integrada WP5–WP8: preparada, sin ejecutar a escala

No enviar estos trabajos desde la conversación ni ejecutar cómputo en el login. El usuario sincroniza el checkout y los ejecuta dentro de asignaciones Slurm. No cerrar aceptación por tener scripts renderizados o pruebas toy aprobadas.

## Launchers existentes revisados

- `validation/tier1/scripts/render_tier1_scripts.py` sirve también para las plantillas de neoplasia/maligno. Añade presupuestos explícitos `--caas-full-perms` y `--rer-perm-batches`, `--asr-cache` para el oráculo congelado, `--keep-work` y `--collect-resources`. Guarda `resolved_project.json`. No cambia el tamaño ni la semilla de los lotes RER, no cambia defaults científicos ni activa POSENRICH.
- `tests/cluster/prepare_wp6_snapshot.sh`: FCS CAAS o RER sobre un nulo guardado; `prepare_wp7_snapshot.sh`: reducción, SCORING, FCS y dos informes sobre el marino guardado. Ambos usan `git archive HEAD`, rechazan un snapshot ya existente y guardan commit, traza, informe y timeline. No son launchers de replay completo ni sustituyen las corridas de producción.
- `adhoc_scripts/{fcs_scaling,fcs_calibrated_ranking_control,fcs_rawmax_ranking_control,fcs_null_as_observed,scoring_perms_scaling}.sbatch` ya aceptan los bloques WP7. El ensamblado global de sus controles se mide aparte del camino productivo. No se crean otros launchers equivalentes.

## Preparación de una corrida de producción

Ejemplo en el cluster después de sincronizar la rama. Primero 1000 ciclos; repetir neoplasia y maligno, con strict activado en sus plantillas. Usar un directorio nuevo por ejecución. ASR_FROZEN debe señalar la caché correspondiente a esa referencia: codeml no es reproducible entre optimizaciones nuevas.

```bash
CHECKOUT="$HOME/scratch/0.Phylophere"
ACCEPT_ROOT="$HOME/scratch/0.Phylophere/.tmp/wp8_accept"
TRAIT_TEMPLATE=neoplasia_def
ACCEPT_N=1000
ASR_FROZEN="/ruta/existente/a/la/cache/ASR/congelada"
ACCEPT_RUN="$ACCEPT_ROOT/${TRAIT_TEMPLATE}_${ACCEPT_N}"
[ ! -e "$ACCEPT_RUN" ] || exit 2
mkdir -p "$ACCEPT_RUN/code" "$ACCEPT_RUN/tmp"
git -C "$CHECKOUT" archive HEAD | tar -x -C "$ACCEPT_RUN/code"
git -C "$CHECKOUT" rev-parse HEAD > "$ACCEPT_RUN/commit.txt"
export TMPDIR="$ACCEPT_RUN/tmp"
python "$ACCEPT_RUN/code/validation/tier1/scripts/render_tier1_scripts.py" \
  --template "$ACCEPT_RUN/code/gui/templates/${TRAIT_TEMPLATE}.json" \
  --outdir "$ACCEPT_RUN/code" --asr-cache "$ASR_FROZEN" \
  --caas-full-perms "$ACCEPT_N" --rer-perm-batches 10 \
  --keep-work --collect-resources
bash -n "$ACCEPT_RUN/code/SBATCH_run_${TRAIT_TEMPLATE}_complete.sh"
bash -n "$ACCEPT_RUN/code/run_${TRAIT_TEMPLATE}_single_complete.sh"
```

El usuario revisa `resolved_project.json` (rutas, catálogo de fenotipos, asignación Slurm, módulos activos y presupuestos) y lanza el SBATCH generado. Los scripts deben permanecer en la raíz del snapshot: su directorio determina `REPO_DIR`. Guardar ID de trabajo, `Rscript -e 'sessionInfo()'`, versiones de Python/Nextflow y los manifests de entradas/caché. No lanzar dos pruebas desde el mismo snapshot/work.

Para 4000 y 10000: directorios nuevos, `--caas-full-perms 4000/10000` y `--rer-perm-batches 40/100`, conservando `rer_perms_per_batch=100` y la semilla guardada. Comprobar que la cosecha RESAMPLE contiene N ciclos reales y suficientes hipótesis; no reciclar columnas como aceptación real. El nulo RER histórico sin semillas por lote no es un oráculo de identidad de una regeneración nueva.

## Referencias y comparación

Referencias locales: `PhyloPhere_validation/golden_581b999/`, `golden_marine_wp4/`, `wp5_check/`; evidencia dirigida nueva de WP6 en `.tmp/wp6_integration/` y WP7 en `.tmp/wp7/`. Referencias cluster congeladas: `~/scratch/2.Primates/2.Primates_results/CAAS_RESULTS/{neoplasia,malignant}_def/*_prevalence_complete`, con manifests en `~/scratch/0.Phylophere/.tmp/golden/`. Verificar las huellas antes de usarlas. Nunca eliminar ni sobrescribir referencias o work existentes.

Comparador existente (en la copia sincronizada de PhyloPhere_validation):

```bash
python "$VALIDATION/unification/compare_null.py" --a "$REFERENCE" --b "$NEW_RUN" --tol 1e-12 > "$ACCEPT_RUN/compare_complete.log"
python "$VALIDATION/unification/compare_null.py" --a "$REFERENCE/strict" --b "$NEW_RUN/strict" --tol 1e-12 > "$ACCEPT_RUN/compare_strict.log"
```

Es diagnóstico, no un criterio automático global: registra también diferencias de formato y estadísticas que cambiaron deliberadamente en WP1/WP6. Criterios específicos:

| Contrato | Criterio |
|---|---|
| Matrices CAAS completo/strict frente al mismo replay y filtros | Mismos genes, ciclos, orden, stat y parts; reconstrucción explícita frente al RDS histórico a 1e-12. No exigir que exista RDS único con exportación false. |
| Roster, removed_units y contadores | Igualdad exacta de claves, ciclos, contadores enteros y universo de posiciones; validar índice con `scoring_caas_perms.R --validate-only`. |
| SCORING y acumulación | Scores/p a 1e-12, claves y universo idénticos. La acumulación usa posiciones positivas agregadas, no requiere cada fila posición × ciclo. |
| FCS calibrado | Oráculo de memoria completa sobre los mismos ciclos y ranking gene_p; AUC centrado, enteros de excedencia y N válidas exactos, p/q a 1e-12, NA y reglas iguales. Los antiguos p analíticos permulados o ranking crudo no son su oráculo. |
| RER | Mismos ajustes/semillas/partición para reproducibilidad nueva; validar observado contra correlación nativa y FCS sobre los mismos bloques. No restringir el nulo a genes con observado disponible. |
| Informes y opcionales | Informes 11/12 sin RDS único; ausencia explícita de 14; DOMINO sin aristas termina con vacíos; exportación true se prueba aparte. |

PEPC guardado (un gen × 1000) sirve para matrices; no acredita rutas positivas con min_g=10. El marino de 121 genes y cáncer de 45 son dirigidos, no escala real. Comparar 1000/4000/10000 como ejecuciones reales independientes con las mismas entradas y semillas, no contra columnas sintéticas repetidas.

## Recursos y limpieza

Conservar `pipeline_info/{trace.tsv,execution.html,timeline.html}`, `.nextflow.log`, params y scripts resueltos; `execution_trace.txt` es el nombre por defecto cuando no se pide la recogida explícita. Guardar `sacct -j JOB --format=JobID,State,ExitCode,Elapsed,AllocCPUS,MaxRSS,MaxDiskRead,MaxDiskWrite` y `du -s -B1` por resultados, work, run y NXF_HOME. Informar picos por familia, intentos y bloque; separar bosque RER, universo de posiciones, GENE_P y exportación global. La memoria constante del pipeline no está demostrada.

Solo para una corrida propia destinada a probar limpieza: guardar SHA-256 de todos los ficheros públicos antes/después y tamaños asignados de las tres carpetas propias; activar `remove_work` en su launcher. Nunca limpiar referencias o work del usuario. WP8 ya midió un caso mínimo real con 1 MiB de carga y 1.093.632 bytes asignados recuperados; esto no valida el ahorro de una corrida completa.
