# WP8. Limpieza, documentación y sincronía

Estado: IMPLEMENTADO Y VALIDADO DE FORMA DIRIGIDA (2026-10-23). Aceptación integrada y recursos reales a escala pendientes. Partida `remove_caap`, `ced33666`. WP5/WP6/WP7 ya implementados; no se repite su trabajo ni se declaran aceptados a escala.

## Limpieza efectiva auditada

| Capa | Comportamiento efectivo |
|---|---|
| Higiene de tareas | Incondicional: SUBSET retira temporales de selección; CORE retira `caas_perms_out` tras extraer sus tablas; MERGE convierte `matrices/` a `matrices.rds` y borra solo ese traspaso crudo. FADE/SELECTION retiran temporales de alineación; copias de fuentes retiran pyc/cache. No dependían ya de remove_intermediates. |
| Publicación | Nextflow copia outputs declarados: bloques completo/strict y b0, roster/índice, contadores y tablas reducidas. El RDS global solo con caas_export_single=true. Las teselas null_tables permanecen en work y no se publican como un almacén global por posición. |
| `remove_intermediates` | Helper sin llamadas tras WP5; retirado de config, GUI, modelo, scripts, ayuda y plantillas. Migración elimina la clave con aviso. No se redefine como poda ni se pierde reutilización. |
| `remove_work` | Solo launcher generado, después de Nextflow con código 0, incluido skip exitoso por pocos contrastes. Elimina WORK_DIR, RUN_DIR y NXF_HOME privado; copia el log a resultados. set -e conserva directorios tras fallo. Ahora valida todos los destinos antes de borrar: bajo WORK_BASE y sin contener resultados/código. No elimina el destino de enlaces de plugins. Nextflow directo no ejecuta esta limpieza. |

Ninguna limpieza elimina bloques/roster/índice públicos. `-resume` requiere conservar work y su caché de ejecución. La exportación global es opcional; su ausencia no bloquea SCORING/FCS/informes.

## GUI, plantillas y pruebas obsoletas

- Retirados los selectores de procesos de permulación ya inexistentes de las dos plantillas LQ; migración elimina esos overrides antiguos sin alterar los actuales.
- Ocho referencias externas a `cancer_multi_final.json` sustituidas por `neoplasia_def.json`. No se recrea la plantilla retirada. La suite externa de limpieza se actualiza al contrato vigente. Los cambios externos están en PhyloPhere_validation/wiring y no viajan en el commit de este repositorio; `tests/test_wp8.py` cubre los contratos dentro del repositorio.
- Tier1 PEPC queda idéntico a su pareja GUI: remove_work, tamaño de bloque, gene_stat, límites 10/0 y parámetros retirados armonizados. Marino completo armoniza esos campos. En toy marino se conservan rutas, módulos activados, filtros, caché y presupuestos particulares; solo se actualizan campos de contrato/defaults ya decididos. La plantilla toy local de validation está ignorada por git y se conserva como tal.
- Traducción de la casilla remove_work en los seis idiomas; tooltip describe también run/NXF_HOME. Nuevos campos WP6/WP7 ya estaban conectados y se mantienen.

## Puertas

- Redundancia: análisis con el código productivo en `.tmp/wp8/redundancy.R`, salida TSV. Solo maligno completo histórico conserva leading edges locales. BH recalculado por ranking/base: tres conjuntos globales; Jaccard 0,25/0,5/0,75 → 3/3/3 grupos; overlap → 2/2/3. p, q y sig permanecen idénticos por aserción. Evidencia limitada, anterior al ranking calibrado actual; se conserva Jaccard 0,5 y agrupación exclusivamente de presentación. No se elige por maximizar significativos.
- POSENRICH: efectivamente omitido desde WP5. Auditoría confirma la cola charset sin difusión en el módulo histórico: solo alimenta el primer lote y explica la pérdida de fuentes registrada como 5/12. No se ha repetido aquella corrida ni se afirma cobertura corregida. El usuario eligió preparar AUC por posición en fase propia; después abrió un fork para evaluar también recuentos como acumulación. `WP10_posenrich_auc.md` registra esa puerta sin activar estadística nueva.
- DOMINO: verificado el RuntimeError por cero aristas; sustituido por vacíos con mensaje de skip. No se llama slicer con red vacía; proceso y CLI de módulos omiten limpiamente. El informe 13 devuelve grafo vacío antes de leer SIF/scores vacíos y sus listas se omiten por falta de genes mapeados.

## Registro y evidencia

- Evidencia: `/home/miguel/IBE-UPF/PhD/PhyloPhere_validation/.tmp/wp8/`. Temporales en IBE, entorno phylophere. Sin trabajos ni cambios en el cluster; sin corrida completa ni 10.000 ciclos.
- `tests/test_wp8.py` + `test_gui_params_json.py`, `test_gui_core_batch.py`, `test_cleanup_options.py`: **76 PASS** (`pytest_final.log`). Carga/guardado/recarga de todas las plantillas GUI y tier, render y bash -n; migración; guardas; launcher de aceptación. Prueba DOMINO ampliada para incluir CLI directo: **1 PASS** dirigido (`domino_final.log`), con 10 casos no seleccionados.
- Prueba remove_work: tres corridas secuenciales reales de un proceso Nextflow desde fuentes propias estables: éxito+borrado, fallo+conservación, éxito+flag false. Directorios propios; carga 1.048.576 bytes, disco asignado recuperado **1.093.632 bytes** en el caso de éxito. Huellas SHA-256 públicas preservadas (incluida salida creada por el proceso); log copiado. Medición mínima, no corrida grande. También se rechaza borrar un work que contiene resultados.
- DOMINO Nextflow desde snapshot estable `domino/code`: dos procesos COMPLETED, cero módulos, `EMPTY_DOMINO_PASS`, traza en `domino/trace_retry.tsv`; STRING local mínimo con score 600 y corte 700, sin descargas. Primer intento del harness pasó dos entradas del mismo basename y falló al staging; corregido el fixture, sin cambiar el pipeline para satisfacerlo.
- GUI inicial + suite general de cableado: **69 PASS / 7 FAIL**, `gui_initial.log`; los siete nombres coinciden exactamente con WP5. Causas comprobadas en aserciones: ruta antigua de acumulación, informe ASR no activado, discovery con configuración incompatible, tres casos de exports retirados y parser estático que no reconoce publishDir por patrón/prefijo. Los 59 casos de las dos suites GUI pasan, frente a 49/10 en WP7. No se declara verde cableado general.
- La ayuda compila desde snapshot estable; salida 1 deliberada de HELP (no fallo de compilación). La pasada general inicial no activa todos los módulos; se comprueba también ayuda de enriquecimiento.
- Se preservan las medidas y pruebas aprobadas WP5–WP7; no se repiten los 49 casos dirigidos ni las seis exportaciones sin cambios que lo justifiquen. Comentarios del core actualizados, sin cambios en cálculo ni formato.
- Launchers existentes adaptados para presupuestos explícitos, ASR congelado, recursos y snapshots no reutilizados. Comandos, referencias, tolerancias y recogida en `WP8_aceptacion_integrada.md`. Renderizado probado; preparación no equivale a ejecución.

## Entrega y pendientes

Implementación: `1416f3ec`. Documentación en el commit temático siguiente. Commits temáticos y push autorizados por el usuario. Preservar modificaciones previas en validation/, WP4, su fila de tabla y prompts; incluir solo deltas WP8 en archivos compartidos. La memoria del proyecto se actualiza tras revisar pruebas.

Pendientes: aceptación integrada WP5–WP8 a 1000/4000/10000 y recursos reales; decisiones WP9; revisión de la comparación recuentos/AUC del fork antes de implementar POSENRICH; adaptación futura de los siete tests generales de cableado. WP8 cierra implementación y documentación dirigida, no estas aceptaciones.
