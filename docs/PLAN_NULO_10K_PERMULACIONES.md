# Escalado del nulo de permulación a 10.000 ciclos y estado del FCS de CAAS

## Contrato vigente tras WP8 (2026-10-23)

Esta sección sustituye las propuestas y descripciones históricas de implementación de este documento. WP5, WP6 y WP7 están implementados y validados de forma dirigida; no quedan por construir sus bloques ni su FCS productivo. La aceptación integrada y las medidas reales a escala siguen pendientes.

- Nulo CAAS canónico: `caas_permulation/blocks/blk_*/matrices.rds`, `cycle_roster.txt`, `null_metadata.json` y contadores/tablas WP5. Cada bloque conserva genes, ciclos, unidades, `gene_stat`, `stat` y `parts` (max/fmax/n). Completo y strict usan el mismo contrato. No hay tabla pública de cada posición × ciclo.
- `caas_export_single=false`: el camino normal no ensambla `caas_perms.rds`. La exportación explícita reconstruye el objeto histórico con coste global deliberado. Se conserva la importación del RDS histórico; el almacén canónico completo permite reutilizar el nulo sin él.
- Ranking CAAS productivo: −log10(p de permulación del gen), con observado y leave-one-out en `gene_p.R`. `gene_stat` sigue parametrizado; WP9 decide entre máximo, F(max) y F(max)^n. El ranking crudo es un control, no el fallback productivo.
- Rutas: `stat = AUC − 0,5`, cola superior CAAS y bilateral RER; empates inclusivos, `(k+1)/(N_válidas+1)`. `sig_perm` es BH sobre `p.perm`, por ranking/base. `redlich_level` es independiente y `sig_top` está retirado. El observado interno usa RDS para conservar empates y recuentos; el TSV es público.
- RER continuo/binario genera bloques con semillas por lote. La disponibilidad nula se evalúa por ciclo y no se restringe a genes disponibles en el observado. Cambiar tamaño de lote cambia las simulaciones; conservar ajustes, semilla y partición al comparar.
- No se afirma memoria constante del pipeline: GENE_P retiene todos los ciclos de un lote de genes; el bosque RER y el universo de posiciones tienen costes separados. Las medidas sintéticas de WP6 no sustituyen recursos reales.
- `remove_intermediates` retirado: era inerte tras WP5. La higiene de tareas es incondicional; ningún bloque, roster o índice público se poda. Los proyectos antiguos conservan solo su elección de `remove_work`. Este último pertenece al launcher: tras éxito elimina work, directorio de ejecución y NXF_HOME privado, con guardas de rutas; no actúa por invocar Nextflow directamente.
- POSENRICH se omite con aviso en el flujo productivo tras WP5. Decisión del usuario en WP8: preparar una fase propia de AUC por posición (opción A), sin activarla aquí; `docs/fases/WP10_posenrich_auc.md`.

Controles: máximo crudo frente a ranking calibrado, nulo como observado, tres estadísticos base y p/BH por gen se distinguen de las pruebas por posición (`p.emp`, `p.emp_fact`) y del AUC por conjunto. Los cuatro `adhoc_scripts/fcs*.sbatch` son oráculos de memoria completa; `tests/test_wp6_blocks.R`, `tests/test_fcs_auc_contract.R` y `tests/test_wp7.py` cubren contratos dirigidos. Los seis RDS PEPC/marine/cáncer completo/strict y los 399 resultados marinos de WP7 están registrados en su fase, no se han vuelto a ejecutar en WP8.

La sensibilidad de agrupación de WP8 usa una tabla histórica con leading edges, no recalcula inferencia productiva ni calibra potencia: Jaccard 0,25/0,5/0,75 conserva los tres grupos BH de maligno completo; overlap da 2/2/3. Se conserva Jaccard 0,5 y la agrupación solo de presentación. Criterios y comandos para aceptación: `docs/fases/WP8_aceptacion_integrada.md`.


Estado: diagnóstico y medidas históricas; implementación vigente y pendientes en el contrato anterior.
Complementa a `docs/FCS_REVISION_RERCONVERGE.md` (puntos P1 a P15 de la revisión del FCS frente a RERconverge).

Convención: **medido** es una cifra leída de una traza de Nextflow, de `sacct` o de una tabla de resultados; **leído en el código** es una
afirmación sobre estructura sin cronometrar; **supuesto** es una extrapolación o una hipótesis sin comprobar.

## 1. Objetivo y decisiones

Objetivo: 10.000 permulaciones por corrida, en CAAS y en RER, con recursos contenidos (RAM, E/S, disco, CPU).

| ID | Decisión | Estado |
|---|---|---|
| D1 | Lachenbruch se elimina del FCS y de POSENRICH. POSENRICH queda fuera de alcance por ahora | Decidido |
| D2 | Criterio de significación del FCS de CAAS: Wilcoxon con BH sobre el p de permulación, por ranking y base de datos (q < 0,15) | Implementado: BH sobre p.perm; Redlich independiente |
| D3 | La reutilización de exports de perm-discovery (`permDiscFiles`) se elimina | Decidido |
| D4 | La acumulación usa el marco hipergeométrico; del nulo solo necesita el conjunto de posiciones con score positivo | Decidido |
| D5 | El coste del FCS a 10.000 ciclos es asumible si también lo es el de RER | Decidido |
| D6 | `remove_intermediates` retirado en WP8; higiene incondicional y `remove_work` solo en launchers tras éxito | Implementado, prueba mínima real; escala pendiente |
| D7 | Tamaño de bloque del replay: fijo (por ciclos) o calculado por memoria | Implementado fijo: 500 por defecto; evaluación integrada pendiente |

## 2. Qué dice el FCS de CAAS (resultados de los controles)

Todos los controles usan "ciclo nulo como observado": un ciclo nulo hace de observado frente a los otros 999. Se comparan el observado y los ciclos
nulos con nivel (genes puntuados) de al menos 0,8 veces el del observado. Conjuntos de al menos 10 genes.

### 2.1 Regla histórica (BH analítico + puerta de p de permulación + dos votos)

Número de conjuntos "Supported" por ciclo nulo de nivel parecido, frente al observado:

| Corrida | Observado | Ciclo nulo (media) | Proporción |
|---|---|---|---|
| Neoplasia estricto | 86 | 35,4 | 41 % |
| Maligno estricto | 19 | 47,6 | 251 % |
| Neoplasia completo | 45 | 36,7 | 82 % |
| Maligno completo | 72 | 43,2 | 60 % |

Con el ranking calibrado por gen (sección 2.2) las proporciones son 34 %, 121 %, 46 % y 47 % (observado 101, 32, 73 y 83). La regla vigente no discrimina
el observado de un ciclo nulo del mismo nivel: sus "Supported" son información de priorización, no enriquecimiento.

### 2.2 Ranking calibrado por gen y regla propuesta (D2)

El ranking es −log10 del p de permulación del gen contra los otros ciclos (el análogo del FCS de RER), con el score de gen `F(max)^n`.
Conjuntos con BH(p de permulación) < 0,15 con el Wilcoxon, sumados sobre los tres rankings (global, top, bottom):

| Corrida | Observado | Ciclos nulos de nivel parecido (mediana, máx.) | Ciclos nulos con al menos el observado |
|---|---|---|---|
| Neoplasia estricto | 1 | 2, 4 | 80 % de 5 |
| Maligno estricto | 0 | 0, 4 | 100 % de 5 |
| Neoplasia completo | 9 | 3, 26 | 15 % de 13 |
| Maligno completo | 15 | 3, 22 | 10 % de 10 |

Conjuntos con p de permulación ≤ 0,001 (esperados por azar: ~30): estricto 10 y 0; completo 47 (nulo medio 25,8) y 48 (nulo medio 23,8; 0 de 10 ciclos nulos
con al menos el observado). En las corridas completas hay un exceso moderado; en las estrictas, ninguno. Con 10 a 13 ciclos nulos por corrida y una semilla, es
sugerente, no significativo.

Lachenbruch: el p de permulación depende del nivel en las corridas estrictas (Spearman 0,81 y 0,93 con el nivel del ciclo) y, con el ranking calibrado, el
observado queda por encima de lo nulo sin respaldo del Wilcoxon. Se elimina (D1).

### 2.3 Ranking con el máximo crudo de posición

Ranking construido igual, pero con el máximo crudo de posición del gen en lugar de `F(max)^n` (neoplasia estricto, 19 ciclos nulos y el observado):
BH(p de permulación) < 0,15 en 1 conjunto frente a una mediana nula de 3 (máx. 4; 80 % de ciclos nulos con al menos el observado); p ≤ 0,001 en 10 conjuntos
(nulo medio 11,2). Es idéntico al resultado con `F(max)^n` (1 frente a 2; 10 frente a 11,6). El cambio de score de gen no modifica el cuadro a nivel de
conjunto, aunque sí lo hace a nivel de gen (sección 2.5).

### 2.4 Score de gen `F(max)^n`

`F` es la ECDF de los scores de posición del mismo etiquetado y dirección; `F(max)` es la fracción de posiciones puntuadas con score ≤ el máximo del gen; `^n`
(n = posiciones puntuadas del gen) lo convierte en P(máximo de n posiciones del pool ≤ el observado), es decir, 1 menos un p de Šidák de la mejor posición.
Medido en 1000 ciclos nulos de neoplasia estricto (dirección `all`):

| Score | Media | Spearman con n | Spearman de la media con la propensión |
|---|---|---|---|
| `F(max)^n` | 0,501 (KS 0,017) | −0,03 | 0,19 |
| `F(max)` | 0,535 | 0,23 | 0,38 |
| Máximo crudo | n/a | 0,23 | 0,41 |

El 80 % de los genes puntuados tienen n = 1 (donde `^n` no actúa); con n > 8 (0,1 % de los genes) la media es 0,35 (conservador). Es una normalización dentro de
cada etiquetado: por construcción casi uniforme, haya muchas o pocas posiciones fuertes. Por eso elimina la magnitud global del observado (sección 2.5).
La fracción de varianza entre genes (ICC) es ~0,10 para los tres scores.

### 2.5 Magnitud del score de posición y p de gen (neoplasia estricto)

- **Pool de scores de posición** (observado frente a ciclos nulos con 0,8 a 1,5 veces las posiciones del observado): media 0,078 frente a 0,057, q90 0,175 frente a
  0,123, q99 0,351 frente a 0,266; entre el 17 y el 21 % de los ciclos nulos parecidos es igual de fuerte. Entre ciclos nulos, más detecciones se asocian a scores
  medios más bajos (Spearman −0,16).
- **`p.emp` de posición**: la recomputación reproduce exactamente la tabla publicada (3638 filas, diferencia 0). Está limitado por el suelo (1/(N + 1)): BH < 0,05 da
  0 posiciones en el observado y en todos los ciclos nulos. En el suelo hay 238 posiciones frente a una mediana nula de 72 (16 % de ciclos nulos con al menos el
  observado); con p ≤ 0,01, 1570 (20 %).
- **p de gen contra los otros ciclos** (2480 genes puntuados): con `F(max)^n`, p ≤ 0,01 en 216 genes frente a una media nula de 215 (38 % de ciclos nulos con al menos
  el observado); con el máximo crudo, 309 frente a 194 (19 %); con p ≤ 0,05, 894 frente a 857 y 1129 frente a 844 (15 %).
- **BH sobre el p de gen** (familia = universo de CAAS estricto, 12.473 genes, no puntuados con p = 1): 0 aciertos en el observado con ambos scores. Para q < 0,05 harían
  falta unos 240 genes en el suelo (≈ 80 para q < 0,15) con 1000 ciclos; el observado tiene 34 con el máximo crudo. Con 10.000 ciclos el suelo baja a 1e-4 y bastarían
  unos 24. Con el máximo crudo, 10 de 1000 ciclos nulos tienen algún acierto con q < 0,15 (máx. 429 en uno): un ciclo globalmente fuerte desplaza todos sus p a la vez.

### 2.6 Lectura y límites

El `caas_score` de posición lleva un exceso de magnitud moderado en estricto (15 a 20 % de ciclos nulos de nivel parecido igual de fuertes), y `F(max)^n` lo elimina por
construcción. Ni con el ranking calibrado ni con el máximo crudo aparece enriquecimiento de conjuntos en las corridas estrictas; en las completas hay un exceso moderado
con el ranking calibrado. No hay evidencia de que la construcción del score esté mal, ni de que el exceso sea señal. Límites: una semilla, entre 5 y 13 ciclos nulos por nivel,
conjuntos que comparten genes (conteos no independientes), emparejado por nivel basado solo en el número de genes o posiciones.

Sin ejecutar: el control posicional y el ranking de máximo crudo en las corridas completas y en maligno; más ciclos nulos (otra semilla) para las completas.

## 3. Recursos medidos

Trazas de Nextflow de las pasadas finales (`pipeline_info/execution_trace.txt`), 1000 ciclos, malignant completo:

| Proceso | Tareas | Tiempo (mediana, máx.) | RSS pico | Notas |
|---|---|---|---|---|
| `CAAS_CORE_BATCHED` | 83 | 13,7 y 21,4 min | mediana 16,6 GB, máx. 30,4 GB (petición 20 GB) | 2 fallos; CPU ~255 % de 400 %; 273 GB leídos y 579 GB escritos (suma de E/S) |
| `CAAS_CORE_MERGE` | 1 | 32,5 min | 4,1 GB | 1,5 GB de E/S |
| `CAAS_CORE_MERGE_STRICT` | 1 | 4,7 min | 2,9 GB | |
| `FCS_COMPUTE` | 4 + 4 | 2 a 4 min | 3 a 4 GB | |
| `SCORING_COMPUTE` | 1 | 2 min | 2,4 GB | |
| `RESAMPLE` | 1 | 0,5 min | 3,5 GB | |

Neoplasia completo: `CAAS_CORE_BATCHED` con RSS máximo de 46,8 GB y 6 fallos. Una agregación de `sacct` de varias corridas mezcla pasadas anteriores a la corrección del merge:
no sirve para dimensionarlo. CPU de `CAAS_CORE_BATCHED`: unas 50 CPU-horas por 1000 ciclos y corrida (83 tareas × 13,7 min × ~2,6 CPUs). Los controles del FCS
(72 tareas): 3,7 a 4,8 GB y 8 a 14 min por tarea (una tarea es un FCS completo con 1000 ciclos).

## 4. Qué hace RERconverge

`getPermsContinuous` (paquete instalado, deserializado en `~/.cc/rer_installed/dump/`): por cada permulación, un paso de álgebra de matrices sobre una matriz RER ya hecha
(correlación de todos los genes con el rasgo permulado). Persiste solo matrices genes × permulaciones (P, Rho y estadístico) y conjuntos × permulaciones (P y
estadístico de enriquecimiento). No guarda nada por rama ni por posición. Su coste por permulación no es comparable al de CAAS (replay de ASR sobre todas las posiciones),
pero la idea de conservar solo el estadístico suficiente de cada consumidor sí se traslada.

## 5. Diagnóstico del escalado con N (leído en el código, salvo lo indicado)

1. **Memoria por filas del gen.** `_parse_discovery_entries` (`gene_wrapper.py:181`) carga como objetos todas las filas de los ciclos del trozo; `score_labelings`
   (`driver.py:185`) devuelve registros de todas las etiquetas; `pool_labelings` agrupa por ciclo base después. El propio código modela ~4400 MB por millón de filas
   (`_REPLAY_MB_PER_MROW`), con RP1 en 14 GB por worker. El agrupado ya es por ciclo base, por lo que procesar de ciclo base en ciclo base debería dar los mismos números.
2. **Término fijo por tarea y por gen proporcional a N × hipótesis.**
   - Cada tarea del pool recibe el diccionario completo de etiquetados y el diseño FOP (`gene_wrapper.py:1151`).
   - `SUBSET_RESAMPLE_PERMS` (`caas_permulation.nf:72`) escribe un único `resample_perms.tab` con todas las hipótesis de los N ciclos (~100 por ciclo) y `fop_pairs.tsv`,
     que reciben todas las tareas de `CAAS_CORE_BATCHED`.
   - `run_ct_perm_replay_batch.sh` lanza un `ct perm-replay` por gen, y `ct` (línea 248) llama a `simtrait_revive` sobre el fichero entero cada vez.
   - **Medido** (40 especies, `resample_perms.tab` real de 11,7 MB): parseo de 0,6, 0,9 y 2,1 s y RSS de 89, 169 y 328 MB para 250, 500 y 1000 ciclos; las máscaras del kernel
     añaden 0,15 a 0,34 s y 28 MB. Son ~2,4 s de CPU por gen (~11 CPU-horas sobre ~16.100 genes, ~20 % de `CAAS_CORE_BATCHED`) a 1000 ciclos; a 10.000, ~3,3 GB por proceso y
     ~24 s por gen (extrapolación lineal), con 4 workers por tarea.
3. **Proceso padre.** `process_all_genes_perms` acumula todos los trozos de un gen (`pending_pooled`) y ejecuta `_perms_worker_finalize` en serie. Su parte en los picos de 30 a 47 GB
   no está medida.
4. **E/S y almacenamiento.**
   - El descubrimiento (`perm_replay.py`, Python) recorre trozos de etiquetados × posiciones, mantiene en listas todos los aciertos del gen y escribe el export en orden de
     posición, como texto (13 a 16 MB por gen a K = 3). El replay lo relee y lo reagrupa por ciclo.
   - El paso B (`_finalize_perm_scores`, `gene_wrapper.py:753`) lee el detalle dos veces (pools por ciclo y genes) y reescribe `perm_pos_cycle_caas.tsv.gz` con todas las filas
     (80,9 millones a 1000 ciclos en cáncer; el 28 % con score de US).
   - `scoring_caas_perms.R` lee entero `gene_cycle_scores.tsv` con `read_tsv` (1,1 GB, ~15 millones de filas a 1000 ciclos en neoplasia completo). **Medido** sobre esa tabla
     recortada (K <= 1000) o reciclada (K > 1000; los recursos valen, la salida no; trabajo 7567372):

     | Ciclos | RSS máx. | Tiempo |
     |---|---|---|
     | 250 | 1,0 GB | 0,7 min |
     | 500 | 1,8 GB | 1,1 min |
     | 1000 | 3,4 GB | 1,8 min |
     | 2000 | 6,8 GB | 2,8 min |
     | 4000 | 13,5 GB | 4,8 min |

     Lineal, ~3,4 MB por ciclo; extrapolado, ~34 GB y ~11 min a 10.000 ciclos. Con la petición actual del merge (24 GB × intento) no cabe al primer intento.
     De las seis matrices que construye, las tres del eje ASR (`global_asr`, `top_asr`, `bottom_asr`, el cuantil 90 de los scores de posición) no tienen consumidor en el repositorio:
     `fcs_enrich.R:660` las concatena con las de CAAS y no las selecciona nunca (los nombres no coinciden con los rankings), `scoring_compute.R:520` solo lee los nombres de columna de la
     matriz global de CAAS, y ningún informe lee el fichero de permulaciones (el comentario de `fcs_enrich.R:663`, que atribuye su uso al informe 11, no se corresponde con el código
     del informe). Falta descartar lectores externos al repositorio.
5. **FCS.** Las matrices de entrada son pequeñas (14.913 genes × N × 8 B = 120 MB por dirección y por 1000 ciclos), pero se construyen copias densas conjuntos × N
   (`fcs_null_enrichstat_vectorized`, `fcs_enrich.R:355`). **Medido** (FCS de neoplasia completo con el ranking calibrado; ciclos reales hasta 1000 y reciclados por encima, de modo que
   los recursos valen y los p-valores no; trabajo 7567366):

   | Ciclos | RSS máx. | Tiempo |
   |---|---|---|
   | 250 | 2,8 GB | 7,5 min |
   | 500 | 3,3 GB | 7,9 min |
   | 1000 | 5,2 GB | 12,0 min |
   | 2000 | 7,9 GB | 16,5 min |
   | 4000 | 14,4 GB | 27,1 min |
   | 10.000 | 32,1 GB | 58,5 min |

   Por encima de 1000 ciclos crece ~3 MB y ~0,3 s por ciclo, y a 10.000 ciclos ocupa 32 GB y casi 1 h (la extrapolación desde 2000 ciclos daba ~30 GB y ~1 h). El p de permulación es
   un recuento por conjunto, aditivo entre bloques de columnas, de modo que el procesamiento por bloques fijaría la memoria.
6. **RESAMPLE** (`permulations.R`), **medido** con los argumentos de la plantilla de neoplasia (`traitname` `neoplasia_prevalence`, `n_trait` `adult_necropsy_count`, `c_trait`
   `neoplasia_necropsy`; trabajo 7567399). A N = 1000 reproduce la cosecha de la corrida: 1184 sorteos, 920 de los 1000 primeros candidatos con 100 hipótesis y 1009 de 1100.

   | N | Tiempo | RSS máx. | `fop_pairs.tsv` | Sorteo en serie | Recuento FOP | Espejo FOP |
   |---|---|---|---|---|---|---|
   | 1000 | 1,1 min | 3,2 GB | 19,4 MB | 0,3 min | 0,3 min | 0,4 min |
   | 4000 | 4,2 min | 3,7 GB | 78,5 MB | 1,0 min | 1,2 min | 1,8 min |
   | 10.000 | 10,8 min | 4,4 GB | 197 MB | 2,4 min | 3,1 min | 4,6 min |

   Tiempo lineal (~0,066 s por ciclo), memoria de ~3,0 GB fijos más ~0,14 MB por ciclo, 100 hipótesis por ciclo (`fop_labelings.tab` de 118 MB a 10.000 ciclos). La cosecha FOP se calcula
   dos veces (recuento y espejo) y suma ~76 % del tiempo; el sorteo en serie, ~22 %. Con estos números RESAMPLE no limita el paso a 10.000 ciclos.
   Otras ineficiencias: el bucle de sorteos es serie; `reject_reasons <- c(...)` es cuadrático en rechazos (irrelevante con ~85 % de aceptación).

## 6. Diseño propuesto

### 6.1 Principio

El nulo se genera por gen (un worker tiene un gen y todos los ciclos) y se consume por ciclo (columnas del FCS, pool de `F(max)^n`) y por posición (`p.emp`). El merge es la
transposición entre ambos ejes. Si cada bloque emite estadísticos suficientes en lugar de filas, esa transposición se hace sobre tablas compactas. Los ciclos son
independientes, así que troceo por ciclos lleva a un merge de concatenación y unión.

### 6.2 Bloques de genes × ciclos y fusión de descubrimiento y replay

- Cada tarea recibe solo los etiquetados y el diseño FOP de su bloque de ciclos, los carga una vez y fusiona descubrimiento y replay por ciclos base, sin export en disco.
  La unidad natural es el trozo de etiquetados del kernel (`_default_b_chunk`), alineado con ciclos base completos.
- La RAM queda acotada por el bloque y no por N ni por el gen; desaparece el troceo por memoria del gen y las recargas de ASR por trozo.
- RESAMPLE escribe `fop_labelings` y `fop_pairs` en bloques de ciclos, para que cada tarea lea solo los suyos.
- Se eliminan el canal `permDiscFiles`, los exports de 13 a 16 MB por gen y la lógica de `prune_exports` (D3).
- `n_detected` solo rellena una columna del detalle y no tiene consumidor; se puede descartar o acumular en un contador.

### 6.3 Salida por bloque (estadísticos suficientes) y consumidores

| Consumidor | Necesita del nulo | Representación propuesta |
|---|---|---|
| FCS (`fcs_enrich.R`) | Matriz genes × N por dirección | Por (gen, ciclo): recuento de posiciones detectadas, bandera de cluster y scores puntuados por dirección (binario); el merge del bloque aplica el filtro de genes, forma los pools exactos por ciclo y deriva `F(max)^n` y el máximo crudo |
| FCS, estadístico nulo | Recuento por conjunto de ciclos con estadístico ≥ observado | Acumulación por bloques de columnas; memoria constante en N. Los controles "ciclo nulo como observado" conservan la matriz completa en un subconjunto de ciclos |
| SCORING (`p.emp`, `p.emp_fact`, `scoring_compute.R:431`) | Por posición: ciclos que la puntúan (`nd`), ciclos con score ≥ observado, clase de propensión, primer ciclo en que aparece | Almacén binario de (posición, ciclo, score) solo con score positivo; contadores por posición |
| Acumulación (`accumulation_test.py`) | Conjunto de posiciones con score positivo en algún ciclo | Indicador por posición (D4) |
| Informe de bioquímica y POSENRICH | Recuentos por (ciclo, esquema) y (posición, ciclo, score) | Fuera de alcance (D1) |
| Informe 11 | Cuantiles y muestra por ciclo (ya compactos) | Sin cambio |

El filtro de genes del nulo (`core.postproc.gene_removal`) se calibra dentro de cada ciclo, con el recuento de posiciones detectadas por (gen, ciclo) y la bandera de cluster, y se decide con todos los
genes del ciclo. Por eso las tareas no pueden emitir un pool ya agregado: emiten esos recuentos y los scores por (gen, ciclo), y el pool de cada ciclo se forma en el merge, sin las unidades
eliminadas. Se prescinde de las filas sin score de US y del formato de texto, no de los scores por (gen, ciclo).

El pool por ciclo se conserva exacto, no en histograma: `F` se eleva a la potencia n y el nulo está muy empatado, de modo que un error pequeño de la ECDF se amplifica n veces
(comentario de `_build_cycle_score_pools`). A 10.000 ciclos el conjunto de pools son unos 2,3 · 10⁸ valores (~1,8 GB en doble precisión; supuesto a partir de las ~2,3 · 10⁷ filas con
score de cáncer a 1000 ciclos).

### 6.4 Cambios por proceso

- `CAAS_CORE_BATCHED`: bloques de genes × ciclos; fusión; emisión de estadísticos; sin exports.
- `CAAS_CORE_MERGE`: unión de pools y concatenación de tablas por (gen, ciclo); sin pasadas sobre decenas de millones de filas. `scoring_caas_perms.R` deja de leer la tabla ancha
  entera o la lee por trozos (medido: ~3,4 MB por ciclo, ~34 GB extrapolados a 10.000). Dejar de construir el eje ASR (3 de las 6 matrices y el cuantil 90 por gen y ciclo) una vez
  descartados lectores externos.
- `RESAMPLE`: salida por bloques; caché de la cosecha FOP del recuento para el espejo (el recuento es ~29 % del tiempo a 10.000 ciclos) y sorteos en paralelo con semillas por sorteo (como el
  espejo, sin cambiar el resultado con el número de workers; el sorteo en serie es ~22 %). Prioridad baja: 10,8 min y 4,4 GB a 10.000 ciclos.
- `FCS_COMPUTE`: procesamiento por bloques de columnas con contadores (D5 lo da por asumible; medido: 32 GB y 58 min a 10.000 ciclos; es el cambio más sencillo si se
  quiere memoria constante).
- Limpieza (D6): `remove_intermediates` (borrar dentro de las tareas lo que ningún proceso posterior lee; en la permulación, con D3, no queda nada que podar) y `remove_work`
  (borrar `work/` al terminar con éxito) como opciones independientes. Falta decidir si `remove_intermediates` sigue condicionando qué se publica en resultados.

### 6.5 Lo que la reducción no resuelve

La CPU sigue siendo lineal en N: unas 50 CPU-horas por 1000 ciclos y corrida en `CAAS_CORE_BATCHED`, unas 500 a 10.000 ciclos. Solo bajaría con ciclos adicionales solo para los genes
candidatos (muestreo secuencial; referencia sin verificar: Besag y Clifford, Biometrika 1991, doi 10.1093/biomet/78.2.301). El FCS no lo admite, porque su nulo necesita todos los genes en
cada ciclo.

### 6.6 Fragmentación

Eje de partición: los ciclos. El nulo es aditivo entre ciclos: el pool, la puntuación de gen y el estadístico del FCS de un ciclo se calculan con los datos de ese ciclo, y los
consumidores globales (p de permulación del FCS, `p.emp`) son sumas de contadores. Un bloque de B ciclos es una unidad independiente de extremo a extremo y B no cambia el resultado.

Cadena fragmentada (propuesta):

1. `RESAMPLE`: etiquetados e hipótesis FOP escritos por bloques de B ciclos base, más el diseño de `b_0` aparte.
2. `CAAS_CORE` (lote de genes × bloque): fusiona descubrimiento y replay por ciclo base y emite, por (gen, ciclo), el recuento de posiciones detectadas, la bandera de cluster y los scores
   puntuados por dirección (binario).
3. `CAAS_CORE_MERGE` por bloque (dispersión sobre bloques): une los lotes de genes del bloque, aplica el filtro de genes (calibrado dentro de cada ciclo), forma los pools por ciclo, calcula los
   scores de gen y escribe la matriz genes × B de los tres rankings. Una reducción final pequeña fija el roster de ciclos y el universo.
4. FCS: dispersión sobre (lote de GMT × bloque). Hoy `FCS_COMPUTE_BATCHED` ya reparte las bases de datos, pero cada tarea carga el RDS completo; con bloques cada tarea lee solo la matriz de
   su bloque y devuelve, por conjunto, el número de ciclos con estadístico igual o más extremo que el observado y el número de ciclos válidos. La reducción suma, calcula el p de permulación
   `(recuento + 1) / (válidos + 1)` y aplica BH y la evidencia. Exige separar en `fcs_run_all` el estadístico observado, el recuento nulo por bloque y el cierre.
5. SCORING (`p.emp`): contadores por posición y bloque (ciclos que la puntúan, ciclos con score al menos igual al observado, primer ciclo en que aparece), sumados en la reducción.
   `p.emp_fact` necesita la propensión global antes de ajustar sus clases, de modo que se hace en dos etapas: suma de `nd`, definición de clases y acumulación de scores por clase.

Memoria por tarea: depende de B y no de N. Con B = 500: matriz genes × 500 de ~60 MB por ranking en el FCS, ~1,7 GB en el cierre del merge (a ~3,4 MB por ciclo, medido) y 50.000
etiquetados por proceso en el replay (~170 MB, medido en el coste de parseo). A 10.000 ciclos son 20 bloques.

Tamaño de bloque: fijo. Una vez el replay procesa por ciclo base, su memoria no depende de B ni de las filas del gen, así que B solo equilibra el número de tareas frente a la sobrecarga
fija por (gen, bloque) (carga del contexto de ASR, pendiente de medir). Un B calculado por memoria solo haría falta si el replay siguiera reteniendo las filas.

Verificado en el código:

- `gene_removal` agrupa las unidades por ciclo (`pools` por `u.labeling`) y calcula los umbrales (cuartiles del recuento por gen y percentil de la densidad de aciertos por longitud) dentro de
  cada ciclo: es local al bloque, siempre que el bloque reúna todos los genes del ciclo (de ahí que el merge sea por bloque y no por lote de genes).
- `p.emp_fact` (`scoring_compute.R:129`) se reduce de forma exacta. Su p es `min(1, (nd + 1) / (N + 1) × S)` con `S = (1 + nº de detecciones de la clase con score > obs − TIE_TOL) / (1 + tamaño de
  la clase)`, así que basta, por posición observada, ese recuento y el tamaño de cada clase, aditivos entre bloques; no hace falta conservar los pools. Las clases (`.fact_breaks`, percentiles de
  las detecciones del nulo) dependen solo del histograma de `nd` por posición, que sale de sumar los contadores de la primera etapa. El filtro de universo (posiciones probadas y observadas) es
  estático por fila, y la curva de saturación sale del mínimo, entre bloques, del primer ciclo en que cada posición aparece.

Por verificar: que el roster, el denominador N, la familia de BH (posiciones del universo) y el resultado final se reconstruyan a partir de las sumas con los mismos números que la versión actual
(prueba de equivalencia de la sección 7).

## 7. Exactitud y validación

Debe conservarse: el pool exacto y la tolerancia `TIE_TOL` (1e-12, comparaciones por `>=` sobre valores casi empatados); el orden de visita de las hipótesis (`labelings.py`: el orden entra
en sumas en coma flotante); las filas contiguas por gen; la semántica del roster de ciclos (N cuenta ciclos sin detección); y que `b_0` coincida con el observado.
Validación: comparación antiguo frente a nuevo en el toy y en un subconjunto real (estadístico de gen idéntico, `p.emp` idéntico, recuentos del FCS idénticos), con el mismo estilo de
las pruebas de streaming existentes en `validation/unification/`.

## 8. Mediciones

Realizadas:

| Qué | Dónde |
|---|---|
| Controles del FCS "ciclo nulo como observado", 4 corridas (72 tareas) | `~/scratch/0.Phylophere/.tmp/fcs_calib_ctrl/` y `fcs_null_ctrl/`; copias en `~/.cc/calib_ctrl/`, `~/.cc/strict_ctrl/` |
| Ranking con máximo crudo, neoplasia estricto (trabajos 7566445, 7566446, 7566449) | `.tmp/fcs_rawmax_ctrl/`; copia en `~/.cc/calib_ctrl/rawmax/` |
| Tres pruebas de posición, gen y BH (scripts `three_tests.R`, `gene_bh.R`) | `~/.cc/calib_ctrl/` |
| RESAMPLE con los argumentos de la plantilla de neoplasia a N = 1000, 4000 y 10.000 (7567399) | `.tmp/resample_scaling/real_N<k>/` |
| Coste fijo del fichero de etiquetados por proceso `ct perm-replay` (7566496) | `.tmp/resample_scaling/labelings_cost/` |
| `scoring_caas_perms.R` a 250 a 4000 ciclos (7567372; tabla recortada o reciclada) | `.tmp/scoring_perms_scaling/` |
| FCS de neoplasia completo con el ranking calibrado a 250 a 10.000 ciclos (7567366; por encima de 1000, columnas recicladas) | `.tmp/fcs_scaling/` |
| Trazas de las pasadas finales de neoplasia y maligno | `~/.cc/resource_scale/` |

La medida de recursos del FCS vale con ciclos reciclados; los p-valores no. Memoria máxima con
`sacct -j <trabajo> --format=JobID,Elapsed,MaxRSS` en los pasos `.batch`.

Pendientes: efecto en el pipeline completo (diseño de la medición por decidir); perfil de memoria de un gen gigante (parte del proceso padre frente a los workers, composición de los
~4,4 KB por fila); tiempo de carga del contexto de ASR por bloque; control posicional y ranking de máximo crudo en las corridas completas.

## 9. Fases propuestas

| Fase | Contenido | Efecto esperado | Riesgo |
|---|---|---|---|
| 1 | Etiquetas cargadas una vez por tarea y bloques por ciclos en `CAAS_CORE_BATCHED`; sin `permDiscFiles` | Elimina el término N × H por gen y por tarea | Medio: cambia la unidad de reparto |
| 2 | Procesamiento por ciclo base dentro del worker (sin objetos de todos los ciclos) y fusión del descubrimiento | RAM del replay acotada y sin export en disco | Medio: debe dar los mismos números |
| 3 | Estadísticos suficientes en lugar de filas y nuevo merge | Almacenamiento y merge lineales y pequeños | Alto: toca SCORING y el informe 11 |
| 4 | FCS por bloques de columnas | Memoria constante en N | Bajo |
| 5 | RESAMPLE: salida por bloques y caché de la cosecha FOP | ~30 % del tiempo de RESAMPLE (10,8 min a 10.000 ciclos) | Bajo |

Cada fase se cierra con la comparación antiguo frente a nuevo de la sección 7.

## 10. Decisiones y preguntas abiertas

1. Regla unificada del FCS para RER y CAAS (P4 de la revisión), con las corridas estrictas de CAAS marcadas como no concluyentes.
2. Tamaño de bloque: fijo o por memoria (D7).
3. Qué publica cada opción de limpieza (D6).
4. Si 10.000 ciclos se aplican a las cuatro corridas o solo a neoplasia estricto, donde está la pregunta científica.
5. Si se conservan las filas sin score de US del nulo (necesarias solo para el informe de bioquímica).
