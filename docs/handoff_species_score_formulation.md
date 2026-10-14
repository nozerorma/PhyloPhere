# Plan de traspaso: nueva formulación del score de posición (sin hipótesis, con especies)

Fecha de redacción: 2026-10-13. Rama de trabajo: `scoring_v2`. Autor del plan: sesión de Claude Code con Miguel Ramón.
Este documento es autocontenido: una conversación nueva puede ejecutarlo sin ver la historia de la anterior.

## 1. Objetivo

Sustituir el pooling del score de posición de `pool_domains` (media ponderada por PSS sobre hipótesis, con shrinkage `M/(M+2)`)
por una formulación **sin hipótesis, sin PSS y con un factor de especies**, decidida por Miguel tras comparar variantes
en PEPC y en un lote de cancer_def. Es una decisión tomada ("me parecen razonables los cambios"); el trabajo pendiente es
**implementarla limpia, probarla y rehacer el nulo**, no volver a discutirla.

## 2. La fórmula acordada

Para una posición, un lado (top o bottom) y los `K` dominios del diseño (clados disjuntos, uno por par de cada hipótesis):

1. **Dominio cambiado:** un dominio ha cambiado al residuo derivado `r` si alguna de sus especies lleva `r` donde el estado
   ancestral era otro (se usan las especies de la unión de las hipótesis; cada especie cuenta una vez, no importa en cuántas
   hipótesis aparezca).
2. **Independencia por pareja de dominios:** `c(a,b) = 1 - P(el derivado ya estaba en el LCA de a y b)`. Es la contribución
   que ya calcula `score_domains_side` (lista `pair_lca` por hipótesis). Entre varias parejas de especies de los mismos dos
   dominios se toma la **media** de las `c` de las combinaciones distintas (decisión de Miguel: media, no máximo).
3. **Score de un dominio:** noisy-OR de las `c` con sus parejas: `s(d) = 1 - Π (1 - c)`; 0 si no tiene pareja.
4. **Fuerza:** media simple de los `K` dominios: `fuerza = Σ s(d) / K`.
5. **Amplitud:** `n_sp / S_lado`. `n_sp` = especies distintas que llevan el derivado **en ese lado** (la especie del par que está
   en el grupo de ese lado) dentro de un dominio con score positivo en alguna hipótesis. **`S_lado` es por lado:** especies
   distintas del grupo de ese lado en el diseño del ciclo, es decir `species1` (grupo top) para el lado top y `species2` (grupo
   bottom) para el lado bottom, unión sobre todos los pares de todas las hipótesis del ciclo. `n_sp` y `S_lado` son por lado;
   no se divide por la suma de los dos grupos (con la suma el score solo llegaría a ~0.5 y se trataría distinto a un lado con
   el grupo más grande).
6. **Score:** `fuerza × n_sp / S_lado`, por lado. Entre 0 y 1 por construcción (`fuerza <= 1`, `n_sp <= S_lado`); un recorte final
   `max(0, min(1, x))` es solo una red de seguridad y no cambia ningún valor. No hay PSS, ni `M`, ni `H`, ni `Meff`, ni
   shrinkage `M/(M+2)`.

Ejemplo de control (K = 3, lado top con S_top = 10): A y B convergen con `c = 0.9`, el residuo está en 2 especies de A y 1 de B
(del grupo top), C no cambia -> `s = (0.9, 0.9, 0)`, fuerza 0.6, amplitud 3/10 = 0.3, score 0.18.

Por qué se divide por S (y por qué por lado): `n_sp` depende de cuántas especies tiene el diseño, y los diseños de los ciclos
nulos difieren del observado (cancer_def: 19 especies de mediana en el nulo frente a 17 observadas, rango 15 a 27; PEPC: de 22
a 73 frente a 48 y 49). Sin dividir, un ciclo nulo con más especies puntúa más solo por tamaño: en cancer_def salen 62
posiciones con p <= 0.05 sin dividir y 175 dividiendo; en PEPC genotípico 23 y 34; en el fenotípico no cambia (24 y 24). Dentro
del observado el orden del score es idéntico (Spearman 1.000); el efecto es sobre la comparación con el nulo.

Consecuencias que Miguel conoce y acepta: el número de dominios deja de mandar de forma estricta (el factor de especies puede
dar la vuelta a un dominio de diferencia en 1 a 6 % de las comparaciones observadas, menos que con el pooling actual, 7 a 13 %;
en cancer_def solo afecta a 10 de las 21 posiciones de 3 dominios, todas con una especie por dominio); orden de significancia
casi igual al de la línea base (Spearman 0.96). Con `S_lado` el score sí ocupa [0, 1] (máximos observados 0.77, 1.00 y 0.85); con
`S` de ambos grupos quedaba en torno a 0.5.

Comparación `S` de ambos grupos frente a `S_lado` (misma réplica, solo cambia el denominador): cancer_def Spearman de -log p 0.98,
p <= 0.05 en 175 y 152 posiciones; PEPC AUC del score 0.825 y 0.825 (genotípico), 0.816 y 0.829 (fenotípico). Rendimiento
equivalente; se elige `S_lado` por fidelidad al significado y por cubrir [0, 1].

## 3. Evidencia que respalda la decisión (no repetir estos experimentos)

PEPC (un gen, 9 sitios de referencia de `validation/truthsets/tier1/pepc_c4.sites.tsv`; AUC del score, genotípico / fenotípico; todas las filas con `n_sp/S` se midieron con `S` de ambos grupos sumados; con `S_lado` el AUC es 0.825 / 0.829 para la fórmula final):

| Formulación | AUC |
|---|---|
| Línea base (hipótesis + PSS acotado + factor de especies) | 0.82 / 0.84 |
| **Sin hipótesis ni PSS + `n_sp/S`** | **0.83 / 0.82** |
| Sin hipótesis, PSS acotado en la media + `n_sp/S` | 0.81 / 0.81 |
| Sin hipótesis, PSS dentro del noisy-OR + `n_sp/S` | 0.81 / 0.82 |
| Sin hipótesis ni PSS, **sin** factor de especies | 0.69 / 0.64 |
| Orígenes independientes (parsimonia) en vez de especies | 0.79 / 0.73 (descartado) |

cancer_def (lote de 200 genes, 463 posiciones observadas, 1000 ciclos nulos, sin sitios de referencia): con factor de especies
todos los enfoques dan el mismo orden de significancia (Spearman 0.90 a 0.99; top 25 de A y B idéntico; esos números se midieron con `S` de ambos grupos, la versión por lado da Spearman 0.98 con ella); sin factor de especies
sale otro orden, con muchas posiciones en el suelo del p y p anticorrelacionado con el soporte. El PSS dentro del noisy-OR es
inestable (casi todo significativo en PEPC, la mitad en cancer_def).

Por qué se quita el PSS: sin ganancia medible; el selector y el nulo (`match_pss`) ya lo usan; los pesos acotados
(`±1/(2K)`) son casi planos; la única versión con efecto real depende de normalizar por el máximo de cada diseño.
Límite honesto: PEPC tiene un rango de PSS estrecho (1.3 a 3.4); no se ha probado con rango amplio y sitios de referencia.

## 3b. Dónde se consume el score (revisado: el cambio de escala no rompe nada)

No hay ningún umbral fijo sobre el valor del score en el código ni en la configuración (búsqueda de comparaciones de
`CAAS_score`/`asr_path_score` con constantes: ninguna; el único `thr` de score es el de red STRING de DOMINO). Todo consumo es por
orden o por puntuación positiva:

- **p.emp, p.emp_fact, p.adj_bh:** comparan cada posición con su propio nulo; la escala no importa, pero sí que observado y nulo
  estén normalizados igual (lo que hace `S_lado`).
- **Score de gen (`size_adj_max`) y FCS:** ECDF del propio pool, por rango.
- **POSENRICH:** la parte hipergeométrica usa solo qué posiciones son positivas; la de rangos, solo el orden.
- **Acumulación:** cuenta posiciones con score positivo, no su magnitud.
- **Informes:** clases por percentil. Los textos que muestren el valor tal cual deben explicar que es una fracción (medianas
  observadas de 0.08 a 0.10 con `S` de ambos grupos; algo mayores con `S_lado`).

## 4. Estado actual del repositorio (al escribir esto)

- Rama `scoring_v2`, HEAD `fe2648c`. Commits relevantes ya hechos: `8cf026f` (FADE), `fede4d3` (chunking por memoria, PSS leído
  del mapa, `n_hyp` en el nulo), `ae6e580` (VEP modo `auto`), `e6bea38` (nulo unificado), `69f2584` (parámetros/GUI),
  `518e3c1` (`n_hypotheses_eff`/`n_hyp_eff`), `fe2648c` (`n_species`/`n_sp`).
- **9 ficheros modificados sin commitear en el árbol de trabajo** (`fop_pool.py` con el interruptor `PHYLOPHERE_POOL_VARIANT`,
  los orígenes cableados, `origins.py`, `tests/fop_pool/test_origins.py`, y toques en `disambiguate_single.py`, `driver.py`,
  `master.py`, `pooling.py`, `models.py`, `gene_wrapper.py`). **No deben commitearse.** Están guardados en la rama
  `exp/pooling-variants` (commit `dfd674f`, ya empujada al remoto) y en `~/.cc/exp_variants.patch`. Primer paso: restaurar el
  árbol con `git checkout -- <esos ficheros>` y borrar `origins.py` y `test_origins.py` (los dos archivos nuevos).
- El cluster (correfoc) tiene el repo en `scoring_v2` limpio y la rama `exp/pooling-variants` traída con `git fetch` (sin checkout).
- Nada de lo nuevo se ha ejecutado de extremo a extremo con Nextflow; los `.nf`, los informes completos y los parquet reales
  siguen sin ejecutar.

## 5. Qué implementar (en este orden)

1. **Restaurar el árbol** (ver arriba). Comprobar `git status` limpio y los tests existentes.
2. **Calcular `S_lado` por ciclo** (dos números: grupo top y grupo bottom). Regla: nunca la suma de los dos grupos.
   - Observado: especies distintas de `species1` (top) y de `species2` (bottom) de las filas de `contrast_hypotheses_pairs.tsv`
     (o de las filas `b_0` de `b0_pairs.tsv`).
   - Nulo: especies distintas de `species1` y `species2` de las filas de `fop_pairs.tsv` de cada ciclo base.
   - Sin FOP (diseño de una hipótesis): especies de los pares de contraste de esa etiqueta, separadas por grupo.
   - Llevarlo a `pool_domains` como argumento (por ejemplo `design_species_by_side: Optional[Dict[str, int]]` con claves
     `"top"` y `"bottom"`); dentro de `_agg(side)` se usa `S_lado = design_species_by_side[side]`. Puntos de cableado:
     `core/labelings.py` (`read_pss` ~279, `observed_pss` ~345: añadir una lectura de especies por ciclo), `core/driver.py`
     (`pool_labelings`, cada llamada a `pool_domains`), `utils/gene_wrapper.py` (carga de `fop_pairs`),
     `convergence/disambiguate_single.py` (`_emit_pooled_side_rows`, `analyze_gene_disambiguation`), `core/observed.py` y
     `observed_b0_main.py`. Si `S_lado` no está disponible, usar factor 1 y avisar en el log (no fallar).
3. **Reescribir `pool_domains`** (`subworkflows/CT_DISAMBIGUATION/local/src/convergence/fop_pool.py`) con la fórmula de la
   sección 2. Quitar: pesos PSS y `delta`, `w_bar`, `f_fop`/`k_0`, `effective_number` y `n_eff`, el interruptor
   `PHYLOPHERE_POOL_VARIANT` y todas sus variantes. Mantener: `agree_num/agree_den`, `convergence_type`, `n_hypotheses`, `n_sp`,
   `participating_hyps`, `domain_scores` (informativo). Referencia de implementación: la rama `union` del experimento
   (`git show exp/pooling-variants:subworkflows/CT_DISAMBIGUATION/local/src/convergence/fop_pool.py`, bloque
   `if _VARIANT.startswith("union")`): construye `comb[(a, b, sp_a, sp_b)] = c` desde `pair_lca` y `domain_der_species`, agrupa por
   pareja de dominios con la media, hace noisy-OR por dominio y media sobre `rep.values()` (los `K` dominios). Añadir el factor
   `n_sp / S_lado` y el recorte de seguridad a [0, 1]. Nota: en el experimento el `core` aún se multiplicaba por `f_fop` y el
   análisis lo deshacía fuera (`s * (M + 2) / M`); en la implementación final **no hay `f_fop`**, así que `asr_path_score` ya es
   el score final.
4. **Quitar `n_hypotheses_eff`/`n_hyp_eff`** de `models.py`, `pooling.py`, `driver.py`, `disambiguate_single.py`, `master.py`
   y `gene_wrapper.py`. **Mantener `n_species`/`n_sp` y `n_hypotheses`/`n_hyp`.**
5. **Subir `SCORE_RULE`** (módulo `core/scores.py`; el guard de `scoring_compute.R` compara el sello entre observado y nulo) para
   que un nulo viejo se rechace. Una nota de memoria anterior dice que el sello no se subió en un cambio previo: revisarlo.
6. **Tests** (carpeta `tests/fop_pool/`): reescribir `test_fop_pool_consistency.py` (hoy prueba el shrinkage `M/(M+2)` y
   `effective_number`, que desaparecen). Casos: ejemplo de control de la sección 2 (0.09); especies contadas una vez aunque
   aparezcan en varias hipótesis; una especie de un dominio con score 0 no cuenta; sin `pair_lca` el score es 0; media (no máximo)
   entre combinaciones; observado y nulo con los mismos registros y la misma `S_lado` dan el mismo score; `S_lado` ausente da
   factor 1; `S_top` y `S_bottom` distintos dan amplitudes distintas para el mismo `n_sp` en cada lado; el score nunca supera 1.
   Mantener los tests de chunking, nulo unificado y VEP existentes (`tests/ct_replay_chunking`, `tests/null_universe`, `tests/vep_side`).
7. **Documentación:** `docs/scoring_v3_core.md` (contrato del score), docstring de `fop_pool.py`, textos del informe 11 y cualquier
   mención de `M/(M+2)`, pesos PSS acotados o `Meff` (`grep -rn "M + 2\|M/(M\|k_0\|delta\|Meff" docs subworkflows README.md`).
   Comentarios solo descriptivos del comportamiento actual (regla del proyecto: sin historia, sin fechas).
8. **Commits** pequeños y por tema (Miguel hace el push). Atribución: línea `Co-Authored-By` del recordatorio de la sesión.

## 6. Verificación (con los arneses ya construidos)

Todo en `~/.cc` (nunca en `/tmp`); si ya no existen, reconstruir con los comandos de la sección 8.

- **PEPC, local:** `~/.cc/pepcval/runv.sh` (réplica del nulo y del observado con el código actual para los dos fenotipos desde los
  directorios de trabajo de `validation/tier1/output/pepc/work/...`) y `~/.cc/pepcval/analyze6.py` (AUC, p empírico, sitios de
  referencia). Esperado con la formulación final: AUC del score cercano a 0.83 (genotípico) y 0.82 (fenotípico), 9 y 8 sitios de
  referencia con p <= 0.05. **Cuidado con la RAM:** el proceso principal del replay llega a 7.5 GB; con 15 GB de máquina hay que
  cerrar el navegador y correr con `WORKERS=1` o 2 y el guardián `~/.cc/pepcval/guard.sh` (para las réplicas si la RAM disponible
  baja de 1.2 GB). Una réplica a la vez.
- **cancer_def, cluster por sbatch:** `~/.cc/poolcmp/poolcmp.sbatch` (réplica del lote de 200 genes con el código de una rama vía
  `git archive`; hay que cambiar `origin/exp/pooling-variants` por la rama o el commit con el código final y que Miguel lo empuje y
  haga `git fetch` en el cluster) y `~/.cc/poolcmp/analyze.py`. Esperado: unas 175 posiciones con p <= 0.05, 26 con p <= 0.01,
  rangos casi iguales a la versión con hipótesis (Spearman 0.96 por p).
- Comprobar que el score observado y el nulo salen en [0, 1] y que ningún p queda en el suelo en masa (señal de sesgo de diseño).

## 7. Restricciones y preferencias (importantes)

- **Cluster (correfoc):** nunca ejecutar cómputo en el nodo de entrada; solo `sbatch` o sesión interactiva; lecturas ligeras sí
  (`ls`, `sacct`, `cat`). Directorio temporal: `~/scratch/0.Phylophere/.tmp` (nunca `/tmp`). **Nunca editar código en el cluster**;
  el código llega por git (Miguel hace el push; en el cluster basta `git fetch` y `git archive` a `.tmp`). `sbatch` se invoca como
  `ssh correfoc "bash -l -c 'sbatch'" < script.sbatch`; en nodos de cómputo no está `micromamba` en el PATH (usar
  `export PATH=/homes/users/mramon/micromamba/envs/phylophere/bin:$PATH`) ni `/usr/bin/time`.
- **Local:** nunca usar `/tmp` (usar `~/.cc`); la máquina es justa de RAM (ver arriba).
- **Estilo:** sin rayas largas (em dashes) en ningún texto; chat en español natural; documentos y comentarios descriptivos,
  sin historia. No ser complaciente: decir si una conclusión es sólida o un límite de los datos. Preferencia explícita de Miguel:
  fórmulas simples; "no complicar la fórmula en exceso".
- **Git:** commitear solo cuando Miguel lo pide (aquí ya lo ha pedido para el trabajo de esta línea); nunca empujar.

## 8. Artefactos y limpieza

- Local: `~/.cc/pepcval*`, `~/.cc/poolcmp` (datos copiados del cluster: 100 MB), `~/.cc/exp_variants.patch`, `~/.cc/carnsel`. Borrar cuando
  se cierre el trabajo (los resultados no son reproducibles solo con el repo: guardar las tablas si se quieren citar).
- Cluster: `~/scratch/0.Phylophere/.tmp/poolcmp` (4 variantes, unos 100 MB cada una). Borrar al terminar.
- Rama remota `exp/pooling-variants`: dejarla como registro de los experimentos o borrarla con el visto bueno de Miguel.

## 9. Pendiente fuera de esta línea (para no perderlo)

- Relanzar Carn completo con el código final (Carn actual está "pre muchas cosas": su scoring es anterior a `abe8b43`).
- VEP: modo `auto` por defecto implementado y probado con la herramienta (70 % de identificadores de TP53 resueltos), pero **no** probado
  con el script del pipeline ni en modo offline; `vep_ensembl` sigue desactivado por defecto. Miguel decidió dejarlo.
- Decisión abierta registrada: el modo `permulation` de la acumulación se mantiene (el pool de `naive` y `cons_decile` es ahora U).
- Ningún `.nf` modificado en esta línea de trabajo se ha ejecutado con Nextflow.
