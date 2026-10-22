# WP1 FCS/RER: estado actual

**Estado: en revisión; no cerrado.** Este registro resume las decisiones vigentes y el trabajo abierto a 2026-10-22. `WP1_fcs_rer.md` y las notas anteriores conservan el historial y decisiones que fueron reemplazadas; para el estado actual, aplicar esta página. Los commits de implementación no equivalen a aceptación integrada.

## Decisiones confirmadas

| Tema | Estado vigente |
|---|---|
| Ranking CAAS | `-log10(p)` confirmado por el usuario; su integración productiva sigue ligada a WP6, no se activa aquí por un interruptor crudo/calibrado. El estadístico que produce ese p sigue abierto: máximo crudo, `F(max)` o `F(max)^n`; corresponde a WP9. |
| Comparador de rutas | `stat = AUC - 0.5` observado y nulo; CAAS unilateral (`null >= observado`), buscando evidencia alta de convergencia; RER bilateral (`abs(null) >= abs(observado)`). Ambos conservan empates inclusivos y pseudo-cuenta. La lateralidad distinta responde a hipótesis distintas y se declara como adaptación frente al paquete. |
| Regla RER/FCS principal | `sig_perm = BH(p.perm) < fdr_thr`, BH por ranking y base de datos. |
| Regla analítica de Redlich | `redlich_level` informa el menor nivel analítico de `0.05, 0.10, 0.15` que cumple `p.adj < nivel` y `p.perm < 0.025`; niveles configurables. Es independiente de `sig_perm`. |
| `sig_top` | Retirado del esquema y de los informes. No debe usarse como nombre alternativo ni conservarse desde entradas antiguas. |
| Signo del efecto | No hay filtro por signo en ninguna regla. Para rankings RER se fuerza prueba bilateral, de forma que `abs(stat)` conserva enriquecimiento en ambas direcciones; el signo se mantiene como descripción del efecto. |
| RER observado y nulo | Corregidos el tratamiento de NA (excluir score no evaluable y conservar ceros medidos), la lateralidad y el uso compartido de `min.sp`, `min.pos`, winsorización y `clade` entre observado/nulo. |
| Parámetros analíticos | El fallback de `analytic_levels` ya está corregido en el cambio local principal. |

## Pendientes y límites

| Tema | Pendiente vigente |
|---|---|
| Estadístico base CAAS | La decisión entre máximo, `F(max)` y `F(max)^n` sigue abierta en WP9; falta el control de ciclo nulo como observado y la evaluación a 10.000 permulaciones. |
| Permulaciones RER | El marino de 1.000 permulaciones existente no se ha regenerado con los ajustes observado/nulo armonizados. Quedan sin ejecutar la corrida a 10.000 y la decisión de defaults. La esquina completa de RER_CONT tampoco está medida. |
| Configuración | Revisar y armonizar valores por defecto, GUI/configuración y plantillas de validación; no asumir que todos quedaron actualizados. |
| Comparador CAAS | Decisión confirmada e implementada: usar `AUC - 0.5`, unilateral en CAAS. Validado en 12 ciclos marinos repartidos entre los 100 guardados, tres rankings y modos crudo/gene_p: 9.576 AUC contra la función nativa y conteos independientes de p, diferencias máximas 0. Falta recalcular las salidas completas afectadas. |
| Empates y pseudo-cuenta | Justificación publicada para `>=` y `(1+k)/(M+1)`, condicionada a los supuestos del muestreo nulo. Se conservan como adaptación explícita; no prueban por sí solos calibración biológica. |
| Duplicados y p extremo | Corregidos localmente: observado RER usa P sin suelo, como el nulo; ±Inf se conservan como extremos y NA/NaN se excluyen. Se rechazan IDs exactos repetidos y colisiones de prefijo antes del punto en RER, sin elegir arbitrariamente un score. Cada ID válido conserva su sufijo. |
| Validación externa | Sincronizados `test_fcs_no_null.py`, `check_fcs_blocks_nf.R` y `check_fcs_blocks_real.R` con la retirada de `sig_top`, tras autorización de escritura fuera del repositorio. No se ha ejecutado la suite externa completa con sus fallos históricos. |

## Verificación mecánica registrada

- `Rscript tests/test_fcs_rules.R` con `/home/miguel/micromamba/envs/phylophere/bin/Rscript`: **PASS**. Comprueba independencia de reglas, niveles custom, input antiguo con `sig_top`, efectos negativos, ausencia de nulo y esquema vacío.
- Render HTML completo de informe 12 con entrada marina precalculada antigua: **PASS**. El HTML resultante se generó con el esquema nuevo y el informe elimina la columna legacy.
- Pytest con entorno `phylophere`, heredando `PATH`, temporales en IBE y sin caché: **17 passed** (`test_fcs_no_null.py` adaptado, `test_fcs_cluster.py`, `test_rer_statistic.py`). Los tres originales externos indicados arriba ya están sincronizados; las copias/diffs de preparación se conservan como artefactos.

## Tablas y artefactos de diagnóstico

Los siguientes TSV conservan las comparaciones de trabajo; su interpretación actual se resume después:

- [`caas_holdout_summary.tsv`](/home/miguel/.codex/visualizations/2026/10/22/01a16449-db13-7541-adac-2fa17f0242e3/wp1_rer_fix/caas_holdout_summary.tsv)
- [`getstat_real_scan.tsv`](/home/miguel/.codex/visualizations/2026/10/22/01a16449-db13-7541-adac-2fa17f0242e3/wp1_rer_fix/getstat_real_scan.tsv)
- [`getstat_synthetic_edges.tsv`](/home/miguel/.codex/visualizations/2026/10/22/01a16449-db13-7541-adac-2fa17f0242e3/wp1_rer_fix/getstat_synthetic_edges.tsv)
- [HTML del informe 12](/home/miguel/.codex/visualizations/2026/10/22/01a16449-db13-7541-adac-2fa17f0242e3/wp1_rer_fix/render/marine_rer_fcs.html)
- [Copias y diffs de tests/checks externos](/home/miguel/.codex/visualizations/2026/10/22/01a16449-db13-7541-adac-2fa17f0242e3/wp1_rer_fix/external_test/)

## Comparador CAAS: evidencia previa y decisión

### Control holdout previo a la decisión

Se usó el nulo real guardado de marino: 100 ciclos, estadístico base `size_adj_max`, mínimo 10 y máximo sin límite. Se hicieron tres particiones de 50 ciclos de referencia y 50 de evaluación, intercambiando ambos lados. Cada ciclo evaluado queda fuera de su referencia en ese fold; el ranking `gene_p` de evaluación se calcula contra la referencia, sin usar una transformación ajustada al conjunto completo. Se mantienen `>=` y la pseudo-cuenta en ambos comparadores. Las cifras de la tabla son proporciones de p de conjunto que cumplen p≤0,05, agrupando rankings y bases:

| Ranking de gen | Comparar AUC | Comparar p analítico |
|---|---:|---:|
| Crudo, como referencia de validación | 3,826 % | 3,883 % |
| −log10(p de gen), decisión de uso | 3,642 % | 3,698 % |

A nivel p≤0,15, las proporciones son 13,458 %/13,577 % para crudo y 13,058 %/13,205 % para `gene_p`. En `gene_p`, entre el 8,9 y el 9,4 % de los p evaluados cambian según el ranking; diferencia máxima absoluta 0,0980. No aparece un exceso global en estos cortes, pero eso no demuestra calibración por base, por gen o en datos reales observados. Las bases pequeñas tienen tasas descriptivas variables; por ejemplo KEGG alcanza 18 % de p≤0,05 en un ranking, con un solo conjunto evaluado, por lo que el promedio global tampoco certifica cada colección.

Límites: son solo 100 ciclos únicos reutilizados; los folds, conjuntos y genes no son réplicas independientes. El tamaño de referencia 50 limita la resolución a 1/51. Hay 119.700 evaluaciones válidas por método y modo de ranking, de 121.500 posibles; faltan 1.800 evaluaciones de HPO que no satisfacen los requisitos para estimar el p en el fold. Se excluyen de las tasas y no se cuentan como resultados negativos. Los recuentos de filas con q<α no se interpretan como FDR. Tampoco es una comparación de potencia. No se ha comparado aquí el estadístico base máximo/F(max)/F(max)^n; eso sigue siendo WP9.

Conclusión operativa: este control mostró que ambas cantidades producen resultados distintos y no demostró equivalencia. La decisión posterior del usuario es adoptar AUC por correspondencia con RERconverge. El control es diagnóstico previo, no validación del código recién cambiado.

Scripts, resultados por base y resumen: [check_caas_tail_choices.R](/home/miguel/.codex/visualizations/2026/10/22/01a16449-db13-7541-adac-2fa17f0242e3/wp1_rer_fix/check_caas_tail_choices.R) y [caas_holdout_results.md](/home/miguel/.codex/visualizations/2026/10/22/01a16449-db13-7541-adac-2fa17f0242e3/wp1_rer_fix/caas_holdout_results.md).

Decisión vigente posterior al control: CAAS compara el AUC centrado (`stat = AUC - 0.5`) con el nulo, usando el mismo estadístico que RERconverge y una cola superior para su hipótesis de evidencia alta. RER conserva el contraste bilateral. El control holdout anterior no valida la implementación nueva; la comprobación independiente sobre el fixture marino se registra al final de esta página. Falta recalcular las salidas completas cuyo `p.perm`, `q.perm` o clasificación cambien.

## Empates, signo e identificadores (auditoría previa)

### Justificación de empates y pseudo-cuenta

[Phipson y Smyth (2010)](https://pubmed.ncbi.nlm.nih.gov/21044043/) explican el problema de devolver cero al muestrear un número finito de permutaciones. [Hemerik y Goeman (2018), secciones 3.3–3.4](https://doi.org/10.1007/s11749-017-0571-1) fundamentan el conteo inclusivo y la inclusión del observado bajo supuestos explícitos de invariancia/muestreo; los empates pueden hacer conservador el procedimiento. No se puede trasladar una garantía de exactitud a cualquier generador filogenético sin verificar esos supuestos.

Hay además una comprobación extrema reproducida con `RERconverge::permpvalenrich`: si los 100 estadísticos nulos son idénticos al observado, el conteo estricto sin pseudo-cuenta devuelve p=0; el inclusivo con pseudo-cuenta devuelve p=1. No hay evidencia de extremidad en ese ejemplo. Mantener `>=` evita ese resultado del conteo estricto; mantener +1 evita ceros cuando el observado supera todos los ciclos muestreados. Esto justifica la adaptación de cálculo, no certifica que el modelo nulo de CAAS sea adecuado.

### Signo

`permpvalenrich` del commit de referencia compara `abs(stat)` y no exige `stat>0`. El AUC ya está centrado en 0,5 al escribir `stat`. En rankings RER firmados, el sentido negativo representa concentración en el extremo negativo y se conserva. En CAAS se mantiene una cola superior sin valor absoluto, porque el objetivo son rutas concentradas entre p de gen pequeños. Esto no añade un filtro posterior `stat > 0`; la elección de cola define el contraste. No se añade un filtro de signo a las reglas. La interpretación de un conjunto negativo de CAAS —depleción de score de convergencia— debe distinguirse de una correlación RER negativa; los nombres `top`/`bottom` de CAAS ya codifican otro eje de dirección. Un filtro posterior, si se quisiera para presentación, sería una selección adicional, no el procedimiento bilateral de referencia.

### Identificadores y p extremos de RER

En marino observado y nulo: 121 genes, ningún ID duplicado ni prefijo duplicado al quitar sufijos tras punto; ningún P=0 ni 0<P<1e-300 (118 P observados y 118.874 celdas nulas finitas). No hay impacto de esos casos en esta referencia.

Los sintéticos sí muestran diferencias: P=1e-310 produce score 310 nativo frente a 300 con el suelo actual; P=0 produce infinito nativo frente a 300; Rho=0/P=0 se elimina como NaN en `getStat`, mientras el suelo puede convertirlo en cero. El observado se limita en `rer_gene_lists.nf` y `scoring_compute.R`, pero el nulo conserva el P crudo. Ese diagnóstico precedía a la corrección actual: se ha retirado el suelo en los productores observados de listas, SCORING e informe RER, usando la misma fórmula que el nulo. P=0 con Rho distinto de cero produce ±Inf; Rho=0/P=0 produce NaN y se excluye del ranking. Los límites usados solo para dibujar evidencia o presentar q no son parte de este estadístico.

El código nativo de `getStat` tampoco selecciona de manera general el score de mayor magnitud entre nombres con el mismo prefijo: en el sintético positivo conserva el primero aunque sea menor; en el negativo puede conservar ambos. No se copiará este comportamiento sin una decisión. La política aplicada rechaza IDs exactos repetidos (también si uno de sus scores es NA), y en RER rechaza además varios IDs con el mismo prefijo anterior al punto. Así el usuario debe resolver el mapeo ambiguo antes del análisis, sin reproducir la selección inconsistente del paquete. Los sufijos de IDs únicos se conservan. Ver [getstat_crosscheck.md](/home/miguel/.codex/visualizations/2026/10/22/01a16449-db13-7541-adac-2fa17f0242e3/wp1_rer_fix/getstat_crosscheck.md).

## Registro histórico a consultar

- [`WP1_fcs_rer.md`](WP1_fcs_rer.md): plan y adendas de implementación, corrección de NA, auditoría de lateralidad y ajustes observados/nulos. Sus descripciones iniciales de `sig_top`, inferencia de signo, ranking CAAS y ejecución pendiente preceden a decisiones posteriores; no tratarlas como vigente si difieren de esta página.
- [`WP9_estadistico_gen.md`](WP9_estadistico_gen.md): diseño y resultados de A1; la elección de estadístico base y controles pendientes siguen abiertos.

## Revisión del comparador AUC e informes (2026-10-22)

- Se corrigió el texto del informe 12 que todavía describía la permulación del p analítico, y el del informe 15 que describía un segundo umbral y una puerta conjunta. El informe 15 puntúa genes mediante `sig_perm`, grupos y `q.perm`; Redlich se informa aparte.
- Se detectó un caso degenerado real: HPO `HP:0000001` y `HP:0000118` abarcan todo el fondo anotado marino. La función observada nativa usa entonces un fallback fuera del GMT y devuelve `stat=-0.5`, mientras el nulo vectorizado exige al menos tres genes anotados fuera del conjunto. El observado aplica ahora el mismo requisito que el nulo. Es una desviación explícita del fallback del paquete, no una reproducción de ese caso degenerado. También se excluye `TCGA_PanCancer_all_genes` por el requisito del grupo de comparación. Las tablas RER marinas pasan de 136 a 133 conjuntos observados por esta corrección.
- Validación nueva: `tests/test_fcs_auc_contract.R` pasó los casos de p extremos, IDs ambiguos y AUC con infinitos, y 9.576 comparaciones nativas marinas sin diferencias en AUC ni p; comprobó invariancia de bloques 5 frente a los 12 ciclos de la muestra. `test_fcs_rules.R` pasó. `test_fcs_rer_missing.R` pasó contra el marino con la expectativa corregida y guardada de 133 conjuntos. La suite externa `test_fcs_null_statistic.py`, con el nuevo validador incluido en la extracción de funciones, pasó 10/10 desde su ubicación real. El fallo de permisos anterior está resuelto; estos controles se han repetido sobre los archivos guardados. Estos controles no regeneran el nulo ni cierran WP6/WP9.

## Entrega y validación integrada (decisión del usuario, 2026-10-22)

Las correcciones y pruebas dirigidas de esta revisión se entregan en un commit temático de WP1. No se exige una corrida completa por WP. La regeneración del nulo RER con los ajustes armonizados, las nuevas tablas completas y las medidas de escala se validarán conjuntamente con la integración de WP6/WP7. Los defaults y las plantillas pendientes deben revisarse antes de esa aceptación integrada. El commit no acredita esas corridas ni cierra las decisiones de WP9. Las funciones compartidas de observado, contadores por bloque y cierre forman parte del cálculo AUC corregido y quedan disponibles para WP6.
