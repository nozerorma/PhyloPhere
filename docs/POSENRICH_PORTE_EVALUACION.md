# POSENRICH sin Lachenbruch: evaluación del porte

Estado: documento de trabajo. POSENRICH queda desactivado, con aviso, mientras no se decida lo que sigue. No hay cambios de código aplicados.
Contexto: `docs/PLAN_NULO_10K_PERMULACIONES.md` (decisión D1: Lachenbruch sale del FCS y de POSENRICH) y plan de implementación WP2.

## 1. Qué hace hoy (verificado en el código)

- **Términos:** conjuntos de posiciones (`gen:posición`), construidos por `build_position_gmt.py` a partir de las bases de datos y de las capas de caracterización.
- **Fondo:** posiciones probadas por CAAStools (`background.output`) en los genes del universo; con nulo, restringido a las posiciones con score > 0 en algún ciclo nulo más las observadas con score > 0.
- **Estadístico** (`posenrich_enrich.py`, `lachenbruch_columns`): prueba de Lachenbruch de dos partes, la misma definición que `fcs_lach_*` de `fcs_enrich.R`: parte 1 (prevalencia), cola superior de la hipergeométrica de las posiciones puntuadas del término entre las posiciones puntuadas de la columna; parte 2 (magnitud), suma de rangos de los scores de las posiciones puntuadas del término frente a las demás, con aproximación normal, corrección de continuidad y término de empates (se exige al menos 2 posiciones puntuadas dentro y fuera). Cada p se acota en 1e-15, se convierte a χ² con 1 g.l. y se suman (χ² con 2 g.l.).
- **Nulo:** cada ciclo de `perm_pos_cycle_caas.tsv.gz` pasa por la misma función como una columna más (`caas_null_matrix`: matriz dispersa posiciones de fondo × ciclos), y `lach_p.perm = (1 + nº de ciclos con χ² ≥ observado) / (N + 1)`.
- **Criterio:** BH sobre `lach_pval` dentro de cada (ranking, base de datos) y `lach_p.perm < fcs_pperm_thr`, con `fcs_fdr_lachenbruch` (cae a `fcs_fdr` si no existe).
- **Entrada del nulo:** por dirección (global, top, bottom), un score por (posición, ciclo); la dirección global toma el máximo sobre lados. `posenrich_prep_caas_null.py` lo reduce una vez para todos los lotes.

## 2. Qué cambia con el plan

- El fichero `perm_pos_cycle_caas.tsv.gz` (80,9 millones de filas a 1000 ciclos en cáncer; ~800 millones a 10.000) desaparece del nulo (WP5): SCORING pasa a contadores por bloque y la acumulación a un indicador por posición.
- Lachenbruch sale del FCS. POSENRICH tiene su propia copia del test, de modo que podría conservarlo, pero ya no tendría el nulo por posición que lo alimenta.
- Mientras tanto: módulo omitido con aviso explícito (informe 14 incluido).

## 3. Opciones

| Opción | Qué implica | Coste y riesgos |
|---|---|---|
| A. Reutilizar la maquinaria del FCS a nivel de posición | El ranking es el score de posición por dirección (cero donde no hay detección); estadístico de Wilcoxon con término de empates sobre el fondo (los ceros quedan empatados en la cola inferior, de modo que la prevalencia entra de forma implícita); p de permulación por bloques de ciclos con contadores aditivos y cierre con BH, como en la fase 4 | Reaprovecha el código y la validación del FCS. Se pierde la separación explícita entre prevalencia y magnitud. Exige un almacén de scores por (posición, ciclo) por bloque, con solo scores positivos (~2,3 · 10⁷ filas a 1000 ciclos, ~2,3 · 10⁸ a 10.000, en binario) |
| B. Conservar Lachenbruch solo en POSENRICH | Lee el almacén binario por bloques y acumula χ² por columna | Mantiene dos pruebas distintas en el repositorio y el almacén por posición lineal en N, que es lo que el plan quiere evitar; la calibración de Lachenbruch dependía del nivel en las corridas estrictas (ver el documento principal, sección 2.2) |
| C. No portar | El análisis por posición queda cubierto por el FCS por gen y por `p.emp` de posición | Menor coste; se pierde la lectura por posición de los conjuntos |

## 4. Criterios para decidir si compensa

1. **Aporte frente al FCS por gen.** En los resultados de cáncer, POSENRICH da unos pocos términos frente a cientos del FCS (nota del 20/10/2026 en la memoria de proyecto), y en la corrida del 19/10 solo se probaron 5 de 12 fuentes de posiciones (anomalía abierta, causa sin investigar). Antes de decidir hay que cerrar esa anomalía.
2. **Potencia a 10.000 ciclos.** Con el suelo de p en 1e-4, la limitación de resolución desaparece en bases grandes; conviene saber si POSENRICH cambia con más ciclos.
3. **Coste.** El almacén por posición es la única salida del nulo que seguiría siendo lineal en N.
4. **Validación.** Cualquier porte necesita igualdad exacta del recuento de ciclos por bloque frente a la versión actual en toy y a escala real (las comprobaciones del plan, sección 7).

## 5. Sin verificar

- Cuántos términos distintos cambian de significación entre Lachenbruch y la opción A con los mismos datos (no se ha calculado).
- Si el conteo de posiciones puntuadas del fondo coincide entre la opción A y la actual al sustituir el nulo por el almacén binario.
- El motivo de las 5 de 12 fuentes probadas.
