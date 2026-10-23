# POSENRICH modular — contrato de entrada, inventario de productores y especificación del estimando

Documento de trabajo WP10 (P0). Fecha: 23/10/2026. Estado del árbol: rama `remove_caap`, HEAD `6e13ba23`, con cambios concurrentes sin commit (migración percentiles → evidencia nominal/BH y anotaciones posicionales). Este contrato fija las entradas que el módulo consume; no modifica productores ajenos.

## 1. Inventario de productores vigentes (auditado 23/10/2026)

| Artefacto | Productor | Contrato vigente |
|---|---|---|
| Universo de genes | `cleaned_background` (background.txt) | Un gen por línea. |
| Scores de gen | `observed_core_scores.py` → `gene_scores.tsv`, `fcs_stats.tsv` | Una fila por gen; `gene_caas_score{,_top,_bottom}`; `size_adj_max` por defecto. |
| Evidencia de gen | `SCORING_EVIDENCE` (`scoring_evidence.R`) → `gene_evidence.tsv`, `gene_evidence_metadata.json` | Una fila por (gene, ranking ∈ {global, top, bottom}): statistic, gene_stat, p.gen, p.gen.adj, nominal, adjusted, N_roster, N_valid, family_id (`CAAS:<rank>`), contract_version. BH sobre el universo completo, nunca por batch. |
| Scores/evidencia de posición | `SCORING_COMPUTE` (`scoring_compute.R` + `null_position_p.R`) → `position_scores.tsv`; enriquecido por `scoring_evidence.R` | Una fila por (Gene, Position, side) **con CAAS_score > 0**. Columnas: Gene, Position, biochem_groupings, is_biochem, descriptores, CAAS_score, side, caas, ancestral_aa, derived_aa, p.emp, p.adj_bh, p.emp_fact, p.adj_bh_fact, nd, null_unseen, caas_nominal, caas_adjusted. **No contiene todas las posiciones evaluadas.** |
| Posiciones detectadas (evaluadas por CT) | `filtered_discovery.tsv` (CT_DISAMBIGUATION) | Una fila por (Gene, Position, side) detectada; `asr_path_score` puede ser ≤ 0 (evaluada sin detección positiva). |
| Universo posicional del null | `null_universe_summary.tsv` (métricas) + contadores del null | N_cycles, null_positive_positions_in_universe (|U|), observed_positive_positions, observed_unseen_positions, universe_size = |U| + observed unseen. Familia BH posicional: U + posiciones observadas con p. |
| Coordenadas | `caas_map_dir` (`<GENE>*.map.tsv`) | Columnas `ori_codon_col`, `status`, `prot_ali_col`, `hg38_aa_pos`, `hg38_nt_coord`; solo filas `selected`; Position = prot_ali_col − 1 (CAAS 0-based). Un directorio con dos MAP del mismo gen se rechaza. |
| Colecciones | `POSITION_ANNOTATIONS` (`build_position_gmt.py`) | GMT: pfam_domains, pfam_clans, genomic_locations (ventanas 1 Mb), custom_features; capas globales: characterization_layers.tsv; tabla plana: position_characterization.tsv. |
| FADE por sitio | `fade_sites_{top,bottom}.csv` (FADE_TABLES) | gene, position (CAAS), max_bf, target_aa; selección upstream, sin umbral adicional. |
| VEP/COSMIC/PrimateAI | módulo VEP → `vep/ensembl_vep_mapped.tsv` etc. | Anotaciones descriptivas; join por gen/posición/dirección; MANE Select preferido. |
| PP2A (piloto) | `bin/prepare_pp2a.py` → `pp2a_prepared/` | source_records.tsv.gz, site_candidates.tsv, unresolved.tsv, manifest.json. **No son GMT listos para CAAS**: falta validar accession/isoforma/péptido contra la referencia productiva y su MAP.** |

Reglas transversales vigentes (migración a p.emp, `docs/plan_migracion_cortes_a_p_emp.md`):

- p.emp posicional = máximo entre lados por (Gene, Position); las filas top/bottom comparten ese p. No son dos pruebas direccionales independientes. Score direccional (CAAS_score por side) y p compartido se mantienen separados.
- Sin percentiles: no slice_25/10/5/1, no flags pct_*, no score por percentiles. Cortes: p.emp < 0,05 nominal; p.adj_bh < 0,05 resaltado.
- FCS CAAS conserva el ranking continuo −log10(p.gen) sobre todos los genes del universo.

## 2. Especificación corta del estimando (antes de implementar inferencia)

- **Unidad:** una posición biológica `(Gene, Position)` con Position 0-based en el alineamiento recortado de CAAS. Una unidad por posición en cada ranking; no se multiplica por hipótesis, transcritos, lados ni filas de evidencia.
- **Dirección:** `top` / `bottom` (score CAAS direccional); `global` = pool de lados (deduplicado por posición). El p.emp es compartido entre lados (máximo); no existe p direccional independiente.
- **Score de ranking (propuesta a validar):** CAAS_score continuo por dirección sobre posiciones evaluables; **no se seleccionan hits antes de calcular AUC**. −log10(p.emp) es una sensibilidad explícita del ranking compartido entre lados, no un score direccional ni una sustitución silenciosa.
- **Estimando P3 (localización dentro de proteínas, prioritario):** AUC de los miembros del término frente a los no miembros **elegibles de la misma proteína**, con medio crédito a empates. Solo informan genes con ambos grupos. Agregación principal propuesta: media de AUC por gen con igual peso (la longitud/carga de una proteína no decide por sí sola); ponderación por pares como sensibilidad. Se fija la principal antes de datos confirmatorios.
- **Control P3 (propuesta, pendiente de puerta de validación):** redistribución condicionada dentro de proteína — permutar etiquetas de membresía dentro de cada proteína, preservando proteína, perfil de scores y cobertura. No se baraja uniformemente entre proteínas ni se aplican desplazamientos circulares. p Monte Carlo unilateral para mayor AUC: excedencias inclusivas, (k+1)/(B_válidas+1), semilla fija. **Ningún control queda validado por este documento.**
- **Estimando P4 (asociación global, vía separada):** AUC del conjunto frente al fondo global compatible, con composición por genes y sensibilidad a retirar cada proteína. Sin calibración global reproducible, solo efecto descriptivo con estado explícito.
- **Estimando P5 (proximidad):** distancia en residuos de referencia al sitio/región más próximo de la misma proteína (0 dentro de región); genómica (ventanas 1 Mb / distancia nt) en contrato separado. La cercanía proteica no se mide en columnas MSA con gaps.

## 3. Estados y denominadores (registro común)

Cada posición del registro clasifica, sin reclasificar a nadie como negativo por ausencia:

- `detected_positive`: en `position_scores.tsv` (CAAS_score > 0).
- `evaluated_no_detection`: en `filtered_discovery.tsv` pero sin score positivo (evaluada sin detección).
- `map_only_not_in_discovery`: columna del MAP sin fila de descubrimiento (evaluación desconocida; no es negativo ni cero).
- `no_map`: en scores/discovery pero sin MAP (coordenadas no proyectables; no evaluable para proximidad).
- `no_annotation`: sin membresía en ninguna colección (ausencia de anotación, no ausencia biológica).

Denominadores explícitos en `coverage.tsv`: posiciones anotables (MAP), posiciones evaluadas (discovery), posiciones con score positivo, posiciones del universo del null (|U| + unseen), genes del universo. Nunca se cambia de denominador al filtrar hits.

## 4. Fixtures

Los fixtures pequeños reproducibles viven en los tests (`tests/test_posenrich_registry.py` y sucesivos): MAP con primer residuo ≠ 1, columna eliminada y gap humano; Pfam con huecos/X; UCR en columna cruda; FUBAR por residuo de referencia; FADE directo; custom validado; colección retirada (eggNOG/COSMIC/PrimateAI) ignorada con aviso. No se usan datos reales de corrida para las pruebas dirigidas.

## 5. Línea base de pruebas (23/10/2026, rama con cambios concurrentes)

PASS: `test_position_collections.py` + `test_runner_reuse.py` (10), `test_position_report.R`, `test_scoring_evidence.R`, `test_fcs_rules.R`, `test_fcs_rer_missing.R`, `test_wp6_blocks.R`, `test_wp6_fcs_cli.R`, `test_wp6_rer_cli.R`, `test_wp6_rer_batches.R`, `test_wp9.R`, pytest WP4/WP6/WP7/WP8/WP9/empirical (46), `test_wp5_null_blocks_cli.py` (19).

Pendiente (migración concurrente, no POSENRICH): `test_fcs_auc_contract.R` falla en una aserción textual obsoleta sobre el informe 5 migrado (busca `rer_stat      = sign(Rho) * -log10(P),`; el contrato vigente es `rer_df$rer_stat <- sign(rer_df$Rho) * -log10(rer_df$P)` sin percentiles). `test_rer_permutations.R` no es ejecutable en este momento (el archivo aparece en el listado pero no se puede abrir; el árbol está en transición por otra conversación). `test_wp6_rer_real.R` requiere datos reales (argumento de corrida).
