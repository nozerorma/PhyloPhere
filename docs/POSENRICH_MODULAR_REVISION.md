# WP10 — correcciones tras revisión de POSENRICH modular

Fecha: 23/10/2026. Base observada: `bb89c6de`, rama `remove_caap`, con cambios
concurrentes preservados. Este registro prevalece sobre el estado de cierre P0–P6
indicado en la entrega inicial.

## Contrato científico corregido

La revisión inicial interpretó un `asr_path_score=0` como un score CAAS cero que
había que incorporar al AUC. Al leer `CT_DISAMBIGUATION/local/src/core/scores.py`
se comprobó que `position_score()` solo devuelve score para valores positivos;
para cero/negativo/faltante devuelve None. Por tanto el ejemplo AUC=0,5 de la
revisión describe otro ranking hipotético, no un oráculo del ranking canónico.
No se altera este productor, FCS, acumulación, gene_stat ni el nulo canónico.

La carencia real es que el módulo describía el ranking de detecciones positivas
como si incluyera todo el universo evaluable y publicaba inferencia con un control
no validado. Las salidas ahora explicitan el condicionamiento y retienen la
inferencia: `status=descriptive_pending_control`, `p_mc=NA`, `q=NA`, `B_valid=0`.
Los parámetros de permutación quedan reservados como solicitud; Nextflow no
lanza permutaciones de localización. El prototipo puro `run_term` solo se mantiene
para pruebas y diseño de controles. No se habilita una bandera para eludir la puerta.

## Cambios implementados

- Registro: conserva filas de discovery sin score, valor crudo, `score_status` y
  estado de evaluación; no imputa score cero ni p=1. Rechaza duplicados por lado,
  p/q posicionales inconsistentes entre lados y desacuerdos numéricos positivos
  entre discovery y la tabla canónica. Propaga el p compartido a lados adicionales.
- AUC: publica vistas `global/top/bottom` con `rank_direction`, sin filtrar por
  significación; sigue condicionado a scores canónicos disponibles. El fondo de
  FUBAR se restringe a `fubar_tested`; otras coberturas experimentales se declaran
  desconocidas. Colecciones incompatibles con asociación global quedan fuera de
  ese componente. No se fabrica un p global ni se reutiliza el de localización.
- Prototipo: semilla por identidad estable del término/ranking/dirección; independiente
  del inventario y del orden de otros términos. Esto no certifica partición futura
  de ciclos ni validación del control.
- Proximidad: score global por posición reproducible mediante el máximo entre
  lados; población positiva declarada; sin p de distancias.
- Propiedad de tablas: con módulo encendido, informe 15 consume la caracterización
  modular; con módulo apagado se conserva POSITION_ANNOTATIONS. No corren ambos
  productores para el mismo informe. GUI prefiere la tabla modular publicada.
- Informes: rutas correctas del manifiesto de caracterización y del resumen de
  universo; parámetros efectivos leídos desde manifiestos publicados, también al
  regenerar desde GUI; comprobación de ranking coherente entre componentes.
- Nextflow: huellas SHA-256 de las dependencias Python/Rmd en cada tarea modular
  para que `-resume` invalide resultados de implementaciones anteriores.
- Ayuda/GUI/config/README y estado WP10: explican el carácter descriptivo y la
  retención de inferencia. No se añaden parámetros productivos nuevos.

## Pruebas y límites

Regresiones nuevas: discovery con cero y NA diferenciados; p compartido
inconsistente; direcciones con efectos opuestos; inferencia siempre retenida en
CLI; restricción FUBAR al fondo testado; flujo aleatorio independiente por término;
preferencia de la tabla modular en regeneración; ausencia de defaults inventados;
HTML con estado y universo correctos. Se conservan los oráculos numéricos del
prototipo (AUC, empates y control iid) sin convertirlos en validación científica.

El smoke incluye el módulo aislado y el workflow ENRICHMENT completo con módulo
encendido/apagado, verifica el enlace real de la tabla hacia el informe 15 y las
tres direcciones de localización. Los primeros intentos del harness integrado
fallaron por defaults de configuración que activaban FCS y por sentinelas con
nombres repetidos; el harness corregido los desactiva explícitamente y usa nombres
únicos. No fueron fallos del motor POSENRICH.

Las evidencias finales y comandos se registran en la entrega de correcciones y
bajo `PhyloPhere_validation/.tmp/posenrich_modular/review/`. No se lanzan trabajos
a escala ni se acepta PP2A sin referencias productivas.

## Puertas que siguen pendientes

1. Universo completo evaluable y ranking para no detecciones: requiere un contrato
   científico explícito; no se obtiene rellenando el MAP o discovery con ceros.
2. Control condicionado: simulaciones independientes con detectabilidad, sesgos de
   anotación, conservación y clusters; aún no validado. El ensayo iid histórico no
   supera esta puerta y aumentar B no la resuelve.
3. Familia BH y tabla inferencial definitiva, validación por lotes/ciclos y aceptación
   integrada a escala después de superar las puertas anteriores.
4. PP2A: correspondencia accession/isoforma/péptido con referencia/MAP productivos.

Estado: componentes descriptivos corregidos y sometidos a pruebas dirigidas;
P3 inferencial y P6 inferencial pendientes. No se declara P0–P6 científicamente cerrado.


## Resultado final de esta corrección

- `python -m pytest tests/ -q`: **111 passed**, 190,34 s; incluye los tres
  smokes Nextflow (módulo aislado y ENRICHMENT con módulo encendido/apagado).
- Seis controles R: **PASS** en `test_scoring_evidence.R`, `test_fcs_rules.R`,
  `test_fcs_auc_contract.R`, `test_wp6_blocks.R`, `test_posenrich_modular_report.R`
  y `test_position_report.R`. Este último conserva la preferencia MANE, los joins
  por dirección y las anotaciones Pfam/FADE del informe 15.
- `git diff --check`: sin errores de formato.
- Verificación contra el diff inicial: once archivos ajenos no editados conservan
  exactamente sus diffs. En `gui/generation/report_registry.py`, compartido,
  se conserva el cambio previo de metadatos FCS y se modifican solo las secciones
  modulares y la preferencia de caracterización del informe 15.

Comando Python usado (entorno phylophere en PATH y TMPDIR bajo la carpeta autorizada):

```
python -m pytest tests/ -q --basetemp=/home/miguel/IBE-UPF/PhD/PhyloPhere_validation/.tmp/posenrich_modular/review/fix_full
```

Para los renders finales:

```
Rscript tests/test_posenrich_modular_report.R <repo> <scratch>/fix_render_final
Rscript tests/test_position_report.R <repo> <scratch>/fix_comparison_final
```

`<scratch>` es `PhyloPhere_validation/.tmp/posenrich_modular/review/`.
Los controles R restantes tienen sus logs bajo `fix_r/`. El aviso de Pandoc
sobre `--mathjax` es de deprecación; los renders terminan correctamente.

No se ha hecho commit, no se han lanzado trabajos grandes y no se declara aceptada
la inferencia científica ni el piloto PP2A. Los cambios están en el árbol de trabajo.
El registro canónico se guarda fuera de `docs/fases/`, que el repositorio ignora,
para que pueda versionarse con el código.
