# PhyloPhere — Claude Code instructions

## Cluster execution policy

Nunca ejecutar cómputo directamente en el nodo de entrada (login node) de un
cluster (Marvin2, correfoc). Formas admitidas de correr algo:

1. Sugerir el comando al usuario para que lo ejecute él mismo.
2. Lanzar una sesión interactiva (`salloc` / `srun --pty ...`) y correr ahí.
3. Enviar el trabajo vía `sbatch`.

Comandos puramente informativos y de bajo coste (squeue, sacct, cat/tail de
logs, ls, df, etc.) sí pueden ejecutarse directamente sobre el login node.

## Directorios temporales en cluster

Nunca usar `/tmp` en los clusters para runs temporales (archivos intermedios,
work dirs de prueba, etc.). Usar en su lugar un directorio `.tmp` dentro del
scratch de cada cluster:

- Marvin2: `/scratch/lab_anavarro/mramon/.tmp`
- correfoc: `~/scratch/0.Phylophere/.tmp`

Crear el subdirectorio si no existe antes de usarlo.

## Gestión de versiones en el cluster

Nunca modificar archivos de código/config directamente en el cluster (Marvin2,
correfoc) por ssh. La gestión de versiones (commits, push, pull, merges) la
hace el usuario manualmente en su propio flujo de git. Cualquier cambio de
código debe hacerse en el repo local y llegar al cluster a través de su propio
proceso de sincronización/git, no editando archivos remotos directamente.
