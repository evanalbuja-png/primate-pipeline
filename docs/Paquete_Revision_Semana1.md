# Paquete de revisión — Semana 1

**Proyecto:** Predicción de susceptibilidad viral en primates neotropicales mediante genómica comparativa de genes de inmunidad innata: estudio de caso del mono araña ecuatoriano (*Ateles fusciceps*)  
**Institución:** INABIO – Laboratorio de Secuenciamiento  
**Responsable:** Daniel Albuja  
**Fecha de cierre Semana 1:** 18 de octubre de 2026  
**Repositorio:** https://github.com/evanalbuja-png/primate-pipeline  
**Etiqueta:** v0.1-semana1  

---

## a) Diagnóstico del sistema (Paso 1)
Ubuntu 24.04.4 LTS (noble)
nproc: 4
Memoria: 9.7 GiB total (≈7.8 GiB available)
Disco en ~: 1007G total, 43G used, 914G available (5%)
conda: 26.3.2  (/home/evan/miniforge3/bin/conda)
mamba: 2.5.0
git: 2.43.0
gh: 2.45.0

**Ruta de trabajo adoptada:** `/mnt/d/Documentos/Proyectos/ateles_immune_genomics`  
(Se eligió esta ruta por limitación de almacenamiento en el sistema de archivos nativo de WSL2.)

---

## b) Entornos conda

### `mamba env list` (entornos del proyecto)
base                 *       /home/evan/miniforge3
ateles_align_filo            /home/evan/miniforge3/envs/ateles_align_filo
ateles_busco                 /home/evan/miniforge3/envs/ateles_busco
ateles_datos                 /home/evan/miniforge3/envs/ateles_datos
ateles_r                     /home/evan/miniforge3/envs/ateles_r
ateles_seleccion             /home/evan/miniforge3/envs/ateles_seleccion


### Archivos en `envs/`

- `ateles_datos.yml` + `ateles_datos.lock.yml`
- `ateles_align_filo.yml` + `ateles_align_filo.lock.yml`
- `ateles_seleccion.yml` + `ateles_seleccion.lock.yml`
- `ateles_r.yml` + `ateles_r.lock.yml`
- `ateles_busco.yml` + `ateles_busco.lock.yml`

Los archivos `.lock.yml` se generaron con `conda env export --no-builds`.

---

## c) Tabla de verificación de herramientas (Paso 4)

| Herramienta       | Entorno            | Versión / Salida                                      | Estado |
|-------------------|--------------------|-------------------------------------------------------|--------|
| datasets          | ateles_datos       | datasets version: 18.38.0                             | OK     |
| blastp            | ateles_datos       | blastp: 2.17.0+                                       | OK     |
| tblastn           | ateles_datos       | tblastn: 2.17.0+                                      | OK     |
| diamond           | ateles_datos       | diamond version 2.2.8                                 | OK     |
| Biopython         | ateles_datos       | 1.88                                                  | OK     |
| gffread           | ateles_datos       | 0.12.9                                                | OK     |
| seqkit            | ateles_datos       | seqkit v2.14.0                                        | OK     |
| mafft             | ateles_align_filo  | v7.526 (2024/Apr/26)                                  | OK     |
| trimal            | ateles_align_filo  | trimAl v1.5.rev1                                      | OK     |
| pal2nal.pl        | ateles_align_filo  | presente                                              | OK     |
| iqtree            | ateles_align_filo  | IQ-TREE version 3.1.4                                 | OK     |
| modeltest-ng      | ateles_align_filo  | ModelTest-NG v0.1.7                                   | OK     |
| mb (MrBayes)      | ateles_align_filo  | MrBayes 3.2.8                                         | OK     |
| orthofinder       | ateles_align_filo  | presente                                              | OK     |
| snakemake         | ateles_align_filo  | 9.27.0                                                | OK     |
| codeml (PAML)     | ateles_seleccion   | presente                                              | OK     |
| HYPHY             | ateles_seleccion   | HYPHY 2.5.103                                         | OK     |
| busco             | ateles_busco       | BUSCO 6.1.0                                           | OK     |
| ggtree            | ateles_r           | 4.0.5                                                 | OK     |
| ComplexHeatmap    | ateles_r           | 2.26.1                                                | OK     |
| ape               | ateles_r           | 5.8.1                                                 | OK     |
| phytools          | ateles_r           | 2.5.2                                                 | OK     |

**Todas las herramientas: OK.**

---

## d) Archivos de Snakemake

### `config/config.yaml`
Contiene: `threads: 4`, lista de 9 especies y 19 genes diana (TLR1–9, TRIM5, APOBEC3G/H, MX1/2, OAS1/2, BST2, ACE2, CD4, DDX58, MAVS).

### `workflow/profiles/default/config.yaml`
cores: 4
use-conda: true
printshellcmds: true
conda-frontend: mamba

### `workflow/Snakefile`
Regla `all` + regla `verificar_entorno` que usa el entorno `ateles_align_filo.yml` y escribe versiones de Snakemake, Python, MAFFT, IQ-TREE y OrthoFinder en `results/verificacion_entorno.txt`.

---

## e) Ejecución de Snakemake

### Dry-run (`snakemake -n`)
DAG construido correctamente. Reglas: `verificar_entorno` + `all`.

### Ejecución real
Primera ejecución falló por `SIGPIPE` (exit 141) causado por `| head` + `pipefail`.  
Se corrigió desactivando `pipefail` en el shell de la regla.  
Segunda ejecución: **éxito**.

### Contenido de `results/verificacion_entorno.txt`
=== Verificación de entorno Snakemake ===
Fecha: Fri Oct  9 21:42:16 -05 2026
Snakemake version: 9.27.0
Python version: Python 3.14.0
MAFFT version: v7.526 (2024/Apr/26)
IQ-TREE version: IQ-TREE version 3.1.4 ...
OrthoFinder: .../bin/orthofinder
Entorno activado correctamente por Snakemake.

---

## f) Estructura del repositorio y `.gitignore`

### `tree -L 3 -a -I '.git'` (resumen)
.
├── .gitignore
├── LICENSE
├── README.md
├── config/ (config.yaml + .gitkeep)
├── data/ (raw/, processed/ — no versionados)
├── docs/ (bitacora.md, verificacion_herramientas.md, logs/, Paquete_Revision_Semana1.md)
├── envs/ (5 YAML + 5 lock.yml)
├── logs/
├── results/ (verificacion_entorno.txt)
└── workflow/ (Snakefile, profiles/, rules/, scripts/)

`.snakemake/` está presente localmente pero **excluido** del repositorio por `.gitignore`.

### `.gitignore`
Excluye: `data/`, `results/` (con excepciones para `.txt/.md/.tsv/.csv/.html`), `.snakemake/`, `__pycache__/`, `*.log`, `*.tmp`, etc.  
Incluye explícitamente los `.lock.yml` y `.yml` de `envs/`.

---

## g) Git y GitHub
git log --oneline --decorate:
bad49df (HEAD -> main, tag: v0.1-semana1, origin/main) docs: README, LICENSE MIT, bitácora y verificación de herramientas (Semana 1)
1f16ca6 feat: Snakefile base con regla verificar_entorno y perfil de ejecución
f47599d chore: entornos conda (YAML + locks) para datos, alineamiento, selección, R y BUSCO
79f1e31 chore: estructura inicial del proyecto y .gitignore
git remote -v:
origin  https://github.com/evanalbuja-png/primate-pipeline.git (fetch)
origin  https://github.com/evanalbuja-png/primate-pipeline.git (push)
git status: On branch main, up to date with origin/main, nothing to commit, working tree clean

- **URL del repositorio:** https://github.com/evanalbuja-png/primate-pipeline  
- **Etiqueta:** https://github.com/evanalbuja-png/primate-pipeline/releases/tag/v0.1-semana1  

---

## h) Errores, advertencias y decisiones tomadas

| Ítem | Descripción | Decisión / Resolución |
|------|-------------|------------------------|
| Ruta de trabajo | Limitación de almacenamiento en filesystem nativo de WSL | Se adoptó `/mnt/d/Documentos/Proyectos/ateles_immune_genomics` |
| IQ-TREE versión | Propuesta mencionaba 2.x; bioconda ofrece 3.1.4 | Se aceptó 3.1.4 (versión actual y funcional) |
| BUSCO versión | Propuesta mencionaba v5; se instaló 6.1.0 | Se aceptó 6.1.0 (versión actual) |
| Snakemake `conda-frontend: mamba` | Advertencia de deprecación en Snakemake 9.x | Se ignoró (no afecta funcionalidad; Snakemake usa libmamba internamente) |
| Exit status 141 (SIGPIPE) | `mafft --version \| head` + `set -o pipefail` | Se corrigió con `set +o pipefail` y redirección directa |
| Archivos basura | `busco_*.log`, `rst`, `rst1`, `rub` | Eliminados antes del primer commit |
| LICENSE | Se creó inicialmente como `.txt` | Renombrado a `LICENSE` (sin extensión) |

---

**Estado de la Semana 1:** COMPLETADA  
Listo para revisión externa y arranque de Semana 2 (revisión bibliográfica + descarga/evaluación del genoma de *Ateles fusciceps*).
