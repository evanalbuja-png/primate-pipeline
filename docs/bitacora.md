# Bitácora del proyecto — Ateles Immune Genomics

**Responsable:** Daniel Albuja  
**Institución:** INABIO – Laboratorio de Secuenciamiento  
**Repositorio:** https://github.com/evanalbuja-png/primate-pipeline  

---

## Semana 1 (9–21 de octubre de 2026)

### Actividades realizadas
- Diagnóstico del sistema (Ubuntu 24.04.4 LTS).
- Creación de estructura de directorios y `.gitignore`.
- Cinco entornos conda aislados (`ateles_datos`, `ateles_align_filo`, `ateles_seleccion`, `ateles_r`, `ateles_busco`) + locks sin `prefix:`.
- Verificación de herramientas clave (todas OK).
- Snakemake 9.27.0 con regla de prueba `verificar_entorno` (conda apunta a `.lock.yml`).
- Documentación: README, LICENSE MIT (Daniel Albuja), bitácora, verificación de herramientas.
- Traslado de la distro WSL Ubuntu-24.04 a `D:\WSL\Ubuntu-24.04` (export/import; `--move` falló con E_ACCESSDENIED / ERROR_FILE_NOT_FOUND).
- Reconfiguración de `.wslconfig`: memory=11GB, processors=6, swap=16GB (swapfile en D:).
- Repo activo en filesystem nativo: `~/ateles_immune_genomics`.
- BUSCO recreado en 5.8.3 (antes 6.1.0).

### Versiones relevantes (cierre Semana 1)
- conda 26.3.2 / mamba 2.5.0
- Snakemake 9.27.0
- IQ-TREE 3.1.4
- OrthoFinder 3.1.5
- MAFFT 7.526
- BUSCO **5.8.3**
- PAML 4.10.10 / HYPHY 2.5.103
- pal2nal.pl v14
- R: ggtree 4.0.5, ComplexHeatmap 2.26.1, ape 5.8.1, phytools 2.5.2

### Desviaciones respecto a la propuesta

| Fecha | Desviación | Detalle / justificación |
|-------|------------|-------------------------|
| 2026-10-09 | IQ-TREE 3.1.4 en lugar de 2.x | Versión actual en bioconda. Compatible con `-m TEST -B 1000 --alrt 1000`. Citar IQ-TREE 3 además de Minh et al. (2020). Prueba de humo OK (soportes SH-aLRT/UFBoot). |
| 2026-10-20 | BUSCO 5.8.3 (recreado; antes 6.1.0) | Propuesta: BUSCO 5.x + linaje `primates_odb10`. Instalado 5.8.3. Dataset disponible: `primates_odb12.2` (no aparece `odb10` en `--list-datasets`). Usar `primates_odb12.2` en análisis. |
| 2026-10-20 | Hardware real vs. propuesto | Propuesta: ≥16 CPU / ≥64 GB RAM. Real: 12 hilos lógicos del host; WSL limitado a **6 processors** y **11 GB RAM** (+ 16 GB swap) vía `.wslconfig`. Suficiente para Semana 1–2; reevaluar para OrthoFinder/PAML a escala. |
| 2026-10-20 | Ruta de trabajo | Inicialmente `/mnt/d/...` por falta de espacio en C:. Tras mover la distro a D:, el repo vive en **`~/ateles_immune_genomics`** (filesystem nativo). Copia antigua: `/mnt/d/Documentos/Proyectos/ateles_immune_genomics_respaldo`. |
| 2026-10-20 | `wsl --manage --move` no disponible | Falló con E_ACCESSDENIED y ERROR_FILE_NOT_FOUND. Se usó export (43.68 GB, 11.1 min) → unregister → import (~19 min). BasePath final: `D:\WSL\Ubuntu-24.04`. |
| 2026-10-09 / 20 | Snakemake `set +o pipefail` | Solo en la regla de prueba `verificar_entorno` (evita SIGPIPE/exit 141). Reglas de análisis futuras mantienen el modo estricto por defecto. |
| 2026-10-20 | `conda-frontend: conda` | Sustituye `mamba` en el perfil para eliminar la advertencia de deprecación en Snakemake 9.x. |

### Próximos pasos (Semana 2)
- Revisión bibliográfica sistematizada.
- Descarga y evaluación de calidad del genoma de *Ateles fusciceps* (BUSCO 5.8.3, linaje `primates_odb12.2`).

### Actualización 2026-10-21 (inicio Semana 2 — Paso 0)

Desviaciones adicionales registradas:

| Fecha      | Desviación                                      | Detalle / justificación |
|------------|-------------------------------------------------|-------------------------|
| 2026-10-21 | OrthoFinder 3.1.5 en lugar de v2                | Versión actual en bioconda (aceptada). Citar OrthoFinder 3 además de Emms & Kelly (2019) en el manuscrito; verificar la referencia oficial de la v3. Comando verificado: `-S diamond` es el default. |
| 2026-10-21 | Binario IQ-TREE = `iqtree3`                     | La instalación de IQ-TREE 3.1.4 expone el ejecutable `iqtree3` (no `iqtree2`). Las reglas futuras de Snakemake deben invocar `iqtree3`. |
| 2026-10-21 | BUSCO: pendiente verificación `primates_odb10`  | Propuesta original: linaje `primates_odb10`. Con BUSCO 5.8.3 el dataset listado es `primates_odb12.2`. **No se acepta aún odb12.2.** Pendiente Semana 3: verificar si se puede descargar/usar `primates_odb10` o documentar la justificación formal del cambio a odb12.2 ante el revisor. |


## Semana 2: decisiones y desviaciones (2026-10-27)
- Canales conda: el archivo ~/.condarc tenia defaults primero (canal comercial de Anaconda) con prioridad estricta, y el lock de ateles_datos no se podia reconstruir (ca-certificates 2026.10.20 solo existe en conda-forge). Se corrigio ~/.condarc (conda-forge, bioconda) y los envs/*.yml pasan a nodefaults.
- envs/ateles_pubmed.yml (python 3.15.2, biopython 1.88, pyyaml 6.0.3): entorno minimo con versiones fijas para las reglas de PubMed.
- PubMed se consulta con Bio.Entrez y no con Entrez Direct: asi se guardan QueryTranslation, ErrorList y WarningList, y la cadena de consulta no pasa por el shell.
- Fase A: la pasada del 22-oct se descarto porque no se pudo verificar que cadenas uso; la oficial es la del 27-oct, con cadenas validadas. Los conteos coinciden salvo +1 en las dos de ACE2.
- Fase B: 15 consultas, 3489 PMIDs sumados, 2909 unicos.
- Pendiente de verificar contra PubMed: PMID 32768453 (TLR7/TLR8 en Alouatta).

## Semana 2: correccion de entornos (2026-10-27)
- El entorno ateles_datos contenia 3 paquetes del canal comercial pkgs/main y su lock no se podia reconstruir sin defaults. Se recreo desde el yml con conda-forge y bioconda, se exporto un lock nuevo sin prefix, y se reconstruyo el entorno definitivo desde ese lock como prueba de reproducibilidad.
- Se marcaron como vigentes (snakemake --touch) las salidas de PubMed tras el cambio de definicion de entornos, para no repetir consultas a NCBI sin necesidad. verificar_entorno se re-ejecuto.
- El commit 452e08b de esta semana incluia un lock de ateles_datos que no se reconstruia; se corrige en el commit siguiente (sin reescribir historia).
