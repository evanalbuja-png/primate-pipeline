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
