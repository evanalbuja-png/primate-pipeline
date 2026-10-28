# Paquete de revisión Semana 1 — cierre

**Proyecto:** Predicción de susceptibilidad viral en primates neotropicales (Ateles fusciceps)  
**Responsable:** Daniel Albuja  
**Repo:** https://github.com/evanalbuja-png/primate-pipeline  
**Etiquetas:** v0.1-semana1 · **v0.1.1-semana1**  
**Fecha de cierre:** 21 de octubre de 2026  

---

## A-bis — Evidencia pendiente (literal)

### OrthoFinder

```
OrthoFinder: v3.1.5
(conda list: orthofinder 3.1.5 hdfd78af_0 bioconda)
```

### config/config.yaml — species y genes

```
species:
  - Ateles_fusciceps
  - Callithrix_jacchus
  - Saimiri_boliviensis
  - Aotus_nancymaae
  - Cebus_imitator
  - Sapajus_apella
  - Macaca_mulatta
  - Pan_troglodytes
  - Homo_sapiens

genes:
  - TLR1
  - TLR2
  - TLR3
  - TLR4
  - TLR7
  - TLR8
  - TLR9
  - TRIM5
  - APOBEC3G
  - APOBEC3H
  - MX1
  - MX2
  - OAS1
  - OAS2
  - BST2
  - ACE2
  - CD4
  - DDX58
  - MAVS
```

- Total TLR: **7**
- Total genes: **19**

### Prueba de humo IQ-TREE 3

- Comando:

```
iqtree3 -s humo.fa -m TEST -B 1000 --alrt 1000 -T AUTO --prefix humo_iqtree
```

- Wall-clock: 8.612 sec
- `.treefile` con soportes SH-aLRT/UFBoot (ej. `83.4/65`, `65.2/42`)
- Log final: Analysis results written; Date and Time: Mon Oct 19 19:00:01 2026

### Identidad

- LICENSE / git: Daniel Albuja <evanalbuja@gmail.com> (sin cambios)

---

## B — Movimiento de la distro y reconfiguración (literal)

### B.1 git status (antes de mover)

```
On branch main
Your branch is up to date with 'origin/main'.
nothing to commit, working tree clean
a9d5515 (HEAD -> main, origin/main) docs: Paquete de revisión Semana 1 para auditor externo
```

### B.2 BasePath original

```
Ubuntu-24.04     C:\\Users\\Evan\\AppData\\Local\\wsl\\{aca7fec7-f300-448a-be83-939b57c2cd82}
```

### B.3 Hilos lógicos

```
12
```

### B.4 Export (respaldo)

```
Inicio export: 10/20/2026 22:23:11
Exportación en curso. (44725 MB)
La operación se completó correctamente.
Fin export: 10/20/2026 22:34:17
Duración: 11.1 minutos
SizeGB: 43,68
D:\\WSL\\backup\\ubuntu-24.04-antes-de-mover.tar
```

### B.5 .wslconfig

```
[wsl2]
memory=11GB
processors=6
swap=16GB
swapFile=D:\\WSL\\swap.vhdx
```

### B.6 wsl --manage --move

```
E_ACCESSDENIED (intento 1)
ERROR_FILE_NOT_FOUND (intentos 2–3, con y sin Admin)
→ Se abandonó --move; se usó export/import.
```

### B.6-alt Unregister + import

```
Unregister: OK
Inicio import: 10/20/2026 22:52:40
Fin import: 10/20/2026 23:11:29
BasePath final: D:\\WSL\\Ubuntu-24.04
```

### B.7 Verificación post-import

```
whoami: evan
uid=1000(evan)
Mem: 10Gi total / Swap: 16Gi
nproc: 6
df -h ~: /dev/sdd 1007G, 912G avail
mamba env list: ateles_* presentes
ext4.vhdx en /mnt/d/WSL/Ubuntu-24.04/ (~45.7 GB)
```

### B.8 Clone nativo

```
pwd: /home/evan/ateles_immune_genomics
tree: envs/, workflow/, config/, README, LICENSE presentes
git ls-remote origin HEAD: a9d5515...
```

---

## C — Correcciones al repositorio (literal)

### C.1 prefix: eliminado de los 5 lock.yml

```
Antes: 5 líneas prefix: (una por lock)
Después: OK: no queda ninguna línea prefix:
conteo por archivo: 0
```

### C.2 Snakefile

```
conda:
  "../envs/ateles_align_filo.lock.yml"
set +o pipefail  → solo línea 17 (regla verificar_entorno)
```

### C.3 Perfil

```
conda-frontend: conda
```

(Advertencia de mamba ausente en la ejecución forceall.)

### C.4 snakemake --forceall

```
Creating conda environment envs/ateles_align_filo.lock.yml...
2 of 2 steps (100%) done
Elapsed time: 0:01:45

results/verificacion_entorno.txt:
Snakemake 9.27.0
Python 3.14.0
MAFFT v7.526
IQ-TREE 3.1.4
OrthoFinder: .../bin/orthofinder
```

### C.5 BUSCO 5.8.3

```
BUSCO 5.8.3
primates_odb12.2 [13188]
lock sin prefix:; 258 líneas
```

### C.6 Bitácora

Sección "Desviaciones respecto a la propuesta" con fechas 2026-10-09 y 2026-10-20 (IQ-TREE 3, BUSCO 5.8.3 + odb12.2, hardware 6 CPU/11 GB, WSL en D:, export/import, pipefail solo en prueba, conda-frontend conda).

### C.7 Git

```
3a42ca9 (HEAD -> main, tag: v0.1.1-semana1, origin/main) docs: bitácora Semana1 con desviaciones...
f1c5385 fix: conda apunta a lock.yml; conda-frontend=conda
3baf5f4 fix: quitar prefix de locks; recrear ateles_busco con BUSCO 5.8.3
refs/tags/v0.1.1-semana1 en origin
```

---

## Lista de decisiones tomadas

1. IQ-TREE 3.1.4 aceptado (prueba de humo con -m/-B/--alrt OK).
2. BUSCO recreado en 5.8.3; linaje operativo `primates_odb12.2` (no odb10 en list-datasets). **Pendiente Semana 3:** verificar descarga de `primates_odb10` o justificar formalmente el uso de odb12.2 (ver actualización 2026-10-21 en bitácora).
3. Distro movida a D: vía export/import (no --move).
4. .wslconfig: 11GB RAM, 6 processors, 16GB swap en D:.
5. Trabajo en `~/ateles_immune_genomics` (nativo); `/mnt/d/..._respaldo` conservado.
6. Locks sin `prefix:`; Snakemake usa `.lock.yml`; `conda-frontend: conda`.
7. `set +o pipefail` solo en regla de prueba.
8. Nombre LICENSE/commits: Daniel Albuja (sin reescritura de historia).
9. Etiqueta de cierre: **v0.1.1-semana1**.
10. OrthoFinder 3.1.5 aceptado (default `-S diamond`); citar OrthoFinder 3 además de Emms & Kelly 2019.
11. Binario de IQ-TREE: `iqtree3` (reglas futuras deben usarlo).

**Estado Semana 1:** CERRADA — lista para auditoría externa y arranque de Semana 2.
