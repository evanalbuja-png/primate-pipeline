# Primate Immune Genomics Pipeline

**Proyecto:** Predicción de susceptibilidad viral en primates neotropicales mediante genómica comparativa de genes de inmunidad innata: estudio de caso del mono araña ecuatoriano (*Ateles fusciceps*)

**Institución:** Instituto Nacional de Biodiversidad (INABIO) – Dirección de Innovación / Laboratorio de Secuenciamiento  
**Responsable:** Daniel Albuja  
**Supervisor:** Pablo Jarrín  
**Duración:** 24 semanas  

## Estado actual

**Semana 1 completada** — Entorno de análisis reproducible configurado.

## Estructura del repositorio
ateles_immune_genomics/
├── config/                 # Configuración del pipeline
├── envs/                   # Definiciones de entornos conda (YAML + locks)
├── workflow/
│   ├── Snakefile           # Flujo principal
│   ├── rules/              # Reglas adicionales (futuras)
│   ├── scripts/            # Scripts auxiliares
│   └── profiles/default/   # Perfil de ejecución
├── data/                   # Datos (no versionados)
│   ├── raw/
│   └── processed/
├── results/                # Resultados (solo tablas/pequeños archivos versionados)
├── logs/                   # Logs de ejecución
└── docs/                   # Documentación y bitácora

## Entornos conda

Se utilizan cinco entornos aislados:

| Entorno            | Propósito                              |
|--------------------|----------------------------------------|
| `ateles_datos`     | Descarga y manipulación de datos       |
| `ateles_align_filo`| Alineamiento, filogenia, OrthoFinder, Snakemake |
| `ateles_seleccion` | PAML y HyPhy                           |
| `ateles_r`         | Visualización y análisis en R          |
| `ateles_busco`     | Evaluación de completitud genómica     |

### Recrear los entornos

```bash
mamba env create -f envs/ateles_datos.yml
mamba env create -f envs/ateles_align_filo.yml
mamba env create -f envs/ateles_seleccion.yml
mamba env create -f envs/ateles_r.yml
mamba env create -f envs/ateles_busco.yml

Los archivos *.lock.yml registran las versiones exactas instaladas.

Ejecutar el flujo
Desde la raíz del proyecto:
Bashmamba run -n ateles_align_filo snakemake \
  --snakefile workflow/Snakefile \
  --directory . \
  --profile workflow/profiles/default \
  --cores 4 \
  --use-conda

Software gráfico (fuera de conda)

AliView: https://ormbunkar.se/aliview/
UCSF ChimeraX: https://www.cgl.ucsf.edu/chimerax/download.html
ColabFold: se ejecuta en Google Colab (no se instala localmente)

Licencia
MIT License © 2026 Daniel Albuja