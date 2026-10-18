Bitácora del proyecto — Ateles Immune Genomics
Semana 1 (9 de octubre de 2026)
Actividades realizadas

Diagnóstico del sistema (Ubuntu 24.04.4 LTS, 4 CPUs, ~9.7 GiB RAM).
Creación de la estructura de directorios del proyecto.
Configuración de .gitignore (excluye datos pesados y .snakemake/).
Creación de cinco entornos conda aislados:
ateles_datos
ateles_align_filo
ateles_seleccion
ateles_r
ateles_busco

Exportación de archivos .lock.yml para cada entorno.
Verificación de todas las herramientas clave (todas OK).
Configuración de Snakemake 9.27.0 con regla de prueba verificar_entorno.
Documentación inicial (README, LICENSE MIT, bitácora).

Versiones relevantes

conda 26.3.2 / mamba 2.5.0
Snakemake 9.27.0
IQ-TREE 3.1.4
MAFFT 7.526
BUSCO 6.1.0
PAML (codeml) y HYPHY 2.5.103
R packages: ggtree 4.0.5, ComplexHeatmap 2.26.1, ape 5.8.1, phytools 2.5.2

Notas / Decisiones

Se trabaja en /mnt/d/Documentos/Proyectos/ateles_immune_genomics por limitación de almacenamiento en el sistema de archivos nativo de WSL.
IQ-TREE y BUSCO se instalaron en versiones superiores a las mencionadas en la propuesta original (3.x y 6.x respectivamente); son las versiones actuales disponibles en bioconda y funcionales.
Se corrigió un problema de SIGPIPE (exit 141) en la regla de prueba de Snakemake desactivando temporalmente pipefail.

Próximos pasos (Semana 2)

Revisión bibliográfica sistematizada.
Descarga y evaluación de calidad del genoma de Ateles fusciceps (BUSCO).