# workflow/rules/pubmed.smk
# Semana 2 — PubMed via Bio.Entrez (script workflow/scripts/pubmed_search.py)

OUTDIR = config["pubmed"]["outdir"]
QUERIES = config["pubmed_queries"]
QUERY_IDS = list(QUERIES.keys())
MAX_FETCH = config["pubmed"].get("max_fetch", 1500)
SCRIPT = "workflow/scripts/pubmed_search.py"
CONFIG = "config/config.yaml"

rule conteo_pubmed:
    output:
        query = f"{OUTDIR}/{{qid}}/query.txt",
        count_file = f"{OUTDIR}/{{qid}}/count.txt",
        esearch = f"{OUTDIR}/{{qid}}/esearch.json",
    conda:
        "../../envs/ateles_pubmed.yml"
    resources:
        ncbi=1
    log:
        f"logs/pubmed/conteo_{{qid}}.log"
    shell:
        """
        python {SCRIPT} count \
            --config {CONFIG} \
            --qid {wildcards.qid} \
            --outdir {OUTDIR}/{wildcards.qid} \
            > {log} 2>&1
        """

rule resumen_busquedas:
    input:
        expand(f"{OUTDIR}/{{qid}}/esearch.json", qid=QUERY_IDS),
    output:
        tsv = f"{OUTDIR}/resumen_busquedas.tsv"
    conda:
        "../../envs/ateles_pubmed.yml"
    params:
        ids = ",".join(QUERY_IDS),
    shell:
        """
        python {SCRIPT} resumen \
            --root {OUTDIR} \
            --ids {params.ids} \
            --output {output.tsv}
        """

rule busqueda_pubmed:
    input:
        count_file = f"{OUTDIR}/{{qid}}/count.txt",
        esearch = f"{OUTDIR}/{{qid}}/esearch.json",
    output:
        pmids = f"{OUTDIR}/{{qid}}/pmids.txt",
        records = f"{OUTDIR}/{{qid}}/records.medline",
    conda:
        "../../envs/ateles_pubmed.yml"
    resources:
        ncbi=1
    log:
        f"logs/pubmed/fetch_{{qid}}.log"
    shell:
        """
        python {SCRIPT} fetch \
            --config {CONFIG} \
            --qid {wildcards.qid} \
            --outdir {OUTDIR}/{wildcards.qid} \
            --max-fetch {MAX_FETCH} \
            > {log} 2>&1
        """
