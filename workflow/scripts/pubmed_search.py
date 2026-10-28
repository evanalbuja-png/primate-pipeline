#!/usr/bin/env python3
"""PubMed search via Bio.Entrez — Semana 2 Ateles immune genomics.

Modes:
  count   : esearch retmax=0 → query.txt, count.txt, esearch.json
  fetch   : download PMIDs + MEDLINE
  resumen : build resumen_busquedas.tsv
  anclas  : resolve DOI→PMID for anchor papers

El término se lee SIEMPRE desde config.yaml (nunca por argv) para evitar
interpretación del shell.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import yaml
from Bio import Entrez
import Bio

TOOL = "primate-pipeline"
SLEEP_SEC = 0.34  # ~3 req/s


def load_config(path: str) -> dict:
    with open(path, encoding="utf-8") as f:
        return yaml.safe_load(f)


def setup_entrez() -> None:
    email = os.environ.get("NCBI_EMAIL")
    if not email:
        sys.exit("ERROR: NCBI_EMAIL no definida. export NCBI_EMAIL=tu@email")
    Entrez.email = email
    Entrez.tool = TOOL
    api_key = os.environ.get("NCBI_API_KEY")
    if api_key:
        Entrez.api_key = api_key


def sha256_text(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


def esearch_raw(term: str, retmax: int = 0, usehistory: str = "n", retries: int = 4) -> dict:
    last_err = None
    for attempt in range(retries):
        try:
            time.sleep(SLEEP_SEC)
            handle = Entrez.esearch(
                db="pubmed", term=term, retmax=retmax, usehistory=usehistory
            )
            record = Entrez.read(handle)
            handle.close()
            return record
        except Exception as e:
            last_err = e
            wait = min(2 ** attempt, 30)
            print(f"  retry {attempt+1}/{retries}: {e} (sleep {wait}s)", file=sys.stderr)
            time.sleep(wait)
    raise RuntimeError(f"esearch failed after {retries} retries: {last_err}")


def cmd_count(args: argparse.Namespace) -> None:
    setup_entrez()
    cfg = load_config(args.config)
    q = cfg["pubmed_queries"][args.qid]
    term = q["term"]
    verrs = validate_term(args.qid, term)
    if verrs:
        sys.exit(f"VALIDACION FALLA {args.qid}: {verrs}")
    tipo = q["tipo"]
    assert isinstance(term, str) and len(term) > 0

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    record = esearch_raw(term, retmax=0)
    count = int(record.get("Count", 0))
    qtrans = str(record.get("QueryTranslation", "") or "")
    error_list = record.get("ErrorList") or {}
    warning_list = record.get("WarningList") or {}
    # normalizar a dict serializable
    if not isinstance(error_list, dict):
        error_list = {"_raw": str(error_list)}
    if not isinstance(warning_list, dict):
        warning_list = {"_raw": str(warning_list)}

    sha = sha256_text(term)
    query_path = outdir / "query.txt"
    query_path.write_text(term + "\n", encoding="utf-8")
    # assert byte-a-byte
    saved = query_path.read_text(encoding="utf-8").rstrip("\n")
    assert saved == term, "query.txt no coincide byte a byte con config term"

    fecha = datetime.now(timezone.utc).isoformat()
    meta = {
        "id": args.qid,
        "tipo": tipo,
        "fecha_utc": fecha,
        "biopython_version": Bio.__version__,
        "sha256": sha,
        "term": term,
        "Count": count,
        "QueryTranslation": qtrans,
        "ErrorList": {k: list(v) if not isinstance(v, str) else v for k, v in error_list.items()},
        "WarningList": {k: list(v) if not isinstance(v, str) else v for k, v in warning_list.items()},
    }
    (outdir / "esearch.json").write_text(
        json.dumps(meta, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    (outdir / "count.txt").write_text(str(count) + "\n", encoding="utf-8")
    print(f"{args.qid}\tcount={count}\tsha={sha[:12]}")



def efetch_raw(**kwargs) -> str:
    """efetch con reintentos y espera."""
    last_err = None
    for attempt in range(4):
        try:
            time.sleep(SLEEP_SEC)
            h = Entrez.efetch(**kwargs)
            data = h.read()
            h.close()
            if isinstance(data, bytes):
                data = data.decode("utf-8", errors="replace")
            return data
        except Exception as e:
            last_err = e
            wait = min(2 ** attempt, 30)
            print(f"  efetch retry {attempt+1}/4: {e} (sleep {wait}s)", file=sys.stderr)
            time.sleep(wait)
    raise RuntimeError(f"efetch failed after retries: {last_err}")

def cmd_fetch(args: argparse.Namespace) -> None:
    setup_entrez()
    cfg = load_config(args.config)
    term = cfg["pubmed_queries"][args.qid]["term"]
    outdir = Path(args.outdir)
    count = int((outdir / "count.txt").read_text().strip() or "0")
    max_fetch = int(args.max_fetch)
    if count > max_fetch:
        sys.exit(f"ERROR: {args.qid} count={count} > max_fetch={max_fetch}")

    pmids_path = outdir / "pmids.txt"
    medline_path = outdir / "records.medline"

    if count == 0:
        pmids_path.write_text("", encoding="utf-8")
        medline_path.write_text("", encoding="utf-8")
        print(f"{args.qid}\tfetch=0")
        return

    record = esearch_raw(term, retmax=count, usehistory="y")
    id_list = list(record.get("IdList", []))
    webenv = record.get("WebEnv")
    query_key = record.get("QueryKey")
    pmids_path.write_text("\n".join(id_list) + ("\n" if id_list else ""), encoding="utf-8")

    batch = 200
    chunks = []
    for start in range(0, len(id_list), batch):
        if webenv and query_key:
            data = efetch_raw(
                db="pubmed",
                rettype="medline",
                retmode="text",
                retstart=start,
                retmax=batch,
                webenv=webenv,
                query_key=query_key,
            )
        else:
            batch_ids = id_list[start : start + batch]
            data = efetch_raw(
                db="pubmed",
                id=",".join(batch_ids),
                rettype="medline",
                retmode="text",
            )
        chunks.append(data)
    medline_path.write_text("".join(chunks), encoding="utf-8")

    n_pmid = sum(1 for line in pmids_path.read_text().splitlines() if line.strip().isdigit())
    n_rec = sum(1 for line in medline_path.read_text().splitlines() if line.startswith("PMID- "))
    if n_pmid != count or n_rec != count:
        sys.exit(
            f"ERROR conteo inconsistente {args.qid}: count={count} pmids={n_pmid} records={n_rec}"
        )
    print(f"{args.qid}\tfetch={n_pmid}")


def cmd_resumen(args: argparse.Namespace) -> None:
    root = Path(args.root)
    ids = args.ids.split(",")
    out_tsv = Path(args.output)
    rows = []
    for qid in ids:
        meta = json.loads((root / qid / "esearch.json").read_text(encoding="utf-8"))
        rows.append({
            "id": qid,
            "tipo": meta.get("tipo", ""),
            "cadena_exacta": meta.get("term", ""),
            "sha256": meta.get("sha256", ""),
            "query_translation": meta.get("QueryTranslation", ""),
            "error_list": json.dumps(meta.get("ErrorList") or {}, ensure_ascii=False),
            "warning_list": json.dumps(meta.get("WarningList") or {}, ensure_ascii=False),
            "fecha_utc": meta.get("fecha_utc", ""),
            "n_resultados": meta.get("Count", ""),
        })
    fields = [
        "id", "tipo", "cadena_exacta", "sha256", "query_translation",
        "error_list", "warning_list", "fecha_utc", "n_resultados",
    ]
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
        w.writeheader()
        w.writerows(rows)
    print(f"Escrito {out_tsv} ({len(rows)} filas)")


def cmd_anclas(args: argparse.Namespace) -> None:
    setup_entrez()
    anchors = [
        ("Sawyer_2005", "10.1073/pnas.0409853102"),
        ("Wlasiuk_Nachman_2010", "10.1093/molbev/msq104"),
        ("Damas_2020", "10.1073/pnas.2010146117"),
        ("Sironi_2015", "10.1038/nrg3905"),
        ("Enard_2016", "10.7554/eLife.12469"),
    ]
    results = []
    for name, doi in anchors:
        term = f"{doi}[aid]"
        rec = esearch_raw(term, retmax=5)
        ids = list(rec.get("IdList", []))
        results.append({"name": name, "doi": doi, "pmids": ids, "count": int(rec.get("Count", 0))})
        print(f"{name}\t{doi}\tpmids={ids}")
    Path(args.output).write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")



import re

def validate_term(qid: str, term: str) -> list[str]:
    """Return list of error messages (empty = OK)."""
    errs = []
    # 1. AND/OR must have spaces on both sides
    if re.search(r"(?<!\s)\b(AND|OR)\b|\b(AND|OR)\b(?!\s)", term):
        errs.append("AND/OR sin espacio a ambos lados")
    # 2. balanced parentheses and even number of double quotes
    if term.count("(") != term.count(")"):
        errs.append(f"paréntesis desbalanceados (={term.count('(')} )={term.count(')')})")
    if term.count('"') % 2 != 0:
        errs.append(f"número impar de comillas dobles ({term.count(chr(34))})")
    # 3a. [tiab] preceded by word/wildcard/closing quote (not space)
    for m in re.finditer(r"\[tiab\]", term):
        i = m.start()
        if i == 0 or term[i-1].isspace():
            errs.append("[tiab] precedido de espacio o inicio")
        # 3b. after [tiab] next char must be space, ) or end
        j = m.end()
        if j < len(term) and term[j] not in (" ", ")"):
            errs.append(f"tras [tiab] caracter inválido {term[j]!r}")
    # 3c. quote parity context
    in_quote = False
    for i, ch in enumerate(term):
        if ch == '"':
            if not in_quote:
                # opening: prev must be space, ( or start
                if i > 0 and term[i-1] not in (" ", "("):
                    errs.append(f"comilla apertura mal precedida en pos {i}")
                in_quote = True
            else:
                # closing: next must be [, space, ) or end
                nxt = term[i+1] if i+1 < len(term) else ""
                if nxt and nxt not in ("[", " ", ")"):
                    errs.append(f"comilla cierre mal seguida en pos {i} por {nxt!r}")
                in_quote = False
    # 4. no double spaces, no non-ASCII
    if "  " in term:
        errs.append("espacios dobles")
    non_ascii = [c for c in term if ord(c) > 127]
    if non_ascii:
        errs.append(f"caracteres no ASCII: {set(non_ascii)}")
    return errs


def cmd_validar(args: argparse.Namespace) -> None:
    cfg = load_config(args.config)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    all_ok = True
    for qid, q in cfg["pubmed_queries"].items():
        term = q["term"]
        errs = validate_term(qid, term)
        status = "ok" if not errs else "falla"
        if errs:
            all_ok = False
        rows.append({
            "id": qid,
            "status": status,
            "n_errores": len(errs),
            "detalle": "; ".join(errs),
            "len": len(term),
            "sha256": sha256_text(term),
        })
        print(f"{qid}\t{status}\t{errs or '-'}")
    fields = ["id", "status", "n_errores", "detalle", "len", "sha256"]
    with out.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"Escrito {out}")
    if not all_ok:
        sys.exit(2)

def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("count")
    p.add_argument("--config", required=True)
    p.add_argument("--qid", required=True)
    p.add_argument("--outdir", required=True)
    p.set_defaults(func=cmd_count)

    p = sub.add_parser("fetch")
    p.add_argument("--config", required=True)
    p.add_argument("--qid", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--max-fetch", type=int, default=1500)
    p.set_defaults(func=cmd_fetch)

    p = sub.add_parser("resumen")
    p.add_argument("--root", required=True)
    p.add_argument("--ids", required=True)
    p.add_argument("--output", required=True)
    p.set_defaults(func=cmd_resumen)

    p = sub.add_parser("anclas")
    p.add_argument("--output", required=True)
    p.set_defaults(func=cmd_anclas)

    p = sub.add_parser("validar")
    p.add_argument("--config", required=True)
    p.add_argument("--output", required=True)
    p.set_defaults(func=cmd_validar)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
