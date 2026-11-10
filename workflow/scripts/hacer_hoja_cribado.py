#!/usr/bin/env python3
"""Genera cribado_manual.xlsx a partir de cribado_base.tsv (v2.1)."""
from __future__ import annotations
import csv
import re
from pathlib import Path
from collections import defaultdict
from openpyxl import Workbook
from openpyxl.worksheet.datavalidation import DataValidation
from openpyxl.utils import get_column_letter
from openpyxl.styles import Font, Alignment, PatternFill
from Bio import Medline

GENES = [
    "TLR1","TLR2","TLR3","TLR4","TLR7","TLR8","TLR9","TRIM5","APOBEC3G","APOBEC3H",
    "MX1","MX2","OAS1","OAS2","BST2","ACE2","CD4","DDX58","MAVS",
]
RE_GENES = {g: re.compile(rf"\b{re.escape(g.lower())}\b", re.I) for g in GENES}
RE_GENES["TRIM5"] = re.compile(r"\btrim5\w*", re.I)
RE_GENES["DDX58"] = re.compile(r"\bddx58\b|\brig-?i\b", re.I)
RE_GENES["BST2"] = re.compile(r"\bbst-?2\b|tetherin", re.I)
RE_GENES["ACE2"] = re.compile(r"\bace-?2\b|angiotensin-converting enzyme 2", re.I)

def load_abstracts(root: Path, fetch: list) -> dict:
    abs_map = {}
    for qid in fetch:
        med = root / qid / "records.medline"
        if not med.exists():
            continue
        with med.open(encoding="utf-8", errors="replace") as fh:
            for rec in Medline.parse(fh):
                pmid = rec.get("PMID")
                if pmid and pmid not in abs_map:
                    abs_map[pmid] = rec.get("AB") or ""
    return abs_map

def genes_in(text: str) -> str:
    t = text.lower()
    hits = [g for g, rx in RE_GENES.items() if rx.search(t)]
    return ";".join(hits)

def main():
    import yaml
    root = Path("results/bibliografia")
    cfg = yaml.safe_load(Path("config/config.yaml").read_text())
    fetch = cfg["pubmed"]["fetch"]
    abs_map = load_abstracts(root, fetch)
    rows = list(csv.DictReader((root / "cribado_base.tsv").open(encoding="utf-8"), delimiter="\t"))

    # dedup: keep first by pmid order, note dup
    by_pmid = {}
    for r in rows:
        if r["pmid"] not in by_pmid:
            by_pmid[r["pmid"]] = r
        else:
            pass
    # mark dups from column
    for r in by_pmid.values():
        r["dup_con"] = r.get("pmid_duplicado") or ""

    def enrich(r):
        ab = abs_map.get(r["pmid"], "")
        title = r.get("title") or ""
        text = f"{title} {ab}"
        return {
            **r,
            "resumen_completo": ab,
            "enlace": f"https://pubmed.ncbi.nlm.nih.gov/{r['pmid']}/",
            "genes_detectados": genes_in(text),
            "terminos_primate": "",  # optional; leave blank or simple
            "dup_con": r.get("pmid_duplicado") or "",
        }

    sheets_map = {
        "P1a": [enrich(r) for r in by_pmid.values() if r["prioridad"] == "P1a"],
        "P2a": [enrich(r) for r in by_pmid.values() if r["prioridad"] == "P2a"],
        "P3a": [enrich(r) for r in by_pmid.values() if r["prioridad"] == "P3a"],
        "P3c_ACE2": [enrich(r) for r in by_pmid.values() if r["prioridad"] == "P3c_ACE2"],
        "P3d_ATELES": [enrich(r) for r in by_pmid.values() if r["prioridad"] == "P3d_ATELES"],
        "Excluidos_P1b_P2b": [enrich(r) for r in by_pmid.values() if r["prioridad"] in ("P1b", "P2b")],
    }
    # sort P1a
    sheets_map["P1a"].sort(
        key=lambda r: (
            0 if r.get("nivel_sugerido") == "A" else 1,
            0 if str(r.get("hit_seleccion")).lower() == "true" else 1,
            -int(r.get("n_consultas") or 0),
            -int(r.get("year") or 0),
        )
    )

    cols = [
        "pmid", "enlace", "year", "title", "journal", "doi", "resumen_completo",
        "prioridad", "nivel_sugerido", "foco_ateles", "consultas", "genes_detectados",
        "terminos_primate", "ancla", "ancla_nombre", "dup_con",
        "paso1_titulo", "paso2_resumen", "motivo_exclusion", "nivel_final",
        "genes_diana", "especies", "notas",
    ]

    wb = Workbook()
    wb.remove(wb.active)
    header_font = Font(bold=True)
    header_fill = PatternFill("solid", fgColor="D9E1F2")

    dv1 = DataValidation(type="list", formula1='"Incluir,Dudoso,Excluir"', allow_blank=True)
    dv2 = DataValidation(type="list", formula1='"Incluir,Excluir"', allow_blank=True)
    dv3 = DataValidation(
        type="list",
        formula1='"sin_primate,sin_gen_diana,modelo_sin_componente_comparativo,humano_clinico,vacuna_terapia,virus_o_celulas,duplicado,otro"',
        allow_blank=True,
    )
    dv4 = DataValidation(type="list", formula1='"A,B,C,foco"', allow_blank=True)

    for name, data in sheets_map.items():
        ws = wb.create_sheet(name)
        for c, h in enumerate(cols, 1):
            cell = ws.cell(1, c, h)
            cell.font = header_font
            cell.fill = header_fill
        for i, r in enumerate(data, 2):
            for c, h in enumerate(cols, 1):
                val = r.get(h, "")
                ws.cell(i, c, val)
        ws.freeze_panes = "A2"
        ws.auto_filter.ref = f"A1:{get_column_letter(len(cols))}{max(len(data)+1, 2)}"
        ws.column_dimensions["D"].width = 50  # title
        ws.column_dimensions["G"].width = 60  # resumen
        ws.column_dimensions["B"].width = 40
        # validations on decision cols
        for dv, col in [(dv1, 17), (dv2, 18), (dv3, 19), (dv4, 20)]:
            dv.add(f"{get_column_letter(col)}2:{get_column_letter(col)}{max(len(data)+1, 2)}")
        ws.add_data_validation(dv1)
        ws.add_data_validation(dv2)
        ws.add_data_validation(dv3)
        ws.add_data_validation(dv4)
        # mark excluidos motivo default
        if name == "Excluidos_P1b_P2b":
            for i in range(2, len(data) + 2):
                ws.cell(i, 19, "virus_o_celulas")
                ws.cell(i, 17, "Excluir")


    # Calibracion_P1a_50
    import random as _random
    audit_path = root / "auditoria_P1a_50.tsv"
    if not audit_path.exists():
        audit_path = Path("docs/bibliografia/auditoria_P1a_50.tsv")
    audit_pmids = []
    if audit_path.exists():
        with audit_path.open(encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                audit_pmids.append(row["pmid"])
    p1a_by_pmid = {r["pmid"]: r for r in sheets_map["P1a"]}
    cal = [p1a_by_pmid[p] for p in audit_pmids if p in p1a_by_pmid]
    _rng = _random.Random(20261110)
    _rng.shuffle(cal)

    # write Calibracion sheet
    name = "Calibracion_P1a_50"
    ws = wb.create_sheet(name)
    for c, h in enumerate(cols, 1):
        cell = ws.cell(1, c, h)
        cell.font = header_font
        cell.fill = header_fill
    for i, r in enumerate(cal, 2):
        for c, h in enumerate(cols, 1):
            ws.cell(i, c, r.get(h, ""))
    ws.freeze_panes = "A2"
    ws.auto_filter.ref = f"A1:{get_column_letter(len(cols))}{max(len(cal)+1, 2)}"
    ws.column_dimensions["D"].width = 50
    ws.column_dimensions["G"].width = 60
    for dv, col in [(dv1, 17), (dv2, 18), (dv3, 19), (dv4, 20)]:
        dv.add(f"{get_column_letter(col)}2:{get_column_letter(col)}{max(len(cal)+1, 2)}")
    ws.add_data_validation(dv1)
    ws.add_data_validation(dv2)
    ws.add_data_validation(dv3)
    ws.add_data_validation(dv4)

    # Criterios
    ws = wb.create_sheet("Criterios")
    criterios_txt = (
        "PASO 1 (título + resumen visible):\n"
        "Incluir = el estudio trata primates (cualquiera; prioridad platirrinos) Y uno de los 19 genes diana "
        "(TLR1,2,3,4,7,8,9; TRIM5; APOBEC3G/H; MX1/2; OAS1/2; BST2; ACE2; CD4 como receptor de entrada; DDX58; MAVS) "
        "o su familia cercana (TLR5/6/10, otros APOBEC3, factores de restricción) Y un componente comparativo: "
        "selección positiva, divergencia o polimorfismo entre especies, restricción o susceptibilidad específica "
        "de especie, o estructura-función comparada.\n"
        "Dudoso = cumple dos de las tres condiciones, o el título es ambiguo, o es un gen fuera de la lista con "
        "selección positiva en primates (motivo_exclusion = gen_fuera_lista).\n"
        "Excluir = primates solo como modelo experimental (vacunas, adyuvantes, terapias, farmacología); "
        "expresión de un gen en un tejido sin comparación entre especies; estudios solo en humanos (clínicos, GWAS); "
        "virus sin componente del hospedador; genes fuera de lista sin selección/evolución.\n\n"
        "PASO 2 (resumen completo): solo para Incluir y Dudoso. Confirmar o excluir; asignar nivel_final, genes_diana y especies.\n\n"
        "NIVELES: A = platirrinos + gen diana (comparativo); B = primates no platirrinos + gen diana (comparativo); "
        "C = mamíferos/vertebrados comparativos o metodológicos (p. ej. ACE2 entre especies, PAML/HyPhy); "
        "foco = P3d_ATELES (patógenos/inmunidad/enfermedad en Ateles y Atelidae). Excluidos de P3d: comportamiento, "
        "endocrinología, reproducción.\n\n"
        "P3a, P3c y P3d: solo por título; Incluir = útil como contexto comparativo o del taxón focal."
    )
    for i, line in enumerate(criterios_txt.split("\n"), 1):
        ws.cell(i, 1, line)
    ws.column_dimensions["A"].width = 120

    out = root / "cribado_manual.xlsx"
    wb.save(out)
    print(f"escrito {out} size={out.stat().st_size}")
    print("filas Calibracion_P1a_50:", len(cal))


    # excluidos P3b tsv
    with open(root / "excluidos_P3b.tsv", "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["pmid", "year", "title", "motivo"], delimiter="\t")
        w.writeheader()
        for r in by_pmid.values():
            if r["prioridad"] == "P3b":
                w.writerow({"pmid": r["pmid"], "year": r["year"], "title": r["title"], "motivo": "sin_gen_ni_primate"})
    print("excluidos_P3b.tsv ok")

if __name__ == "__main__":
    main()
