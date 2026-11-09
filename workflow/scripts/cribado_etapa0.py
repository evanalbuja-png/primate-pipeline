#!/usr/bin/env python3
"""Etapa 0 de cribado bibliográfico — prefiltro de PMIDs de las 18 consultas PubMed.

Uso:
  python workflow/scripts/cribado_etapa0.py --test
  python workflow/scripts/cribado_etapa0.py \\
      --config config/config.yaml \\
      --root results/bibliografia \\
      --output results/bibliografia/cribado_base.tsv \\
      --audit-p3 results/bibliografia/auditoria_P3b_100.tsv \\
      --audit-p1 results/bibliografia/auditoria_P1_50.tsv
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

from Bio import Medline

SEED = 20261027
ANCHORS = {
    "15689398": "Sawyer_2005",
    "20410160": "Wlasiuk_Nachman_2010",
    "32826334": "Damas_2020",
    "32583896": "Torosin_2020_TLR7_Alouatta",
    "32768453": "Torosin_2020_TLR7_TLR8_primates",
}

# --- regex (re.I) ---
RE_PLATIRRINO = re.compile(
    r"\b(alouatta|aotus|ateles|callithrix|cebus|sapajus|saimiri|saguinus|"
    r"callicebus|plecturocebus|lagothrix|brachyteles|pithecia|cacajao|"
    r"chiropotes|callimico|leontopithecus|mico)\b"
    r"|platyrrhin\w*"
    r"|\bnew[- ]world (monkey|monkeys|primate|primates)\b"
    r"|\b(marmosets?|tamarins?|capuchins?|howler monkeys?|owl monkeys?|"
    r"squirrel monkeys?|spider monkeys?|woolly monkeys?|titi monkeys?|"
    r"saki monkeys?|uakaris?)\b"
    r"|neotropical primate\w*",
    re.I,
)

RE_PRIMATE_EXTRA = re.compile(
    r"primate\w*"
    r"|\bmonkeys?\b"
    r"|\bsimian\w*"
    r"|catarrhin\w*"
    r"|\bold[- ]world (monkey|monkeys|primate|primates)\b"
    r"|\b(macaques?|macaca|rhesus|cynomolgus|baboons?|papio|chlorocebus|"
    r"african green monkeys?|vervets?|mangabeys?|cercopithec\w*|"
    r"chimpanzees?|bonobos?|gorillas?|orangutans?|gibbons?|great apes?|"
    r"pan troglodytes|lemurs?|galagos?|tarsiers?)\b",
    re.I,
)

RE_GEN_DIANA = re.compile(
    r"\btlr[1-4]\b|\btlr[7-9]\b|\btrim5\w*|\bapobec3[gh]\b|\bmx1\b|\bmx2\b|\bmx[ab]\b"
    r"|myxovirus resistance|\boas[12]\b|oligoadenylate synthetase|\bbst-?2\b|tetherin"
    r"|\bcd317\b|\bace-?2\b|angiotensin-converting enzyme 2|\bddx58\b|\brig-?i\b"
    r"|retinoic acid[- ]inducible gene|\bmavs\b|\bcardif\b|\bips-?1\b",
    re.I,
)

RE_CD4 = re.compile(r"\bcd4\b", re.I)
RE_CD4_CTX = re.compile(
    r"receptor|entry|co-?receptor|gp120|\benv\b|tropism|binding|host range|"
    r"species-specific|cross-species",
    re.I,
)

RE_GEN_FAMILIA = re.compile(
    r"\btlr(5|6|10)\b|toll[- ]like receptor\w*|\bapobec3\w*|\bapobec\b|"
    r"restriction factor\w*|\bmx proteins?\b|\boas\b",
    re.I,
)

RE_SEL = re.compile(
    r"positive(ly)? select\w*|adaptive evolution|\bd[nN]\s?/\s?d[sS]\b|"
    r"\bomega\b(?= *(>|=|ratio))|selective pressure\w*|\bpurifying selection\b|"
    r"\bevolutionar\w* (rate|conflict|arms race)|arms race",
    re.I,
)

RE_VIR = re.compile(
    r"virus|viral|retrovir\w*|lentivir\w*|\bhiv\b|\bsiv\b|coronavir\w*|\bsars\b|"
    r"yellow fever|flavivir\w*|herpesvir\w*|antiviral",
    re.I,
)


def norm_match(text: str) -> str:
    """Minúsculas, espacios colapsados; conserva guiones/barras/símbolos."""
    t = (text or "").lower()
    t = re.sub(r"\s+", " ", t).strip()
    return t


def norm_titulo(text: str) -> str:
    """Solo letras y dígitos en minúscula (sin espacios ni signos)."""
    return re.sub(r"[^a-z0-9]", "", (text or "").lower())


def flags_for(text: str) -> dict:
    tm = norm_match(text)
    hit_plat = bool(RE_PLATIRRINO.search(tm))
    hit_prim = hit_plat or bool(RE_PRIMATE_EXTRA.search(tm))
    hit_gen = bool(RE_GEN_DIANA.search(tm))
    cd4_amb = False
    if RE_CD4.search(tm):
        if RE_CD4_CTX.search(tm):
            hit_gen = True
        else:
            cd4_amb = True
    hit_fam = bool(RE_GEN_FAMILIA.search(tm))
    hit_sel = bool(RE_SEL.search(tm))
    hit_vir = bool(RE_VIR.search(tm))
    return dict(
        hit_platirrino=hit_plat,
        hit_primate=hit_prim,
        hit_gen_diana=hit_gen,
        hit_gen_familia=hit_fam,
        hit_seleccion=hit_sel,
        hit_virus=hit_vir,
        cd4_ambiguo=cd4_amb,
    )


def run_tests() -> int:
    cases = [
        ("marmosets", "plat", True),
        ("howler monkeys", "plat", True),
        ("squirrel monkeys", "plat", True),
        ("New World monkeys", "plat", True),
        ("platyrrhines", "plat", True),
        ("Ateles fusciceps", "plat", True),
        ("tamarind fruit", "plat", False),
        ("mice and rats", "plat", False),
        ("TRIM5alpha", "gen", True),
        ("TRIM5α", "gen", True),
        ("TLR7", "gen", True),
        ("RIG-I", "gen", True),
        ("ACE2", "gen", True),
        ("MxB", "gen", True),
        ("BST-2", "gen", True),
        ("CD4 receptor binding", "gen", True),
        ("CD4+ T cells proliferate", "gen_amb", True),
        ("APOBEC3F", "fam", True),
        ("Toll-like receptor 5", "fam", True),
        ("TLR10", "fam", True),
        ("dN/dS", "sel", True),
        ("positive selection", "sel", True),
        ("positively selected", "sel", True),
    ]
    fails = 0
    for text, kind, expect in cases:
        f = flags_for(text)
        if kind == "plat":
            got = f["hit_platirrino"]
        elif kind == "gen":
            got = f["hit_gen_diana"]
        elif kind == "gen_amb":
            got = (not f["hit_gen_diana"]) and f["cd4_ambiguo"]
            expect = True
        elif kind == "fam":
            got = f["hit_gen_familia"]
        elif kind == "sel":
            got = f["hit_seleccion"]
        else:
            got = None
        ok = got == expect
        if not ok:
            fails += 1
        print(f"{'OK' if ok else 'FAIL'}\t{kind}\texpect={expect}\tgot={got}\t{text!r}\t{f}")
    print(f"fails={fails}")
    return 1 if fails else 0


def extract_doi(rec: dict) -> str:
    for field in ("AID", "LID"):
        vals = rec.get(field)
        if not vals:
            continue
        if isinstance(vals, str):
            vals = [vals]
        for v in vals:
            if "[doi]" in v.lower():
                return v.lower().replace("[doi]", "").strip()
            if re.match(r"10\.\d+/", v):
                return v.strip()
    return ""


def load_records(root: Path, fetch_ids: list, queries_meta: dict):
    pmid_data = {}  # pmid -> aggregate
    for qid in fetch_ids:
        tipo = queries_meta.get(qid, {}).get("tipo", "")
        med = root / qid / "records.medline"
        if not med.exists():
            continue
        with med.open(encoding="utf-8", errors="replace") as fh:
            for rec in Medline.parse(fh):
                pmid = rec.get("PMID")
                if not pmid:
                    continue
                if pmid not in pmid_data:
                    pmid_data[pmid] = {
                        "rec": rec,
                        "consultas": set(),
                        "tipos": set(),
                    }
                pmid_data[pmid]["consultas"].add(qid)
                pmid_data[pmid]["tipos"].add(tipo)
    return pmid_data


def assign_priority(hit: dict, tipos: set) -> tuple[str, str]:
    """Return (prioridad, nivel_sugerido)."""
    non_lit = tipos - {"literal"}
    has_non_lit = bool(non_lit)
    gen = hit["hit_gen_diana"] or hit["hit_gen_familia"]
    if hit["hit_primate"] and gen and has_non_lit:
        pri = "P1"
    elif hit["hit_primate"] and gen and not has_non_lit:
        pri = "P2"
    elif gen and hit["hit_seleccion"] and not hit["hit_primate"]:
        pri = "P3a"
    else:
        pri = "P3b"
    if hit["hit_platirrino"] and gen:
        nivel = "A"
    elif hit["hit_primate"] and gen and not hit["hit_platirrino"]:
        nivel = "B"
    elif pri == "P3a":
        nivel = "C"
    else:
        nivel = ""
    return pri, nivel


def run_pipeline(args: argparse.Namespace) -> None:
    import random
    import yaml

    cfg = yaml.safe_load(Path(args.config).read_text(encoding="utf-8"))
    fetch = list(cfg["pubmed"]["fetch"])
    qmeta = cfg["pubmed_queries"]
    root = Path(args.root)
    pmid_data = load_records(root, fetch, qmeta)

    rows = []
    for pmid, blob in pmid_data.items():
        rec = blob["rec"]
        title = rec.get("TI") or ""
        ab = rec.get("AB") or ""
        ot = rec.get("OT") or []
        if isinstance(ot, str):
            ot = [ot]
        mh = rec.get("MH") or []
        if isinstance(mh, str):
            mh = [mh]
        texto = " ".join([title, ab] + list(ot) + list(mh))
        hit = flags_for(texto)
        pri, nivel = assign_priority(hit, blob["tipos"])
        dp = rec.get("DP") or ""
        m = re.search(r"(\d{4})", dp)
        year = m.group(1) if m else ""
        pt = rec.get("PT") or []
        if isinstance(pt, str):
            pt = [pt]
        rows.append({
            "pmid": pmid,
            "year": year,
            "title": title,
            "journal": rec.get("JT") or "",
            "doi": extract_doi(rec),
            "titulo_norm": norm_titulo(title),
            "consultas": ";".join(sorted(blob["consultas"])),
            "n_consultas": len(blob["consultas"]),
            "tipos": ";".join(sorted(blob["tipos"])),
            "solo_literal": str(blob["tipos"] <= {"literal"}).lower(),
            "pt": ";".join(pt),
            "tiene_resumen": str(bool(ab)).lower(),
            "hit_platirrino": str(hit["hit_platirrino"]).lower(),
            "hit_primate": str(hit["hit_primate"]).lower(),
            "hit_gen_diana": str(hit["hit_gen_diana"]).lower(),
            "hit_gen_familia": str(hit["hit_gen_familia"]).lower(),
            "hit_seleccion": str(hit["hit_seleccion"]).lower(),
            "hit_virus": str(hit["hit_virus"]).lower(),
            "cd4_ambiguo": str(hit["cd4_ambiguo"]).lower(),
            "prioridad": pri,
            "nivel_sugerido": nivel,
            "ancla": str(pmid in ANCHORS).lower(),
            "ancla_nombre": ANCHORS.get(pmid, ""),
            "dup_titulo": "false",
            "pmid_duplicado": "",
            "ab_preview": (ab or "")[:300],
        })

    # duplicados por doi o titulo_norm
    by_key = defaultdict(list)
    for r in rows:
        if r["doi"]:
            by_key[("doi", r["doi"])].append(r["pmid"])
        if r["titulo_norm"]:
            by_key[("ti", r["titulo_norm"])].append(r["pmid"])
    pmid_to_dup = {}
    for key, plist in by_key.items():
        uniq = sorted(set(plist))
        if len(uniq) > 1:
            for p in uniq:
                others = [x for x in uniq if x != p]
                pmid_to_dup[p] = others[0]
    for r in rows:
        if r["pmid"] in pmid_to_dup:
            r["dup_titulo"] = "true"
            r["pmid_duplicado"] = pmid_to_dup[r["pmid"]]

    fields = list(rows[0].keys()) if rows else []
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda x: x["pmid"]))

    # auditorías
    rng = random.Random(SEED)
    p3b = [r for r in rows if r["prioridad"] == "P3b"]
    p1 = [r for r in rows if r["prioridad"] == "P1"]
    def write_audit(path, pool, n):
        sample = pool if len(pool) <= n else rng.sample(pool, n)
        with open(path, "w", encoding="utf-8", newline="") as fh:
            cols = ["pmid", "year", "title", "prioridad", "nivel_sugerido", "ab_preview"]
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(sample)
        print(f"audit {path}: {len(sample)} de {len(pool)}")

    write_audit(args.audit_p3, p3b, 100)
    write_audit(args.audit_p1, p1, 50)
    print(f"cribado_base: {len(rows)} filas -> {args.output}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--test", action="store_true")
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--root", default="results/bibliografia")
    ap.add_argument("--output", default="results/bibliografia/cribado_base.tsv")
    ap.add_argument("--audit-p3", default="results/bibliografia/auditoria_P3b_100.tsv")
    ap.add_argument("--audit-p1", default="results/bibliografia/auditoria_P1_50.tsv")
    args = ap.parse_args()
    if args.test:
        sys.exit(run_tests())
    run_pipeline(args)


if __name__ == "__main__":
    main()
