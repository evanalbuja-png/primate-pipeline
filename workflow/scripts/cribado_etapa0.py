#!/usr/bin/env python3
"""Etapa 0 del cribado bibliográfico — prefiltro PubMed (v2.1).

v2.1: filtra nombres de virus (incl. herpesvirus ateles); foco_ateles; P3d_ATELES y líneas celulares antes de hit_primate_real /
hit_platirrino_real; prioridades P1a/P1b/P2a/P2b/P3a/P3c_ACE2/P3b.

Uso:
  python workflow/scripts/cribado_etapa0.py --test
  python workflow/scripts/cribado_etapa0.py --config config/config.yaml \\
      --root results/bibliografia \\
      --output results/bibliografia/cribado_base.tsv
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

from Bio import Medline

VERSION = "v2.1"
SEED = 20261027
ANCHORS = {
    "15689398": "Sawyer_2005",
    "20410160": "Wlasiuk_Nachman_2010",
    "32826334": "Damas_2020",
    "32583896": "Torosin_2020_TLR7_Alouatta",
    "32768453": "Torosin_2020_TLR7_TLR8_primates",
}

# --- patrones a ELIMINAR antes de primate_real / plat_real ---
RE_STRIP = re.compile(
    r"herpes\w*\s?virus saimiri|saimiriine|\bhvs\b|"
    r"simian immunodeficiency virus\w*|simian[- ]human immunodeficiency|"
    r"simian foamy|simian virus \d+|simian retrovirus|simian varicella|"
    r"simian hemorrhagic|simian herpes\w*|simian sarcoma|"
    r"simian t[- ]lymphotropic|mason-pfizer monkey virus|"
    r"squirrel monkey retrovirus|monkeypox|monkey pox|"
    r"herpesvirus aotus|herpes\w*\s?virus ateles|ateline (alpha|gamma)?herpesvirus|\bhva\b|\bsiv\w*\b|"
    r"\bvero(\s?e6)?\b|\bcos-?7\b|llc-?mk2|\bma-?104\b|\bcv-?1\b|"
    r"\bbs-?c-?1\b|\bb95-?8\b|(african green )?monkey kidney|monkey cells?",
    re.I,
)

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

RE_ACE2 = re.compile(r"\bace-?2\b|angiotensin-converting enzyme 2", re.I)
RE_ACE2_CROSS = re.compile(
    r"cross-species|host range|animal species|species-specific|"
    r"wild (animals|mammals)|domestic and wild|vertebrates|orthologs?|"
    r"across species|different species|multiple species|"
    r"species susceptib\w+|susceptibility of (animals|mammals|species)",
    re.I,
)


def norm_match(text: str) -> str:
    t = (text or "").lower()
    return re.sub(r"\s+", " ", t).strip()


def norm_titulo(text: str) -> str:
    return re.sub(r"[^a-z0-9]", "", (text or "").lower())


def strip_virus_cell(text: str) -> str:
    return RE_STRIP.sub(" ", text)


def flags_for(text: str) -> dict:
    tm = norm_match(text)
    tm_clean = norm_match(strip_virus_cell(tm))

    hit_plat = bool(RE_PLATIRRINO.search(tm))
    hit_prim = hit_plat or bool(RE_PRIMATE_EXTRA.search(tm))
    hit_plat_real = bool(RE_PLATIRRINO.search(tm_clean))
    hit_prim_real = hit_plat_real or bool(RE_PRIMATE_EXTRA.search(tm_clean))

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
    hit_ace2 = bool(RE_ACE2.search(tm))
    hit_ace2_cross = hit_ace2 and bool(RE_ACE2_CROSS.search(tm))

    RE_FOCO_ATELES = re.compile(
        r"\bateles\b|spider monkeys?|atelidae|atelid\b|atelines?\b", re.I
    )
    foco_ateles = bool(RE_FOCO_ATELES.search(tm_clean))

    return dict(
        hit_platirrino=hit_plat,
        hit_primate=hit_prim,
        hit_platirrino_real=hit_plat_real,
        hit_primate_real=hit_prim_real,
        hit_gen_diana=hit_gen,
        hit_gen_familia=hit_fam,
        hit_seleccion=hit_sel,
        hit_virus=hit_vir,
        cd4_ambiguo=cd4_amb,
        hit_ace2=hit_ace2,
        hit_ace2_cross=hit_ace2_cross,
        foco_ateles=foco_ateles,
        texto_match=tm,
        texto_match_clean=tm_clean,
    )


def assign_priority_v2(hit: dict, tipos: set) -> tuple[str, str]:
    non_lit = tipos - {"literal"}
    has_non_lit = bool(non_lit)
    gen = hit["hit_gen_diana"] or hit["hit_gen_familia"]
    # v1 priority first
    if hit["hit_primate"] and gen and has_non_lit:
        pri_v1 = "P1"
    elif hit["hit_primate"] and gen and not has_non_lit:
        pri_v1 = "P2"
    elif gen and hit["hit_seleccion"] and not hit["hit_primate"]:
        pri_v1 = "P3a"
    else:
        pri_v1 = "P3b"

    if gen and hit["hit_primate_real"] and has_non_lit:
        pri = "P1a"
    elif pri_v1 == "P1":
        pri = "P1b"
    elif gen and hit["hit_primate_real"] and not has_non_lit:
        pri = "P2a"
    elif pri_v1 == "P2":
        pri = "P2b"
    elif pri_v1 == "P3a":
        pri = "P3a"
    elif pri_v1 == "P3b" and hit["hit_ace2_cross"]:
        pri = "P3c_ACE2"
    else:
        pri = "P3b"

    # P3d_ATELES: foco Ateles fuera de P1a/P2a
    if hit.get("foco_ateles") and pri not in ("P1a", "P2a", "P1b", "P2b"):
        pri = "P3d_ATELES"

    if hit["hit_platirrino_real"] and gen:
        nivel = "A"
    elif hit["hit_primate_real"] and gen and not hit["hit_platirrino_real"]:
        nivel = "B"
    elif pri in ("P3a", "P3c_ACE2"):
        nivel = "C"
    elif pri == "P3d_ATELES":
        nivel = "foco"
    else:
        nivel = ""
    return pri, nivel, pri_v1


def run_tests() -> int:
    cases = [
        # v1 cases
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
        # v2 virus/cell strip → real flags
        ("herpesvirus saimiri", "plat_real", False),
        ("simian immunodeficiency virus", "prim_real", False),
        ("SIVmac239", "prim_real", False),
        ("Mason-Pfizer monkey virus", "prim_real", False),
        ("simian foamy virus", "prim_real", False),
        ("monkeypox", "prim_real", False),
        ("Vero E6 cells", "prim_real", False),
        ("African green monkey kidney cells", "prim_real", False),
        ("COS-7", "prim_real", False),
        ("rhesus macaques infected with SIV", "prim_real", True),
        ("owl monkey TRIMCyp", "plat_real", True),
        ("Callithrix jacchus", "plat_real", True),
        ("New World monkeys", "plat_real", True),
        (
            "Predicting susceptibility to SARS-CoV-2 based on ACE2 across species",
            "ace2_cross",
            True,
        ),
        ("ACE2 expression in lung tissue of COVID-19 patients", "ace2_cross", False),
        ("Herpesvirus ateles", "foco_false", True),
        ("SARS-CoV-2 infection in Ateles fusciceps", "foco", True),
        ("brown spider monkeys", "foco", True),
        ("Atelidae", "foco", True),
        ("spider mite", "foco", False),
    ]
    fails = 0
    for text, kind, expect in cases:
        f = flags_for(text)
        if kind == "plat":
            got = f["hit_platirrino"]
        elif kind == "plat_real":
            got = f["hit_platirrino_real"]
        elif kind == "prim_real":
            got = f["hit_primate_real"]
        elif kind == "gen":
            got = f["hit_gen_diana"]
        elif kind == "gen_amb":
            got = (not f["hit_gen_diana"]) and f["cd4_ambiguo"]
            expect = True
        elif kind == "fam":
            got = f["hit_gen_familia"]
        elif kind == "sel":
            got = f["hit_seleccion"]
        elif kind == "ace2_cross":
            got = f["hit_ace2_cross"]
        elif kind == "foco":
            got = f["foco_ateles"]
        elif kind == "foco_false":
            got = (not f["foco_ateles"]) and (not f["hit_platirrino_real"])
            expect = True
        else:
            got = None
        ok = got == expect
        if not ok:
            fails += 1
        print(
            f"{'OK' if ok else 'FAIL'}\t{kind}\texpect={expect}\tgot={got}\t{text!r}"
        )
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
    pmid_data = {}
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
                    pmid_data[pmid] = {"rec": rec, "consultas": set(), "tipos": set()}
                pmid_data[pmid]["consultas"].add(qid)
                pmid_data[pmid]["tipos"].add(tipo)
    return pmid_data


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
        pri, nivel, pri_v1 = assign_priority_v2(hit, blob["tipos"])
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
            "hit_platirrino_real": str(hit["hit_platirrino_real"]).lower(),
            "hit_primate_real": str(hit["hit_primate_real"]).lower(),
            "hit_gen_diana": str(hit["hit_gen_diana"]).lower(),
            "hit_gen_familia": str(hit["hit_gen_familia"]).lower(),
            "hit_seleccion": str(hit["hit_seleccion"]).lower(),
            "hit_virus": str(hit["hit_virus"]).lower(),
            "cd4_ambiguo": str(hit["cd4_ambiguo"]).lower(),
            "hit_ace2": str(hit["hit_ace2"]).lower(),
            "hit_ace2_cross": str(hit["hit_ace2_cross"]).lower(),
            "foco_ateles": str(hit["foco_ateles"]).lower(),
            "prioridad": pri,
            "prioridad_v1": pri_v1,
            "nivel_sugerido": nivel,
            "ancla": str(pmid in ANCHORS).lower(),
            "ancla_nombre": ANCHORS.get(pmid, ""),
            "dup_titulo": "false",
            "pmid_duplicado": "",
            "ab_preview": (ab or "")[:300],
            "version": VERSION,
        })

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
                pmid_to_dup[p] = [x for x in uniq if x != p][0]
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

    rng = random.Random(SEED)

    def write_audit(path, pool, n):
        sample = pool if len(pool) <= n else rng.sample(pool, n)
        with open(path, "w", encoding="utf-8", newline="") as fh:
            cols = ["pmid", "year", "title", "prioridad", "nivel_sugerido", "ab_preview"]
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(sample)
        print(f"audit {path}: {len(sample)} de {len(pool)}")

    pools = {
        args.audit_p1a: ([r for r in rows if r["prioridad"] == "P1a"], 50),
        args.audit_p1b: ([r for r in rows if r["prioridad"] == "P1b"], 30),
        args.audit_p3c: ([r for r in rows if r["prioridad"] == "P3c_ACE2"], 30),
        args.audit_p3b: ([r for r in rows if r["prioridad"] == "P3b"], 100),
    }
    for path, (pool, n) in pools.items():
        write_audit(path, pool, n)
    print(f"cribado_base {VERSION}: {len(rows)} filas -> {args.output}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--test", action="store_true")
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--root", default="results/bibliografia")
    ap.add_argument("--output", default="results/bibliografia/cribado_base.tsv")
    ap.add_argument("--audit-p1a", default="results/bibliografia/auditoria_P1a_50.tsv")
    ap.add_argument("--audit-p1b", default="results/bibliografia/auditoria_P1b_30.tsv")
    ap.add_argument("--audit-p3c", default="results/bibliografia/auditoria_P3c_ACE2_30.tsv")
    ap.add_argument("--audit-p3b", default="results/bibliografia/auditoria_P3b_100.tsv")
    args = ap.parse_args()
    if args.test:
        sys.exit(run_tests())
    run_pipeline(args)


if __name__ == "__main__":
    main()
