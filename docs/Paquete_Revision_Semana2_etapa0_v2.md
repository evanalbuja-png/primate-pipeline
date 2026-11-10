# Paquete de revisión Semana 2 — Etapa 0 v2

**Commit v1 archivada:** 55695aa
**Versión script:** v2
**Fecha ejecución:** 2026-11-09 / 2026-11-10 UTC
**Estado:** Etapa 0 v2 ejecutada; 5/5 anclas en P1a; sin Zotero ni cribado manual

---
## a) Script completo

```python
#!/usr/bin/env python3
"""Etapa 0 del cribado bibliográfico — prefiltro PubMed (v2).

v2: filtra nombres de virus y líneas celulares antes de hit_primate_real /
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

VERSION = "v2"
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
    r"herpesvirus aotus|\bsiv\w*\b|"
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

    if hit["hit_platirrino_real"] and gen:
        nivel = "A"
    elif hit["hit_primate_real"] and gen and not hit["hit_platirrino_real"]:
        nivel = "B"
    elif pri in ("P3a", "P3c_ACE2"):
        nivel = "C"
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

```

---
## b) --test completo

```
OK	plat	expect=True	got=True	'marmosets'
OK	plat	expect=True	got=True	'howler monkeys'
OK	plat	expect=True	got=True	'squirrel monkeys'
OK	plat	expect=True	got=True	'New World monkeys'
OK	plat	expect=True	got=True	'platyrrhines'
OK	plat	expect=True	got=True	'Ateles fusciceps'
OK	plat	expect=False	got=False	'tamarind fruit'
OK	plat	expect=False	got=False	'mice and rats'
OK	gen	expect=True	got=True	'TRIM5alpha'
OK	gen	expect=True	got=True	'TRIM5α'
OK	gen	expect=True	got=True	'TLR7'
OK	gen	expect=True	got=True	'RIG-I'
OK	gen	expect=True	got=True	'ACE2'
OK	gen	expect=True	got=True	'MxB'
OK	gen	expect=True	got=True	'BST-2'
OK	gen	expect=True	got=True	'CD4 receptor binding'
OK	gen_amb	expect=True	got=True	'CD4+ T cells proliferate'
OK	fam	expect=True	got=True	'APOBEC3F'
OK	fam	expect=True	got=True	'Toll-like receptor 5'
OK	fam	expect=True	got=True	'TLR10'
OK	sel	expect=True	got=True	'dN/dS'
OK	sel	expect=True	got=True	'positive selection'
OK	sel	expect=True	got=True	'positively selected'
OK	plat_real	expect=False	got=False	'herpesvirus saimiri'
OK	prim_real	expect=False	got=False	'simian immunodeficiency virus'
OK	prim_real	expect=False	got=False	'SIVmac239'
OK	prim_real	expect=False	got=False	'Mason-Pfizer monkey virus'
OK	prim_real	expect=False	got=False	'simian foamy virus'
OK	prim_real	expect=False	got=False	'monkeypox'
OK	prim_real	expect=False	got=False	'Vero E6 cells'
OK	prim_real	expect=False	got=False	'African green monkey kidney cells'
OK	prim_real	expect=False	got=False	'COS-7'
OK	prim_real	expect=True	got=True	'rhesus macaques infected with SIV'
OK	plat_real	expect=True	got=True	'owl monkey TRIMCyp'
OK	plat_real	expect=True	got=True	'Callithrix jacchus'
OK	plat_real	expect=True	got=True	'New World monkeys'
OK	ace2_cross	expect=True	got=True	'Predicting susceptibility to SARS-CoV-2 based on ACE2 across species'
OK	ace2_cross	expect=False	got=False	'ACE2 expression in lung tissue of COVID-19 patients'
fails=0
```

---
## c) cribado_meta.txt

```
version: v2
fecha_utc: 2026-11-10T01:58:55.148102+00:00
script_sha256: cebb02f8c9026d05a5b417d9130086287cee7a42f6b742ca3ea67bb8fc3891a9
medline_sha256	base1_literal	9ec0f321e04c8daadcaf900edbb8f5099f16f2bc158f1e697f0543931f268227	records=837
medline_sha256	base1_refinada	0e6b4d31a7af63894443cdd681364ac5ee5de4ecba430f58d99b43e59aec3add	records=59
medline_sha256	base2_refinada	7f0eb88a827cb66af55be82a69f2331bb910b48b1f24d7f20980bb461cd8983f	records=348
medline_sha256	base3_literal	16b93b87495f8a1fe869b00279b0e5d3ddecee7e47df69a2e99acd700132dc2e	records=1061
medline_sha256	base3_refinada	c00b8acc67f2cd188b192a12dc48da6e0b17213fa78f06a0bd82111b93453eb0	records=353
medline_sha256	base4_literal	ad896e599c7df8559e608ff73927ea17f0c3af7bc23e4964ec4b9b0d601deded	records=13
medline_sha256	base4_refinada	2769d95c05bdfe1125b857ea7dc122e8221733a37718aa67835873f36a2b4c7d	records=143
medline_sha256	grupo_TLR	47ae98942cc49093db364a95fad91b7468566cdbabba2cfec45e3951bd3149d4	records=31
medline_sha256	grupo_APOBEC3	4e7980b27700a1f6692945575fbaf47d5d89cf31b4c47d9772fda48478030e75	records=135
medline_sha256	grupo_MX_OAS	1ac3dba36bb20136eb5e93355ad227d52373e85666685c103e8133d2fdbc6e12	records=27
medline_sha256	grupo_BST2	8f7d81330f48087ef7350b071ec9661aa9ba6410ec7bb270f05f59a4a3abea73	records=66
medline_sha256	grupo_ACE2	d774c7f9cd093091e87d0fb732c9f8a4fbb465b9e6783a8986090fc51a315e3e	records=103
medline_sha256	grupo_CD4	d2990b642a6af15d45db09e7cbcbfb1bbd2b56ac03eca4528280ee2833d6d224	records=143
medline_sha256	grupo_DDX58_MAVS	eee808d0d2643308924da892a779d4e674c53e430c1cb4161c89ca8440078c1d	records=28
medline_sha256	complementaria_Ateles	a4c2066a660217c936a6655c1811fad7407438e3a4fe02c705c938f8535c24bc	records=142
medline_sha256	plat_TLR_sensores	b51010b57dd370724c46395db3c992237f91b728e9a6f89a5abf15b1181fa8b3	records=21
medline_sha256	plat_restriccion	eaf98290cb3f1f9dbc0c7509375872a482e80161e21545428e2e6a2b4ce0d00a	records=60
medline_sha256	plat_receptores	bc53a8c32506b6b025793325ade102dee647de22e494b59404a4343346b957bd	records=164
pmids_unicos_leidos: 3102
```

---
## d) Tablas de resultados

### Prioridades v2 (todos / sin duplicados)
- P1a: 828 / 803
- P1b: 37 / 37
- P2a: 53 / 53
- P2b: 2 / 2
- P3a: 41 / 41
- P3c_ACE2: 222 / 205
- P3b: 1919 / 1917

### Transición v1 → v2
- P1 → P1a: 828
- P1 → P1b: 37
- P2 → P2a: 53
- P2 → P2b: 2
- P3a → P3a: 41
- P3b → P3b: 1919
- P3b → P3c_ACE2: 222

### nivel_sugerido
- A: 139
- B: 742
- C: 263
- (vacío): 1958
- platirrinos reales en P1a: 139

### Anclas (todas P1a)
| pmid | nombre | prioridad | nivel | plat_real | prim_real | gen | fam |
|---|---|---|---|---|---|---|---|
| 15689398 | Sawyer | P1a | B | false | true | true | false |
| 20410160 | Wlasiuk | P1a | B | false | true | false | true |
| 32826334 | Damas | P1a | B | false | true | true | false |
| 32583896 | Torosin | P1a | A | true | true | true | true |
| 32768453 | Torosin2 | P1a | A | true | true | true | true |

### gen × P1a
- gen_diana: 690
- gen_familia: 138

### P1a por consulta (PMID puede estar en varias)
- base2_refinada: 348
- base4_refinada: 143
- grupo_CD4: 135
- grupo_APOBEC3: 135
- base3_refinada: 125
- grupo_ACE2: 102
- grupo_BST2: 66
- plat_restriccion: 57
- plat_receptores: 52
- base3_literal: 41
- base1_refinada: 30
- grupo_TLR: 30
- grupo_DDX58_MAVS: 28
- base1_literal: 28
- grupo_MX_OAS: 27
- plat_TLR_sensores: 17
- base4_literal: 13
- complementaria_Ateles: 5
