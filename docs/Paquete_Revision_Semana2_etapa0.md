# Paquete de revisión Semana 2 — Etapa 0 (prefiltro)

**Proyecto:** Predicción de susceptibilidad viral en primates neotropicales (Ateles fusciceps)
**Responsable:** Daniel Albuja
**Repo:** https://github.com/evanalbuja-png/primate-pipeline
**Commit script:** 69dfdab
**Fecha de ejecución:** 2026-11-08 / 2026-11-09 UTC
**Estado:** Etapa 0 ejecutada; sin Zotero ni cribado manual

---
## a) Script y regla

### workflow/scripts/cribado_etapa0.py

```python
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

```

### Regla cribado_etapa0 (pubmed.smk)

```
rule cribado_etapa0:
    input:
        medlines = expand(f"{OUTDIR}/{{qid}}/records.medline", qid=config["pubmed"]["fetch"]),
        esearch = expand(f"{OUTDIR}/{{qid}}/esearch.json", qid=config["pubmed"]["fetch"]),
    output:
        tsv = f"{OUTDIR}/cribado_base.tsv",
        audit_p3 = f"{OUTDIR}/auditoria_P3b_100.tsv",
        audit_p1 = f"{OUTDIR}/auditoria_P1_50.tsv",
    conda:
        "../../envs/ateles_pubmed.yml"
    log:
        "logs/pubmed/cribado_etapa0.log"
    shell:
        """
        python workflow/scripts/cribado_etapa0.py \
            --config config/config.yaml \
            --root {OUTDIR} \
            --output {output.tsv} \
            --audit-p3 {output.audit_p3} \
            --audit-p1 {output.audit_p1} \
            > {log} 2>&1
        """

```

---
## b) Salida completa de `--test`

```
OK	plat	expect=True	got=True	'marmosets'	{'hit_platirrino': True, 'hit_primate': True, 'hit_gen_diana': False, 'hit_gen_familia': False, 'hit_seleccion': False, 'hit_virus': False, 'cd4_ambiguo': False}
OK	plat	expect=True	got=True	'howler monkeys'	{'hit_platirrino': True, 'hit_primate': True, 'hit_gen_diana': False, 'hit_gen_familia': False, 'hit_seleccion': False, 'hit_virus': False, 'cd4_ambiguo': False}
OK	plat	expect=True	got=True	'squirrel monkeys'	{'hit_platirrino': True, 'hit_primate': True, 'hit_gen_diana': False, 'hit_gen_familia': False, 'hit_seleccion': False, 'hit_virus': False, 'cd4_ambiguo': False}
OK	plat	expect=True	got=True	'New World monkeys'	{'hit_platirrino': True, 'hit_primate': True, 'hit_gen_diana': False, 'hit_gen_familia': False, 'hit_seleccion': False, 'hit_virus': False, 'cd4_ambiguo': False}
OK	plat	expect=True	got=True	'platyrrhines'	{'hit_platirrino': True, 'hit_primate': True, 'hit_gen_diana': False, 'hit_gen_familia': False, 'hit_seleccion': False, 'hit_virus': False, 'cd4_ambiguo': False}
OK	plat	expect=True	got=True	'Ateles fusciceps'	{'hit_platirrino': True, 'hit_primate': True, 'hit_gen_diana': False, 'hit_gen_familia': False, 'hit_seleccion': False, 'hit_virus': False, 'cd4_ambiguo': False}
OK	plat	expect=False	got=False	'tamarind fruit'	{'hit_platirrino': False, 'hit_primate': False, 'hit_gen_diana': False, 'hit_gen_familia': False, 'hit_seleccion': False, 'hit_virus': False, 'cd4_ambiguo': False}
OK	plat	expect=False	got=False	'mice and rats'	{'hit_platirrino': False, 'hit_primate': False, 'hit_gen_diana': False, 'hit_gen_familia': False, 'hit_seleccion': False, 'hit_virus': False, 'cd4_ambiguo': False}
OK	gen	expect=True	got=True	'TRIM5alpha'	{'hit_platirrino': False, 'hit_primate': False, 'hit_gen_diana': True, 'hit_gen_familia': False, 'hit_seleccion': False, 'hit_virus': False, 'cd4_ambiguo': False}
OK	gen	expect=True	got=True	'TRIM5α'	{'hit_platirrino': False, 'hit_primate': False, 'hit_gen_diana': True, 'hit_gen_familia': False, 'hit_seleccion': False, 'hit_virus': False, 'cd4_ambiguo': False}
OK	gen	expect=True	got=True	'TLR7'	{'hit_platirrino': False, 'hit_primate': False, 'hit_gen_diana': True, 'hit_gen_familia': False, 'hit_seleccion': False, 'hit_virus': False, 'cd4_ambiguo': False}
OK	gen	expect=True	got=True	'RIG-I'	{'hit_platirrino': False, 'hit_primate': False, 'hit_gen_diana': True, 'hit_gen_familia': False, 'hit_seleccion': False, 'hit_virus': False, 'cd4_ambiguo': False}
OK	gen	expect=True	got=True	'ACE2'	{'hit_platirrino': False, 'hit_primate': False, 'hit_gen_diana': True, 'hit_gen_familia': False, 'hit_seleccion': False, 'hit_virus': False, 'cd4_ambiguo': False}
OK	gen	expect=True	got=True	'MxB'	{'hit_platirrino': False, 'hit_primate': False, 'hit_gen_diana': True, 'hit_gen_familia': False, 'hit_seleccion': False, 'hit_virus': False, 'cd4_ambiguo': False}
OK	gen	expect=True	got=True	'BST-2'	{'hit_platirrino': False, 'hit_primate': False, 'hit_gen_diana': True, 'hit_gen_familia': False, 'hit_seleccion': False, 'hit_virus': False, 'cd4_ambiguo': False}
OK	gen	expect=True	got=True	'CD4 receptor binding'	{'hit_platirrino': False, 'hit_primate': False, 'hit_gen_diana': True, 'hit_gen_familia': False, 'hit_seleccion': False, 'hit_virus': False, 'cd4_ambiguo': False}
OK	gen_amb	expect=True	got=True	'CD4+ T cells proliferate'	{'hit_platirrino': False, 'hit_primate': False, 'hit_gen_diana': False, 'hit_gen_familia': False, 'hit_seleccion': False, 'hit_virus': False, 'cd4_ambiguo': True}
OK	fam	expect=True	got=True	'APOBEC3F'	{'hit_platirrino': False, 'hit_primate': False, 'hit_gen_diana': False, 'hit_gen_familia': True, 'hit_seleccion': False, 'hit_virus': False, 'cd4_ambiguo': False}
OK	fam	expect=True	got=True	'Toll-like receptor 5'	{'hit_platirrino': False, 'hit_primate': False, 'hit_gen_diana': False, 'hit_gen_familia': True, 'hit_seleccion': False, 'hit_virus': False, 'cd4_ambiguo': False}
OK	fam	expect=True	got=True	'TLR10'	{'hit_platirrino': False, 'hit_primate': False, 'hit_gen_diana': False, 'hit_gen_familia': True, 'hit_seleccion': False, 'hit_virus': False, 'cd4_ambiguo': False}
OK	sel	expect=True	got=True	'dN/dS'	{'hit_platirrino': False, 'hit_primate': False, 'hit_gen_diana': False, 'hit_gen_familia': False, 'hit_seleccion': True, 'hit_virus': False, 'cd4_ambiguo': False}
OK	sel	expect=True	got=True	'positive selection'	{'hit_platirrino': False, 'hit_primate': False, 'hit_gen_diana': False, 'hit_gen_familia': False, 'hit_seleccion': True, 'hit_virus': False, 'cd4_ambiguo': False}
OK	sel	expect=True	got=True	'positively selected'	{'hit_platirrino': False, 'hit_primate': False, 'hit_gen_diana': False, 'hit_gen_familia': False, 'hit_seleccion': True, 'hit_virus': False, 'cd4_ambiguo': False}
fails=0
```

---
## c) snakemake -n y ejecución real

### Dry-run (resumen)
```
Job stats:
job               count
--------------  -------
cribado_etapa0        1
total                 1
input: 18 records.medline + 18 esearch.json
output: cribado_base.tsv, auditoria_P3b_100.tsv, auditoria_P1_50.tsv
```

### Ejecución real
```
[Sun Nov  8 20:29:05 2026] localrule cribado_etapa0
Finished jobid: 0 (Rule: cribado_etapa0)
1 of 1 steps (100%) done
Elapsed time: 0:00:06.250614
```

---
## d) cribado_meta.txt

```
fecha_utc: 2026-11-09T01:29:10.790459+00:00
commit: 69dfdab
script_sha256: cbd2c3a1b4f602dd11d00f6635d4477bec1a27046eff97bde74a85b6d2a3bb71
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

## e) Resultados cuantitativos

PMIDs unicos totales: 3102
PMIDs sin duplicados (dup_titulo=false): 3058

### Por prioridad (todos / sin duplicados)
- P1: 865 / 840
- P2: 55 / 55
- P3a: 41 / 41
- P3b: 2141 / 2122

### nivel_sugerido
- (vacio): 2141
- A: 186
- B: 734
- C: 41

### gen detectado × prioridad (solo P1+P2)
- gen_diana × P1: 725
- gen_diana × P2: 51
- gen_familia × P1: 140
- gen_familia × P2: 4

hit_platirrino: 399
cd4_ambiguo: 172

### PMIDs por consulta de origen
- base1_literal: 837
- base1_refinada: 59
- base2_refinada: 348
- base3_literal: 1061
- base3_refinada: 353
- base4_literal: 13
- base4_refinada: 143
- complementaria_Ateles: 142
- grupo_ACE2: 103
- grupo_APOBEC3: 135
- grupo_BST2: 66
- grupo_CD4: 143
- grupo_DDX58_MAVS: 28
- grupo_MX_OAS: 27
- grupo_TLR: 31
- plat_TLR_sensores: 21
- plat_receptores: 164
- plat_restriccion: 60

### Prioridad de los 193 PMIDs nuevos de plat_*
- P1: 113
- P2: 0
- P3a: 0
- P3b: 80

### f) Anclas
| pmid | prioridad | nivel | plat | primate | gen_diana | gen_familia | seleccion | virus | consultas |
|---|---|---|---|---|---|---|---|---|---|
| 15689398 | P1 | B | false | true | true | false | true | true | base1_refinada;base4_refinada |
| 20410160 | P1 | B | false | true | false | true | true | true | base1_literal;base1_refinada;grupo_TLR |
| 32826334 | P1 | B | false | true | true | false | false | true | base3_refinada;grupo_ACE2 |
| 32583896 | P1 | A | true | true | true | true | true | true | plat_TLR_sensores |
| 32768453 | P1 | A | true | true | true | true | true | true | grupo_TLR;plat_TLR_sensores |

---
## g) Títulos plat_TLR_sensores, plat_restriccion y muestra plat_receptores

Ver salida literal del Bloque 2 del mensaje al revisor (2026-11-08): listados completos de 21 + 60 títulos y 25 de 154 nuevos en plat_receptores (semilla 20261027). Se reproducen abajo.

### plat_TLR_sensores

```
42923898 | 2026 | Dietary Phytoadditives Enhance Growth Performance, Digestive Enzyme Activity,  Antioxidant Status, Immune Gene Expression, and Vibrio parahaemolyticus Resistance in Pacific White Shrimp (Penaeus vannamei). | NUEVO
42860161 | 2026 | Effect of galactose side-chain modification on immunomodulatory activity of  tamarind xyloglucan oligosaccharides. | NUEVO
42540570 | 2025 | Effects of tamarind pulp powder supplementation on growth performance,  haemato-immunological parameters, and immune-related gene expression in pacific white leg shrimp (Litopenaeus vannamei). | NUEVO
38177528 | 2024 | Highly secreted tryptophanyl tRNA synthetase 1 as a potential theranostic target  for hypercytokinemic severe sepsis. | NUEVO
35704759 | 2022 | Circadian protection against bacterial skin infection by epidermal  CXCL14-mediated innate immunity. | NUEVO
34128045 | 2021 | Innate immunity stimulation via CpG oligodeoxynucleotides ameliorates Alzheimer's  disease pathology in aged squirrel monkeys. | NUEVO
32768453 | 2020 | Comparing the selective landscape of TLR7 and TLR8 across primates reveals unique  sites under positive selection in Alouatta.
32583896 | 2020 | Positively selected variants in functionally important regions of TLR7 in  Alouatta guariba clamitans with yellow fever virus exposure in Northern Argentina. | NUEVO
32194391 | 2020 | Class C CpG Oligodeoxynucleotide Immunomodulatory Response in Aged Squirrel  Monkey (Saimiri Boliviensis Boliviensis).
31352229 | 2019 | Construction of a mavs-inactivated MDCK cell line for facilitating the  propagation of canine distemper virus (CDV). | NUEVO
31201676 | 2019 | Effect of gluten diet on blood innate immune gene expressions and stool  consistency in Spix's Saddleback Tamarin (Leontocebus fuscicollis) raised in captivity. | NUEVO
30555535 | 2018 | Mucoadhesive role of tamarind xyloglucan on inflammation attenuates ulcerative  colitis. | NUEVO
25867333 | 2015 | Sequence characterization and phylogenetic analysis of Toll-like receptor (TLR) 4  gene in the Tibetan macaque (Macaca thibetana).
25329054 | 2014 | Protective efficacy of a Plasmodium vivax circumsporozoite protein-based vaccine  in Aotus nancymaae is associated with antibodies to the repeat region. | NUEVO
24189556 | 2013 | [Development and application of next generation SE36 malaria vaccine formulated  with a novel adjuvant: approach to travelers' vaccine]. | NUEVO
23291928 | 2013 | TLR9 adjuvants enhance immunogenicity and protective efficacy of the SE36/AHG  malaria vaccine in nonhuman primate models. | NUEVO
19726537 | 2009 | Evolutionary trajectories of primate genes involved in HIV pathogenesis.
17093192 | 2007 | GB virus B disrupts RIG-I signaling by NS3/4A-mediated cleavage of the adaptor  protein MAVS.
17071591 | 2006 | Phagocytes containing a disease-promoting Toll-like receptor/Nod ligand are  present in the brain during demyelinating disease in primates. | NUEVO
16912902 | 2006 | Evolution of an intronic microsatellite polymorphism in Toll-like receptor 2  among primates.
15900500 | 2005 | Structural and functional characterisation of the Toll like receptor 9 of Aotus  nancymaae, a non-human primate model for malaria vaccine development. | NUEVO
```

### plat_restriccion

```
41803134 | 2026 | CXCR4-tropic HIV-1 infection in an immunocompetent monkey model. | NUEVO
40526717 | 2025 | Sp140L functions as a herpesvirus restriction factor suppressing viral  transcription and activating interferon-stimulated genes.
40261013 | 2025 | Deciphering lentiviral Vpr/Vpx determinants required for HUSH and SAMHD1  antagonism highlights the molecular plasticity of these evolutionary conflicts.
39713285 | 2025 | Sp140L Is a Novel Herpesvirus Restriction Factor.
32849820 | 2020 | Transcriptome of the Southern Muriqui Brachyteles arachnoides  (Primates:Platyrrhini), a Critically Endangered New World Monkey: Evidence of Adaptive Evolution.
30713099 | 2019 | Recurrent Loss-of-Function Mutations Reveal Costs to OAS1 Antiviral Activity in  Primates.
30170890 | 2018 | Activation of NF-kappaB induced by TRIMCyp showing a discrepancy between owl monkey  and northern pig-tailed macaque. | NUEVO
29187540 | 2018 | General Model for Retroviral Capsid Pattern Recognition by TRIM5 Proteins. | NUEVO
27876849 | 2016 | Characterization of two distinct early post-entry blocks to HIV-1 in common  marmoset lymphocytes.
26764007 | 2016 | TRIM5alpha Degradation via Autophagy Is Not Required for Retroviral Restriction. | NUEVO
26557426 | 2015 | Investigation into the effects of antioxidant-rich extract of Tamarindus indica  leaf on antioxidant enzyme activities, oxidative stress and gene expression profiles in HepG2 cells. | NUEVO
26512082 | 2016 | HIV-1 Adapts To Replicate in Cells Expressing Common Marmoset APOBEC3G and BST2. | NUEVO
25822622 | 2015 | The effect of exon 7 deletion during the evolution of TRIMCyp fusion proteins on  viral restriction, cytoplasmic body formation and multimerization.
25502651 | 2015 | Inhibition of microtubules and dynein rescues human immunodeficiency virus type 1  from owl monkey TRIMCyp-mediated restriction in a cellular context-specific fashion. | NUEVO
24162735 | 2014 | Orthologous mammalian APOBEC3A cytidine deaminases hypermutate nuclear DNA. | NUEVO
24158812 | 2014 | The cyclosporin A washout assay to detect HIV-1 uncoating in infected cells. | NUEVO
23505372 | 2013 | Fates of retroviral core components during unrestricted and TRIM5-restricted  infection. | NUEVO
22574228 | 2012 | Co-evolution of primate SAMHD1 and lentivirus Vpx leads to the loss of the vpx  gene in HIV-1 ancestor.
22345017 | 2012 | [The primate TRIMCyp fusion genes and mechanism of restricting retroviruses  replication].
22291694 | 2012 | TRIM5alpha and Species Tropism of HIV/SIV. | NUEVO
22278248 | 2012 | Herpesvirus saimiri antagonizes nuclear domain 10-instituted intrinsic immunity  via an ORF3-mediated selective degradation of cellular protein Sp100. | NUEVO
21789205 | 2011 | Origin and evolution of TRIM proteins: new insights from the complete TRIM  repertoire of zebrafish and pufferfish. | NUEVO
21247355 | 2011 | Recent insights into the mechanism and consequences of TRIM5alpha retroviral  restriction. | NUEVO
20808866 | 2010 | Conformational adaptation of Asian macaque TRIMCyp directs lineage specific  antiviral activity.
20219908 | 2010 | The specificity of TRIM5 alpha-mediated restriction is influenced by its  coiled-coil domain.
20130055 | 2010 | Species-specific inhibition of foamy viruses from South American monkeys by New  World Monkey TRIM5alpha proteins.
19931648 | 2010 | Evolution of TRIM5alpha B30.2 (SPRY) domain in New World primates.
19741300 | 2009 | Potent inhibition of HIV-1 by TRIM5-cyclophilin fusion proteins engineered from  human components. | NUEVO
19726537 | 2009 | Evolutionary trajectories of primate genes involved in HIV pathogenesis.
19553332 | 2009 | A New World primate deficient in tetherin-mediated restriction of human  immunodeficiency virus type 1.
18799582 | 2008 | Identification of postentry restrictions to Mason-Pfizer monkey virus infection  in New World monkey cells.
18613956 | 2008 | Implication of TRIM alpha and TRIMCyp in interferon-induced anti-retroviral  restriction activities.
18497858 | 2008 | Proteasomal degradation of TRIM5alpha during retrovirus restriction. | NUEVO
18389077 | 2008 | Evolution of a TRIM5-CypA splice isoform in old world monkeys.
18385767 | 2008 | Tissue-specific restriction of cyclophilin A-independent HIV-1- and SIV-derived  lentiviral vectors.
18367529 | 2008 | Restriction of foamy viruses by primate Trim5alpha.
18287035 | 2008 | Independent evolution of an antiviral TRIMCyp in rhesus macaques.
18287034 | 2008 | Independent genesis of chimeric TRIM5-cyclophilin proteins in two primate  species.
18287033 | 2008 | TRIMCyp expression in Old World primates Macaca nemestrina and Macaca  fascicularis.
18172386 | 2007 | A novel fusion gene, TRIM5-Cyclophilin A in the pig-tailed macaque determines its  susceptibility to HIV-1 infection. | NUEVO
18077709 | 2008 | Heat shock perturbs TRIM5alpha restriction of human immunodeficiency virus type  1. | NUEVO
17920096 | 2007 | Comparative requirements for the restriction of retrovirus infection by  TRIM5alpha and TRIMCyp.
17609268 | 2007 | Fusion of cyclophilin A to Fv1 enables cyclosporine-sensitive restriction of  human and feline immunodeficiency viruses. | NUEVO
17574642 | 2007 | The ability of multimerized cyclophilin A to restrict retrovirus infection.
16973579 | 2006 | Proteasome inhibition reveals that a functional preintegration complex  intermediate can be generated during restriction by diverse TRIM5 proteins.
16828831 | 2006 | Antiretroviral potential of human tripartite motif-5 and related proteins. | NUEVO
16650449 | 2006 | Requirements for capsid-binding and an effector function in TRIMCyp-mediated  restriction of HIV-1. | NUEVO
16641261 | 2006 | Cyclophilin A renders human immunodeficiency virus type 1 sensitive to Old World  monkey but not human TRIM5 alpha antiviral activity. | NUEVO
16571822 | 2006 | Trim-cyclophilin A fusion proteins can restrict human immunodeficiency virus type  1 infection at two distinct phases in the viral life cycle.
16439561 | 2006 | Arsenic counteracts human immunodeficiency virus type 1 restriction by various  TRIM5 orthologues in a cell type-dependent manner.
16306627 | 2005 | Restriction of human immunodeficiency virus type 1 by TRIM-CypA occurs with rapid  kinetics and independently of cytoplasmic bodies, ubiquitin, and proteasome activity. | NUEVO
16282502 | 2005 | Evolution of cyclophilin A and TRIMCyp retrotransposition in New World primates.
16203999 | 2005 | Cyclophilin A is required for TRIM5alpha-mediated resistance to HIV-1 in Old  World monkey cells. | NUEVO
16199531 | 2005 | Naturally occurring capsid substitutions render HIV-1 cyclophilin A independent  in human cells and TRIM-cyclophilin-resistant in Owl monkey cells. | NUEVO
16183097 | 2005 | TRIM5alpha association with cytoplasmic bodies is not required for antiretroviral  activity.
16140735 | 2005 | Differential restriction of human immunodeficiency virus type 2 and simian  immunodeficiency virus SIVmac by TRIM5alpha alleles.
15919943 | 2005 | Disruption of human TRIM5alpha antiviral activity by nonhuman primate  orthologues.
15767395 | 2005 | Retrovirus restriction by TRIM5alpha variants from Old World and New World  primates.
15326303 | 2004 | A Trim5-cyclophilin A fusion protein found in owl monkey kidney cells can  restrict HIV-1. | NUEVO
15243629 | 2004 | Cyclophilin A retrotransposition into TRIM5 explains owl monkey resistance to  HIV-1.
```

### plat_receptores nuevos (muestra 25, semilla 20261027)

```
nuevos = 154
20448482 | 2010 | Induction of progressive demyelinating autoimmune encephalomyelitis in common  marmoset monkeys using MOG34-56 peptide in incomplete freund adjuvant.
35909394 | 2022 | Chronic upper airway inflammation related to high Th2 cytokines in Mendelian  susceptibility to mycobacterial disease case.
22593679 | 2012 | Antitumor and immunopotentiating activity of polysaccharide PST001 isolated from  the seed kernel of Tamarindus indica: an in vivo study in mice.
17142786 | 2006 | Familial NK cell deficiency associated with impaired IL-2- and IL-15-dependent  survival of lymphocytes.
2403285 | 1990 | Polyhalogenated dibenzo-p-dioxins and dibenzofurans and the immune system. 1.  Effects on peripheral lymphocyte subpopulations of a non-human primate (Callithrix jacchus) after treatment with 2,3,7,8-tetrachlorodibenzo-p-dioxin (TCDD).
9064482 | 1997 | Comparative study on age-dependent development of surface receptors on peripheral  blood lymphocytes in children and young nonhuman primates (marmosets).
26063421 | 2015 | Identification of Owl Monkey CD4 Receptors Broadly Compatible with Early-Stage  HIV-1 Isolates.
15661043 | 2005 | Herpesvirus saimiri immortalization of Aotus T lymphocytes specific for an  immunogenically modified peptide of Plasmodium falciparum merozoite surface antigen 2.
39226323 | 2024 | A nonhuman primate model for genital herpes simplex virus 2 infection that  results in vaginal vesicular lesions, virus shedding, and seroconversion.
38990001 | 2024 | Genomic analysis and identification of a novel superantigen, SargEY, in  Staphylococcus argenteus isolated from atopic dermatitis lesions.
7525717 | 1994 | Polyclonal B cell activation induced by herpesvirus saimiri-transformed human  CD4+ T cell clones. Role for membrane TNF-alpha/TNF-alpha receptors and CD2/CD58 interactions.
8892618 | 1996 | Herpesvirus saimiri-transformed human CD4+ T cells can provide polyclonal B cell  help via the CD40 ligand as well as the TNF-alpha pathway and through release of lymphokines.
1626863 | 1992 | Immunophenotypic characterization of owl monkey peripheral blood mononuclear  cells.
37112786 | 2023 | Immunogenicity and In Vivo Protective Effects of Recombinant Nucleocapsid-Based  SARS-CoV-2 Vaccine Convacell((R)).
30713763 | 2018 | Assessment of the Parameters of Adaptive Cell-Mediated Immunity in Naive Common  Marmosets (Callithrix jacchus).
8676080 | 1996 | T cell interleukin-17 induces stromal cells to produce proinflammatory and  hematopoietic cytokines.
15234557 | 2004 | Characterization of a reduced peptide bond analogue of a promiscuous CD4 T cell  epitope derived from the Plasmodium falciparum malaria vaccine candidate merozoite surface protein 1.
1578139 | 1992 | IL-12 receptor. II. Distribution and regulation of receptor expression.
16540180 | 2006 | The efficacy of combined therapy of arsenic trioxide and alpha interferon in  human T-cell leukemia virus type-1-infected squirrel monkeys (Saimiri sciureus).
12196293 | 2002 | Herpesviral protein targets a cellular WD repeat endosomal protein to  downregulate T lymphocyte receptor expression.
9140094 | 1997 | Interactions between autologous CD4+ and CD8+ T lymphocytes and human squamous  cell carcinoma of the head and neck.
18209026 | 2008 | Fast progression of recombinant human myelin/oligodendrocyte glycoprotein  (MOG)-induced experimental autoimmune encephalomyelitis in marmosets is associated with the activation of MOG34-56-specific cytotoxic T cells.
16797207 | 2006 | Immunogenicity of a recombinant protein containing the Plasmodium vivax vaccine  candidate MSP1(19) and two human CD4+ T-cell epitopes administered to non-human primates (Callithrix jacchus jacchus).
7836781 | 1995 | Selection of monoclonal antibodies for the identification of lymphocyte surface  antigens in the New World primate Saguinus oedipus oedipus (cotton top tamarin).
9743383 | 1998 | Signaling through a CD3 gamma-deficient TCR/CD3 complex in immortalized mature  CD4+ and CD8+ T lymphocytes.
```

---
## h) Registro Torosin et al. 2020 (PMID 32583896)

```
PMID- 32583896
DP  - 2020 Sep
TI  - Positively selected variants in functionally important regions of TLR7 in 
      Alouatta guariba clamitans with yellow fever virus exposure in Northern 
      Argentina.
LID - 10.1002/ajpa.24086 [doi]
JT  - American journal of physical anthropology
AID - 10.1002/ajpa.24086 [doi]
```

---
## i) Decisiones y problemas

- 5/5 anclas caen en P1 de forma natural (sin forzar).
- P1 = 865 (840 sin duplicados): cribado manual en dos pasadas recomendable.
- 193 PMIDs nuevos de plat_*: 113 en P1, 80 en P3b.
- Genoma Ateles (38244218): no está en el conjunto de 18 consultas; PRJNA1009451 no aparece en campo SI del MEDLINE.
- Sironi 2015 y Enard 2016 no están en el conjunto; se añadirán a la matriz con fuente_busqueda=propuesta.
- No se importó a Zotero ni se hizo cribado manual.
