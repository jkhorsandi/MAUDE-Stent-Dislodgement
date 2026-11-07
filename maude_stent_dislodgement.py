#!/usr/bin/env python3
"""
Reproducible extraction and analysis pipeline for:

  Stent Dislodgement With Contemporary Drug-Eluting Stents:
  A Reproducible Extraction Framework and Quantified Misclassification
  in the FDA MAUDE Database

Running this script end to end reproduces every count, proportion, and
confidence interval reported in the manuscript, subject only to openFDA
snapshot drift. The query was executed 15 July 2026 against a snapshot
last updated 21 July 2026; later snapshots accrue reports retrospectively
and will yield slightly larger counts.

Dependencies: Python 3.9+, standard library only.
Usage:        python3 maude_stent_dislodgement.py
"""

import collections
import hashlib
import json
import math
import re
import time
import urllib.request

BASE = "https://api.fda.gov/device/event.json"

# ---------------------------------------------------------------------------
# 1. SEARCH DEFINITION
# ---------------------------------------------------------------------------
# Seven brand-name expressions cover the eight named platforms. The ORSIRO
# expression also returns ORSIRO MISSION records; the RESOLUTE ONYX
# expression also returns RX and OTW delivery-system variants.
BRANDS = [
    "SYNERGY XD", "SYNERGY MEGATRON", "RESOLUTE ONYX", "ONYX FRONTIER",
    "XIENCE SIERRA", "XIENCE SKYPOINT", "ORSIRO",
]
DATE_FILTER = "date_received:[19910101+TO+20260601]"   # cohort entry field
SEARCH = ("(" + "+OR+".join(f'device.brand_name:"{b.replace(" ", "+")}"'
                            for b in BRANDS) + ")" + "+AND+" + DATE_FILTER)

# ---------------------------------------------------------------------------
# 2. CLASSIFICATION RULES (all prespecified)
# ---------------------------------------------------------------------------
DISLODGE = re.compile(r"\bdislodg\w*", re.I)

# Loss terms must fall within 80 characters of "stent" in either direction.
PROXIMITY = re.compile(
    r"\bstent\b[^.]{0,80}?\b(embolis\w*|emboliz\w*|detach\w*|migrat\w*|"
    r"separat\w*|came off|fell off|not present|missing|lost|loss)\b"
    r"|\b(embolis\w*|emboliz\w*|detach\w*|migrat\w*|separat\w*|lost|loss)\b"
    r"[^.]{0,80}?\bstent\b", re.I)

# Exclusions applied only when no independent "dislodg*" match is present.
IFU_BOILERPLATE = re.compile(
    r"listed in the[^.]{0,120}instructions for use"
    r"|known patient effects of coronary", re.I)
NON_STENT_COMPONENT = re.compile(
    r"\b(hub|tip|shaft|hypotube|marker band)\b[^.]{0,25}\bdetach", re.I)

# Post hoc refinement used only as a sensitivity analysis.
MODAL = re.compile(
    r"\b(may|might|could|can|reduces the chance of|risk of)\b[^.]{0,45}"
    r"\b(dislodg|damag\w+ or dislodg)", re.I)
COMPONENT_ANY = re.compile(
    r"\b(hub|shaft|tip|hypotube|balloon|guide ?wire|marker band|inflation lumen)\b"
    r"[^.]{0,30}\b(separat\w*|detach\w*|dislodg\w*)", re.I)

VARIABLES = {
    "anatomic":       r"\btortuos\w*|\btortuous\b|\bcalcifi\w*|\bcalcific\b|\bangulat\w*",
    "calcification":  r"\bcalcifi\w*|\bcalcific\b",
    "tortuosity":     r"\btortuos\w*|\btortuous\b",
    "severe_calc":    r"(?:severe|severely|heavy|heavily)(?:ly)? calcifi",
    "retrieval":      r"\bsnar\w*|\bretriev\w*|\bextract\w*|\bgoose ?neck\b",
    "mi_thrombosis":  r"\bmyocardial infarction\b|\bthromb\w*|\bocclu\w*|\bSTEMI\b|\bNSTEMI\b",
    "perforation":    r"\bperforat\w*|\bdissect\w*|\btamponade\b|\bpericardial effusion\b",
    "surgery":        (r"\bCABG\b|\bbypass graft\w*\b|\bcoronary artery bypass\b"
                       r"|\bemergen\w+ surg\w*|\bsurgical(?:ly)? (?:remov|retriev|extract)\w*"
                       r"|\bsternotomy\b"),
    "ex_vivo":        (r"outside (?:of )?the (?:patient|patients? (?:body|anatomy)|body)"
                       r"|no patient (?:involvement|contact|was involved)"
                       r"|(?:device|stent|system) was (?:therefore )?not used"
                       r"|before introduction into|during unpack\w+"
                       r"|while (?:unpacking|removing the stent protector)"),
    "mech_cross":     r"(?:could not|failed to|unable to) (?:cross|pass|advance|be advanced|be delivered)",
    "mech_balloon":   (r"balloon (?:rupture|burst|ruptured)"
                       r"|(?:came off|slipped off|separated from|fell off)[^.]{0,25}(?:the )?balloon"
                       r"|off the balloon"),
    "mech_prior":     r"previously (?:deployed|implanted|placed) stent|mesh of the previously|entangled[^.]{0,40}stent",
    "mech_guide":     (r"(?:dislodg\w*|came off|stuck|caught|remained)[^.]{0,60}\b(?:in|into|within|inside) "
                       r"the (?:guide|guiding) (?:catheter|extension)"
                       r"|\b(?:guide|guiding) (?:catheter|extension)[^.]{0,60}(?:dislodg\w*|came off|stuck|caught)"),
    "disp_retrieved": (r"(?:was|were|successfully)\s+(?:retrieved|snared|extracted)"
                       r"|retrieved (?:with|using|via|by)|snare[^.]{0,40}(?:remov|retriev|extract)"),
    "disp_crushed":   r"\bcrush\w+|adapted to the vessel wall|apposed[^.]{0,25}vessel wall|\bjailed\b(?! balloon)",
    "disp_left":      (r"(?:stent|it) was not (?:removed|retrieved)"
                       r"|(?:stent|it) (?:was |were )?left in (?:the )?(?:patient|vessel|artery|lesion)"
                       r"|no attempts? (?:was|were) made to remove"),
    "disp_periph":    (r"emboliz\w+[^.]{0,40}\b(?:iliac|femoral|peroneal|profunda|brachial|aorta|renal)"
                       r"|\b(?:iliac|peroneal|profunda|abdominal aorta)\b[^.]{0,40}(?:stent|dislodg|emboliz)"
                       r"|below the renal"),
    "disp_surgical":  r"remov\w+[^.]{0,35}(?:during|at|via) (?:CABG|surgery|sternotomy)|surgical(?:ly)? (?:remov|retriev|extract)\w*",
    "vessel_named":   r"\bleft anterior descending\b|\bLAD\b|\bright coronary\b|\bRCA\b|\bcircumflex\b|\bLCX\b|\bleft main\b|\bSVG\b",
    "left_main":      r"\bleft main\b|\bLMCA\b",
    "predilatation":  r"\bpre-?dilat\w+",
    "calcium_device": r"\brotational atherectomy\b|\brotablat\w*|\bburr\b|\borbital atherectomy\b|\blithotripsy\b|\bshockwave\b",
    "force_used":     r"excessive force was (?:used|applied)|force was applied",
    "force_field":    r"excessive force",
    "no_harm":        (r"no (?:adverse |further )?(?:reported )?(?:adverse )?patient "
                       r"(?:effect|injury|sequela|harm|consequence)|no further patient injury"),
    "death_terms":    r"\bdied\b|\bdeath\b|\bexpired\b|\bfatal\w*|\bdeceased\b|\bpassed away\b",
}

# Negation guards. A record matching the guard is not counted for that variable.
NEGATION = {
    "mech_prior": r"(?:did not|didn't|not) pass through a[^.]{0,20}previously",
    "disp_left":  r"nothing was left in",
}

PLATFORM_ORDER = ["SYNERGY MEGATRON", "SYNERGY", "ONYX FRONTIER", "RESOLUTE ONYX",
                  "XIENCE SIERRA", "XIENCE SKYPOINT", "ORSIRO MISSION", "ORSIRO"]


# ---------------------------------------------------------------------------
# 3. HELPERS
# ---------------------------------------------------------------------------
def wilson(k, n, z=1.96):
    """Wilson score interval, returned as percentages."""
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (round(100 * (centre - half), 1), round(100 * (centre + half), 1))


def prr(a, n_platform, a_all, n_all):
    """Proportional reporting ratio with log-normal 95% CI and Yates chi-square."""
    b = n_platform - a
    c = a_all - a
    d = (n_all - n_platform) - c
    ratio = (a / (a + b)) / (c / (c + d))
    se = math.sqrt(1 / a - 1 / (a + b) + 1 / c - 1 / (c + d))
    total = a + b + c + d
    chi = (total * (abs(a * d - b * c) - total / 2) ** 2
           / ((a + b) * (c + d) * (a + c) * (b + d)))
    return (ratio,
            math.exp(math.log(ratio) - 1.96 * se),
            math.exp(math.log(ratio) + 1.96 * se),
            chi,
            ratio >= 2 and chi >= 4 and a >= 3)          # Evans criteria


def normalize(text):
    """Normalization used for narrative-level deduplication."""
    return re.sub(r"\W+", " ", text.upper()).strip()


def platform_of(brand):
    b = (brand or "").upper()
    for fam in PLATFORM_ORDER:
        if fam in b:
            return fam
    return None


# ---------------------------------------------------------------------------
# 4. DOWNLOAD
# ---------------------------------------------------------------------------
def download():
    records, skip = [], 0
    while skip <= 20000:
        url = f"{BASE}?search={SEARCH}&limit=100&skip={skip}"
        payload = None
        for _ in range(5):
            try:
                with urllib.request.urlopen(url, timeout=120) as resp:
                    payload = json.loads(resp.read().decode())
                break
            except urllib.error.HTTPError as exc:
                if exc.code == 404:
                    payload = {"results": []}
                    break
                time.sleep(4)
            except Exception:
                time.sleep(4)
        results = (payload or {}).get("results", [])
        if not results:
            break
        for rec in results:
            device = (rec.get("device") or [{}])[0]
            records.append({
                "narrative": " ".join(t.get("text", "")
                                      for t in (rec.get("mdr_text") or [])
                                      if t.get("text")),
                "event_type": rec.get("event_type"),
                "date_received": rec.get("date_received"),
                "date_of_event": rec.get("date_of_event"),
                "brand": device.get("brand_name", ""),
                "product_problems": rec.get("product_problems") or [],
                "source_type": rec.get("source_type") or [],
                "country": rec.get("reporter_country_code"),
            })
        skip += 100
    return records


# ---------------------------------------------------------------------------
# 5. PIPELINE
# ---------------------------------------------------------------------------
def deduplicate(records):
    """Collapse records sharing an identical normalized narrative."""
    groups = collections.defaultdict(list)
    for rec in records:
        key = hashlib.md5(normalize(rec["narrative"])[:1500].encode()).hexdigest()
        groups[key].append(rec)
    return [max(grp, key=lambda r: len(r["narrative"])) for grp in groups.values()]


def is_dislodgement(text, refined=False):
    """Primary classification rule. refined=True applies the sensitivity variant."""
    if refined:
        hits = [m for m in DISLODGE.finditer(text)
                if not MODAL.search(text[max(0, m.start() - 70):m.end() + 25])
                and not COMPONENT_ANY.search(text[max(0, m.start() - 45):m.end() + 25])]
        if hits:
            return True
        return bool(PROXIMITY.search(text)) and not COMPONENT_ANY.search(text)

    if not (DISLODGE.search(text) or PROXIMITY.search(text)):
        return False
    if NON_STENT_COMPONENT.search(text) and not DISLODGE.search(text):
        return False
    if IFU_BOILERPLATE.search(text) and not DISLODGE.search(text):
        return False
    return True


def count(cohort, name):
    pattern = re.compile(VARIABLES[name], re.I)
    guard = re.compile(NEGATION[name], re.I) if name in NEGATION else None
    return sum(1 for r in cohort
               if pattern.search(r["narrative"])
               and not (guard and guard.search(r["narrative"])))


def main():
    print("Downloading ...")
    raw = download()
    unique = deduplicate(raw)
    cohort = [r for r in unique if is_dislodgement(r["narrative"])]
    refined = [r for r in unique if is_dislodgement(r["narrative"], refined=True)]
    ex_vivo = re.compile(VARIABLES["ex_vivo"], re.I)
    in_patient = [r for r in cohort if not ex_vivo.search(r["narrative"])]

    n, m = len(cohort), len(in_patient)
    print(f"\nretrieved                 {len(raw):,}")
    print(f"redundant narratives      {len(raw) - len(unique):,} "
          f"({100 * (len(raw) - len(unique)) / len(raw):.1f}%)")
    print(f"unique reports            {len(unique):,}")
    print(f"analytic cohort           {n:,} ({100 * n / len(unique):.1f}%)")
    print(f"refined rule (sensitivity){len(refined):,}")
    print(f"in-patient subgroup       {m:,}\n")

    events = collections.Counter(r["event_type"] for r in cohort)
    for key in ("Malfunction", "Injury", "Death", "Other"):
        print(f"  {key:<14}{events.get(key, 0):>6,} {100 * events.get(key, 0) / n:5.1f}%")

    print()
    for name in VARIABLES:
        if name in ("ex_vivo", "death_terms", "force_field"):
            continue
        k = count(cohort, name)
        print(f"  {name:<16}{k:>6,} {100 * k / n:5.1f}%  95% CI {wilson(k, n)}")

    print("\nProportional reporting ratio (Evans criteria)")
    totals = collections.Counter()
    events_by = collections.Counter()
    for rec in unique:
        fam = platform_of(rec["brand"])
        if fam:
            totals[fam] += 1
    for rec in cohort:
        fam = platform_of(rec["brand"])
        if fam:
            events_by[fam] += 1
    n_all, a_all = sum(totals.values()), sum(events_by.values())
    for fam in sorted(totals, key=lambda f: -events_by[f] / totals[f]):
        ratio, lo, hi, chi, signal = prr(events_by[fam], totals[fam], a_all, n_all)
        print(f"  {fam:<18}{ratio:5.2f} ({lo:.2f}-{hi:.2f})  chi2 {chi:7.1f}  "
              f"{'SIGNAL' if signal else '-'}")

    print("\nReporting life cycle (dislodgement share by device-year)")
    first_seen = {}
    for fam in totals:
        dates = sorted(r["date_received"] for r in unique
                       if platform_of(r["brand"]) == fam and r["date_received"])
        first_seen[fam] = int(dates[0][:4]) + int(dates[0][4:6]) / 12
    buckets = collections.defaultdict(lambda: [0, 0])
    cohort_ids = {id(r) for r in cohort}
    for rec in unique:
        fam = platform_of(rec["brand"])
        if not fam or not rec["date_received"]:
            continue
        age = int((int(rec["date_received"][:4])
                   + int(rec["date_received"][4:6]) / 12) - first_seen[fam])
        bucket = min(age, 6)
        buckets[bucket][1] += 1
        if id(rec) in cohort_ids:
            buckets[bucket][0] += 1
    for bucket in sorted(buckets):
        a, b = buckets[bucket]
        label = f"year {bucket}" if bucket < 6 else "year 6+"
        print(f"  {label:<10}{a:>6,}/{b:<7,} {100 * a / b:5.1f}%")


if __name__ == "__main__":
    main()
