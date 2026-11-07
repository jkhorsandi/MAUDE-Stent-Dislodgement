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



# ===========================================================================
# 6. MANTEL-HAENSZEL ADD-ON (early 0-3 vs late 4+ device-years, stratified by platform)
#    Everything above is your original pipeline, unchanged. Device-year uses your
#    definition: int((year + month/12) - first_seen_year_month_of_that_platform).
# ===========================================================================
import csv
import os

EARLY_MAX = 3
CACHE_FILE = "maude_raw.json"
TARGET = {"ONYX FRONTIER": (440, 852), "ORSIRO": (570, 1643), "XIENCE SKYPOINT": (355, 1497),
          "RESOLUTE ONYX": (1110, 4878), "ORSIRO MISSION": (177, 797), "SYNERGY": (472, 2675),
          "XIENCE SIERRA": (361, 2211), "SYNERGY MEGATRON": (42, 395)}


def download_cached():
    if os.path.exists(CACHE_FILE):
        return json.load(open(CACHE_FILE))
    raw = download()
    json.dump(raw, open(CACHE_FILE, "w"))
    return raw


def chi2_sf(x, df):
    """Upper-tail chi-square probability, standard library only (closed forms, df >= 1)."""
    if x <= 0:
        return 1.0
    if df % 2 == 0:
        term, s = 1.0, 0.0
        for i in range(df // 2):
            if i:
                term *= (x / 2) / i
            s += term
        return math.exp(-x / 2) * s
    acc, term = 0.0, math.sqrt(x)
    for r in range(1, (df - 1) // 2 + 1):
        if r > 1:
            term *= x / (2 * r - 1)
        acc += term
    return math.erfc(math.sqrt(x / 2)) + 2 * math.exp(-x / 2) / math.sqrt(2 * math.pi) * acc


def woolf_or(a, b, c, d):
    if min(a, b, c, d) == 0:
        a, b, c, d = a + .5, b + .5, c + .5, d + .5
    o = a * d / (b * c)
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return o, math.exp(math.log(o) - 1.96 * se), math.exp(math.log(o) + 1.96 * se)


def mantel_haenszel(strata):
    """strata: list of (name, a, b, c, d) with a,b = early event/non-event; c,d = late."""
    u, dropped = [], []
    for p, a, b, c, d in strata:
        if min(a + b, c + d, a + c, b + d) == 0:
            dropped.append(p)
        else:
            u.append((a, b, c, d, a + b + c + d))
    Rs = sum(a * d / n for a, b, c, d, n in u)
    Ss = sum(b * c / n for a, b, c, d, n in u)
    o = Rs / Ss
    v = (sum((a + d) / n * a * d / n for a, b, c, d, n in u) / (2 * Rs ** 2)
         + sum(((a + d) / n * b * c / n + (b + c) / n * a * d / n) for a, b, c, d, n in u) / (2 * Rs * Ss)
         + sum((b + c) / n * b * c / n for a, b, c, d, n in u) / (2 * Ss ** 2))
    lo, hi = math.exp(math.log(o) - 1.96 * math.sqrt(v)), math.exp(math.log(o) + 1.96 * math.sqrt(v))
    E = sum((a + b) * (a + c) / n for a, b, c, d, n in u)
    V = sum((a + b) * (c + d) * (a + c) * (b + d) / (n ** 2 * (n - 1)) for a, b, c, d, n in u)
    chi = (abs(sum(x[0] for x in u) - E) - .5) ** 2 / V
    ea, va = [], []                                   # Breslow-Day with Tarone adjustment
    for a, b, c, d, n in u:
        n1, n0, m1 = a + b, c + d, a + c
        A_, B_, C_ = o - 1, -(o * (n1 + m1) + (n0 - m1)), o * n1 * m1
        if abs(A_) < 1e-12:
            roots = [-C_ / B_]
        else:
            s = math.sqrt(B_ * B_ - 4 * A_ * C_)
            roots = [(-B_ - s) / (2 * A_), (-B_ + s) / (2 * A_)]
        x = [z for z in roots if max(0, m1 - n0) - 1e-9 <= z <= min(n1, m1) + 1e-9][0]
        ea.append(x)
        va.append(1 / (1 / x + 1 / (n1 - x) + 1 / (m1 - x) + 1 / (n0 - m1 + x)))
    bd = (sum((x[0] - e) ** 2 / w for x, e, w in zip(u, ea, va))
          - (sum(x[0] for x in u) - sum(ea)) ** 2 / sum(va))
    k = len(u)
    return dict(o=o, lo=lo, hi=hi, chi=chi, p=chi2_sf(chi, 1), bd=bd,
                bdp=chi2_sf(bd, k - 1) if k > 1 else float("nan"), k=k, dropped=dropped)


def fp(p):
    return "NA" if p != p else ("< .001" if p < .001 else f"= {p:.3f}")


def build_strata(rows, cut):
    st = collections.defaultdict(lambda: [0, 0, 0, 0])
    for fam, age, flag in rows:
        i = (0 if age <= cut else 2) + (0 if flag else 1)
        st[fam][i] += 1
    return [(fam, *v) for fam, v in sorted(st.items())]


def main_mh():
    raw = download_cached()
    unique = deduplicate(raw)

    first_seen = {}
    for fam in {platform_of(r["brand"]) for r in unique} - {None}:
        dates = sorted(r["date_received"] for r in unique
                       if platform_of(r["brand"]) == fam and r["date_received"])
        first_seen[fam] = int(dates[0][:4]) + int(dates[0][4:6]) / 12

    rows, n_cohort = [], 0
    for r in unique:
        fam = platform_of(r["brand"])
        flag = is_dislodgement(r["narrative"])
        n_cohort += flag
        if fam and r["date_received"]:
            age = int((int(r["date_received"][:4]) + int(r["date_received"][4:6]) / 12) - first_seen[fam])
            rows.append((fam, age, flag))
    with open("reports_mh.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["platform", "device_year", "dislodgement"]); w.writerows(rows)

    print(f"=== RECONCILIATION vs manuscript ===\nretrieved {len(raw):,} (16,328) | unique {len(unique):,} "
          f"(14,948) | cohort {n_cohort:,} (3,527)")
    print(f"{'platform':<18}{'here n/N':>14}{'manuscript':>14}")
    for fam, (tn, tN) in TARGET.items():
        sub = [x for x in rows if x[0] == fam]
        print(f"{fam:<18}{f'{sum(x[2] for x in sub)}/{len(sub)}':>14}{f'{tn}/{tN}':>14}")

    st = build_strata(rows, EARLY_MAX)
    A, B, C, D = (sum(s[i] for s in st) for i in (1, 2, 3, 4))
    co, cl, ch = woolf_or(A, B, C, D)
    print(f"\n=== PRIMARY: early (0-{EARLY_MAX}) vs late ({EARLY_MAX + 1}+) device-years ===")
    print(f"Crude pooled: early {A}/{A+B} ({100*A/(A+B):.1f}%) vs late {C}/{C+D} ({100*C/(C+D):.1f}%) "
          f"[manuscript: 2,291/8,809 vs 1,236/6,139]; OR {co:.2f} ({cl:.2f}-{ch:.2f})")
    for p, a, b, c, d in st:
        s = ("not estimable (no late period)" if c + d == 0 else
             "OR %.2f (%.2f-%.2f)" % woolf_or(a, b, c, d))
        print(f"  {p:<18} early {a}/{a+b:<6} late {c}/{c+d:<6} {s}")
    m = mantel_haenszel(st)
    print(f"\nMANTEL-HAENSZEL OR (early vs late), {m['k']} strata: {m['o']:.2f} "
          f"(95% CI {m['lo']:.2f}-{m['hi']:.2f}); CMH chi-square {m['chi']:.1f}, P {fp(m['p'])}")
    print(f"Breslow-Day (Tarone) homogeneity: chi-square {m['bd']:.1f}, df {m['k']-1}, P {fp(m['bdp'])}")
    print("Dropped (no early or no late reports):", ", ".join(m["dropped"]) or "none")

    ex = [s for s in st if s[0] != "ONYX FRONTIER"]
    ea_, eb_ = sum(s[1] for s in ex), sum(s[2] for s in ex)
    print(f"\nPooled early share incl. Onyx Frontier {100*A/(A+B):.1f}% | excl. {100*ea_/(ea_+eb_):.1f}% "
          f"| late {100*C/(C+D):.1f}%")

    big = max(st, key=lambda s: sum(s[1:]))[0]
    mb = mantel_haenszel([s for s in st if s[0] != big])
    print(f"\nSENSITIVITY excluding {big}: MH OR {mb['o']:.2f} ({mb['lo']:.2f}-{mb['hi']:.2f}), P {fp(mb['p'])}")
    for cut in (2, 4):
        mc = mantel_haenszel(build_strata(rows, cut))
        print(f"SENSITIVITY early = 0-{cut}: MH OR {mc['o']:.2f} ({mc['lo']:.2f}-{mc['hi']:.2f}), P {fp(mc['p'])}")
    for p, a, b, c, d in st:
        if p == "RESOLUTE ONYX":
            o, lo, hi = woolf_or(a, b, c, d)
            print(f"RESOLUTE ONYX only: early {a}/{a+b} ({100*a/(a+b):.1f}%) vs late {c}/{c+d} "
                  f"({100*c/(c+d):.1f}%); OR {o:.2f} ({lo:.2f}-{hi:.2f})")
    with open("mh_strata_counts.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["platform", "early_event", "early_nonevent", "late_event", "late_nonevent"])
        w.writerows(st)
    print("\nSaved: maude_raw.json (cache), reports_mh.csv, mh_strata_counts.csv")


if __name__ == "__main__":
    main_mh()
