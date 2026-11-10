#!/usr/bin/env python
"""
1x2 figure: de-novo Jolma-motif recovery WITHOUT vs WITH Jolma-analog filters, per method.
 Left  (A): raw recall/precision sweeps (from denovo_jolma_recovery_sweep.tsv): each method's
            dedup sweep (vsearch/STREME self-seg, STREME anchor) + STREME-only N-motif sweep.
 Right (B): apply Jolma-analog filters (low-complexity removal + prevalence top-N = primary/
            secondary cap) to each DCA arm's dedup0.8 motif set; plot the baseline->top5->top3->top2
            trajectory per method.
"""
import os, re, sys, contextlib, io
from collections import Counter
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
import os, sys
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)                                   # nolp_nussinov.py
sys.path.insert(0, os.path.dirname(_HERE))                  # the raptcouple src/ package
from src.helper_jolma import load_pwms, calc_entropy_on_PWM
from src.plmc import read_params

DS = "revision/query_free/JZM2020"; SELEXD = "SELEXdata/JZM2020"
XLSX = "公開用/example/Jolma2020/data/others/Supplemental_Tables_S1-S8.xlsx"
ORD = ["A", "C", "G", "U"]
THR, MIN_LEN, MAX_GAP, MIN_OVL, DEDUP, MATCH, MAX_PAD = 1.5, 5, 1, 4, 0.8, 0.8, 3
CORE = 6; LOWCOMP_HOMO, LOWCOMP_DINUC = 0.70, 0.80; TOPN = [5, 3, 2]


def norm(p):
    p = p.copy(); s = p.sum(0); s[s == 0] = 1; return p / s
def cosine(p, q):
    d = np.linalg.norm(p, axis=0) * np.linalg.norm(q, axis=0); d[d == 0] = np.nan
    return (p * q).sum(0) / d
def maxcos(a, b, pad=3):
    if a is None or b is None or a.shape[1] == 0 or b.shape[1] == 0: return np.nan
    if a.shape[1] < b.shape[1]: a, b = b, a
    Lp, Lq = a.shape[1], b.shape[1]; mp = min(pad, Lq // 2); ys = []
    for i in range(-mp, Lp - Lq + mp + 1):
        if i < 0:
            ov = Lq + i
            if ov < MIN_OVL: continue
            ys.append(np.nanmean(cosine(a[:, :ov], b[:, -i:])))
        elif i > Lp - Lq:
            cut = i - (Lp - Lq); ov = Lq - cut
            if ov < MIN_OVL: continue
            ys.append(np.nanmean(cosine(a[:, i:], b[:, :ov])))
        else:
            ys.append(np.nanmean(cosine(a[:, i:i + Lq], b)))
    return max(ys) if ys else np.nan
def best_offset(pwm, pf):
    width = pwm.shape[1]; length = pf.shape[1]; mp = min(MAX_PAD, width // 2); best = (-np.inf, None)
    for off in range(-mp, length - width + mp + 1):
        if off < 0:
            ov = width + off
            if ov <= 0: continue
            sc = np.nanmean(cosine(pf[:, :ov], pwm[:, -off:]))
        elif off > length - width:
            cut = off - (length - width); ov = width - cut
            if ov <= 0: continue
            sc = np.nanmean(cosine(pf[:, off:], pwm[:, :ov]))
        else:
            sc = np.nanmean(cosine(pf[:, off:off + width], pwm))
        if sc > best[0]: best = (sc, off)
    return best[1], width
def trim_one(df):
    p = norm(df.reindex(ORD)); ent = np.asarray(calc_entropy_on_PWM(p)); idx = np.where(ent < THR)[0]
    if len(idx) == 0: return None
    s, e = int(idx[0]), int(idx[-1])
    while e - s + 1 < MIN_LEN and (s > 0 or e < p.shape[1] - 1):
        le = ent[s - 1] if s > 0 else np.inf; ri = ent[e + 1] if e < p.shape[1] - 1 else np.inf
        if le <= ri: s -= 1
        else: e += 1
    return p.values[:, s:e + 1] if e - s + 1 >= MIN_LEN else None
def segment(acgu):
    pwm = norm(acgu); ent = np.asarray(calc_entropy_on_PWM(pwm)); info = ent < THR
    L = len(info); bl = []; i = 0
    while i < L:
        if not info[i]: i += 1; continue
        s = e = i; gap = 0; k = i
        while k + 1 < L:
            k += 1
            if info[k]: e = k; gap = 0
            else:
                gap += 1
                if gap > MAX_GAP: break
        if e - s + 1 >= MIN_LEN: bl.append(norm(pwm[:, s:e + 1]))
        i = e + 1
    return bl
def dedup(blocks):
    cl = []
    for b in blocks:
        pl = False
        for c in cl:
            if maxcos(b, c["rep"]) >= DEDUP:
                if b.shape[1] > c["rep"].shape[1]: c["rep"] = b
                pl = True; break
        if not pl: cl.append({"rep": b})
    return [c["rep"] for c in cl]
def consensus(pwm): return "".join("ACGU"[int(np.argmax(col))] for col in pwm.T)
def core_kmer(pwm):
    info = 2 - np.asarray(calc_entropy_on_PWM(pwm)); cons = consensus(pwm)
    if len(cons) <= CORE: return cons
    best = (-1, 0)
    for i in range(len(cons) - CORE + 1):
        sc = info[i:i + CORE].sum()
        if sc > best[0]: best = (sc, i)
    return cons[best[1]:best[1] + CORE]
def is_lowcomp(pwm):
    cons = consensus(pwm); L = len(cons)
    homo = Counter(cons).most_common(1)[0][1] / L
    dinuc = sum(1 for i in range(L - 2) if cons[i] == cons[i + 2]) / max(1, L - 2)
    return homo >= LOWCOMP_HOMO or dinuc >= LOWCOMP_DINUC
def parse_count(h):
    m = re.search(r'-(\d+)-(\d+)-[\d.]+(?:-\d+)?$', h); return int(m.group(2)) if m else 1
def read_fi(mp):
    with contextlib.redirect_stdout(io.StringIO()): p = read_params(mp)
    full = pd.DataFrame(p["fi"], columns=list(p["alphabet"])).T
    matched = [i for i, n in enumerate(p["target_seq"]) if n in "AUGC"]
    return norm(full.reindex(ORD).iloc[:, matched].values)
def parse_streme(txt):
    d = {}; L = open(txt).read().splitlines(); i = 0
    while i < len(L):
        if L[i].startswith("MOTIF"):
            name = L[i].split()[1]; cons = name.split("-", 1)[1] if "-" in name else name
            j = i + 1
            while j < len(L) and not L[j].startswith("letter-probability"): j += 1
            if j >= len(L): break
            w = int(re.search(r'w=\s*(\d+)', L[j]).group(1))
            rows = [[float(x) for x in L[k].split()] for k in range(j + 1, j + 1 + w)]
            df = pd.DataFrame(np.array(rows).T, index=["A", "C", "G", "T"]).rename(index={"T": "U"}).reindex(ORD)
            d[cons] = norm(df.values); i = j + 1 + w
        else: i += 1
    return d

def build_arm(rbp, arm):
    """return list of de-novo motif PWMs (dedup0.8) for an arm."""
    if arm == "vsearch":
        res = f"{DS}/{rbp}/topN/results.tsv"; bl = []
        if not os.path.exists(res): return []
        for _, r in pd.read_csv(res, sep="\t").iterrows():
            mp = f"{DS}/{rbp}/topN/rep{int(r['rep']):02d}.model_params"
            if os.path.exists(mp):
                try: bl += segment(read_fi(mp))
                except Exception: pass
        return dedup(bl)
    if arm == "streme":
        res = f"{DS}/{rbp}/streme_query/results.tsv"; bl = []
        if not os.path.exists(res): return []
        for _, r in pd.read_csv(res, sep="\t").iterrows():
            sid = str(r["seed_id"])
            if sid in ("NOSITE", "MSAfail", "nan", "") or sid.startswith("foldfail"): continue
            mp = f"{DS}/{rbp}/streme_query/rep{int(r['rep']):02d}.model_params"
            if os.path.exists(mp):
                try: bl += segment(read_fi(mp))
                except Exception: pass
        return dedup(bl)
    if arm == "anchor":
        res = f"{DS}/{rbp}/streme_query/results.tsv"; stxt = f"{DS}/{rbp}/streme10k/out/streme.txt"
        if not (os.path.exists(res) and os.path.exists(stxt)): return []
        sp = parse_streme(stxt); bl = []
        for _, r in pd.read_csv(res, sep="\t").iterrows():
            spwm = sp.get(str(r["motif"]))
            if spwm is None: continue
            mp = f"{DS}/{rbp}/streme_query/rep{int(r['rep']):02d}.model_params"
            if not os.path.exists(mp): continue
            try:
                pf = read_fi(mp); off, w = best_offset(spwm, pf)
                if off is None: continue
                lo, hi = max(0, off), min(pf.shape[1], off + w)
                if hi - lo >= MIN_OVL: bl.append(norm(pf[:, lo:hi]))
            except Exception: pass
        return dedup(bl)
    if arm == "streme_only":
        stxt = f"{DS}/{rbp}/streme10k/out/streme.txt"
        if not os.path.exists(stxt): return []
        out = []
        for v in parse_streme(stxt).values():
            t = trim_one(pd.DataFrame(v, index=ORD))
            if t is not None: out.append(t)
        return out

def load_reads(rbp):
    f = f"{SELEXD}/{rbp}/{rbp}.count.ann.all_selex.rmlow.unique.fa"
    if not os.path.exists(f): return []
    reads = []; cnt = 0
    with open(f) as fh:
        for line in fh:
            if line.startswith(">"): cnt = parse_count(line[1:].strip())
            else: reads.append((line.strip().replace("T", "U"), cnt))
    return reads


pwms = load_pwms(XLSX)
genes = {}
for k in pwms: genes.setdefault(k.split("_")[0], []).append(k)
rbps = [d for d in sorted(os.listdir(DS))
        if os.path.isdir(f"{DS}/{d}") and os.path.exists(f"{DS}/{d}/streme_query/results.tsv")]
ARMS = ["vsearch", "streme", "anchor", "streme_only"]

store = {}
for rbp in rbps:
    gene = rbp.split("-construct")[0]
    jol = [trim_one(pwms[k].reindex(ORD)) for k in genes.get(gene, [])]; jol = [x for x in jol if x is not None]
    if not jol: continue
    reads = load_reads(rbp)
    armmot = {}
    for arm in ARMS:
        mot = []
        for pwm in build_arm(rbp, arm):
            core = core_kmer(pwm)
            prev = sum(c for sq, c in reads if core in sq) if reads else 0
            mot.append(dict(pwm=pwm, prevalence=prev, lowcomp=is_lowcomp(pwm)))
        armmot[arm] = mot
    store[rbp] = dict(jol=jol, arms=armmot)
print(f"genes: {len(store)}", flush=True)

def evaluate(arm, lowcomp_filter, topn):
    tj = tjm = td = tdm = 0; ng = 0
    for rbp, s in store.items():
        jol = s["jol"]; mot = s["arms"][arm]
        kept = [m for m in mot if (not m["lowcomp"]) or not lowcomp_filter]
        if topn is not None: kept = sorted(kept, key=lambda m: -m["prevalence"])[:topn]
        dn = [m["pwm"] for m in kept]
        ng += 1; td += len(dn)
        for j in jol:
            tj += 1
            if max([maxcos(j, x) for x in dn], default=0) >= MATCH: tjm += 1
        for x in dn:
            if max([maxcos(x, j) for j in jol], default=0) >= MATCH: tdm += 1
    return dict(recall=tjm / tj, precision=tdm / td if td else 0, mpg=td / ng)

# Panel B trajectories: baseline -> (a)+top5 -> (a)+top3 -> (a)+top2
rows = []
for arm in ARMS:
    rows.append(dict(arm=arm, stage="baseline", **evaluate(arm, False, None)))
    for tn in TOPN:
        rows.append(dict(arm=arm, stage=f"(a)+top{tn}", **evaluate(arm, True, tn)))
fb = pd.DataFrame(rows)
fb.to_csv(f"{DS}/denovo_jolma_filters_byarm.tsv", sep="\t", index=False)
print(fb.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

# ---- 1x2 figure ----
sweep = pd.read_csv(f"{DS}/denovo_jolma_recovery_sweep.tsv", sep="\t")
COL = {"vsearch": "#1f77b4", "streme": "#2ca02c", "anchor": "#d62728", "streme_only": "black"}
LAB = {"vsearch": "vsearch->DCA (self-seg)", "streme": "STREME->DCA self-seg", "anchor": "STREME->DCA anchor", "streme_only": "STREME-only PWM"}
MK = {"vsearch": "-D", "streme": "-D", "anchor": "-D", "streme_only": "-^"}
SW = {"vsearch": "vsearch->DCA (self-seg)", "streme": "STREME->DCA (self-seg)", "anchor": "STREME->DCA (anchor)"}
fig, (axA, axB) = plt.subplots(1, 2, figsize=(14, 6), sharex=True, sharey=True)
# Panel A: raw dedup sweeps (DCA arms) + STREME-only (top-N sweep)
for arm in ["vsearch", "streme", "anchor"]:
    s = sweep[sweep["method"] == SW[arm]].sort_values("dedup")
    axA.plot(s["recall"], s["precision"], "-o", color=COL[arm], lw=1.8, ms=6, label=LAB[arm])
    for _, r in s.iterrows():
        axA.annotate(f"d{r['dedup']}\n{r['mpg']:.1f}/g", (r["recall"], r["precision"]),
                     textcoords="offset points", xytext=(6, -2), fontsize=6.5, color=COL[arm])
so = sweep[sweep["method"] == "STREME-only (PWM)"].sort_values("nmotif")
axA.plot(so["recall"], so["precision"], "-^", color="black", lw=1.4, ms=6, label="STREME-only PWM")
for _, r in so.iterrows():
    axA.annotate(f"N{int(r['nmotif'])}\n{r['mpg']:.1f}/g", (r["recall"], r["precision"]),
                 textcoords="offset points", xytext=(6, -2), fontsize=6.5, color="black")
axA.set_title("(A) WITHOUT Jolma filters  (dedup / top-N sweep)")
axA.set_xlabel("recall"); axA.set_ylabel("precision"); axA.grid(alpha=0.3); axA.legend(loc="center left", fontsize=8)
# Panel B: per-arm filter trajectory
for arm in ARMS:
    t = fb[fb["arm"] == arm]
    t = t.set_index("stage").loc[["baseline", "(a)+top5", "(a)+top3", "(a)+top2"]].reset_index()
    axB.plot(t["recall"], t["precision"], MK[arm], color=COL[arm], lw=2.2, ms=8, label=LAB[arm])
    for _, r in t.iterrows():
        axB.annotate(f"{r['stage'].replace('(a)+','')}\n{r['mpg']:.1f}/g", (r["recall"], r["precision"]),
                     textcoords="offset points", xytext=(6, -2), fontsize=6.5, color=COL[arm])
axB.set_title("(B) WITH Jolma filters  (low-complexity + prevalence top-N)")
axB.set_xlabel("recall"); axB.grid(alpha=0.3); axB.legend(loc="center left", fontsize=8)
fig.suptitle("De-novo Jolma-motif recovery per method: without vs with Jolma-analog filtering", fontsize=12)
plt.tight_layout(); plt.savefig(f"{DS}/denovo_jolma_filters_byarm.png", dpi=150)
print("wrote denovo_jolma_filters_byarm.png")
