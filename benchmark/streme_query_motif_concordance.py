#!/usr/bin/env python
"""
Motif concordance for STREME-query DCA outputs (arm B of the 3-way), set-wise vs Jolma.

PURE de-novo design (matches streme_query_method_comparison.py):
  candidate motif = a STREME motif (1:1 with a streme-query rep). For each rep we anchor ITS OWN
  STREME motif onto that rep's DCA fi (best_offset) and read off the same-width fi window =
  "DCA-refined PWM" (arm B). De-novo motif SET = dedup of those windows (cosine>=DEDUP_THR);
  NO recurrence filter (STREME already separates motifs; dedup merges redundant copies).

  arm A = streme-only raw STREME PWMs (trim_one, no DCA, no recurrence) for the B-vs-A contrast.

Jolma is used only as the reference to score recovery (cosine, recall@RECALL_T); it never
enters localization.

Usage:
  python benchmark/streme_query_motif_concordance.py <dataset_dir> [bq_subdir] [streme_subdir]
"""
import os, re, sys, contextlib, io
import numpy as np, pandas as pd
import os, sys
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)                                   # nolp_nussinov.py
sys.path.insert(0, os.path.dirname(_HERE))                  # the raptcouple src/ package
from src.helper_jolma import load_pwms, calc_entropy_on_PWM
from src.plmc import read_params

DSDIR = sys.argv[1]
BQ_SUB = sys.argv[2] if len(sys.argv) > 2 else "streme_query"
STREME_SUB = sys.argv[3] if len(sys.argv) > 3 else "streme10k"
OUT = os.path.join(DSDIR, f"{BQ_SUB}_motif_concordance.tsv")
XLSX = "公開用/example/Jolma2020/data/others/Supplemental_Tables_S1-S8.xlsx"

ORD = ["A", "C", "G", "U"]
THR, MIN_LEN, MIN_OVL, RECALL_T = 1.5, 5, 4, 0.8
DEDUP_THR = 0.8      # anchored fi-window PWMs with cosine>=this are the same motif
MAX_PADDING = 3


def norm(p):
    p = p.copy(); s = p.sum(axis=0); s[s == 0] = 1.0; return p / s

def cosine(p, q):
    d = np.linalg.norm(p, axis=0) * np.linalg.norm(q, axis=0); d[d == 0] = np.nan
    return (p * q).sum(axis=0) / d

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

def best_offset(pwm, pos_field):
    width = pwm.shape[1]; length = pos_field.shape[1]; mp = min(MAX_PADDING, width // 2)
    best = (-np.inf, None)
    for offset in range(-mp, length - width + mp + 1):
        if offset < 0:
            ov = width + offset
            if ov <= 0: continue
            sc = np.nanmean(cosine(pos_field[:, :ov], pwm[:, -offset:]))
        elif offset > length - width:
            cut = offset - (length - width); ov = width - cut
            if ov <= 0: continue
            sc = np.nanmean(cosine(pos_field[:, offset:], pwm[:, :ov]))
        else:
            sc = np.nanmean(cosine(pos_field[:, offset:offset + width], pwm))
        if sc > best[0]: best = (sc, offset)
    return best[1], width

def trim_one(pwm_df):
    pwm = norm(pwm_df)
    ent = np.asarray(calc_entropy_on_PWM(pwm)); idx = np.where(ent < THR)[0]
    if len(idx) == 0: return None
    start, end = int(idx[0]), int(idx[-1])
    while end - start + 1 < MIN_LEN and (start > 0 or end < pwm.shape[1] - 1):
        le = ent[start - 1] if start > 0 else np.inf; re_ = ent[end + 1] if end < pwm.shape[1] - 1 else np.inf
        if le <= re_: start -= 1
        else: end += 1
    return pwm.values[:, start:end + 1] if end - start + 1 >= MIN_LEN else None

def parse_streme(txt):
    """STREME motifs -> {consensus: normalized ACGU PWM} and (ordered list of trimmed PWMs)."""
    named = {}; lines = open(txt).read().splitlines(); i = 0
    while i < len(lines):
        if lines[i].startswith("MOTIF"):
            name = lines[i].split()[1]; cons = name.split("-", 1)[1] if "-" in name else name
            j = i + 1
            while j < len(lines) and not lines[j].startswith("letter-probability"): j += 1
            if j >= len(lines): break
            w = int(re.search(r'w=\s*(\d+)', lines[j]).group(1))
            rows = [[float(x) for x in lines[k].split()] for k in range(j + 1, j + 1 + w)]
            df = pd.DataFrame(np.array(rows).T, index=["A", "C", "G", "T"]).rename(index={"T": "U"}).reindex(ORD)
            named[cons] = df; i = j + 1 + w
        else:
            i += 1
    return named

def dedup(blocks):
    """cluster anchored fi-window PWMs by cosine>=DEDUP_THR; keep one rep per cluster (widest).
    NO recurrence filter."""
    clusters = []
    for b in blocks:
        placed = False
        for c in clusters:
            if maxcos(b, c["rep"]) >= DEDUP_THR:
                if b.shape[1] > c["rep"].shape[1]: c["rep"] = b
                placed = True; break
        if not placed:
            clusters.append({"rep": b})
    return [c["rep"] for c in clusters]

def streme_query_windows(outdir, streme_named):
    """anchored fi-window PWMs: for each rep, anchor its own STREME motif onto its fi."""
    results = os.path.join(outdir, "results.tsv")
    if not os.path.exists(results): return []
    res = pd.read_csv(results, sep="\t")
    blocks = []
    for _, r in res.iterrows():
        sid = str(r["seed_id"])
        if sid in ("NOSITE", "MSAfail", "nan", "") or sid.startswith("foldfail"): continue
        spwm_df = streme_named.get(str(r["motif"]))
        if spwm_df is None: continue
        mp = os.path.join(outdir, f"rep{int(r['rep']):02d}.model_params")
        if not os.path.exists(mp): continue
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                p = read_params(mp)
            full = pd.DataFrame(p["fi"], columns=list(p["alphabet"])).T
            matched = [i for i, n in enumerate(p["target_seq"]) if n in "AUGC"]
            pos_field = full.reindex(ORD).iloc[:, matched].values
            spwm = norm(spwm_df.values)
            off, w = best_offset(spwm, pos_field)
            if off is None: continue
            lo, hi = max(0, off), min(pos_field.shape[1], off + w)
            if hi - lo < MIN_OVL: continue
            blocks.append(norm(pos_field[:, lo:hi]))
        except Exception:
            pass
    return blocks


print("loading Jolma PWMs ...", flush=True)
pwms = load_pwms(XLSX)
genes = {}
for k in pwms: genes.setdefault(k.split("_")[0], []).append(k)

rows = []
for rbp_dir in sorted(os.path.join(DSDIR, x) for x in os.listdir(DSDIR)):
    if not os.path.isdir(rbp_dir): continue
    rbp = os.path.basename(rbp_dir)
    bq_out = os.path.join(rbp_dir, BQ_SUB)
    if not os.path.exists(os.path.join(bq_out, "results.tsv")): continue
    stxt = os.path.join(rbp_dir, STREME_SUB, "out", "streme.txt")
    streme_named = parse_streme(stxt) if os.path.exists(stxt) else {}
    # --- arm B: streme-query distinct motifs (STREME-anchored fi window + dedup, no recurrence) ---
    bblocks = streme_query_windows(bq_out, streme_named)
    bq = dedup(bblocks)
    # --- arm A: streme-only raw PWMs (trim only, no recurrence) ---
    streme = []
    for m in streme_named.values():
        t = trim_one(m.reindex(ORD))
        if t is not None: streme.append(t)
    # --- jolma reference ---
    gene = rbp.split("-construct")[0]
    jolma = []
    for key in genes.get(gene, []):
        t = trim_one(pwms[key].reindex(ORD))
        if t is not None: jolma.append(t)

    def bestcos(A, B): return max([maxcos(x, y) for x in A for y in B], default=np.nan)
    bj = bestcos(bq, jolma)
    aj = bestcos(streme, jolma)
    ba = bestcos(bq, streme)
    jrec_b = [max([maxcos(j, x) for x in bq], default=np.nan) for j in jolma]
    jrec_a = [max([maxcos(j, x) for x in streme], default=np.nan) for j in jolma]
    rows.append(dict(
        rbp=rbp, n_stremeq=len(bq), n_stremeonly=len(streme), n_jolma=len(jolma),
        stremeq_vs_jolma=round(bj, 3) if bj == bj else np.nan,
        stremeonly_vs_jolma=round(aj, 3) if aj == aj else np.nan,
        stremeq_vs_stremeonly=round(ba, 3) if ba == ba else np.nan,
        jolma_recall_stremeq=np.mean([v >= RECALL_T for v in jrec_b if v == v]) if any(v == v for v in jrec_b) else np.nan,
        jolma_recall_stremeonly=np.mean([v >= RECALL_T for v in jrec_a if v == v]) if any(v == v for v in jrec_a) else np.nan,
    ))

df = pd.DataFrame(rows)
df.to_csv(OUT, sep="\t", index=False)
print(f"wrote {OUT}  (rows={len(df)})")

def med(c):
    v = df[c].dropna()
    return f"{v.median():.3f} (n={len(v)})" if len(v) else "NA"
print("\n=== medians vs Jolma ===")
print(" stremeq_vs_jolma   :", med("stremeq_vs_jolma"))
print(" stremeonly_vs_jolma:", med("stremeonly_vs_jolma"))
print(" jolma_recall_stremeq   (mean):", round(df["jolma_recall_stremeq"].dropna().mean(), 3) if df["jolma_recall_stremeq"].notna().any() else "NA")
print(" jolma_recall_stremeonly(mean):", round(df["jolma_recall_stremeonly"].dropna().mean(), 3) if df["jolma_recall_stremeonly"].notna().any() else "NA")
print(" motif/gene stremeq (mean):", round(df["n_stremeq"].mean(), 2), " stremeonly (mean):", round(df["n_stremeonly"].mean(), 2))
pair = df.dropna(subset=["stremeq_vs_jolma", "stremeonly_vs_jolma"])
if len(pair) >= 3:
    from scipy.stats import wilcoxon
    try:
        w = wilcoxon(pair["stremeq_vs_jolma"], pair["stremeonly_vs_jolma"])
        win = int((pair["stremeq_vs_jolma"] > pair["stremeonly_vs_jolma"]).sum())
        print(f"\n B vs A (cos vs Jolma) paired n={len(pair)}: B>A {win}/{len(pair)}  Wilcoxon p={w.pvalue:.2e}")
    except Exception as e:
        print(" wilcoxon failed:", e)
