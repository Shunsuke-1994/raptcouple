#!/usr/bin/env python
"""
STREME-query motif-local STRUCTURE comparison across DCA / MIp / RNAalifold / R-scape,
set-wise vs Jolma. PURE de-novo design (no Jolma in localization):

  candidate motif = a STREME motif (1:1 with a streme-query rep). For each rep we anchor
  ITS OWN STREME motif onto that rep's DCA fi (best_offset) -> window [off-K, off+w+K] in
  matched coords. The structure window is therefore defined ENTIRELY by STREME->seed->DCA;
  Jolma is never used to place a window. De-novo motifs are obtained by deduping the anchored
  fi-window PWMs (cosine>=DEDUP_THR); NO recurrence filter (STREME already separates motifs,
  so we merge redundant copies rather than drop low-recurrence ones).

  A de-novo motif is structured(method) if a motif-local base pair lies in its STREME-anchored
  window, aggregated over the cluster's member reps by MAJORITY vote (single-rep motif = that
  rep's call; ties -> structured). Jolma is used ONLY for answer-checking: each Jolma motif is
  matched to its best de-novo motif by PWM cosine (cos>=MATCH_T), and on matched motifs the
  structured call is compared to the Jolma motif_type.

R-scape input: <DSDIR>/<SUBDIR>_rscape/<RBP>__repNN_1.cov
STREME motifs : <DSDIR>/<RBP>/<STREME_SUB>/out/streme.txt  (joined to results.tsv `motif` col)

Usage:
  python benchmark/streme_query_method_comparison.py <dataset_dir> [subdir] [rbp_filter_csv]
"""
import collections, contextlib, io, os, re, subprocess, sys, tempfile, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from statsmodels.stats import contingency_tables

import os, sys
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)                                   # nolp_nussinov.py
sys.path.insert(0, os.path.dirname(_HERE))                  # the raptcouple src/ package
from nolp_nussinov import nussinov_nolp
from src.helper_jolma import calc_entropy_on_PWM, load_metadata, load_pwms
from src.plmc import detect_coupling, read_params

DSDIR = sys.argv[1]
SUBDIR = sys.argv[2] if len(sys.argv) > 2 else "streme_query"
RBP_FILTER = set(sys.argv[3].split(",")) if len(sys.argv) > 3 else None
STREME_SUB = "streme10k"
RSCAPE_DIR = os.path.join(DSDIR, f"{SUBDIR}_rscape")
# candidate-extraction mode (seeding stays 1:1 STREME motif -> rep regardless):
#   streme  = one window per rep, anchored at the rep's own STREME motif (current default)
#   segment = self-segment each rep's fi into entropy blocks (decoupled from STREME localization)
#   union   = both
CAND_MODE = os.environ.get("CAND_MODE", "streme")
# Optional Jolma-analog filtering of the de-novo motif SET before structure eval:
#   FILTER_TOPN=N -> drop low-complexity motifs, then keep top-N by prevalence (primary/secondary
#   cap). Lets us compare linear/structured detection with vs without Jolma-like curation.
FILTER_TOPN = int(os.environ["FILTER_TOPN"]) if os.environ.get("FILTER_TOPN") else None
SELEXD = "SELEXdata/JZM2020"
CORE_K, LOWCOMP_HOMO, LOWCOMP_DINUC = 6, 0.70, 0.80
_suf = SUBDIR if CAND_MODE == "streme" else f"{SUBDIR}_{CAND_MODE}"
OUT_PREFIX = os.path.join(DSDIR, _suf + (f"_top{FILTER_TOPN}" if FILTER_TOPN else ""))
XLSX = "公開用/example/Jolma2020/data/others/Supplemental_Tables_S1-S8.xlsx"
RNAALIFOLD = os.environ.get("RNAALIFOLD", "RNAalifold")

ORD = ["A", "C", "G", "U"]; ALPH = "AUGC."; Q = 5; LAMBDA = 0.5
THR, MIN_LEN, MIN_OVL, MAX_GAP = 1.5, 5, 4, 1
DEDUP_THR = 0.8          # anchored fi-window PWMs with cosine>=this are the same motif
ZT, ML, MATCH_T, WINDOW_K, MAX_PADDING = 3.0, 3, 0.8, 6, 3
RNAALIFOLD_MAX_ROWS = 5000
METHODS = ["dca", "mip", "rnaalifold", "rscape"]
# "streme" = pure de-novo headline window (STREME motif anchored on fi).
# "jolma_anchored" = DIAGNOSTIC ONLY: localize the matched Jolma motif on the SAME reps' fi
#   (same window definition as the jolma-query benchmark) -> isolates the de-novo-window effect
#   from the de-novo-seed effect when comparing streme-query vs jolma-query. Jolma is used here
#   only to place an evaluation window for the comparison, never in detection.
WINS = ["streme", "jolma_anchored"]


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
    """localize query pwm on pos_field (matched coords, 4xL); return (offset, width).
    offset may be negative or run past the end (partial overlap, padded up to MAX_PADDING)."""
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

def segment_pos(acgu):
    """self-segment matched-coord fi (4xL) into contiguous high-information blocks -> [(s,e), ...].
    Blocks shorter than MIN_LEN are DROPPED (not extended): unlike trim_one() for a single motif,
    extending a short fi run here only appends low-information columns -> tested to give no recall
    gain (matched stays 138) and slightly lower structure kappa (0.640->0.618), so drop is kept.
    Real short motifs (e.g. ELAVL3 GUUUUUG, ~7nt) already exceed MIN_LEN and survive."""
    pwm = norm(acgu); ent = np.asarray(calc_entropy_on_PWM(pwm)); info = ent < THR
    L = len(info); blocks = []; i = 0
    while i < L:
        if not info[i]: i += 1; continue
        s = e = i; gap = 0; k = i
        while k + 1 < L:
            k += 1
            if info[k]: e = k; gap = 0
            else:
                gap += 1
                if gap > MAX_GAP: break
        if e - s + 1 >= MIN_LEN: blocks.append((s, e))
        i = e + 1
    return blocks

def consensus(pwm):
    return "".join("ACGU"[int(np.argmax(col))] for col in pwm.T)

def core_kmer(pwm):
    info = 2 - np.asarray(calc_entropy_on_PWM(pwm)); cons = consensus(pwm)
    if len(cons) <= CORE_K: return cons
    best = (-1, 0)
    for i in range(len(cons) - CORE_K + 1):
        sc = info[i:i + CORE_K].sum()
        if sc > best[0]: best = (sc, i)
    return cons[best[1]:best[1] + CORE_K]

def is_lowcomp(pwm):
    from collections import Counter as _C
    cons = consensus(pwm); L = len(cons)
    homo = _C(cons).most_common(1)[0][1] / L
    dinuc = sum(1 for i in range(L - 2) if cons[i] == cons[i + 2]) / max(1, L - 2)
    return homo >= LOWCOMP_HOMO or dinuc >= LOWCOMP_DINUC

def parse_count(h):
    m = re.search(r'-(\d+)-(\d+)-[\d.]+(?:-\d+)?$', h); return int(m.group(2)) if m else 1

def load_reads(rbp):
    f = f"{SELEXD}/{rbp}/{rbp}.count.ann.all_selex.rmlow.unique.fa"
    if not os.path.exists(f): return []
    reads = []; cnt = 0
    with open(f) as fh:
        for line in fh:
            if line.startswith(">"): cnt = parse_count(line[1:].strip())
            else: reads.append((line.strip().replace("T", "U"), cnt))
    return reads

def trim_one(pwm_df):
    pwm = norm(pwm_df.reindex(ORD)); ent = np.asarray(calc_entropy_on_PWM(pwm)); idx = np.where(ent < THR)[0]
    if len(idx) == 0: return None
    s, e = int(idx[0]), int(idx[-1])
    while e - s + 1 < MIN_LEN and (s > 0 or e < pwm.shape[1] - 1):
        le = ent[s - 1] if s > 0 else np.inf; re_ = ent[e + 1] if e < pwm.shape[1] - 1 else np.inf
        if le <= re_: s -= 1
        else: e += 1
    return pwm.values[:, s:e + 1] if e - s + 1 >= MIN_LEN else None

def parse_streme(txt):
    """STREME motifs -> {consensus: normalized ACGU PWM (4 x w)} (join key = results.tsv `motif`)."""
    d = {}; lines = open(txt).read().splitlines(); i = 0
    while i < len(lines):
        if lines[i].startswith("MOTIF"):
            name = lines[i].split()[1]; cons = name.split("-", 1)[1] if "-" in name else name
            j = i + 1
            while j < len(lines) and not lines[j].startswith("letter-probability"): j += 1
            if j >= len(lines): break
            w = int(re.search(r'w=\s*(\d+)', lines[j]).group(1))
            rows = [[float(x) for x in lines[k].split()] for k in range(j + 1, j + 1 + w)]
            df = pd.DataFrame(np.array(rows).T, index=["A", "C", "G", "T"]).rename(index={"T": "U"}).reindex(ORD)
            d[cons] = norm(df.values); i = j + 1 + w
        else:
            i += 1
    return d

def db_pairs(db):
    stk = []; pr = set()
    for k, ch in enumerate(db):
        if ch in "([{<": stk.append(k)
        elif ch in ")]}>": pr.add((stk.pop(), k))
    return pr

def load_msa(path, length):
    rows = []; cur = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if cur: rows.append("".join(cur)); cur = []
            elif line: cur.append(line)
    if cur: rows.append("".join(cur))
    return [r.upper().replace("T", "U") for r in rows if len(r) == length]

def mi_msa(alignment, length):
    idx = {c: p for p, c in enumerate(ALPH)}; nseq = len(alignment)
    x = np.zeros((nseq, length * Q), dtype=np.float32)
    for n, row in enumerate(alignment):
        for i, ch in enumerate(row):
            x[n, i * Q + idx.get(ch, 4)] = 1.0
    xtx = (x.T @ x) / nseq; fi = x.mean(0).reshape(length, Q); mi = np.zeros((length, length))
    for i in range(length - 1):
        for j in range(i + 1, length):
            f2 = (1 - LAMBDA) * xtx[i * Q:(i + 1) * Q, j * Q:(j + 1) * Q] + LAMBDA / Q / Q
            fo = np.outer((1 - LAMBDA) * fi[i] + LAMBDA / Q, (1 - LAMBDA) * fi[j] + LAMBDA / Q)
            mi[i, j] = mi[j, i] = np.sum(f2 * np.log(f2 / fo))
    mi_i = mi.sum(1) / (length - 1); mi_bar = mi[np.triu_indices(length, 1)].mean()
    mip = mi - np.outer(mi_i, mi_i) / mi_bar; np.fill_diagonal(mip, 0); return mip

def fold_pairs(params, score_matrix, matched, pos2idx, rna_seq):
    tmp = dict(params); tmp["FN_apc"] = score_matrix
    pairs = detect_coupling(tmp, threshold=ZT, min_dist=ML + 1, sanity_check=False, use_mask_min_dist=True)
    st = np.zeros((len(matched), len(matched)))
    for _nuc, pos in pairs:
        a, b = pos
        if a in pos2idx and b in pos2idx: st[pos2idx[a], pos2idx[b]] = score_matrix[a, b]
    _s, db = nussinov_nolp(rna_seq, ML, st, sanity_check=False)
    return db_pairs(db)

def run_rnaalifold(msa_path, length):
    # RNAalifold output depends ONLY on the MSA (not on candidate mode / window), so cache the
    # consensus dot-bracket in a sidecar per rep MSA -> any later run (streme/segment/union) is instant.
    cache = msa_path + ".rnaalifold"
    if os.path.exists(cache):
        try:
            s = open(cache).read().strip()
            return None if s == "NONE" else (s if len(s) == length else None)
        except Exception:
            pass
    rows = []; cur = []; hdr = None
    with open(msa_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur: rows.append((hdr, "".join(cur))); cur = []
                hdr = line
            elif line: cur.append(line)
    if cur: rows.append((hdr, "".join(cur)))
    rows = [(h, s.replace(".", "-").upper().replace("T", "U")) for h, s in rows if len(s) == length][:RNAALIFOLD_MAX_ROWS]
    result = None
    if rows:
        fasta = "\n".join(f"{h}\n{s}" for h, s in rows) + "\n"
        with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as tmp:
            tmp.write(fasta); path = tmp.name
        try:
            proc = subprocess.run([RNAALIFOLD, "--noPS", "-r", "-f", "F", path], capture_output=True, text=True, timeout=300, cwd=tempfile.gettempdir())
        finally:
            os.unlink(path)
        for line in proc.stdout.splitlines():
            tok = line.strip().split(" ")[0]
            if set(tok) <= set("().,[]{}<>~_") and len(tok) == length:
                result = tok.replace(",", ".").replace("~", ".").replace("_", "."); break
    try:
        with open(cache, "w") as fh: fh.write(result if result else "NONE")
    except Exception:
        pass
    return result

def rnaalifold_pairs(msa_path, length, pos2idx):
    struct = run_rnaalifold(msa_path, length)
    if not struct: return set()
    raw = sorted((pos2idx[a], pos2idx[b]) for a, b in db_pairs(struct) if a in pos2idx and b in pos2idx)
    rs = set(raw); kept = set()
    for a, b in raw:
        if b - a < ML + 1: continue
        if (a + 1, b - 1) in rs or (a - 1, b + 1) in rs: kept.add((a, b))
    return kept

def rscape_pairs(cov_path, pos2idx):
    if not os.path.exists(cov_path): return set(), "missing"
    pairs = set()
    with open(cov_path) as fh:
        for line in fh:
            if not (line.startswith("~") or (line and line[0].isdigit())): continue
            parts = line.split()
            try: a, b = int(parts[1]) - 1, int(parts[2]) - 1
            except Exception: continue
            if a in pos2idx and b in pos2idx and pos2idx[b] - pos2idx[a] >= ML + 1:
                pairs.add((pos2idx[a], pos2idx[b]))
    return pairs, ""

def count_local(pairs, lo, hi):
    return sum(1 for a, b in pairs if lo <= a <= hi and lo <= b <= hi)


# ---------- Jolma reference (gene -> [(seed, type, trimmed pwm)]) ----------
print("loading Jolma PWMs + metadata ...", flush=True)
pwms = load_pwms(XLSX); md = load_metadata(XLSX).reset_index(drop=True)
jolma_by_gene = collections.defaultdict(list)
for _, r in md.iterrows():
    gene = str(r["Gene Symbol"]); seed = str(r["seed"]).upper().replace("T", "U")
    mt = "structured" if r["motif_type"] == "structured" else "linear"
    pwm = None
    for key, m in pwms.items():
        parts = key.split("_")
        if parts[0] == gene and seed in [p.upper() for p in parts]:
            pwm = trim_one(m); break
    if pwm is not None: jolma_by_gene[gene].append((seed, mt, pwm))

def rep_methods(model_path, msa_path, cov_path):
    """return (bp{method}, usable{method}, pos_field) in matched coords (raw fi). No segmentation."""
    with contextlib.redirect_stdout(io.StringIO()):
        p = read_params(model_path)
    ts = p["target_seq"]; matched = [i for i, n in enumerate(ts) if n in "AUGC"]
    pos2idx = {pp: i for i, pp in enumerate(matched)}; rna_seq = "".join(ts[i] for i in matched)
    full = pd.DataFrame(p["fi"], columns=list(p["alphabet"])).T
    pos_field = full.reindex(ORD).iloc[:, matched].values            # 4 x L_matched (raw fi)
    bp = {}; usable = {}
    bp["dca"] = fold_pairs(p, p["FN_apc"], matched, pos2idx, rna_seq); usable["dca"] = True
    ali = load_msa(msa_path, len(ts))
    bp["mip"] = fold_pairs(p, mi_msa(ali, len(ts)), matched, pos2idx, rna_seq); usable["mip"] = True
    if os.environ.get("SKIP_RNAALIFOLD"):
        bp["rnaalifold"] = set(); usable["rnaalifold"] = False   # fast path: DCA/MIp/Rscape only
    else:
        try:
            bp["rnaalifold"] = rnaalifold_pairs(msa_path, len(ts), pos2idx); usable["rnaalifold"] = True
        except Exception:
            bp["rnaalifold"] = set(); usable["rnaalifold"] = False
    rs, err = rscape_pairs(cov_path, pos2idx); bp["rscape"] = rs; usable["rscape"] = (err == "")
    return bp, usable, pos_field


# ---------- main ----------
rows = []; funnel = []
for rbp in sorted(os.listdir(DSDIR)):
    rbp_dir = os.path.join(DSDIR, rbp)
    if not os.path.isdir(rbp_dir): continue
    if RBP_FILTER and rbp not in RBP_FILTER: continue
    gene = rbp.split("-construct")[0]
    if gene not in jolma_by_gene: continue
    outdir = os.path.join(rbp_dir, SUBDIR); res = os.path.join(outdir, "results.tsv")
    if not os.path.exists(res): continue
    stxt = os.path.join(rbp_dir, STREME_SUB, "out", "streme.txt")
    streme_pwms = parse_streme(stxt) if os.path.exists(stxt) else {}
    rdf = pd.read_csv(res, sep="\t")
    units = []   # candidate units: each is a (off,w,fi-window) on some rep, tagged with that rep's bp/usable/pf
    for _, rr in rdf.iterrows():
        sid = str(rr["seed_id"])
        if sid in ("NOSITE", "MSAfail", "nan", "") or sid.startswith("foldfail"): continue
        rep = int(rr["rep"]); mp = os.path.join(outdir, f"rep{rep:02d}.model_params")
        msa = os.path.join(outdir, f"rep{rep:02d}.msa")
        cov = os.path.join(RSCAPE_DIR, f"{rbp}__rep{rep:02d}_1.cov")
        if not (os.path.exists(mp) and os.path.exists(msa)): continue
        try:
            bp, usable, pos_field = rep_methods(mp, msa, cov)
        except Exception:
            continue
        Lm = pos_field.shape[1]; cands = []   # (off, w) windows in matched coords
        if CAND_MODE in ("streme", "union"):
            spwm = streme_pwms.get(str(rr["motif"]))
            if spwm is not None:
                off, w = best_offset(spwm, pos_field)
                if off is not None: cands.append((off, w))
        if CAND_MODE in ("segment", "union"):
            for (s, e) in segment_pos(pos_field): cands.append((s, e - s + 1))
        for (off, w) in cands:
            lo0, hi0 = max(0, off), min(Lm, off + w)
            if hi0 - lo0 < MIN_OVL: continue
            units.append(dict(sid=sid, off=off, w=w, bp=bp, usable=usable,
                              win=norm(pos_field[:, lo0:hi0]), pf=pos_field))
    if not units: continue
    # dedup candidate fi-window PWMs -> de-novo motifs (NO recurrence filter)
    clusters = []
    for rp in units:
        placed = False
        for c in clusters:
            if maxcos(rp["win"], c["rep"]) >= DEDUP_THR:
                c["members"].append(rp); c["seeds"].add(rp["sid"])
                if rp["win"].shape[1] > c["rep"].shape[1]: c["rep"] = rp["win"]
                placed = True; break
        if not placed:
            clusters.append(dict(rep=rp["win"], members=[rp], seeds={rp["sid"]}))
    motifs = clusters   # keep all (dedup only, no recurrence)
    if FILTER_TOPN is not None:   # Jolma-analog filter: drop low-complexity, keep top-N by prevalence
        reads = load_reads(rbp)
        for c in motifs:
            core = core_kmer(c["rep"])
            c["prev"] = sum(cnt for sq, cnt in reads if core in sq) if reads else 0
            c["lowcomp"] = is_lowcomp(c["rep"])
        motifs = sorted([c for c in motifs if not c["lowcomp"]], key=lambda c: -c["prev"])[:FILTER_TOPN]
        if not motifs: continue
    funnel.append(dict(rbp=rbp, n_units=len(units), n_denovo_motifs=len(motifs),
                       n_distinct_seeds=len({rp["sid"] for rp in units})))

    def cluster_struct(cluster, win_fn):
        """per-method (call, usable) by majority vote over member reps; win_fn(rep)->(lo,hi) or None."""
        out = {}
        for m in METHODS:
            votes = []
            for rp in cluster["members"]:
                if not rp["usable"][m]: continue
                win = win_fn(rp)
                if win is None: continue
                votes.append(count_local(rp["bp"][m], win[0], win[1]) > 0)
            if not votes:
                out[m] = (False, False); continue
            out[m] = (sum(votes) * 2 >= len(votes), True)   # majority; ties -> structured
        return out
    def streme_win(rp): return (rp["off"] - WINDOW_K, rp["off"] + rp["w"] + WINDOW_K)
    streme_calls = [cluster_struct(c, streme_win) for c in motifs]

    # answer-checking: each Jolma motif -> best de-novo motif by PWM cosine (fi-window)
    for (seed, mt, jpwm) in jolma_by_gene[gene]:
        cos = [maxcos(jpwm, c["rep"]) for c in motifs]
        bi = int(np.nanargmax([v if v == v else -1 for v in cos]))
        bc = cos[bi]; matched = (bc == bc and bc >= MATCH_T)
        scall = streme_calls[bi]
        def jwin(rp, _j=jpwm):
            off, w = best_offset(_j, rp["pf"])
            return None if off is None else (off - WINDOW_K, off + w + WINDOW_K)
        jcall = cluster_struct(motifs[bi], jwin)   # diagnostic: Jolma-anchored window on same reps
        row = dict(rbp=rbp, gene=gene, seed=seed, jolma_type=mt,
                   match_cos=round(bc, 3) if bc == bc else np.nan, matched=matched,
                   n_denovo_motifs=len(motifs), n_reps=len(motifs[bi]["members"]),
                   n_seeds=len(motifs[bi]["seeds"]))
        for m in METHODS:
            row[f"{m}_streme"], su = scall[m]
            row[f"{m}_streme_usable"] = matched and su
            row[f"{m}_jolma_anchored"], ju = jcall[m]
            row[f"{m}_jolma_anchored_usable"] = matched and ju
        rows.append(row)

df = pd.DataFrame(rows)
fdf = pd.DataFrame(funnel)
df.to_csv(OUT_PREFIX + "_method_comparison.tsv", sep="\t", index=False)
fdf.to_csv(OUT_PREFIX + "_method_funnel.tsv", sep="\t", index=False)

def metrics(df, win):
    out = []
    for method in METHODS:
        u = df[df[f"{method}_{win}_usable"]]
        bs = bl = om = oj = 0
        for _, r in u.iterrows():
            j = r["jolma_type"] == "structured"; call = bool(r[f"{method}_{win}"])
            if call and j: bs += 1
            elif (not call) and (not j): bl += 1
            elif call: om += 1
            else: oj += 1
        n = bs + bl + om + oj; mp_, mn = bs + om, bl + oj; jp, jn = bs + oj, bl + om
        po = (bs + bl) / n if n else np.nan
        pe = (mp_ * jp + mn * jn) / (n * n) if n else np.nan
        kappa = np.nan if not n or np.isclose(1 - pe, 0) else (po - pe) / (1 - pe)
        den = np.sqrt(mp_ * mn * jp * jn); mcc = np.nan if den == 0 else (bs * bl - om * oj) / den
        try: mcn = contingency_tables.mcnemar([[bl, om], [oj, bs]]).pvalue
        except Exception: mcn = np.nan
        out.append(dict(window=win, method=method, usable=n, both_struct=bs, both_lin=bl,
                        only_method=om, only_jolma=oj, recall=bs / jp if jp else np.nan,
                        fpr=om / jn if jn else np.nan, cohen_kappa=kappa, mcc=mcc, mcnemar_p=mcn))
    return pd.DataFrame(out)

met = pd.concat([metrics(df, w) for w in WINS], ignore_index=True)
met.to_csv(OUT_PREFIX + "_method_metrics.tsv", sep="\t", index=False)
print(f"wrote {OUT_PREFIX}_method_{{comparison,funnel,metrics}}.tsv")
print(f"rows(jolma motifs)={len(df)}  rbps={df['rbp'].nunique() if len(df) else 0}  matched={int(df['matched'].sum()) if len(df) else 0}")
if len(df):
    for w in WINS:
        print(f"\n=== window={w} ===")
        print(met[met.window == w].drop(columns=["window"]).to_string(index=False, float_format=lambda x: f"{x:.3g}"))
