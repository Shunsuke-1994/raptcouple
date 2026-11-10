#!/usr/bin/env python
"""Does the query-free candidate set contain a motif corresponding to each PDB RNA?

The PDB comparison (pdb_concordance.py) localizes candidates by sequence coverage, so the
selected candidate is simply the one that covers the PDB RNA best. That a selected candidate
does not correspond therefore says nothing about whether a corresponding one exists. This
script tests correspondence directly, with one rule applied to all 14 RBPs and both datasets.

Element. For each PDB RNA (DSSR dot-bracket in pdb_concordance_reselected.tsv), the element is
the hairpin closed by the outermost pair of its longest helix: all DSSR pairs inside that
window, and its unpaired nucleotides, which stand in for the recognized loop or single-stranded
nucleotides. Ties between helices of equal length go to the 5'-most.

Candidates. Every query-free model of the RBP (both arms, streme_query and topN), taken as the
full model before motif segmentation: the positional nucleotide distribution of its matched
columns and its MSA. Whether the best mapping lies inside a segmented motif is reported, so a
stem arm lost by segmentation can be told apart from a candidate-generation failure.

Mapping and scores. The element is placed on each model at every ungapped offset at which it
fits entirely. Two components, each in [0, 1]:
  sequence - mean probability, under the model's positional distribution, of the PDB nucleotide
             at the element's unpaired positions. Substitutions lower the score; they are not
             rejected outright.
  pairing  - fraction of the element's DSSR pairs whose mapped columns are pairable (Watson-Crick
             or G-U) in at least half of the MSA sequences that have no gap at either column.
             Pairability is a property of the SELEX family, so a stem with a different sequence,
             as for RC3H2, still scores; and it does not use RaptCouple's coupling call.
Correspondence requires both, so the joint score is the smaller of the two after the sequence
score is rescaled so that 0 is chance (0.25, the mean probability of any nucleotide) and 1 is a
perfect match. A loop without the PDB stem, or a stem without the loop, therefore scores near 0.
Correspondence is assessed on the maximum joint score over all models and offsets.

Null. Searching many models and offsets finds a good-looking mapping by chance, so the null
repeats the same search after permuting each model's matched columns: one permutation is applied
to the positional distribution and to the pairability matrix, which keeps every column's
composition and every column pair's pairability but places them at random. The p-value is
(1 + number of permutations whose maximum reaches the observed one) / (1 + permutations).
Correspondence is called at p < 0.05. The sequence component is tested alone in the same way and
reported, so that a recovered loop without the PDB stem is visible.

For the best mapping, RaptCouple's predicted pairs (z >= 3, noLP, on the full model) inside the
window are counted, with how many coincide with mapped DSSR pairs. Outcome:
  corresponding, stem predicted     - structure prediction can be evaluated
  corresponding, no stem predicted  - a corresponding candidate exists and its stem was not
                                      predicted (a miss)
  no corresponding candidate        - workflow: not recovered; structure evaluation: not
                                      evaluable

LARP7. All rows mentioning LARP7 in the Jolma2020 supplementary spreadsheet are written out, so
the SELEX construct can be compared with the C-terminal xRRM in 6D12.

Run from the repository root where revision/query_free/{JZM2020,Laverty2023}/ and the Jolma2020
spreadsheet are available (the external drive):
    python revision/structure/pdb_validation/pdb_correspondence.py
Outputs in revision/structure/pdb_validation/: pdb_correspondence.tsv, pdb_correspondence.md,
larp7_construct.tsv.
    python revision/structure/pdb_validation/pdb_correspondence.py --selftest
checks the scoring and the null on synthetic families and needs no data.

Reproducing without the external drive. The test uses only, per model, the positional
distribution of the matched columns, the pairability matrix computed from the MSA, the predicted
pairs and the motif segments. `--export-inputs` writes exactly these, for every RBP and dataset,
to correspondence_inputs/ (a few MB, tracked in git). When the models are not available, the
script reads that folder instead and gives the same result; the LARP7 spreadsheet step then
keeps the committed larp7_construct.tsv.
    python revision/structure/pdb_validation/pdb_correspondence.py --export-inputs
"""
import argparse
import csv
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
INPUTS = HERE / "correspondence_inputs"
RAW_PROBE = HERE.parents[1] / "query_free/JZM2020/RC3H2/topN"     # present only with the external drive
ORD = "ACGU"
IDX = {b: k for k, b in enumerate(ORD)}
PAIRS = [("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"), ("G", "U"), ("U", "G")]
PAIRABLE = 0.5
CHANCE = 0.25                     # mean probability of any nucleotide under a positional distribution
ALPHA = 0.05
PERMUTATIONS = 1000
SEED = 20261001
DATASETS = ["JZM2020", "Laverty2023"]
DATASET_LABEL = {"JZM2020": "Jolma2020", "Laverty2023": "Laverty2022"}


# ---------- structure of the PDB RNA ----------
def pairs_from_dbn(dbn):
    stack, pairs = [], []
    for k, c in enumerate(dbn):
        if c == "(":
            stack.append(k)
        elif c == ")":
            pairs.append((stack.pop(), k))
    return sorted(pairs)


def element(dbn):
    """The hairpin closed by the outermost pair of the longest helix (0-based, inclusive)."""
    pairs = pairs_from_dbn(dbn)
    if not pairs:
        return None
    pset, best = set(pairs), None
    for i, j in pairs:
        if (i - 1, j + 1) in pset:          # not the outermost pair of its stack
            continue
        n, a, b = 0, i, j
        while (a, b) in pset:
            n, a, b = n + 1, a + 1, b - 1
        if best is None or n > best[0]:
            best = (n, i, j)
    n, start, end = best
    inside = [(i, j) for i, j in pairs if start <= i and j <= end]
    paired = {k for pair in inside for k in pair}
    return dict(start=start, end=end, stack=n, pairs=inside,
                unpaired=[k for k in range(start, end + 1) if k not in paired])


# ---------- the SELEX family ----------
def pairability(rows, columns):
    """W x W fraction of sequences pairable (WC or G-U) at each column pair, gaps excluded."""
    length = len(rows[0])
    if any(len(r) != length for r in rows):
        raise ValueError("MSA rows differ in length")
    arr = np.frombuffer("".join(rows).encode(), dtype="S1").reshape(len(rows), length)[:, columns]
    onehot = {b: (arr == b.encode()).astype(np.float32) for b in ORD}
    valid = sum(onehot.values())
    pairable = sum(onehot[a].T @ onehot[b] for a, b in PAIRS)
    counted = valid.T @ valid
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(counted > 0, pairable / counted, np.nan)


# ---------- scoring ----------
def prepare(el, pdb_seq):
    """Element positions relative to its start, keeping only positions with A, C, G or U."""
    seq = pdb_seq.upper().replace("T", "U")
    un = [k for k in el["unpaired"] if seq[k] in IDX]
    pr = [(i, j) for i, j in el["pairs"] if seq[i] in IDX and seq[j] in IDX]
    return dict(span=el["end"] - el["start"] + 1,
                un=np.array([k - el["start"] for k in un], int),
                base=np.array([IDX[seq[k]] for k in un], int),
                pi=np.array([i - el["start"] for i, _ in pr], int),
                pj=np.array([j - el["start"] for _, j in pr], int))


def scores(pwm, pair_matrix, ref):
    """Sequence, pairing and joint scores at every offset where the element fits."""
    width = pwm.shape[1]
    offsets = np.arange(0, width - ref["span"] + 1)
    if len(offsets) == 0:
        return offsets, np.array([]), np.array([]), np.array([])
    if len(ref["un"]):
        sequence = pwm[ref["base"][None, :], offsets[:, None] + ref["un"][None, :]].mean(1)
    else:
        sequence = np.full(len(offsets), np.nan)
    if len(ref["pi"]):
        values = pair_matrix[offsets[:, None] + ref["pi"][None, :], offsets[:, None] + ref["pj"][None, :]]
        pairing = np.nan_to_num(values, nan=0.0).__ge__(PAIRABLE).mean(1)
    else:
        pairing = np.full(len(offsets), np.nan)
    rescaled = np.clip((sequence - CHANCE) / (1 - CHANCE), 0, None)
    joint = np.nanmin(np.vstack([rescaled, pairing]), axis=0)
    return offsets, sequence, pairing, joint


def best_over(models, ref, key):
    """(score, model index, offset index) maximizing `key` ('joint' or 'sequence')."""
    best = (-np.inf, None, None)
    for m, (pwm, pair_matrix) in enumerate(models):
        offsets, sequence, pairing, joint = scores(pwm, pair_matrix, ref)
        values = joint if key == "joint" else sequence
        if len(values) and np.isfinite(values).any():
            k = int(np.nanargmax(values))
            if values[k] > best[0]:
                best = (float(values[k]), m, k)
    return best


def null_maxima(models, ref, rng, permutations=PERMUTATIONS):
    """Maxima of the joint and sequence scores after permuting each model's matched columns."""
    joint, sequence = np.empty(permutations), np.empty(permutations)
    for b in range(permutations):
        shuffled = []
        for pwm, pair_matrix in models:
            order = rng.permutation(pwm.shape[1])
            shuffled.append((pwm[:, order], pair_matrix[np.ix_(order, order)]))
        joint[b] = best_over(shuffled, ref, "joint")[0]
        sequence[b] = best_over(shuffled, ref, "sequence")[0]
    return joint, sequence


def p_value(observed, null):
    return (1 + int(np.sum(null >= observed - 1e-12))) / (1 + len(null))


# ---------- data (needs the external drive) ----------
def load_models(pc, gene, dataset):
    """Full query-free models of `gene` in `dataset`: positional distribution, MSA pairability,
    predicted pairs and motif segments, all in matched-column coordinates."""
    out = []
    for directory in pc.denovo_dirs(gene):
        if directory.split("/")[2] != dataset:
            continue
        for arm in ["streme_query", "topN"]:
            for path in sorted(Path(directory, arm).glob("rep*.model_params")):
                msa = path.with_suffix(".msa")
                if not msa.exists():
                    continue
                folded = pc.rep_fold(str(path))
                if folded is None:
                    continue
                params, _rna, _energy, ss = folded
                pwm, matched = pc.fi_acgu(params)
                rows, current = [], ""
                for line in open(msa):
                    if line.startswith(">"):
                        if current:
                            rows.append(current)
                        current = ""
                    else:
                        current += line.strip()
                if current:
                    rows.append(current)
                rows = [r.upper().replace("T", "U") for r in rows]
                try:
                    pair_matrix = pairability(rows, matched)
                except ValueError as error:
                    print(f"skip {path}: {error}", file=sys.stderr)
                    continue
                out.append(dict(path=os.path.relpath(path, "revision/query_free"), pwm=pwm,
                                pair_matrix=pair_matrix, predicted=set(pairs_from_dbn(ss)),
                                segments=pc.segment(pwm)))
    return out


def save_inputs(records, gene, dataset):
    """Write the per-model inputs of the test (see the module docstring) for one RBP and dataset."""
    arrays = {"paths": np.array([r["path"] for r in records])}
    for k, r in enumerate(records):
        arrays[f"pwm_{k}"] = r["pwm"]
        arrays[f"pair_{k}"] = r["pair_matrix"]
        arrays[f"predicted_{k}"] = np.array(sorted(r["predicted"]), int).reshape(-1, 2)
        arrays[f"segments_{k}"] = np.array(sorted(r["segments"]), int).reshape(-1, 2)
    (INPUTS / dataset).mkdir(parents=True, exist_ok=True)
    np.savez_compressed(INPUTS / dataset / f"{gene}.npz", **arrays)


def load_inputs(gene, dataset):
    path = INPUTS / dataset / f"{gene}.npz"
    if not path.exists():
        return []
    data = np.load(path)
    return [dict(path=str(p), pwm=data[f"pwm_{k}"], pair_matrix=data[f"pair_{k}"],
                 predicted={tuple(x) for x in data[f"predicted_{k}"].tolist()},
                 segments=[tuple(x) for x in data[f"segments_{k}"].tolist()])
            for k, p in enumerate(data["paths"])]


def describe(record, ref, el, offset):
    """Predicted and shared pairs inside the mapped window, and whether a segment contains it."""
    lo, hi = offset, offset + ref["span"] - 1
    inside = {(i, j) for i, j in record["predicted"] if lo <= i and j <= hi}
    mapped = {(offset + i - el["start"], offset + j - el["start"]) for i, j in el["pairs"]}
    within = any(s <= lo and hi <= e for s, e in record["segments"])
    return len(inside), len(inside & mapped), within


def run(permutations=PERMUTATIONS, export=False):
    sys.path.insert(0, str(HERE))
    os.chdir(HERE.parents[2])                       # repository root, as pdb_concordance expects
    raw = RAW_PROBE.exists()
    if export and not raw:
        raise SystemExit("--export-inputs needs the query-free models (external drive)")
    if raw:
        import pdb_concordance as pc                # loads the Jolma2020 spreadsheet
        print("inputs: query-free models and MSAs (external drive)")
    else:
        pc = None
        print(f"inputs: {INPUTS.relative_to(HERE.parents[2])} (models not available or --from-inputs)")

    rng = np.random.default_rng(SEED)
    rows = list(csv.DictReader(open(HERE / "pdb_concordance_reselected.tsv"), delimiter="\t"))
    results = []
    for row in rows:
        el = element(row["pdb_dbn"])
        ref = prepare(el, row["pdb_seq"])
        for dataset in DATASETS:
            records = load_models(pc, row["rbp"], dataset) if raw else load_inputs(row["rbp"], dataset)
            if export:
                save_inputs(records, row["rbp"], dataset)
            base = dict(rbp=row["rbp"], pdb=row["pdb"], dataset=DATASET_LABEL[dataset],
                        element=f"{el['start'] + 1}-{el['end'] + 1}", element_pairs=len(ref["pi"]),
                        element_unpaired=len(ref["un"]), n_models=len(records))
            if not records:
                results.append({**base, "outcome": "not assayed"})
                continue
            models = [(r["pwm"], r["pair_matrix"]) for r in records]
            joint, m, k = best_over(models, ref, "joint")
            if m is None:
                results.append({**base, "outcome": "element longer than every model"})
                continue
            seq_best, seq_m, seq_k = best_over(models, ref, "sequence")
            null_joint, null_seq = null_maxima(models, ref, rng, permutations)
            offsets, sequence, pairing, _ = scores(*models[m], ref)
            offset = int(offsets[k])
            predicted, shared, within = describe(records[m], ref, el, offset)
            p_joint = p_value(joint, null_joint)
            if p_joint >= ALPHA:
                outcome = "no corresponding candidate"
            elif predicted >= 2:
                outcome = "corresponding, stem predicted"
            else:
                outcome = "corresponding, no stem predicted"
            results.append({**base, "best_model": records[m]["path"], "best_offset": offset + 1,
                            "sequence_score": round(float(sequence[k]), 3),
                            "pairing_score": round(float(pairing[k]), 3),
                            "joint_score": round(joint, 3), "p_joint": round(p_joint, 4),
                            "best_sequence_model": records[seq_m]["path"],
                            "best_sequence_score": round(seq_best, 3),
                            "p_sequence": round(p_value(seq_best, null_seq), 4),
                            "within_segment": within, "predicted_pairs_in_window": predicted,
                            "shared_with_pdb": shared, "outcome": outcome})
            print(f"{row['rbp']:8s} {DATASET_LABEL[dataset]:11s} models={len(records):3d} "
                  f"joint={joint:.3f} p={p_joint:.4f} -> {outcome}")

    fields = ["rbp", "pdb", "dataset", "element", "element_pairs", "element_unpaired", "n_models",
              "best_model", "best_offset", "sequence_score", "pairing_score", "joint_score", "p_joint",
              "best_sequence_model", "best_sequence_score", "p_sequence", "within_segment",
              "predicted_pairs_in_window", "shared_with_pdb", "outcome"]
    with open(HERE / "pdb_correspondence.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", restval="")
        writer.writeheader()
        writer.writerows(results)
    write_summary(results, permutations)
    if raw:
        write_larp7(pc.XLSX)
    else:
        print("LARP7: spreadsheet not available; keeping the committed larp7_construct.tsv")


def write_summary(results, permutations):
    lines = ["# PDB correspondence: all query-free candidates", "",
             f"Rule and null as in `pdb_correspondence.py` ({permutations} permutations, seed {SEED}, "
             f"alpha {ALPHA}, pairable if >= {PAIRABLE} of sequences). Element positions are 1-based "
             "on the PDB RNA; offsets are 1-based matched columns of the model.", "",
             "| RBP | PDB | dataset | models | joint (seq / pair) | p | p (seq only) | in a segment | "
             "predicted / shared pairs | outcome |", "|---|---|---|---|---|---|---|---|---|---|"]
    for r in results:
        if "joint_score" not in r:
            lines.append(f"| {r['rbp']} | {r['pdb']} | {r['dataset']} | {r['n_models']} | | | | | | {r['outcome']} |")
            continue
        lines.append(f"| {r['rbp']} | {r['pdb']} | {r['dataset']} | {r['n_models']} | "
                     f"{r['joint_score']} ({r['sequence_score']} / {r['pairing_score']}) | {r['p_joint']} | "
                     f"{r['p_sequence']} | {'yes' if r['within_segment'] else 'no'} | "
                     f"{r['predicted_pairs_in_window']} / {r['shared_with_pdb']} | {r['outcome']} |")
    (HERE / "pdb_correspondence.md").write_text("\n".join(lines) + "\n")


def write_larp7(xlsx):
    import pandas as pd
    found = []
    for sheet, frame in pd.read_excel(xlsx, sheet_name=None, header=None).items():
        for index, values in frame.iterrows():
            cells = [str(v) for v in values.tolist() if str(v) != "nan"]
            if any("LARP7" in c.upper() for c in cells):
                found.append([sheet, index] + cells)
    with open(HERE / "larp7_construct.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sheet", "row", "cells"])
        writer.writerows(found)
    print(f"LARP7: {len(found)} spreadsheet rows written to larp7_construct.tsv")


# ---------- self-test on synthetic families ----------
def synthetic(rng, n=400, width=40, stem=None, loop=None, start=10):
    """Aligned family of random sequences, optionally with a covarying stem and a fixed loop."""
    seqs = rng.choice(list(ORD), size=(n, width))
    if stem:
        complement = {"A": "U", "U": "A", "G": "C", "C": "G"}
        for s in seqs:
            arm = rng.choice(list(ORD), size=stem)
            s[start:start + stem] = arm
            s[start + stem + len(loop):start + 2 * stem + len(loop)] = [complement[b] for b in arm[::-1]]
    if loop and "N" not in loop:
        seqs[:, start + (stem or 0):start + (stem or 0) + len(loop)] = list(loop)
    rows = ["".join(s) for s in seqs]
    counts = np.array([[sum(r[j] == b for r in rows) for j in range(width)] for b in ORD], float)
    return counts / counts.sum(0), pairability(rows, list(range(width)))


def selftest():
    rng = np.random.default_rng(1)
    # the PDB element: a 6-bp stem with a UGCAU loop, stem sequence unrelated to the family's
    pdb_seq = "GGAUCCUGCAUGGAUCC"
    dbn = "((((((.....))))))"
    el = element(dbn)
    assert (el["start"], el["end"], el["stack"], len(el["pairs"]), len(el["unpaired"])) == (0, 16, 6, 6, 5)
    ref = prepare(el, pdb_seq)
    cases = {
        "hairpin with the loop": [synthetic(rng, stem=6, loop="UGCAU")],
        "loop without a stem": [synthetic(rng, loop="UGCAU")],
        "stem without the loop": [synthetic(rng, stem=6, loop="NNNNN")],
        "random family": [synthetic(rng)],
    }
    results = {}
    for name, models in cases.items():
        joint = best_over(models, ref, "joint")[0]
        seq = best_over(models, ref, "sequence")[0]
        null_joint, null_seq = null_maxima(models, ref, np.random.default_rng(2), permutations=200)
        results[name] = (p_value(joint, null_joint), p_value(seq, null_seq))
        print(f"{name:24s} joint={joint:.3f} p={results[name][0]:.3f}  sequence={seq:.3f} p={results[name][1]:.3f}")
    assert results["hairpin with the loop"][0] < ALPHA, "a corresponding hairpin must be detected"
    assert results["loop without a stem"][1] < ALPHA, "a recovered loop must show in the sequence score"
    assert results["loop without a stem"][0] >= ALPHA, "a loop alone is not a corresponding stem-loop"
    assert results["stem without the loop"][0] >= ALPHA, "a stem alone is not a corresponding stem-loop"
    assert results["random family"][0] >= ALPHA, "a random family must not correspond"
    # nested structure and ties: the 5'-most of two equal helices
    el2 = element("((((....))))..((((....))))")
    assert (el2["start"], el2["end"]) == (0, 11)
    print("self-test passed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--selftest", action="store_true", help="check the method on synthetic data")
    parser.add_argument("--permutations", type=int, default=PERMUTATIONS)
    parser.add_argument("--export-inputs", action="store_true",
                        help="also write the per-model inputs to correspondence_inputs/ (needs the models)")
    parser.add_argument("--from-inputs", action="store_true",
                        help="read correspondence_inputs/ even when the models are available")
    args = parser.parse_args()
    if args.from_inputs:
        RAW_PROBE = Path("/nonexistent")
    selftest() if args.selftest else run(args.permutations, export=args.export_inputs)
