"""Query-free RaptCouple: discover motifs and structures in a SELEX pool without a query.

The target-guided pipeline (run_jackhmmer.py -> train_potts.py) starts from a query the
user supplies. This script generates the queries from the pool itself and then runs the
same steps for each of them:

    pool -> seeds (STREME or vsearch) -> jackhmmer -> Potts model -> motif candidates -> fold

Usage:
    python scripts/run_query_free.py --config <config.yaml> --arm streme --out <dir>
    python scripts/run_query_free.py --config <config.yaml> --arm vsearch --out <dir>

The config is the same YAML used by the other scripts; `MSA_parameters` and
`Potts_parameters` are read from it, and `MSA_parameters.target_id` is ignored because the
seeds replace it. `Query_free_parameters` is optional:

    Query_free_parameters:
      n_seeds: 10           # STREME motifs, or vsearch clusters, to use
      identity: 0.7         # vsearch clustering identity
      abundance_round: 7    # count reads from this round only (default: the latest round)
      top_n: 3              # motif candidates kept per seed, by prevalence
      streme_dir: ...       # reuse an existing STREME output directory
      uc: ...               # reuse an existing vsearch .uc file

Needs STREME (MEME Suite) for the STREME arm and vsearch for the vsearch arm.
"""
import argparse
import csv
import json
import os
import sys
from collections import Counter

import yaml

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from src import motif, seeds
from src.hmmer import jackhmmer, load_fasta, save_msa
from src.plmc import fit, read_params
from src.structure import fold


def parse_args():
    parser = argparse.ArgumentParser(description="Query-free motif and structure discovery")
    parser.add_argument("--config", required=True, help="Path to the config YAML")
    parser.add_argument("--arm", choices=["streme", "vsearch"], required=True,
                        help="How to generate seeds")
    parser.add_argument("--out", required=True, help="Output directory")
    parser.add_argument("--pool", help="Pool FASTA (default: MSA_parameters.all_fasta)")
    return parser.parse_args()


def read_msa(path):
    """[(id, sequence)] of an aligned FASTA."""
    records, name, sequence = [], None, ""
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    records.append((name, sequence))
                name, sequence = line[1:].strip(), ""
            else:
                sequence += line.strip()
    if name is not None:
        records.append((name, sequence))
    return records


def consensus_of(sequences):
    if not sequences:
        return ""
    width = len(sequences[0])
    return "".join(Counter(s[i] for s in sequences).most_common(1)[0][0] for i in range(width))


def build_model(seed, pool, msa_params, potts_params, out_dir, index):
    """jackhmmer + plmc for one seed. Returns (model params path, MSA records)."""
    prefix = os.path.join(out_dir, f"seed{index:02d}")
    msa_path, param_path = prefix + ".msa", prefix + ".model_params"
    if not os.path.exists(param_path):
        iterations = jackhmmer(
            selex_data=pool, sequence=seed.sequence,
            max_iters=msa_params.get("iters", 10),
            T=msa_params.get("T", 5), domT=msa_params.get("domT", 5),
            incT=msa_params.get("incT", 5), incdomT=msa_params.get("incdomT", 5),
            F1=msa_params.get("F1", 0.02), F2=msa_params.get("F2", 1e-3),
            F3=msa_params.get("F3", 1e-4), print_result=False)
        save_msa(iterations[-1], msa_path)
        records = read_msa(msa_path)
        # jackhmmer renames the seed; fall back to the first row if it is not found
        target = next((name for name, _ in records if seed.id in name), records[0][0])
        fit(fasta_file=msa_path, target=target, param_file=param_path,
            coupling_file=prefix + ".coupling",
            vocab=potts_params.get("vocab", "AUGC."),
            threshold=potts_params.get("sim_threshold", 0.05), print_result=False)
    return param_path, read_msa(msa_path)


def main():
    args = parse_args()
    config = yaml.safe_load(open(args.config))
    msa_params = config.get("MSA_parameters", {})
    potts_params = config.get("Potts_parameters", {})
    options = config.get("Query_free_parameters", {}) or {}
    pool_fasta = args.pool or msa_params["all_fasta"]
    os.makedirs(args.out, exist_ok=True)

    n_seeds = options.get("n_seeds", 10)
    abundance_round = options.get("abundance_round")

    # 1. Seeds, without using any reference motif.
    if args.arm == "vsearch":
        uc = options.get("uc") or os.path.join(args.out, "clusters.uc")
        if not os.path.exists(uc):
            print(f"clustering {pool_fasta} with vsearch ...", flush=True)
            seeds.run_vsearch(pool_fasta, uc, identity=options.get("identity", 0.7))
        found = seeds.cluster_seeds(pool_fasta, uc, top_n=n_seeds, round_number=abundance_round)
    else:
        streme_dir = options.get("streme_dir") or os.path.join(args.out, "streme")
        if not os.path.exists(os.path.join(streme_dir, "streme.txt")):
            print(f"running STREME on {pool_fasta} ...", flush=True)
            seeds.run_streme(pool_fasta, streme_dir, n_motifs=n_seeds)
        found = seeds.motif_seeds(pool_fasta, streme_dir, n_motifs=n_seeds)
    found = seeds.unique(found)
    print(f"{len(found)} seed(s) from the {args.arm} arm", flush=True)

    # 2. The standard pipeline for each seed, then its motif candidates and structure.
    print(f"loading search pool {pool_fasta} ...", flush=True)
    pool = load_fasta(pool_fasta)
    reads = motif.load_reads(pool_fasta)
    rows = []
    for index, seed in enumerate(found, 1):
        param_path, records = build_model(seed, pool, msa_params, potts_params, args.out, index)
        params = read_params(param_path)
        _, energy, structure, _ = fold(params, min_loop_length=3, threshold=3)
        candidates = motif.candidates(params, reads=reads, top_n=options.get("top_n", 3))
        for candidate in candidates:
            rows.append({
                "seed": index, "seed_id": seed.id, "seed_label": seed.label,
                "msa_depth": len(records), "motif": candidate["consensus"],
                "core_kmer": candidate["core_kmer"], "prevalence": candidate["prevalence"],
                "structure": structure, "n_base_pairs": structure.count("("),
                "fold_score": round(energy, 3),
                "consensus": consensus_of([s for _, s in records]).replace("T", "U"),
            })
        print(f"seed{index:02d}: depth={len(records)} motifs={len(candidates)} "
              f"bp={structure.count('(')}", flush=True)

    results = os.path.join(args.out, "results.tsv")
    fields = ["seed", "seed_id", "seed_label", "msa_depth", "motif", "core_kmer", "prevalence",
              "structure", "n_base_pairs", "fold_score", "consensus"]
    with open(results, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    with open(os.path.join(args.out, "settings.json"), "w") as handle:
        json.dump({"config": args.config, "arm": args.arm, "pool": pool_fasta,
                   "n_seeds": n_seeds, "abundance_round": abundance_round,
                   "top_n": options.get("top_n", 3)}, handle, indent=4)
    print(f"\nwrote {results} ({len(rows)} motif candidates)")


if __name__ == "__main__":
    main()
