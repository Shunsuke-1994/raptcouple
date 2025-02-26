import argparse
import os
import json
from src.structure import fold
from src.coupling import read_params

def main():
    parser = argparse.ArgumentParser(description="Fold RNA using coupling scores")
    parser.add_argument("--coupling", required=True, help="Path to the coupling file")
    parser.add_argument("--min_loop_len", type=int, default=3, help="Minimum loop length")
    parser.add_argument("--z_threshold", type=float, default=2, help="Z-score threshold for coupling")
    parser.add_argument("--output", required=True, help="Output file to save the results")
    args = parser.parse_args()

    # Read the coupling parameters
    params = read_params(args.coupling)

    # Perform RNA folding
    rna, energy, ss, score_table = fold(params, min_loop_length=args.min_loop_len, threshold=args.z_threshold)

    # Save the results
    result = {
        "rna": rna,
        "energy": energy,
        "secondary_structure": ss,
        "model_params": args.coupling,
        "min_loop_len": args.min_loop_len,
        "z_threshold": args.z_threshold
    }
    with open(args.output, "w") as f:
        json.dump(result, f, indent=4)

if __name__ == "__main__":
    main()
