import argparse
import json
import os,sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from src.structure import fold
from src.plmc import read_params

def main():
    parser = argparse.ArgumentParser(description="Fold RNA using coupling scores")
    parser.add_argument("--coupling", required=True, help="Path to the coupling file")
    parser.add_argument("--min_loop_len", type=int, default=3, help="Minimum loop length")
    parser.add_argument("--z_threshold", type=float, default=3, help="Z-score threshold for coupling")
    parser.add_argument("--no-nolp", dest="nolp", action="store_false",
                        help="Disable no-lonely-pair constraint (allow isolated base pairs)")
    parser.add_argument("--canonical-only", dest="canonical_only", action="store_true",
                        help="Keep only A-U, G-C and G-U pairs of the focus sequence "
                             "(by default pairs are not restricted, as in the paper)")
    parser.add_argument("--output", required=True, help="Output file to save the results")
    args = parser.parse_args()

    params = read_params(args.coupling)

    rna, energy, ss, score_table = fold(params, min_loop_length=args.min_loop_len,
                                        threshold=args.z_threshold, nolp=args.nolp,
                                        sanity_check=args.canonical_only)

    result = {
        "rna": rna,
        "energy": energy,
        "secondary_structure": ss,
        "model_params": args.coupling,
        "min_loop_len": args.min_loop_len,
        "z_threshold": args.z_threshold,
        "nolp": args.nolp,
        "canonical_only": args.canonical_only
    }
    with open(args.output, "w") as f:
        json.dump(result, f, indent=4)

if __name__ == "__main__":
    main()
