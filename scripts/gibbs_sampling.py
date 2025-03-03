import argparse
import numpy as np
import sys 
sys.path.append(".")
from src.potts import PottsModel
from src.util import onehot2seq  # assuming onehot2seq converts a 2D one-hot array to a sequence string

def parse_args():
    parser = argparse.ArgumentParser(description="Generate sequences via Gibbs sampling and write to FASTA format.")
    parser.add_argument("--param_file", type=str, required=True, help="Path to parameter file.")
    parser.add_argument("--out", type=str, required=True, help="Output FASTA file path.")
    parser.add_argument("--n_samples", type=int, default=10, help="Number of samples to generate.")
    parser.add_argument("--random_seed", type=int, default=42, help="Random seed for reproducibility.")
    return parser.parse_args()

def main():
    args = parse_args()
    # Build model from parameter file; build_from_file handles is_dna/is_gapped correctly.
    model = PottsModel.build_from_file(args.param_file)
    
    # Generate samples using Gibbs sampling.
    samples = model.gibbs_sampling(n_samples=args.n_samples, random_seed=args.random_seed)
    
    # Write to FASTA file with energy info.
    with open(args.out, "w") as f:
        for i, sample in enumerate(samples):
            energy = model.eval_energy(sample)
            seq = onehot2seq(sample)
            f.write(f">sample_{i} energy:{energy:.4f}\n{seq}\n")
    print(f"FASTA file saved to {args.out}")

if __name__ == "__main__":
    main()
