import argparse
import numpy as np
import sys, os
sys.path.append("./")
from src.potts import PottsModel
from src.plmc import read_params
from src.util import onehot2seq, seq2onehot

def parse_args():
    parser = argparse.ArgumentParser(description="Predict mutation effects using PottsModel.")
    parser.add_argument("--param_file", type=str, required=True, help="Path to parameter file.")
    parser.add_argument("--mutations", type=str, required=False, help="Comma-separated list of mutations (e.g., A15G,C22U).")
    parser.add_argument("--mutations_file", type=str, required=False, help="File with one mutation per line.")
    return parser.parse_args()

def parse_mutation(mutation_str, alphabet):
    """Parse a mutation string like 'A15G' or 'A21.' into (from_idx, pos, to_idx).

    Indices follow the state order of the model (``alphabet``). T is accepted
    for U and U for T, so mutations can be written in either nucleic acid.
    """
    def state(nuc):
        if nuc not in alphabet:
            if nuc == "T" and "U" in alphabet:
                nuc = "U"
            elif nuc == "U" and "T" in alphabet:
                nuc = "T"
            else:
                raise ValueError(f"Nucleotide {nuc!r} in mutation {mutation_str} not in alphabet {alphabet}")
        return alphabet.index(nuc)

    from_idx = state(mutation_str[0])
    to_idx = state(mutation_str[-1])
    pos = int(mutation_str[1:-1]) - 1  # Convert to 0-indexed
    return from_idx, pos, to_idx

def get_alphabet(params):
    """State order of the Potts model as stored by plmc (e.g. 'AUGC.').

    The model's states are ordered as in this string, so mutation indices must
    be taken from it rather than from an alphabetical order.
    """
    alphabet = params["alphabet"]
    if len(alphabet) not in (4, 5):
        raise ValueError(f"Unsupported alphabet: {alphabet}")
    return alphabet

def main():
    args = parse_args()
    
    # Build model from parameter file
    model = PottsModel.build_from_file(args.param_file)
    
    # Determine alphabet (state order of the model)
    alphabet = get_alphabet(read_params(args.param_file))
    print(f"Using alphabet: {alphabet}")
    
    # Get target sequence in the same alphabet
    target_seq = onehot2seq(model.spins, is_dna=("T" in alphabet), is_gapped=("." in alphabet))
    print(f"Target sequence: {target_seq}")
    
    # Parse mutations
    mutations = []
    if args.mutations:
        mutations.extend(args.mutations.split(','))
    if args.mutations_file:
        with open(args.mutations_file, 'r') as f:
            mutations.extend([line.strip() for line in f if line.strip()])
    
    if not mutations:
        print("No mutations specified. Use --mutations or --mutations_file.")
        sys.exit(1)
    
    # Predict effect for each mutation
    print("\nPredicting mutation effects:")
    print("Mutation\tEnergy Change")
    print("-" * 50)
    
    for mut_str in mutations:
        try:
            from_idx, pos, to_idx = parse_mutation(mut_str, alphabet)
            
            # Verify the original nucleotide matches the target sequence
            if target_seq[pos] != alphabet[from_idx]:
                print(f"Warning: Original nucleotide in mutation {mut_str} doesn't match target sequence ({target_seq[pos]} at position {pos+1})")
                continue
            
            # Compute energy change
            delta_energy = model.compute_delta_energy([(from_idx, pos, to_idx)])
                        
            print(f"{mut_str}\t{delta_energy:.4f}")
            
        except (ValueError, IndexError) as e:
            print(f"Error processing mutation {mut_str}: {e}")
    
if __name__ == "__main__":
    main()
