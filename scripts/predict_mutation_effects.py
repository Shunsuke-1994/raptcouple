import argparse
import numpy as np
import sys, os
sys.path.append("./")
from src.potts import PottsModel
from src.util import onehot2seq, seq2onehot

def parse_args():
    parser = argparse.ArgumentParser(description="Predict mutation effects using PottsModel.")
    parser.add_argument("--param_file", type=str, required=True, help="Path to parameter file.")
    parser.add_argument("--mutations", type=str, required=False, help="Comma-separated list of mutations (e.g., A15G,C22U).")
    parser.add_argument("--mutations_file", type=str, required=False, help="File with one mutation per line.")
    return parser.parse_args()

def parse_mutation(mutation_str, alphabet):
    """Parse a mutation string like 'A15G'"""
    from_nuc = mutation_str[0]
    to_nuc = mutation_str[-1]
    pos = int(mutation_str[1:-1]) - 1  # Convert to 0-indexed
    
    # Validate nucleotides are in the alphabet
    if from_nuc not in alphabet or to_nuc not in alphabet:
        raise ValueError(f"Nucleotides in mutation {mutation_str} not found in alphabet {alphabet}")
    
    # Get the indices in the alphabet
    from_idx = alphabet.index(from_nuc)
    to_idx = alphabet.index(to_nuc)
    
    return from_idx, pos, to_idx

def get_alphabet_from_model(model):
    """Determine the alphabet based on the model parameters"""
    if model.num_states == 4:
        # RNA or DNA
        target_seq = onehot2seq(model.spins)
        if 'T' in target_seq:
            return "ACGT"
        else:
            return "ACGU"
    elif model.num_states == 5:
        # RNA or DNA with gaps
        target_seq = onehot2seq(model.spins)
        if 'T' in target_seq:
            return "ACGT."
        else:
            return "ACGU."
    else:
        raise ValueError(f"Unsupported number of states: {model.num_states}")

def main():
    args = parse_args()
    
    # Build model from parameter file
    model = PottsModel.build_from_file(args.param_file)
    
    # Determine alphabet
    alphabet = get_alphabet_from_model(model)
    print(f"Using alphabet: {alphabet}")
    
    # Get target sequence
    target_seq = onehot2seq(model.spins)
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
    print("Mutation\tEnergy Change\tPrediction")
    print("-" * 50)
    
    for mut_str in mutations:
        try:
            from_idx, pos, to_idx = parse_mutation(mut_str, alphabet)
            
            # Verify the original nucleotide matches the target sequence
            if target_seq[pos] != mut_str[0]:
                print(f"Warning: Original nucleotide in mutation {mut_str} doesn't match target sequence ({target_seq[pos]} at position {pos+1})")
                continue
            
            # Compute energy change
            delta_energy = model.compute_delta_energy([(from_idx, pos, to_idx)])
            
            # Determine effect based on energy change
            if delta_energy < 0:
                effect = "Beneficial"
            elif delta_energy > 0:
                effect = "Deleterious"
            else:
                effect = "Neutral"
            
            print(f"{mut_str}\t{delta_energy:.4f}\t{effect}")
            
        except (ValueError, IndexError) as e:
            print(f"Error processing mutation {mut_str}: {e}")
    
if __name__ == "__main__":
    main()
