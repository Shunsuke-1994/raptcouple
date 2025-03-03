import argparse
import numpy as np
import sys, os
sys.path.append("./")
from src.potts import PottsModel
from src.util import onehot2seq  # to convert one-hot state to sequence string if needed

def parse_args():
    parser = argparse.ArgumentParser(description="Simulated annealing using PottsModel.")
    parser.add_argument("--param_file", type=str, required=True, help="Path to parameter file.")
    parser.add_argument("--T_init", type=float, default=10.0, help="Initial temperature.")
    parser.add_argument("--T_end", type=float, default=0.1, help="Final temperature.")
    parser.add_argument("--max_steps", type=int, default=1000, help="Maximum number of steps.")
    parser.add_argument("--random_init", action="store_true", help="Randomly initialize spins.")
    parser.add_argument("--random_seed", type=int, default=42, help="Random seed for reproducibility.")
    parser.add_argument("--print_log", action="store_true", help="Print log during annealing.")
    return parser.parse_args()

def main():
    args = parse_args()
    model = PottsModel.build_from_file(args.param_file)
    print("Initial Energy:", model.compute_energy())
    print("Initial Sequence:", onehot2seq(model.spins))

    # Perform simulated annealing.
    model.sim_anneal(
        T_init=args.T_init,
        T_end=args.T_end,
        max_steps=args.max_steps,
        random_init_state=args.random_init,
        random_seed=args.random_seed,
        print_log=args.print_log
    )
    
    # After annealing, compute the energy and obtain sequence representation.
    final_energy = model.compute_energy()
    final_sequence = onehot2seq(model.spins)
    
    print("Final Energy:", final_energy)
    print("Final Sequence:", final_sequence)

if __name__ == "__main__":
    main()
