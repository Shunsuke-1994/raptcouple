# nhmmer
import sys
sys.path.append(".")
from src.hmmer import load_fasta, nhmmer, save_hmm, save_msa
from plmc import plmc
import subprocess
import matplotlib.pyplot as plt

def fetch_seq_by_id(fasta_file, target_id):
    """
    get target sequence
    """
    cmd = "grep -A 1 " + target_id + " " + fasta_file + " | tail -n 1"
    res = subprocess.run(cmd, shell = True, capture_output = True)
    sequence = res.stdout.decode("utf-8").strip()
    return sequence


def main(args):
    """
    run nhmmer
    """
    sequence = fetch_seq_by_id(args.fasta, args.target)
    selex_data = load_fasta(args.fasta)
    hits = nhmmer(
        selex_data = selex_data,
        sequence = sequence,
        T = args.T,
        domT = args.domT,
        incT = args.incT,
        incdomT = args.incdomT,
        F1 = args.F1, 
        F2 = args.F2,
        F3 = args.F3
        )
    
    # if args.plot:
    plt.figure()
    plt.hist([hit.score for hit in hits], bins = 50)
    plt.grid(alpha = 0.3)
    plt.title(f"F1={args.F1}, F2={args.F2}, F3={args.F1}")
    plt.xlabel("bit score")
    plt.xlim(0, 100)
    plt.show()

    return 

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run jackhmmer and plmc on a fasta file")
    parser.add_argument("--fasta", help="fasta file", required=True)
    parser.add_argument("--target", help="target id of fasta file", required=True)
    parser.add_argument("--F1", help="threshold F1 for jackhmmer (default: 0.02)", default = 0.02, type = float)
    parser.add_argument("--F2", help="threshold F2 for jackhmmer (default: 1e-3)", default = 1e-3, type = float)
    parser.add_argument("--F3", help="threshold F3 for jackhmmer (default: 1e-4)", default = 1e-4, type = float)

    parser.add_argument("--T", help="threshold T for jackhmmer (default: 5)", default =5, type = float)
    parser.add_argument("--domT", help="threshold domT for jackhmmer (default: 5)", default =5, type = float)
    parser.add_argument("--incT", help="threshold incT for jackhmmer (default: 5)", default =5, type = float)
    parser.add_argument("--incdomT", help="threshold incdomT for jackhmmer (default: 5)", default =5, type = float)

    # parser.add_argument("--save_dir", help="directory to save results", required=True)
    # parser.add_argument("--plot", help="flag for plotting the result", action="store_true", required=True)
    args = parser.parse_args()
    main(args)
