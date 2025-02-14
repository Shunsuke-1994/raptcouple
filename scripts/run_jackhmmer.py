# run jackhmmer from a query ID in a fasta file
import sys
sys.path.append(".")
from src.hmmer import load_fasta, jackhmmer, save_hmm, save_msa
from src.coupling import plmc
import subprocess
import os

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
    run jackhmmer.
    """
    sequence = fetch_seq_by_id(args.fasta, args.target)
    print("target\t:", args.target)
    print("sequence:", sequence)
    selex_data = load_fasta(args.fasta)
    iterations = jackhmmer(
        selex_data = selex_data,
        sequence = sequence,
        max_iters = args.iters,
        T = args.T,
        domT = args.domT,
        incT = args.incT,
        incdomT = args.incdomT,
        F1 = args.F1, 
        F2 = args.F2,
        F3 = args.F3
        )
    # try:
    #     # for unkwown reason, likely to fail
    #     hmm_file = save_hmm(iterations[-1], os.path.join(args.save_dir, args.prefix+args.target + ".hmm"), name=args.target)
    #     print("HMM saved\t:", hmm_file)
    # except:
    #     print("Failed to save HMM file")
    #     pass
    try:
        msa_file = save_msa(iterations[-1], os.path.join(args.save_dir, args.prefix+args.target + ".msa"))
        print("MSA saved\t:", msa_file)
    except:
        print("Failed to save MSA file")
        pass
    # return hmm_file, msa_file
    # return msa_file, hmm_file


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run jackhmmer and plmc on a fasta file")
    parser.add_argument("--fasta", help="fasta file", required=True)
    parser.add_argument("--target", help="target id of fasta file", required=True)
    parser.add_argument("--save_dir", help="directory to save results", required=True)
    parser.add_argument("--prefix", help="prefix of output files", default = "")


    parser.add_argument("--iters", help="number of iterations for jackhmmer (default: 10)", default = 10, type = int)

    parser.add_argument("--F1", help="threshold F1 for jackhmmer (default: 0.02)", default = 0.02, type = float)
    parser.add_argument("--F2", help="threshold F2 for jackhmmer (default: 1e-3)", default = 1e-3, type = float)
    parser.add_argument("--F3", help="threshold F3 for jackhmmer (default: 1e-4)", default = 1e-4, type = float)

    parser.add_argument("--T", help="threshold T for jackhmmer (default: 5)", default =5, type = float)
    parser.add_argument("--domT", help="threshold domT for jackhmmer (default: 5)", default =5, type = float)
    parser.add_argument("--incT", help="threshold incT for jackhmmer (default: 5)", default =5, type = float)
    parser.add_argument("--incdomT", help="threshold incdomT for jackhmmer (default: 5)", default =5, type = float)

    parser.add_argument("--print_result", action="store_true")
    args = parser.parse_args()
    main(args)
