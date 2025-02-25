# run jackhmmer from a query ID in a fasta file
import sys
sys.path.append(".")
from src.hmmer import load_fasta, jackhmmer, save_hmm, save_msa
from src.coupling import plmc
import subprocess
import os, yaml

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
    config = yaml.safe_load(open(args.config, "r"))["MSA_parameters"]
    sequence = fetch_seq_by_id(config["all_fasta"], config["target_id"])
    print("target\t:", config["target_id"])
    print("sequence:", sequence)
    selex_data = load_fasta(config["all_fasta"])
    iterations = jackhmmer(
        selex_data = selex_data,
        sequence = sequence,
        max_iters = config["iters"],
        T = config["T"],
        domT = config["domT"],
        incT = config["incT"],
        incdomT = config["incdomT"],
        F1 = config["F1"],
        F2 = config["F2"],
        F3 = config["F3"]
        )
    try:
        msa_file = save_msa(iterations[-1], os.path.join(config["save_dir"], config["prefix"]+config["target_id"] + ".msa"))
        print("MSA saved\t:", msa_file)
    except:
        print("Failed to save MSA file")
        pass


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run jackhmmer and plmc on a fasta file")
    parser.add_argument("--config", help="config file", required=True)
    args = parser.parse_args()
    main(args)
