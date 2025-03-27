import sys
import os 
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(current_dir))
from src.plmc import fit
import subprocess, yaml

def main(args):
    config = yaml.safe_load(open(args.config, "r"))["Potts_parameters"]
    base, _ = os.path.splitext(config["input_fasta"])
    param_file = base + config["suffix"] + ".model_params"
    coupling_file = base + config["suffix"] + ".coupling"
    target_id = config["target_id"]

    # cmd_get_id = f"""head -n1 {config["input_fasta"]}"""
    cmd_get_complete_id = f"""grep {target_id} {config["input_fasta"]}"""
    res_id = subprocess.run(cmd_get_complete_id, shell = True, capture_output =True)
    target_complete_id = res_id.stdout.decode().strip().replace(">", "")

    cmd = fit(
            fasta_file=config["input_fasta"],
            target=target_complete_id,
            param_file=param_file,
            coupling_file=coupling_file,
            vocab = config["vocab"],
            threshold=config["sim_threshold"],
            print_result=config["print_result"]
        )
    res = subprocess.run(cmd, shell=True, capture_output=True)
    print(res.stdout.decode())
    print(res.stderr.decode())
    return cmd

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", help="config file", required=True)
    args = parser.parse_args()
    print(main(args))
