import subprocess
import os
import yaml
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from src.fastaptamer import fastaptamer_count, fastaptamer_enrich

def cutadapt_of_all_fasta(config, quiet=False):
    adapter_5 = config["adapter_5"]
    adapter_3 = config["adapter_3"]
    N_random = config["N_random"]

    option = ""
    if adapter_3 != "": option += f"-a {adapter_3} "
    if adapter_5 != "": option += f"-g {adapter_5} "

    for fasta_file in config["fasta_annotation"].keys():

        fasta_trim_file = os.path.join(config["data_dir"], fasta_file.replace(".fa", ".trim.fa"))
        
        cmd = f"""cutadapt {option} -l {str(N_random)} {os.path.join(config["data_dir"], fasta_file)} -o {fasta_trim_file}"""
        res = subprocess.run(cmd, shell=True, capture_output=True)
        if not quiet:
            print(res.stdout.decode("utf-8"))
            print(res.stderr.decode("utf-8"))
    return [os.path.join(config["data_dir"], fasta_file.replace(".fa", ".trim.fa")) for fasta_file in config["fasta_annotation"].keys()]

def fastaptamer_count_of_all_trimfasta(config, quiet=False):
    for fasta in config["fasta_annotation"].keys():
        fasta_trim = os.path.join(config["data_dir"], fasta.replace(".fa", ".trim.fa"))

        res = fastaptamer_count(fasta_trim, fasta_trim.replace(".fa", ".count.fa"))
        if not quiet:
            print(res.stdout.decode("utf-8"))
            print(res.stderr.decode("utf-8"))
    return [os.path.join(config["data_dir"], fasta.replace(".fa", ".trim.count.fa")) for fasta in config["fasta_annotation"].keys()]

def main(args):
    config = yaml.safe_load(open(args.config, "r"))
    cutadapt_of_all_fasta(config, quiet=args.quiet)
    fastas_trim_count = fastaptamer_count_of_all_trimfasta(config, quiet=args.quiet)   
    fastaptamer_enrich(args.outtxt, *fastas_trim_count, quiet=args.quiet)
    return 

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=str, help="selex config file path")
    parser.add_argument("--outtxt", required=True, type=str, help="output file of fastaptamer_enrich")
    parser.add_argument("--quiet", "-q", action="store_true", help="quiet mode")
    args = parser.parse_args()

    main(args)