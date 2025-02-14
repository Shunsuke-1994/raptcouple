import sys
import os 
# sys.path.append("/Users/sumishunsuke/Desktop/RNA/VariationalNeuralPotts")
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(current_dir))
from src.coupling import plmc
import subprocess

def main(args):
    base, _ = os.path.splitext(args.input_fasta)
    param_file = base + args.suffix + ".model_params"
    coupling_file = base + args.suffix + ".coupling"

    # not essential, but plmc requires target_id.
    # get the first id from input_fasta
    if args.target is not None:
        target_id = args.target
    else:
        cmd_get_id = f"""head -n1 {args.input_fasta}"""
        res_id = subprocess.run(cmd_get_id, shell = True, capture_output =True)
        target_id = res_id.stdout.decode().strip().replace(">", "")
    # print(target_id)

    cmd = plmc(
            fasta_file=args.input_fasta,
            target=target_id,
            param_file=param_file,
            coupling_file=coupling_file,
            vocab = args.vocab,
            threshold=args.threshold,
            print_result=args.print_result
        )
    res = subprocess.run(cmd, shell=True, capture_output=True)
    print(res.stdout.decode())
    print(res.stderr.decode())
    return cmd

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_fasta", type=str, required=True, help="input fasta file (=alignment)")
    parser.add_argument("--target", type=str, required=False, help="specified target", default=None)    
    parser.add_argument("--threshold", default = 0.05, type=float, help="similarity threshold for plmc (default: 0.05)")
    parser.add_argument("--vocab", default="AUGC.", type = str, help="vocabulary for plmc (default: AUGC.)")
    parser.add_argument("--iters", default=200, type=int, help="number of iterations for plmc (default: 200)")
    parser.add_argument("--suffix", default="", type=str, help="suffix for output files (default: '')")
    parser.add_argument("--print_result", action="store_true", help="print result of plmc")
    args = parser.parse_args()
    print(main(args))
