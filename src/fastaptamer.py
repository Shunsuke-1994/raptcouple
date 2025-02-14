# wrapper of fastaptamer
import subprocess


def fastaptamer_count(fasta_file, fasta_count_file):
    cmd = f"""fastaptamer_count -i {fasta_file} -o {fasta_count_file}"""
    res = subprocess.run(cmd, shell=True, capture_output=True)
    return res

def fastaptamer_compare(outtxt, fasta1, fasta2, quiet = False, output_all = False):
    q = "-q" if quiet else ""
    a = "-a" if output_all else ""
    cmd = f"""fastaptamer_compare -x {fasta1} -y {fasta2} -o {outtxt} {q} {a}"""
    res = subprocess.run(cmd, shell=True, capture_output=True)
    return res

def fastaptamer_enrich(outtxt, *fastas, quiet = False, verbose = False):
    q = "-q" if quiet else ""
    v = "-v" if verbose else ""
    if len(fastas) == 2:
        in_options = f"-x {fastas[0]} -y {fastas[1]}"
    elif len(fastas) == 3:
        in_options = f"-x {fastas[0]} -y {fastas[1]} -z {fastas[2]}"
    else: 
        raise ValueError("fastaptamer_enrich takes 2 or 3 fasta files")
    
    cmd = f"""fastaptamer_enrich {in_options} -o {outtxt} {q} {v}"""
    res = subprocess.run(cmd, shell=True, capture_output=True)
    return res
