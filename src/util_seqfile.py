import subprocess
from Bio import SeqIO
import gzip
from Bio import AlignIO

def get_nth_id(fasta, nth = 1):
    cmd_get_id = f"""head -n{2*nth-1} {fasta}"""
    res_id = subprocess.run(cmd_get_id, shell = True, capture_output =True)
    nth_id = res_id.stdout.decode().split("\n")[-2].strip().replace(">", "")
    return nth_id

def fastq2fasta(fastq, is_gzip = False):
    if is_gzip:
        fasta = fastq.replace(".fastq.gz", ".fa")
        records = SeqIO.parse(gzip.open(fastq, "rt") , "fastq")
    else:
        fasta = fastq.replace(".fastq", ".fa")
        records = SeqIO.parse(fastq, "fastq")

    # with open(fasta, "w") as f:
    #     for record in records:
    #         f.write(">" + record.id + "\n")
    #         f.write(str(record.seq) + "\n")
    SeqIO.write(records, fasta, "fasta")
    return fasta

def convert_revcom_fasta(file):
    records_revcom = []
    for record in SeqIO.parse(file, "fasta"):
        record.seq = record.seq.reverse_complement()
        records_revcom.append(record)
    SeqIO.write(records_revcom, file.replace(".fa", ".revcom.fa"),"fasta")
    return file.replace(".fa", ".revcom.fa")


import re 
def add_gap_to_ss_from_gapseq(ss, seq_gap):
    """
    a function for adding gaps to secondary structure from gapped sequence
    """
    ss_gap = ""
    count_nongap = 0
    for n_gap in seq_gap:
        if n_gap != ".":
            # print(count_nongap)
            ss_gap += ss[count_nongap] 
            count_nongap += 1

        else:
            ss_gap += "."

    return ss_gap

def write_sto_withSS_from_msa(msa_file, out_sto, ss):
    align = AlignIO.read(msa_file, "fasta")
    # target_aln_record = [aln for aln in align if target_id in aln.id][0]
    # target_aln_seq = str(target_aln_record.seq)

    # positions_dots = [m.start() for m in re.finditer("\.", target_aln_seq)]
    # ss_list = list(ss)
    # for i in positions_dots:
    #     ss_list.insert(i, ".")

    # aln_ss = "".join(ss_list).replace("(", "<").replace(")", ">")
    # align.column_annotations["secondary_structure"] = aln_ss
    align.column_annotations["secondary_structure"] = ss

    AlignIO.write(align, out_sto, "stockholm")
    return out_sto