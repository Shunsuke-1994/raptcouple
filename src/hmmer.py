# wrapper of pyhmmer    
import os
import subprocess
import pyhmmer
from pyhmmer import easel
from pyhmmer.hmmer import hmmalign
from math import ceil, floor

# see the user manual of HMMER3 for the details. 
# http://eddylab.org/software/hmmer3/3.1b2/Userguide.pdf

def load_fasta(fasta_file):
    with easel.SequenceFile(fasta_file, digital = True, alphabet = easel.Alphabet.rna()) as fasta:
        selex_data = fasta.read_block()
    return selex_data

def nhmmer(selex_data, sequence, T, domT, incT, incdomT, F1 = 0.02, F2 = 1e-3, F3 = 1e-5):
    original = easel.TextSequence(sequence = sequence)
    original_d = original.digitize(easel.Alphabet.rna())
    pli = pyhmmer.plan7.Pipeline(easel.Alphabet.rna(), T = T, domT = domT, incT = incT, incdomT = incdomT, F1=F1, F2=F2, F3=F3)
    hits = pli.search_seq(original_d, selex_data)
    return hits

def jackhmmer(selex_data, sequence, max_iters, T, domT, incT, incdomT, F1 = 0.02, F2 = 1e-3, F3 = 1e-5, print_result = True):
    original = easel.TextSequence(sequence = sequence)
    original_d = original.digitize(easel.Alphabet.rna())
    pli = pyhmmer.plan7.Pipeline(easel.Alphabet.rna(), T = T, domT = domT, incT = incT, incdomT = incdomT, F1=F1, F2=F2, F3=F3)
    search = pli.iterate_seq(original_d, selex_data)
    
    iterations = []
    for i in range(max_iters):
        iteration = next(search)
        iterations.append(iteration)
        # worst_score = min(hit.score for hit in iteration.hits if hit.included)
        if print_result:
            print(f"Iteration={iteration.iteration} Hits={len(iteration.hits):3} Included={len(iteration.hits.included):3} Converged={search.converged}")

        if search.converged:
            break

    return iterations

def save_hmm(iteration_res, hmm_file, name, description = None, accession = None):
    hmm = iteration_res.hmm
    hmm.name = bytes(name, 'utf-8')
    hmm.description = bytes(description, 'utf-8') if description else None
    hmm.accession = bytes(accession, 'utf-8') if accession else None

    # set cutoffs
    try:
        noise_score = max([hit.score for hit in iteration_res.hits if not hit.included])
        noise_domscore = max([hit.best_domain.score for hit in iteration_res.hits if not hit.included])
        hmm.cutoffs.noise = ceil(noise_score), ceil(noise_domscore)
    except:
        hmm.cutoffs.noise = 0, 0

    gathering_score = min([hit.score for hit in iteration_res.hits if hit.included])
    gathering_domscore = min([hit.best_domain.score for hit in iteration_res.hits if hit.included])
    hmm.cutoffs.gathering = floor(gathering_score), floor(gathering_domscore)


    with open(hmm_file, "wb") as f_hmm:
        hmm.write(f_hmm)
        
    return hmm_file

def save_msa(iteration_res, msa_file, file_format = "fasta"):
    assert file_format in ["fasta", "stockholm"], "file_format must be fasta or stockholm"


    if file_format == "fasta":
        with open(msa_file, "w") as f:
            textmsa = iteration_res.hits.to_msa(easel.Alphabet.rna())
            for name, aligned in zip(textmsa.names, textmsa.alignment):
                name_tmp = name.decode() 
                f.write(">" + name_tmp + "\n")
                f.write(aligned.replace("-", ".") + "\n")

    if file_format == "stockholm":
        fasta_file = os.path.splitext(msa_file)[0] + ".fasta"
        
        with open(fasta_file, "w") as f:
            textmsa = iteration_res.hits.to_msa(easel.Alphabet.rna())
            for name, aligned in zip(textmsa.names, textmsa.alignment):
                name_tmp = name.decode() 
                f.write(">" + name_tmp + "\n")
                f.write(aligned.replace("-", ".") + "\n")

        cmd = f"esl-reformat stockholm {fasta_file} > {msa_file}"
        res = subprocess.run(cmd, shell = True, capture_output=True)
        print(res.stdout.decode())
        print(res.stderr.decode())

    return msa_file

def align_to_hmm(hmm_file, fasta_file, output_file, all_consensus_cols = True):
    hmm = pyhmmer.plan7.HMMFile(hmm_file).read()
    with easel.SequenceFile(fasta_file, digital=True, alphabet=easel.Alphabet.rna()) as sequence_file:
        sequences = sequence_file.read_block()
        msa = hmmalign(hmm, sequences, all_consensus_cols = all_consensus_cols)
    
    with open(output_file, "w") as f:
        for name, aligned in zip(msa.names, msa.alignment):
            f.write(f">{name.decode()}\n")
            f.write(f"{aligned}\n")
    return msa

def load_hmm(hmm_file):
    return pyhmmer.plan7.HMMFile(hmm_file).read()