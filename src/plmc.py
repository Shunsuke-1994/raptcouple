# wrapper of plmc
import pandas as pd
import subprocess
import numpy as np
from Bio import SeqIO, AlignIO
from collections import Counter
import matplotlib.pyplot as plt 
import seaborn as sns

# you need to specify the path to plmc
PATH_TO_PLMC = "~/desktop/RNA/plmc/bin/plmc" # "path/to/plmc"


def fit(fasta_file, target, param_file, coupling_file, vocab, threshold = 0.01, print_result = True):
    alignments = list(SeqIO.parse(fasta_file, "fasta"))
    le = round(0.2*(len(alignments[0])-1), 1)
    cmd = f"""{PATH_TO_PLMC} -o {param_file} -c {coupling_file} -a {vocab} -f {target} -le {str(le)} -lh 0.01 -m 200 -t {threshold} {fasta_file}"""
    res = subprocess.run(cmd, shell = True, capture_output = True)
    if print_result:
        print(res.stdout.decode("utf-8"))
        print(res.stderr.decode("utf-8"))
    return cmd

def read_coupling_score(coupling_file):
    """
    The average Frobenius norm of DCA with average product correction (APC).
    """
    with open(coupling_file, "r") as c:
        df_c = pd.read_csv(c, sep = " ", header=None)[[0, 2, 5]]
    L = df_c[2].max()
    np_c = np.zeros([L, L])
    np_c[:,:] = np.nan
    for x,y,v in np.array(df_c):
        np_c[int(x)-1, int(y)-1] = v
    return np_c

def read_params(paramfile):
    """
    This function reads a binary file of parameters for a pairwise maximum entropy
    model (undirected graphical model). The model describes
    the distribution of sequences of length L drawn from an alphabet with q
    characters.
    ref: https://github.com/debbiemarkslab/plmc/blob/18c9e55e3bd2f14f4968be19a807b401996c929a/scripts/read_params.m

    Parameters
    ----------
    paramfile : str
        The path to the binary file containing the parameters.

    Returns
    -------
    params : dict
        A dictionary containing the parameters of the model.

    """

    # Define the precision
    PRECISION = np.float32
    
    # Open the binary file
    with open(paramfile, 'rb') as f:
        params = {}
        L = np.fromfile(f, dtype=np.int32, count=1)[0]
        q = np.fromfile(f, dtype=np.int32, count=1)[0]
        params["numSeqs"] = np.fromfile(f, dtype=np.int32, count=1)[0]
        params["numInvalidSeqs"] = np.fromfile(f, dtype=np.int32, count=1)[0]
        params["numIter"] = np.fromfile(f, dtype=np.int32, count=1)[0]
        params["theta"] = np.fromfile(f, dtype=PRECISION, count=1)[0]
        params["lambda_h"] = np.fromfile(f, dtype=PRECISION, count=1)[0]
        params["lambda_J"] = np.fromfile(f, dtype=PRECISION, count=1)[0]
        params["lambda_group"] = np.fromfile(f, dtype=PRECISION, count=1)[0]
        params["nEff"] = np.fromfile(f, dtype=PRECISION, count=1)[0]
        params["alphabet"] = "".join(map(chr, np.fromfile(f, dtype=np.int8, count=q)))
        num_weights = params["numSeqs"] + params["numInvalidSeqs"]
        params["weights"] = np.fromfile(f, dtype=PRECISION, count=num_weights)
        params["target_seq"] = "".join(map(chr, np.fromfile(f, dtype=np.int8, count=L)))
        params["offset_map"] = np.fromfile(f, dtype=np.int32, count=L)
        params["fi"] = np.fromfile(f, dtype=PRECISION, count=q*L).reshape((L, q))
        params["hi"] = np.fromfile(f, dtype=PRECISION, count=q*L).reshape((L, q))
        
        # Reading marginals and couplings
        params["fij"] = np.zeros((L, L, q, q))
        params["Jij"] = np.zeros((L, L, q, q))
        block_size = L * (L-1) * q * q // 2
        for matrix in ("fij", "Jij"):
            block = np.fromfile(f, dtype=PRECISION, count=block_size)
            offset = 0
            for i in range(L-1):
                for j in range(i+1, L):
                    slice_ = block[offset:offset+q*q].reshape((q, q))
                    params[matrix][j, i, :, :] = slice_
                    params[matrix][i, j, :, :] = slice_.T
                    offset += q*q

    # calc Frobenius norm with APC
    # same as read_params
    def calculate_FN(params, slice=None):
        """
        Calculate the Frobenius norm of the Jij parameter.

        Parameters
        ----------
        params : dict
            The parameters dictionary.
        slice : slice, optional
            The slice of the Jij parameter to use.

        Returns
        -------
        FN : ndarray
            The Frobenius norms.
        APC : ndarray
            The Average Product Correction.
        CN : ndarray
            The corrected norms.

        """
        if slice is None:
            slice = np.s_[:]

        N = params["Jij"].shape[0]
        FN = np.zeros((N, N))

        for i in range(N-1):
            for j in range(i+1, N):
                FN[i, j] = np.linalg.norm(params["Jij"][i, j, slice, slice])
                FN[j, i] = FN[i, j]

        # Average Product Correction
        # FN: Frobenius norm
        # APC: Average Product Correction
        # CN: Corrected Norm
        FN_means = np.mean(FN, axis=0)*N/(N-1)
        FN_means_all = np.mean(FN)*N/(N-1)
        APC = np.outer(FN_means, FN_means) / FN_means_all
        CN = FN - APC
        CN = CN - np.diag(np.diag(CN))

        return FN, APC, CN

    params["FN"], _, params["FN_apc"] = calculate_FN(params)

    return params

def calc_mi_from_msa(fasta_file):
    align = AlignIO.read(fasta_file, "fasta")
    seqlen = len(align[0])
    i_j_mi_ij = []
    mitable = np.zeros([align.get_alignment_length(), align.get_alignment_length()])

    for i in range(seqlen):
        for j in range(seqlen):

            aln_size = len(align)
            if i != j:
                col_i = align[:,i]
                col_j = align[:,j]
                pairs_count = Counter([i+j for i,j in zip(col_i, col_j)])
                total_count = sum(pairs_count.values())
                mi_ij = [(c/total_count)*(np.log2(c*total_count/(col_i.count(p[0])*col_j.count(p[1])))) for p, c in pairs_count.items()]
                mitable[i,j] = sum(mi_ij)
            else:
                # entrpy
                col_i = align[:, i]
                nuc_count = Counter(col_i)
                mi_diag = [(c/aln_size)*(np.log2(c) - np.log2(aln_size) + 2) for nuc, c in nuc_count.items() if nuc != "."]
                # mi_diag = [(c/aln_size)*(np.log2(c) - np.log2(aln_size) + 2) for nuc, c in nuc_count.items()]
                mitable[i,j] = sum(mi_diag)

    return mitable
 
def cleanup_pairs(detected_pairs, min_dist, sanity_check = True):
    """
    check if detected pairs are valid as normal RNA.
    min_dist: minimum distance between two positions.
    """
    clean_pairs = []
    for pair in detected_pairs:
        # check if pair is distance over than min_dist
        if abs(pair[1][0]-pair[1][1]) >= min_dist:

            # check if pair is (A,U) |(U,A) | (G,C) | (C,G) | (G,U) | (U,G)
            if sanity_check:
                if set(pair[0]) in [{"A","U"}, {"G","C"}, {"G","U"}]:
                    clean_pairs.append(pair)
            else:
                clean_pairs.append(pair)
            # print(pair)
    return clean_pairs

def detect_coupling(
        params, 
        threshold=3, 
        min_dist=3, 
        sanity_check = True, 
        use_mask_min_dist = True
        ):
    """
    detect coupling over than threshold
    """
    # looplen = 0 (min_dist = 1) --> k = 1
    # looplen = 1 (min_dist = 2) --> k = 2
    FN_apc = np.triu(params["FN_apc"], k = min_dist if use_mask_min_dist else 0)
    FN_apc[FN_apc==0] = np.nan
    FN_apc[FN_apc<np.nanmean(FN_apc)+threshold*np.nanstd(FN_apc)] = np.nan
    FN_apc = np.argwhere(~np.isnan(FN_apc))
    detected_pairs = []
    for pos in FN_apc:
        detected_pairs.append(
            (
                (params["target_seq"][pos[0]], params["target_seq"][pos[1]]),
                (pos[0], pos[1])
            )
        )

    return cleanup_pairs(detected_pairs, min_dist=min_dist, sanity_check=sanity_check)

# def fold_soft_constraint(file_model_params):
#     params = read_params(file_model_params)
#     target_seq = params["target_seq"]
#     sc = np.zeros([len(params["target_seq"])+1, len(params["target_seq"])+1])
#     sc[1:, 1:] = params["FN_apc"]

#     fc = RNA.fold_compound(params["target_seq"])
#     fc.sc_set_bp(sc)

#     return target_seq, fc.mfe()

def fold_evalue(file_model_params, evalue = 0.01):
    params = read_params(file_model_params)
    target_seq = params["target_seq"]

    return target_seq

def visualize_coupling(
        file_model_params, 
        cmap = "cividis", 
        figsize = (6,4)
        ):
    params = read_params(file_model_params)
    FN_apc = params["FN_apc"]
    target_seq = params["target_seq"]

    fig, ax = plt.subplots(figsize = figsize)
    sns.heatmap(FN_apc, cmap = cmap, ax = ax, square = True, cbar = True, cbar_kws = {"label":"Frobenius Norm"})
    ax.set_xticks([i+0.5 for i in range(len(target_seq))], target_seq);
    ax.set_yticks([i+0.5 for i in range(len(target_seq))], target_seq);

    return(fig, ax)
