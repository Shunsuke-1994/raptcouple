import pandas as pd 
import numpy as np
import RNA 

def load_metadata(path_to_supptables):
    df_info = pd.read_excel(path_to_supptables, sheet_name=1, skiprows=15)
    df_linear = pd.DataFrame(
        df_info.values[
            [i for i in range(df_info.shape[0]) if i%5 == 0]
            ],
    ).dropna(how="all", axis = 1)
    df_linear.columns = "Gene Symbol	Barcode	Batch	seed	multinomial	cycle	Mono or multimer	Motif type	Experiment	Number of RBDs in the construct".split("\t") + ["Note"]
    df_linear = df_linear.fillna({"Number of RBDs in the construct":0}).astype({"Number of RBDs in the construct":int})
    df_linear["motif_type"] = "linear"

    df_info = pd.read_excel(path_to_supptables, sheet_name=2, skiprows=13)
    df_struct = pd.DataFrame(
        df_info.values[
            [i for i in range(df_info.shape[0]) if i%5 == 0]
            ],
    ).dropna(how="all", axis = 1)
    df_struct.columns = "Gene Symbol	Barcode	Batch	seed	multinomial	cycle	Motif type	Experiment	Number of RBDs in the construct".split("\t") #+ ["None"]
    df_struct = df_struct.fillna({"Number of RBDs in the construct":0}).astype({"Number of RBDs in the construct":int})
    df_struct["motif_type"] = "structured"

    df_metadata = pd.concat([df_linear, df_struct]).reset_index(drop = True)
    df_metadata["seed_regex"] = [ambiguous_rna_to_regex(seed) for seed in df_metadata["seed"]]

    return df_metadata

def load_pwms(path_to_supptables):
    """
    load PWM of each motif from supplementary tables of Jolma et al. 2020. 
    https://genome.cshlp.org/content/suppl/2020/07/23/gr.258848.119.DC1/Supplemental_Tables_S1-S8.xlsx
    Each motif is named as concat of metadata: "gene_barcode_batch_seed_multinominal_cycle_monomer_or_multimer_motif_type_experiment"
    """

    dict_gene2motif_jolma = {}

    df_motif_jolma_linear = pd.read_excel(path_to_supptables, sheet_name="S2-PWMs-linear", skiprows=15)
    df_motif_jolma_structured = pd.read_excel(path_to_supptables, sheet_name="S3-PWMs-structured", skiprows=13)

    # add linear motifs
    for i,row in df_motif_jolma_linear.iterrows():
        if str(row["Base position"]) not in "ACGU":
            gene = str(row["Base position"])
            barcode = row[1]
            batch = row[2]
            seed = row[3]
            multinominal = int(row[4])
            cycle = row[5]
            monomer_or_multimer = row[6]
            motif_type = row[7]
            experiment = row[8]
            tf_key = f"{gene}_{barcode}_{batch}_{seed}_{multinominal}_{cycle}_{monomer_or_multimer}_{motif_type}_{experiment}"
            motif_matrix = np.zeros((4, len(row)-1))
            motif_matrix[:,:] = np.nan
            motif_matrix = pd.DataFrame(motif_matrix, index = ["A", "C", "G", "U"])
        else:
            motif_matrix.loc[row["Base position"]] = row[1:].values.astype(float)
            if row["Base position"] == "U":
                motif_matrix.dropna(axis=1, how="all", inplace=True)
                for col in motif_matrix.columns:
                    motif_matrix[col] = motif_matrix[col] / motif_matrix[col].sum()
                dict_gene2motif_jolma[tf_key] = motif_matrix

    # add structured motifs
    for i,row in df_motif_jolma_structured.iterrows():
        if str(row["Base position"]) not in "ACGU":
            gene = str(row["Base position"])
            barcode = row[1]
            batch = row[2]
            seed = row[3]
            multinominal = int(row[4])
            cycle = row[5]
            motif_type = row[6]
            experiment = row[7]
            tf_key = f"{gene}_{barcode}_{batch}_{seed}_{multinominal}_{cycle}_{motif_type}_{experiment}"
            motif_matrix = np.zeros((4, len(row)-1))
            motif_matrix[:,:] = np.nan
            motif_matrix = pd.DataFrame(motif_matrix, index = ["A", "C", "G", "U"])
        else:
            motif_matrix.loc[row["Base position"]] = row[1:].values.astype(float)
            if row["Base position"] == "U":
                motif_matrix.dropna(axis=1, how="all", inplace=True)
                for col in motif_matrix.columns:
                    motif_matrix[col] = motif_matrix[col] / motif_matrix[col].sum()

                dict_gene2motif_jolma[tf_key] = motif_matrix


    return dict_gene2motif_jolma

def ambiguous_rna_to_regex(seq):
    conversion_dict = {
        'R': '[AG]',
        'Y': '[UC]',
        'K': '[GU]',
        'M': '[AC]',
        'S': '[GC]',
        'W': '[AU]',
        'B': '[GUC]',
        'D': '[GAU]',
        'H': '[AUC]',
        'V': '[GAC]',
        'N': '[AUGC]'
    }

    regex_seq = ''.join([conversion_dict[base] if base in conversion_dict else base for base in seq])
    return regex_seq

def calc_entropy_on_PWM(df_pwm, eps= 1e-20):
    """
    df_pwm: (pd.DataFrame) with shape (4, n)
    return: entroyp (float)  with shape (n, 1)
    """
    return -np.nansum(df_pwm * np.log2(df_pwm + eps), axis=0)

def trim_pwm_with_entropy(df_pwm, threshold, eps=1e-20):
    """
    df_pwm: (pd.DataFrame) with shape (4, n)
    threshold: (float) 
    return: df_pwm_trimmed (pd.DataFrame) with shape (4, n_trimmed)
    """
    entropy = calc_entropy_on_PWM(df_pwm, eps = eps)
    idx_untrimmed = np.where(entropy < threshold)[0]
    return df_pwm[range(idx_untrimmed[0],idx_untrimmed[-1]+1)]

def remove_isolated_pairs(structure):
    structure_list = list(structure)
    length = len(structure_list)
    
    pt = RNA.ptable(structure)
    dict_pt = {i-1: j-1 for i, j in enumerate(pt) if i>0}

    # if i is not end
    for i in range(length):
        if 1<=i<length-1:
            if "".join(structure_list[i-1:i+2]) == ".(." or "".join(structure_list[i-1:i+2]) == ".).":
                structure_list[i] = "."
                structure_list[dict_pt[i]] = "."

        elif i == 0:
            if "".join(structure_list[0:2]) == "(.":
                structure_list[0] = "."
                structure_list[dict_pt[i]] = "."
        elif i == length-1:
            if "".join(structure_list[-2:]) == ".)":
                structure_list[-1] = "."
                structure_list[dict_pt[i]] = "."

    return ''.join(structure_list)
