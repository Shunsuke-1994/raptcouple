# Purpose: RNA secondary structure prediction using Nussinov algorithm
# TODO: replace with viennaRNA package to use fully customizable scoring table
import numpy as np
from src.plmc import detect_coupling


def nussinov(rna, min_loop_length=3, score_table=None, sanity_check=True):
    def valid_pair(n1, n2):
        if sanity_check:
            return {n1, n2} in [{"A", "U"}, {"C", "G"}, {"G", "U"}]
        else:
            return True

    n = len(rna)
    dp = [[0 for _ in range(n)] for _ in range(n)]
    if score_table is None:
        score_table = np.ones((n, n))

    # dynamic programming
    for l in range(1, n):
        for i in range(n - l):
            j = i + l 

            dp[i][j] = max(dp[i + 1][j], 
                          dp[i][j - 1], 
                          (dp[i + 1][j - 1] + score_table[i][j]) if valid_pair(rna[i], rna[j]) and l > min_loop_length else dp[i + 1][j - 1],
                          max(dp[i][k] + dp[k + 1][j] for k in range(i, j))
                          )
    
    def traceback(i, j, brackets):
        if i < j:
            if dp[i][j] == dp[i + 1][j]:
                traceback(i + 1, j, brackets)
            elif dp[i][j] == dp[i][j - 1]:
                traceback(i, j - 1, brackets)
            elif dp[i][j] == dp[i + 1][j - 1] + score_table[i][j]: # + (1 if valid_pair(rna[i], rna[j]) else 0):
                # brackets.append((i, j))
                brackets[i] = "("
                brackets[j] = ")"
                traceback(i + 1, j - 1, brackets)
            else: # bifurcation
                for k in range(i, j):
                    if dp[i][j] == dp[i][k] + dp[k + 1][j]:
                        traceback(i, k, brackets)
                        traceback(k + 1, j, brackets)
                        break

    brackets = ["." for _ in range(n)]
    traceback(0, n - 1, brackets)
    return dp[0][n - 1], "".join(brackets)

def fold(
        params, 
        min_loop_length=3, 
        threshold = 3, 
        only_match_col = True, 
        sanity_check = True
        ):
    
    if only_match_col:
        pairs = detect_coupling(params, threshold=threshold, min_dist=min_loop_length+1, sanity_check=sanity_check, use_mask_min_dist=True)
        matched_cols = [i for i,n in enumerate(params["target_seq"]) if n in "AUGC"]
        pos2index = {pos:i for i,pos in enumerate(matched_cols)}
        rna = "".join([n for i,n in enumerate(params["target_seq"]) if i in matched_cols])
        score_table = np.zeros([len(matched_cols), len(matched_cols)])
    else:
        pairs = detect_coupling(params, threshold=threshold, min_dist=min_loop_length+1, sanity_check=sanity_check, use_mask_min_dist=False)
        rna = params["target_seq"]
        score_table = np.zeros(params["FN_apc"].shape)

    for pair in pairs:
        nuc, pos = pair
        if only_match_col:
            # pairs can contain non-matching columns, so we need to exclude them
            if (pos[0] in pos2index.keys()) and (pos[1] in pos2index.keys()):
                # if abs(pos2index[pos[1]] - pos2index[pos[0]]) > min_loop_length:
                # remove if. see 08-04-check
                score_table[pos2index[pos[0]], pos2index[pos[1]]] = params["FN_apc"][pos[0], pos[1]]
        else:
            score_table[pos[0], pos[1]] = params["FN_apc"][pos[0], pos[1]]

    energy, ss = nussinov(rna, min_loop_length=min_loop_length, score_table=score_table, sanity_check=sanity_check)
    return rna, energy, ss, score_table



if __name__ == "__main__":
    rna = "GCAAAGCC"
    score, structure = nussinov(rna)
    print(rna)
    print(structure)

    rna = "GGCCAAGGCC"
    score, structure = nussinov(rna)
    print(rna)
    print(structure)

