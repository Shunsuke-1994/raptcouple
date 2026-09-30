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
            elif valid_pair(rna[i], rna[j]) and (j - i) > min_loop_length and dp[i][j] == dp[i + 1][j - 1] + score_table[i][j]:
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

NEG = -1e18

def nussinov_nolp(rna, min_loop_length=2, score_table=None, sanity_check=True):
    """noLP Nussinov: enforces helices of length >= 2 inside the DP (no isolated base pairs).

    Two matrices:
      W[i][j] : best score on [i..j], no lonely pairs anywhere.
      V[i][j] : best score on [i..j] with (i,j) paired as the outer pair of a helix of >= 2
                stacked pairs.
    """
    def _valid(a, b):
        if not sanity_check:
            return True
        return {a, b} in [{"A", "U"}, {"C", "G"}, {"G", "U"}]

    n = len(rna)
    if n < 2:
        return 0.0, "." * n
    if score_table is None:
        score_table = np.ones((n, n))
    W = [[0.0] * n for _ in range(n)]
    V = [[NEG] * n for _ in range(n)]

    def can_pair(i, j):
        return _valid(rna[i], rna[j]) and (j - i) > min_loop_length and score_table[i][j] > 0

    for l in range(1, n):
        for i in range(n - l):
            j = i + l
            best = NEG
            if can_pair(i, j):
                s_ij = score_table[i][j]
                if i + 1 < j - 1 and V[i + 1][j - 1] > NEG / 2:
                    best = max(best, s_ij + V[i + 1][j - 1])
                if can_pair(i + 1, j - 1):
                    inner = W[i + 2][j - 2] if (i + 2 <= j - 2) else 0.0
                    best = max(best, s_ij + score_table[i + 1][j - 1] + inner)
            V[i][j] = best
            w = max(W[i + 1][j] if i + 1 <= j else 0.0,
                    W[i][j - 1] if i <= j - 1 else 0.0)
            if V[i][j] > NEG / 2:
                w = max(w, V[i][j])
            for k in range(i, j):
                w = max(w, W[i][k] + W[k + 1][j])
            W[i][j] = w

    brackets = ["."] * n

    def tb_V(i, j):
        brackets[i] = "("; brackets[j] = ")"
        s_ij = score_table[i][j]
        if i + 1 < j - 1 and V[i + 1][j - 1] > NEG / 2 and abs(V[i][j] - (s_ij + V[i + 1][j - 1])) < 1e-9:
            tb_V(i + 1, j - 1); return
        inner = W[i + 2][j - 2] if (i + 2 <= j - 2) else 0.0
        if can_pair(i + 1, j - 1) and abs(V[i][j] - (s_ij + score_table[i + 1][j - 1] + inner)) < 1e-9:
            brackets[i + 1] = "("; brackets[j - 1] = ")"
            if i + 2 <= j - 2:
                tb_W(i + 2, j - 2)
            return
        raise RuntimeError(f"tb_V: no matching branch at ({i},{j})")

    def tb_W(i, j):
        if i >= j:
            return
        if abs(W[i][j] - W[i + 1][j]) < 1e-9:
            tb_W(i + 1, j); return
        if abs(W[i][j] - W[i][j - 1]) < 1e-9:
            tb_W(i, j - 1); return
        if V[i][j] > NEG / 2 and abs(W[i][j] - V[i][j]) < 1e-9:
            tb_V(i, j); return
        for k in range(i, j):
            if abs(W[i][j] - (W[i][k] + W[k + 1][j])) < 1e-9:
                tb_W(i, k); tb_W(k + 1, j); return

    tb_W(0, n - 1)
    return W[0][n - 1], "".join(brackets)


def fold(
        params,
        min_loop_length=3,
        threshold = 3,
        only_match_col = True,
        sanity_check = True,
        nolp = True
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

    if nolp:
        energy, ss = nussinov_nolp(rna, min_loop_length=min_loop_length, score_table=score_table, sanity_check=sanity_check)
    else:
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

    # noLP Nussinov checks
    # 4-bp helix folds
    rna = "GGGGAAAACCCC"
    st = np.zeros((len(rna), len(rna)))
    for i, j in [(0, 11), (1, 10), (2, 9), (3, 8)]:
        st[i][j] = 1
    _, db = nussinov_nolp(rna, 2, st, sanity_check=True)
    assert db == "((((....))))", db
    # a single isolated pair is not retained
    rna2 = "GAAAAAAC"
    st2 = np.zeros((len(rna2), len(rna2))); st2[0][7] = 5
    _, db2 = nussinov_nolp(rna2, 2, st2, sanity_check=True)
    assert db2 == "........", db2
    # a 2-stack is kept at min_loop 2 and rejected at min_loop 3 (inner loop of two bases)
    rna3 = "GGUUCC"
    st3 = np.zeros((len(rna3), len(rna3))); st3[0][5] = 1; st3[1][4] = 1
    _, db3 = nussinov_nolp(rna3, 2, st3, sanity_check=True)
    assert db3 == "((..))", db3
    _, db4 = nussinov_nolp(rna3, 3, st3, sanity_check=True)
    assert db4 == "......", db4
    print("noLP checks passed")
