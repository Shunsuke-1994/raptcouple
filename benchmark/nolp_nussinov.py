"""
no-lonely-pair (noLP) Nussinov-style folding for RaptCouple coupling scores.

Enforces helices of length >= 2 INSIDE the DP (no isolated base pairs), instead of the
post-hoc `remove_isolated_pairs`. Max-score base pairing with a coupling score_table.

Two matrices:
  W[i][j] : best score on [i..j], no lonely pairs anywhere.
  V[i][j] : best score on [i..j] with (i,j) paired as the OUTER pair of a helix of >=2
            stacked pairs ("stem"). The innermost pair of a stem is supplied by the
            length-2 base case (so it never needs its own V).

V[i][j] = score(i,j) + max(
            V[i+1][j-1],                              # extend a >=2 stem inward (-> total >=3)
            score(i+1,j-1) + W[i+2][j-2]              # base case: exactly a 2-stack, then loop region
          )    (only if valid_pair & loop-length constraints hold for the involved pairs)
W[i][j] = max(W[i+1][j], W[i][j-1], V[i][j], max_k W[i][k]+W[k+1][j])
"""
import numpy as np

NEG = -1e18

def _valid(a, b, sanity_check):
    if not sanity_check:
        return True
    return {a, b} in [{"A", "U"}, {"C", "G"}, {"G", "U"}]

def nussinov_nolp(rna, min_loop_length=2, score_table=None, sanity_check=True):
    n = len(rna)
    if n < 2:
        return 0.0, "." * n
    if score_table is None:
        score_table = np.ones((n, n))
    W = [[0.0] * n for _ in range(n)]
    V = [[NEG] * n for _ in range(n)]

    def can_pair(i, j):
        # only place a base pair where a coupling was actually detected (score>0):
        # the structure consists solely of detected couplings; noLP then requires >=2 stacked.
        return _valid(rna[i], rna[j], sanity_check) and (j - i) > min_loop_length and score_table[i][j] > 0

    for l in range(1, n):
        for i in range(n - l):
            j = i + l
            best = NEG
            if can_pair(i, j):
                s_ij = score_table[i][j]
                # extend an inner stem (>=3 total)
                if i + 1 < j - 1 and V[i + 1][j - 1] > NEG / 2:
                    best = max(best, s_ij + V[i + 1][j - 1])
                # base case: 2-stack (i,j)&(i+1,j-1), then W inside
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
        # (i,j) is the outer pair of a >=2 stem
        brackets[i] = "("; brackets[j] = ")"
        s_ij = score_table[i][j]
        # branch 1: extend an inner >=2 stem
        if i + 1 < j - 1 and V[i + 1][j - 1] > NEG / 2 and abs(V[i][j] - (s_ij + V[i + 1][j - 1])) < 1e-9:
            tb_V(i + 1, j - 1); return
        # branch 2: base-case 2-stack (verify explicitly, per codex review)
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


if __name__ == "__main__":
    rna = "GGGGAAAACCCC"
    scores = np.zeros((len(rna), len(rna)))
    for i, j in [(0, 11), (1, 10), (2, 9), (3, 8)]:
        scores[i, j] = 1
    assert nussinov_nolp(rna, 3, scores) == (4.0, "((((....))))")
    scores = np.zeros((8, 8)); scores[0, 7] = 5
    assert nussinov_nolp("GAAAAAAC", 3, scores) == (0.0, "........")
    scores = np.zeros((6, 6)); scores[0, 5] = scores[1, 4] = 1
    assert nussinov_nolp("GGAACC", 2, scores) == (2.0, "((..))")
    assert nussinov_nolp("GGAACC", 3, scores) == (0.0, "......")
    print("noLP stem, isolated-pair, and loop-length checks passed")
