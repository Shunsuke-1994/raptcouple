import numpy as np

vocab_rna = {'A': 0, 'U': 1, 'G': 2, 'C': 3}
vocab_rna_rev = {0: 'A', 1: 'U', 2: 'G', 3: 'C'}
vocab_dna = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
vocab_dna_rev = {0: 'A', 1: 'T', 2: 'G', 3: 'C'}

vocab_rna_gapped = {'A': 0, 'U': 1, 'G': 2, 'C': 3, '.': 4}
vocab_rna_gapped_rev = {0: 'A', 1: 'U', 2: 'G', 3: 'C', 4: '.'}
vocab_dna_gapped = {'A': 0, 'T': 1, 'G': 2, 'C': 3, '.': 4}
vocab_dna_gapped_rev = {0: 'A', 1: 'T', 2: 'G', 3: 'C', 4: '.'}


def int2nuc(i, is_dna=True, is_gapped=True):
    if is_dna:
        if is_gapped:
            return vocab_dna_gapped_rev[i]
        else:
            return vocab_dna_rev[i]
    else:
        if is_gapped:
            return vocab_rna_gapped_rev[i]
        else:
            return vocab_rna_rev[i]
    
def nuc2int(nuc, is_dna=True, is_gapped=True):
    if is_dna:
        if is_gapped:
            return vocab_dna_gapped[nuc]
        else:
            return vocab_dna[nuc]
    else:
        if is_gapped:
            return vocab_rna_gapped[nuc]
        else:
            return vocab_rna[nuc]

def seq2onehot(seq, max_len, is_dna=True, is_gapped=True):
    """
    if not is_gapped: AUCG -> [[1, 0, 0, 0], [0, 1, 0, 0], ...]
    if is_gapped: AUCG. -> [[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], ...]
    """
    onehot = np.zeros((max_len, 4 if not is_gapped else 5))
    for i in range(len(seq)):
        onehot[i, nuc2int(seq[i], is_dna=is_dna, is_gapped=is_gapped)] = 1
    return onehot

def onehot2seq(onehot, is_dna=True, is_gapped=True):
    """
    if not is_gapped: [[1, 0, 0, 0], [0, 1, 0, 0], ...] -> AUCG
    if is_gapped: [[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], ...] -> AUCG.
    """
    seq = ""
    for i in range(onehot.shape[0]):
        seq += int2nuc(np.argmax(onehot[i]), is_dna=is_dna, is_gapped=is_gapped)
    return seq

def int2onehot(i, num_state):
    """
    if not is_gapped: 0 -> [1, 0, 0, 0]
    if is_gapped: 0 -> [1, 0, 0, 0, 0]
    """
    onehot = np.zeros(num_state)
    onehot[i] = 1
    return onehot

def onehot2int(onehot):
    """
    if not is_gapped: [1, 0, 0, 0] -> 0
    if is_gapped: [1, 0, 0, 0, 0] -> 0
    """
    return np.argmax(onehot)


if __name__ == "__main__":
    print('Testing util.py')
    print(np.eye(4) == seq2onehot('AUGC', 4, is_dna=False, is_gapped=False))
    print(np.eye(5) == seq2onehot('AUGC.', 5, is_dna=False, is_gapped=True))
    print("ATGC" == onehot2seq(seq2onehot('AUGC', 4, is_dna=False)))
    print(int2onehot(0, is_dna=False, is_gapped=False))
    print("A" == int2nuc(nuc2int('A')))
    print("G" == int2nuc(nuc2int('G')))
    print('Testing util.py done.')
