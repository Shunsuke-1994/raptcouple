# Purpose: extract sequence-motif candidates from a trained Potts model and compare motifs.
#
# A Potts model stores positional nucleotide frequencies (`fi`) over the columns of its
# alignment. Informative stretches of those frequencies are the candidate motifs. The
# helpers below implement the procedure used in the RaptCouple paper:
#
#   positional frequencies -> keep the columns the focus sequence occupies, renormalize
#   over A/C/G/U -> cut into blocks of low positional entropy (self-segmentation) ->
#   merge near-duplicate blocks -> drop low-complexity blocks -> rank by prevalence.
#
# Motifs are 4 x L arrays with rows ordered A, C, G, U and columns summing to 1.
import re
from collections import Counter

import numpy as np

ALPHABET = "ACGU"

# Defaults used for the benchmarks in the paper.
ENTROPY_THRESHOLD = 1.5   # bits; a column is informative below this
MIN_BLOCK_LEN = 5         # nucleotides; shorter informative blocks are discarded
MAX_INTERNAL_GAP = 1      # uninformative columns tolerated inside a block
DEDUP_SIMILARITY = 0.8    # cosine similarity at which two blocks are the same motif
MIN_OVERLAP = 4           # columns that must overlap when sliding two motifs
MAX_PADDING = 3           # columns a motif may hang over the end of the other
CORE_KMER = 6             # width of the most informative core used for prevalence
LOWCOMP_HOMOPOLYMER = 0.70
LOWCOMP_DINUCLEOTIDE = 0.80


def normalize(pwm):
    """Column-normalize a 4 x L array; all-zero columns are left unchanged."""
    pwm = np.asarray(pwm, dtype=float).copy()
    total = pwm.sum(0)
    total[total == 0] = 1
    return pwm / total


def entropy(pwm, eps=1e-20):
    """Positional Shannon entropy in bits, one value per column."""
    pwm = np.asarray(pwm, dtype=float)
    return -np.nansum(pwm * np.log2(pwm + eps), axis=0)


def cosine(pwm_a, pwm_b):
    """Column-wise cosine similarity between two equally wide motifs."""
    denominator = np.linalg.norm(pwm_a, axis=0) * np.linalg.norm(pwm_b, axis=0)
    denominator[denominator == 0] = np.nan
    return (pwm_a * pwm_b).sum(0) / denominator


def similarity(pwm_a, pwm_b, max_padding=MAX_PADDING, min_overlap=MIN_OVERLAP):
    """Best mean cosine similarity over all offsets of the shorter motif along the longer.

    The shorter motif may hang over either end by up to ``max_padding`` columns, in which
    case only the overlapping columns are compared; offsets overlapping fewer than
    ``min_overlap`` columns are skipped. Returns NaN if neither motif has any column.
    """
    if pwm_a is None or pwm_b is None or pwm_a.shape[1] == 0 or pwm_b.shape[1] == 0:
        return np.nan
    if pwm_a.shape[1] < pwm_b.shape[1]:
        pwm_a, pwm_b = pwm_b, pwm_a
    long_len, short_len = pwm_a.shape[1], pwm_b.shape[1]
    padding = min(max_padding, short_len // 2)
    scores = []
    for offset in range(-padding, long_len - short_len + padding + 1):
        if offset < 0:
            overlap = short_len + offset
            if overlap < min_overlap:
                continue
            scores.append(np.nanmean(cosine(pwm_a[:, :overlap], pwm_b[:, -offset:])))
        elif offset > long_len - short_len:
            overlap = short_len - (offset - (long_len - short_len))
            if overlap < min_overlap:
                continue
            scores.append(np.nanmean(cosine(pwm_a[:, offset:], pwm_b[:, :overlap])))
        else:
            scores.append(np.nanmean(cosine(pwm_a[:, offset:offset + short_len], pwm_b)))
    return max(scores) if scores else np.nan


def positional_frequencies(params):
    """Motif-shaped positional frequencies of a model, over the focus sequence's columns.

    ``params`` is the dictionary returned by ``src.plmc.read_params``. Alignment columns
    where the focus sequence has a gap are insertions relative to it and are dropped;
    the remaining columns are renormalized over A, C, G and U, so a gap state in the
    model does not contribute.
    """
    alphabet = list(params["alphabet"])
    frequencies = np.asarray(params["fi"], dtype=float).T          # states x columns
    rows = [alphabet.index("U" if base == "U" else base) if base in alphabet else
            alphabet.index("T") for base in ALPHABET]
    occupied = [i for i, base in enumerate(params["target_seq"]) if base in alphabet[:4]]
    return normalize(frequencies[rows][:, occupied])


def segment(pwm, entropy_threshold=ENTROPY_THRESHOLD, min_length=MIN_BLOCK_LEN,
            max_gap=MAX_INTERNAL_GAP):
    """Cut positional frequencies into informative blocks (self-segmentation).

    A column is informative when its entropy is below ``entropy_threshold``. Blocks may
    contain up to ``max_gap`` consecutive uninformative columns and must be at least
    ``min_length`` columns wide.
    """
    pwm = normalize(pwm)
    informative = entropy(pwm) < entropy_threshold
    blocks, i = [], 0
    while i < len(informative):
        if not informative[i]:
            i += 1
            continue
        start = end = i
        gap, k = 0, i
        while k + 1 < len(informative):
            k += 1
            if informative[k]:
                end, gap = k, 0
            else:
                gap += 1
                if gap > max_gap:
                    break
        if end - start + 1 >= min_length:
            blocks.append(normalize(pwm[:, start:end + 1]))
        i = end + 1
    return blocks


def deduplicate(blocks, threshold=DEDUP_SIMILARITY):
    """Merge blocks whose similarity reaches ``threshold``, keeping the widest of each group."""
    kept = []
    for block in blocks:
        for i, other in enumerate(kept):
            if similarity(block, other) >= threshold:
                if block.shape[1] > other.shape[1]:
                    kept[i] = block
                break
        else:
            kept.append(block)
    return kept


def consensus(pwm):
    """Most probable base at each column."""
    return "".join(ALPHABET[int(np.argmax(column))] for column in np.asarray(pwm).T)


def core_kmer(pwm, width=CORE_KMER):
    """The ``width``-nucleotide window of the consensus carrying the most information."""
    sequence = consensus(pwm)
    if len(sequence) <= width:
        return sequence
    information = 2 - entropy(pwm)
    start = max(range(len(sequence) - width + 1), key=lambda i: information[i:i + width].sum())
    return sequence[start:start + width]


def is_low_complexity(pwm, homopolymer=LOWCOMP_HOMOPOLYMER, dinucleotide=LOWCOMP_DINUCLEOTIDE):
    """True for consensus sequences dominated by one base or by a two-base repeat."""
    sequence = consensus(pwm)
    length = len(sequence)
    most_common = Counter(sequence).most_common(1)[0][1] / length
    repeated = sum(1 for i in range(length - 2) if sequence[i] == sequence[i + 2]) / max(1, length - 2)
    return most_common >= homopolymer or repeated >= dinucleotide


def read_count(header):
    """Read count encoded in a FASTAptamer header (``...-rank-count-rpm``); 1 if absent."""
    match = re.search(r"-(\d+)-(\d+)-[\d.]+(?:-\d+)?$", header)
    return int(match.group(2)) if match else 1


def load_reads(fasta_file):
    """[(sequence, read count)] from an annotated FASTA, with T read as U."""
    reads, count = [], 1
    with open(fasta_file) as handle:
        for line in handle:
            if line.startswith(">"):
                count = read_count(line[1:].strip())
            else:
                reads.append((line.strip().upper().replace("T", "U"), count))
    return reads


def prevalence(pwm, reads):
    """Total reads containing the motif's core k-mer."""
    core = core_kmer(pwm)
    return sum(count for sequence, count in reads if core in sequence)


def candidates(params, reads=None, remove_low_complexity=True, top_n=3, **segment_kwargs):
    """Motif candidates of one trained model, ranked by prevalence.

    ``params`` comes from ``src.plmc.read_params``. Supply ``reads`` (see ``load_reads``)
    from the SELEX pool the model was trained on to rank candidates by prevalence and to
    keep the ``top_n`` most prevalent; without reads the candidates are returned in the
    order they occur along the sequence. Returns a list of dictionaries with the motif
    (``pwm``), its ``consensus``, ``core_kmer``, ``prevalence`` and ``low_complexity``.
    """
    blocks = deduplicate(segment(positional_frequencies(params), **segment_kwargs))
    found = [{"pwm": block,
              "consensus": consensus(block),
              "core_kmer": core_kmer(block),
              "prevalence": prevalence(block, reads) if reads else None,
              "low_complexity": is_low_complexity(block)} for block in blocks]
    if remove_low_complexity:
        found = [m for m in found if not m["low_complexity"]]
    if reads and top_n is not None:
        found = sorted(found, key=lambda m: -m["prevalence"])[:top_n]
    return found


if __name__ == "__main__":
    # A motif and a shifted copy of it are the same motif.
    rng = np.random.default_rng(0)
    motif = normalize(rng.random((4, 12)) ** 4)
    assert np.isclose(similarity(motif, motif), 1.0)
    assert similarity(motif, motif[:, 2:]) > DEDUP_SIMILARITY
    assert len(deduplicate([motif, motif[:, 1:-1]])) == 1

    # Self-segmentation keeps an informative block inside uninformative flanks.
    flat = np.full((4, 6), 0.25)
    sharp = np.zeros((4, 7))
    sharp[[0, 1, 2, 3, 0, 1, 2], range(7)] = 1.0
    blocks = segment(np.hstack([flat, sharp, flat]))
    assert len(blocks) == 1 and blocks[0].shape[1] == 7, [b.shape for b in blocks]
    assert consensus(blocks[0]) == "ACGUACG", consensus(blocks[0])

    # Low-complexity consensus sequences are recognised.
    poly = np.zeros((4, 10)); poly[0, :] = 1.0
    assert is_low_complexity(poly)
    repeat = np.zeros((4, 10)); repeat[[0, 1] * 5, range(10)] = 1.0
    assert is_low_complexity(repeat)
    assert not is_low_complexity(sharp)

    # The core k-mer is the most informative window.
    mixed = np.hstack([np.full((4, 3), 0.25), sharp])
    assert core_kmer(mixed) == "ACGUAC", core_kmer(mixed)

    print("motif checks passed")
