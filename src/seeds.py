# Purpose: choose query sequences (seeds) from a SELEX pool without prior motif knowledge.
#
# The target-guided mode of RaptCouple starts from a query the user supplies. In query-free
# mode the query is generated from the pool itself, by one of two independent routes:
#
#   STREME arm   enriched motifs are discovered with STREME, and the pool read carrying
#                each motif is used as the seed.
#   vsearch arm  reads are clustered with vsearch, and the most abundant member of each
#                of the most abundant clusters is used as the seed.
#
# Both routes return the same `Seed` records, so the downstream pipeline (jackhmmer ->
# Potts model -> motif candidates -> folding) is identical.
#
# Reads are expected to carry FASTAptamer annotations, as written by
# scripts/merge_and_cutadapt_all_rounds.py: `<round>-<rank>-<count>-<rpm>`. The round is
# written differently by different studies (`6R`, `Cycle4`, `round7`, or a bare number),
# and a merged read carries one such record per round it was seen in.
import re
import subprocess
from collections import defaultdict, namedtuple

# id:       header of the chosen read
# sequence: its sequence, as RNA
# label:    what selected it (a STREME motif consensus, or the cluster id)
# count:    read count of the chosen read
# support:  reads backing the choice (cluster abundance, or the seed's own count)
Seed = namedtuple("Seed", "id sequence label count support")

RECORD = re.compile(r"(?:^|[_-])(?:Cycle|round)?(\d+)R?-(\d+)-(\d+)-[0-9.]+(?:-\d+)?")


def as_rna(sequence):
    return sequence.strip().upper().replace("T", "U")


def records(header):
    """All FASTAptamer records in a header, as (round, rank, count).

    A merged read carries one record per round it was seen in, so a header can hold
    several. Accepts the plain, `Cycle` and `round` spellings.
    """
    return [(int(m.group(1)), int(m.group(2)), int(m.group(3))) for m in RECORD.finditer(header)]


def best_record(header):
    """The record from the latest round, breaking ties by rank then by read count."""
    found = records(header)
    if not found:
        raise ValueError(f"No FASTAptamer record in header: {header}")
    return sorted(found, key=lambda r: (-r[0], r[1], -r[2]))[0]


def read_count(header, round_number=None):
    """Read count of a header, from ``round_number`` if given, else from the latest round."""
    if round_number is None:
        found = records(header)
        return best_record(header)[2] if found else 0
    for number, _, count in records(header):
        if number == round_number:
            return count
    return 0


def load_pool(fasta_file):
    """[(header, RNA sequence)] from a FASTA file, in file order."""
    pool, header = [], None
    with open(fasta_file) as handle:
        for line in handle:
            if line.startswith(">"):
                header = line[1:].strip()
            elif header is not None:
                pool.append((header, as_rna(line)))
                header = None
    return pool


# ---------------------------------------------------------------- vsearch arm
def run_vsearch(fasta_file, uc_file, identity=0.7, executable="vsearch"):
    """Cluster a pool with vsearch and write the .uc cluster file."""
    subprocess.run([executable, "--cluster_fast", str(fasta_file), "--id", str(identity),
                    "--strand", "plus", "--uc", str(uc_file)], check=True, capture_output=True)
    return uc_file


def read_clusters(uc_file):
    """{cluster representative: [member headers]} from a vsearch .uc file."""
    members = defaultdict(list)
    with open(uc_file) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "S":
                members[fields[8]].append(fields[8])
            elif fields[0] == "H":
                members[fields[9]].append(fields[8])
    return members


def cluster_seeds(fasta_file, uc_file, top_n=10, round_number=None):
    """Seeds from the ``top_n`` most abundant vsearch clusters.

    Clusters are ranked by the summed read counts of their members; the seed of a cluster
    is its most abundant member. ``round_number`` restricts counting to one SELEX round.
    """
    pool = dict(load_pool(fasta_file))
    counts = {header: read_count(header, round_number) for header in pool}
    clusters = [(members, sum(counts.get(m, 0) for m in members))
                for members in read_clusters(uc_file).values()]
    clusters.sort(key=lambda item: -item[1])
    seeds = []
    for members, abundance in clusters[:top_n]:
        best = max(members, key=lambda m: counts.get(m, 0))
        seeds.append(Seed(best, pool[best], f"cluster of {len(members)} reads",
                          counts.get(best, 0), abundance))
    return seeds


# ---------------------------------------------------------------- STREME arm
def run_streme(fasta_file, output_dir, n_motifs=10, executable="streme"):
    """Discover motifs in a pool with STREME and write its output directory."""
    subprocess.run([executable, "--rna", "--p", str(fasta_file), "--oc", str(output_dir),
                    "--nmotifs", str(n_motifs)], check=True, capture_output=True)
    return output_dir


def read_streme_motifs(streme_txt, n_motifs=None):
    """[(motif id, consensus)] from streme.txt, in STREME's own order."""
    motifs = []
    with open(streme_txt) as handle:
        for line in handle:
            if line.startswith("MOTIF "):
                name = line.split()[1]
                rank, consensus = name.split("-", 1)
                motifs.append((int(rank), name, consensus))
    return [(name, consensus) for _, name, consensus in sorted(motifs)[:n_motifs]]


def read_streme_sites(sites_tsv, input_fasta, motif_ids):
    """{motif id: {sequences of the STREME input that carry the motif}}.

    ``sites.tsv`` identifies a site by the header of the sequence STREME was run on, so
    the headers are resolved against ``input_fasta``, the FASTA given to STREME.
    """
    sequences = dict(load_pool(input_fasta))
    wanted = set(motif_ids)
    sites = defaultdict(set)
    with open(sites_tsv) as handle:
        header = next(handle).rstrip("\n").split("\t")
        motif_column = header.index("motif_ID") if "motif_ID" in header else 0
        id_column = header.index("seq_ID") if "seq_ID" in header else 2
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(motif_column, id_column):
                continue
            if fields[motif_column] in wanted and fields[id_column] in sequences:
                sites[fields[motif_column]].add(sequences[fields[id_column]])
    return sites


def motif_seeds(pool_fasta, streme_dir, n_motifs=10, streme_input=None):
    """Seeds from the reads carrying the top STREME motifs.

    ``streme_input`` is the FASTA STREME was run on, which may be a subset of the pool
    (for example its most abundant reads); it defaults to ``pool_fasta``. Each motif is
    mapped back to the whole pool, and among the reads carrying it the seed is the one
    from the latest round, then the best rank, then the highest read count, so the choice
    does not depend on file order.
    """
    import os
    streme_txt = os.path.join(streme_dir, "streme.txt")
    sites_tsv = os.path.join(streme_dir, "sites.tsv")
    input_fasta = streme_input or pool_fasta

    motifs = read_streme_motifs(streme_txt, n_motifs)
    sites = read_streme_sites(sites_tsv, input_fasta, [name for name, _ in motifs])
    carries = defaultdict(list)
    for name, _ in motifs:
        for sequence in sites[name]:
            carries[sequence].append(name)

    def rank_key(record):
        round_number, rank, count, header, _ = record
        return (-round_number, rank, -count, header)   # latest round, best rank, most reads

    best = {}
    for header, sequence in load_pool(pool_fasta):
        if sequence not in carries:
            continue
        candidate = (*best_record(header), header, sequence)
        for name in carries[sequence]:
            if name not in best or rank_key(candidate) < rank_key(best[name]):
                best[name] = candidate

    seeds = []
    for name, consensus in motifs:
        if name not in best:
            continue
        _, _, count, header, sequence = best[name]
        seeds.append(Seed(header, sequence, consensus, count, count))
    return seeds


def unique(seeds):
    """Drop seeds whose sequence already appeared, keeping the first occurrence."""
    seen, kept = set(), []
    for seed in seeds:
        if seed.sequence not in seen:
            seen.add(seed.sequence)
            kept.append(seed)
    return kept


if __name__ == "__main__":
    # The four round spellings seen in the datasets of the paper.
    assert records("Ishida2020-6R-1-2626-55264.43-0") == [(6, 1, 2626)]
    assert records("RC3H1_GA40NCTGATT_AAN_D_Cycle4-1-213-1333.07-0") == [(4, 1, 213)]
    assert records("0_BOLL_TGTTCG40NGAC_EMJ_4-91-44-166.51-3") == [(4, 91, 44)]
    assert records("ACMP_round7-3-10-1.0") == [(7, 3, 10)]
    # A merged read carries one record per round.
    assert records("Ishida2020-3R-1-9-90.52_Ishida2020-4R-1-585-7400.19") == [(3, 1, 9), (4, 1, 585)]
    # The latest round wins, then the better rank, then the higher count.
    assert best_record("a_round2-9-100-1.0_round7-3-10-1.0") == (7, 3, 10)
    assert read_count("a_round2-9-100-1.0_round7-3-10-1.0") == 10
    assert read_count("a_round2-9-100-1.0_round7-3-10-1.0", round_number=2) == 100
    assert read_count("a_round2-9-100-1.0", round_number=5) == 0
    assert as_rna("acgt\n") == "ACGU"
    assert [s.id for s in unique([Seed("a", "ACGU", "", 1, 1), Seed("b", "ACGU", "", 1, 1),
                                  Seed("c", "GGGG", "", 1, 1)])] == ["a", "c"]
    print("seed checks passed")
