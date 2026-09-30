"""Extract sequence-motif candidates from trained models, as in the query-free benchmark.

`scripts/run_query_free.py` runs the whole pipeline from a SELEX pool. This sample starts
one step later, from models that are already in the repository, and shows the candidate
extraction on its own: positional frequencies -> informative blocks -> deduplication ->
low-complexity removal.

    python example/reproduce/03_motif_candidates.py [gene ...]

With no argument it uses a few RNA-binding proteins from the Jolma2020 collection. Runs in
seconds; needs only files in the repository.
"""
import contextlib
import io
import os
import sys
from glob import glob

sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))
from src import motif
from src.plmc import read_params

MODELS = "example/Jolma2020/outputs"
DEFAULT = ["RBM4", "RBFOX1", "ELAVL1", "BOLL"]


def models_of(gene):
    """Model files whose name starts with the gene, e.g. `0_RBM4_...model_params`."""
    return sorted(glob(os.path.join(MODELS, f"*_{gene}_*.model_params")))


genes = sys.argv[1:] or DEFAULT
for gene in genes:
    files = models_of(gene)
    if not files:
        print(f"{gene}: no model in {MODELS}")
        continue
    path = files[0]
    with contextlib.redirect_stdout(io.StringIO()):     # read_params prints a summary
        params = read_params(path)

    frequencies = motif.positional_frequencies(params)
    blocks = motif.segment(frequencies)
    unique = motif.deduplicate(blocks)
    candidates = motif.candidates(params, remove_low_complexity=True, top_n=None)

    print(f"\n{gene}  ({os.path.basename(path)})")
    print(f"  {frequencies.shape[1]} columns occupied by the focus sequence")
    print(f"  {len(blocks)} informative block(s) -> {len(unique)} after deduplication "
          f"-> {len(candidates)} after removing low-complexity")
    for candidate in candidates:
        print(f"    {candidate['consensus']:<26} core {candidate['core_kmer']}")
    dropped = [motif.consensus(block) for block in unique if motif.is_low_complexity(block)]
    for sequence in dropped:
        print(f"    {sequence:<26} dropped as low-complexity")

print("\nTo rank candidates by prevalence and keep the top three, as in the paper's")
print("benchmark, pass the pool reads:")
print("    reads = motif.load_reads('<pool>.fa')")
print("    motif.candidates(params, reads=reads, top_n=3)")
