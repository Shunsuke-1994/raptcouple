"""Reproduce the coupling-supported secondary structure of the Ishida2020 aptamer (Fig. 4a).

Folds the bundled Potts model with the settings stated in the paper's Methods: base pairs
are placed only where the APC-corrected coupling z-score is at least 3, the minimum hairpin
loop is 3 nucleotides, and the no-lonely-pair constraint requires every pair to be stacked
on a neighbour.

    python example/reproduce/01_structure_from_couplings.py

Runs in seconds; needs only the model file in the repository.
"""
import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))
from src.plmc import detect_coupling, read_params
from src.structure import fold

MODEL = "example/Ishida2020/outputs/Ishida2020-6R-1-2626-55264.43.model_params"
EXPECTED = ".(((....((............))...))).........."

params = read_params(MODEL)
sequence, energy, structure, _ = fold(params, min_loop_length=3, threshold=3, nolp=True)

print(f"model     {MODEL}")
print(f"sequence  {sequence}")
print(f"structure {structure}")
print(f"score     {energy:.3f}   base pairs {structure.count('(')}")

# detect_coupling works in alignment columns; the folded sequence keeps only the columns
# the focus sequence occupies, so map the positions before printing them.
columns = [i for i, base in enumerate(params["target_seq"]) if base in "AUGC"]
position = {column: k for k, column in enumerate(columns)}
couplings = detect_coupling(params, threshold=3, min_dist=4)
print(f"\ncouplings with z-score >= 3 and separation > 3: {len(couplings)}")
for (base_i, base_j), (i, j) in couplings:
    paired = structure[position[i]] == "(" and structure[position[j]] == ")"
    print(f"  {base_i}{position[i] + 1:>3}-{base_j}{position[j] + 1:<3}"
          f"  {'in the structure' if paired else 'not retained by noLP folding'}")

# The same model folded without the no-lonely-pair constraint keeps isolated pairs.
_, _, legacy, _ = fold(params, min_loop_length=3, threshold=3, nolp=False)
print(f"\nwithout noLP: {legacy}  ({legacy.count('(')} base pairs)")

assert structure == EXPECTED, f"expected {EXPECTED}, got {structure}"
print("\nmatches the structure reported in the paper")
