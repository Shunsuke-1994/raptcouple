"""Reproduce the deletion-effect correlation for the Ishida2020 aptamer (Fig. 6b, left).

Two steps:

1. Score deletions with the Potts model. Deleting a nucleotide is evaluated as a change to
   the model's gap state, so the effect is an energy difference the model can compute
   directly. This is what `scripts/predict_mutation_effects.py` does.
2. Correlate the predicted effects of the measured deletion variants with their relative
   binding. The measured values and the scores used in the paper are in
   `example/Ishida2020/variants/mutation_effect_predictions.csv`.

    python example/reproduce/02_deletion_effects.py

Runs in seconds; needs only files in the repository. Reported in the paper as
Spearman's rho = 0.710, P = 0.003 over 15 distinct deletion variants.
"""
import csv
import os
import sys

import numpy as np
from scipy.stats import spearmanr

sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))
from src.plmc import read_params
from src.potts import PottsModel
from src.util import onehot2seq

MODEL = "example/Ishida2020/outputs/Ishida2020-6R-1-2626-55264.43.model_params"
MEASURED = "example/Ishida2020/variants/mutation_effect_predictions.csv"

# ---------------------------------------------------------------- 1. scoring deletions
params = read_params(MODEL)
alphabet = params["alphabet"]
model = PottsModel.build_from_file(MODEL)
sequence = onehot2seq(model.spins, is_dna=("T" in alphabet), is_gapped=("." in alphabet))
gap = alphabet.index(".")

print(f"model alphabet {alphabet}; '.' is the gap state used to score a deletion")
print(f"target sequence ({len(sequence)} alignment columns)\n  {sequence}\n")

print("deleting each of the first informative positions:")
print(f"{'position':>9}  {'base':>4}  {'effect (-dH)':>13}")
deleted = [i for i, base in enumerate(sequence) if base != "."][:8]
for column in deleted:
    state = alphabet.index(sequence[column])
    effect = -model.compute_delta_energy([(state, column, gap)])
    print(f"{column + 1:>9}  {sequence[column]:>4}  {effect:>13.3f}")

# Deleting two positions at once is not the sum of the two single deletions: the model
# includes the coupling between them.
first, second = deleted[0], deleted[1]
both = -model.compute_delta_energy([(alphabet.index(sequence[first]), first, gap),
                                    (alphabet.index(sequence[second]), second, gap)])
separately = sum(-model.compute_delta_energy([(alphabet.index(sequence[c]), c, gap)])
                 for c in (first, second))
print(f"\npositions {first + 1} and {second + 1} together: {both:.3f}")
print(f"the two single deletions added:      {separately:.3f}")
print(f"difference (their coupling):         {both - separately:+.3f}")

# ---------------------------------------------------------------- 2. measured variants
rows = list(csv.DictReader(open(MEASURED)))
predicted = np.array([-float(row["delta energy"]) for row in rows])
binding = np.array([float(row["relative binding"]) for row in rows])
rho, p_value = spearmanr(predicted, binding)

print(f"\n{len(rows)} measured deletion variants")
print(f"Spearman's rho = {rho:.3f}, P = {p_value:.3f}   (paper: 0.710, 0.003)")
print(f"\n{'predicted (-dH)':>16}  {'relative binding (%)':>20}")
for k in np.argsort(-predicted):
    print(f"{predicted[k]:>16.3f}  {binding[k]:>20.1f}")

best = int(np.argmax(binding))
rank = int(np.where(np.argsort(-predicted) == best)[0][0]) + 1
print(f"\nthe variant with the highest measured binding ranks {rank} of {len(rows)} "
      f"by predicted effect")

assert abs(rho - 0.710) < 0.001 and p_value < 0.005, (rho, p_value)
print("matches the correlation reported in the paper")
