# Worked examples

Short, self-contained scripts that reproduce results from the RaptCouple paper using only
files in this repository. Each runs in seconds and prints what it checks, so they double as
a starting point for your own analyses.

Run them from the repository root:

```bash
python example/reproduce/01_structure_from_couplings.py
python example/reproduce/02_deletion_effects.py
python example/reproduce/03_motif_candidates.py
```

| Script | What it shows | Paper |
|---|---|---|
| `01_structure_from_couplings.py` | Folding a trained model: couplings with z-score >= 3 become base pairs, and the no-lonely-pair constraint keeps only stacked pairs | Structure of the Ishida2020 2'-F aptamer |
| `02_deletion_effects.py` | Scoring mutations, including deletions (a change to the model's gap state) and the coupling between two simultaneous changes; then the correlation with measured binding | Deletion-effect correlation, Spearman's rho = 0.710, P = 0.003 |
| `03_motif_candidates.py` | Extracting sequence-motif candidates from a trained model: informative blocks, deduplication, low-complexity removal | Candidate extraction of the query-free benchmark |

Scripts 01 and 02 assert the published values, so they fail if the result changes.

## Going further

These start from models that are already trained. To run the whole pipeline on your own
data, see the README in the repository root:

- **target-guided** (you supply the query): `merge_and_cutadapt_all_rounds.py` ->
  `run_jackhmmer.py` -> `train_potts.py` -> `fold_by_coupling.py`
- **query-free** (the queries are generated from the pool): `run_query_free.py`

The datasets under `example/` each carry a config YAML and the trained models, so you can
also start from any of them.
