# Benchmark scripts

The analysis scripts behind the benchmarks of the paper, as they were run. They are kept here
for transparency rather than repackaged as a library: each one is a single file that reads the
outputs of the pipeline described in the root README and writes a table.

Run them from the repository root with the environment of `environment.yaml`.

| Script | Produces | Inputs |
|---|---|---|
| `pdb_correspondence.py` | correspondence between predicted elements and PDB structures (Fig. 3, Supplementary Table 2) | included — runs as-is |
| `streme_query_method_comparison.py` | structure-annotation concordance of DCA, MIp, R-scape and RNAalifold against Jolma2020 (Fig. 2c) | query-free run + R-scape and STREME outputs |
| `streme_query_motif_concordance.py` | sequence-motif recovery, alignment-derived vs STREME-only PWMs (Fig. 2b) | query-free run + STREME outputs |
| `denovo_jolma_filters_byarm.py` | recall and precision across the candidate-filtering steps (Supplementary Fig. 1) | the sweep table written by the two scripts above |
| `compare_deletion_correlations.py` | RaptCouple vs MPBind deletion-effect correlations (Discussion) | the published values, included in the file |

`nolp_nussinov.py` is the folding routine used by the comparison scripts; it is the
reference implementation of the no-lonely-pair constraint of `src/structure.py`.

## PDB correspondence

Runs without any external data:

```bash
python benchmark/pdb_correspondence.py --from-inputs
python benchmark/pdb_correspondence.py --selftest     # four synthetic cases
```

`correspondence_inputs/` holds, for each gene and dataset, the positional probabilities, the
pairable fractions, the predicted pairs and the candidate segments of the twenty models that
were compared, exported from the trained models. `pdb_concordance_reselected.tsv` lists the
PDB elements. The result is written to `pdb_correspondence.tsv` and `pdb_correspondence.md`,
and reproduces the published table exactly.

Without `--from-inputs` the script looks for the trained models themselves and re-exports the
inputs from them; that path needs the SELEX data directory, which is not part of this
repository.

## The Jolma2020 benchmarks

`streme_query_method_comparison.py` and `streme_query_motif_concordance.py` take a dataset
directory holding one subdirectory per RNA-binding protein, each with the models of a
query-free run (`results.tsv` and `repNN.model_params`), the STREME output that seeded it
(`<RBP>/streme10k/out/streme.txt`) and, for the comparison script, R-scape covariation output
(`<subdir>_rscape/<RBP>__repNN_1.cov`). Build that directory with
`scripts/run_query_free.py --arm streme` on the Jolma2020 pool, run R-scape on the resulting
alignments, then:

```bash
python benchmark/streme_query_method_comparison.py <dataset_dir> streme_query
python benchmark/streme_query_motif_concordance.py <dataset_dir> streme_query
```

Both also need the Jolma2020 reference motifs and metadata read by `src/helper_jolma.py`.
The trained models of the Jolma2020 and Laverty2023 runs are in `example/`, so the comparison
can be repeated from them once the STREME and R-scape outputs are regenerated.
