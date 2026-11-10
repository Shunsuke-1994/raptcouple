# PDB correspondence: all query-free candidates

Rule and null as in `pdb_correspondence.py` (1000 permutations, seed 20261001, alpha 0.05, pairable if >= 0.5 of sequences). Element positions are 1-based on the PDB RNA; offsets are 1-based matched columns of the model.

| RBP | PDB | dataset | models | joint (seq / pair) | p | p (seq only) | in a segment | predicted / shared pairs | outcome |
|---|---|---|---|---|---|---|---|---|---|
| RC3H1 | 5F5F | Jolma2020 | 20 | 0.727 (0.795 / 1.0) | 0.009 | 0.004 | no | 5 / 4 | corresponding, stem predicted |
| RC3H1 | 5F5F | Laverty2022 | 20 | 0.666 (0.749 / 0.833) | 0.001 | 0.004 | no | 5 / 4 | corresponding, stem predicted |
| RC3H2 | 4ZLD | Jolma2020 | 20 | 0.667 (0.897 / 0.667) | 0.9451 | 0.9281 | yes | 0 / 0 | no corresponding candidate |
| RC3H2 | 4ZLD | Laverty2022 | 0 | | | | | | not assayed |
| RBFOX1 | 2N82 | Jolma2020 | 60 | 0.414 (0.56 / 0.5) | 0.4086 | 0.6513 | yes | 0 / 0 | no corresponding candidate |
| RBFOX1 | 2N82 | Laverty2022 | 20 | 0.305 (0.479 / 0.5) | 0.2637 | 0.2188 | yes | 0 / 0 | no corresponding candidate |
| SNRPA | 1AUD | Jolma2020 | 20 | 0.5 (0.629 / 0.5) | 0.991 | 0.991 | no | 0 / 0 | no corresponding candidate |
| SNRPA | 1AUD | Laverty2022 | 20 | 0.375 (0.532 / 0.5) | 0.8931 | 1.0 | yes | 0 / 0 | no corresponding candidate |
| SNRNP70 | 4PKD | Jolma2020 | 20 | 0.22 (0.415 / 0.5) | 0.9061 | 0.97 | no | 0 / 0 | no corresponding candidate |
| SNRNP70 | 4PKD | Laverty2022 | 0 | | | | | | not assayed |
| SNRPB2 | 1A9N | Jolma2020 | 20 | 0.385 (0.539 / 0.4) | 0.4835 | 0.6254 | yes | 0 / 0 | no corresponding candidate |
| SNRPB2 | 1A9N | Laverty2022 | 0 | | | | | | not assayed |
| LARP7 | 6D12 | Jolma2020 | 20 | 0.685 (0.764 / 0.714) | 0.1189 | 0.6264 | no | 0 / 0 | no corresponding candidate |
| LARP7 | 6D12 | Laverty2022 | 0 | | | | | | not assayed |
| NOVA2 | 1EC6 | Jolma2020 | 20 | 0.349 (0.512 / 0.5) | 0.8062 | 0.996 | yes | 0 / 0 | no corresponding candidate |
| NOVA2 | 1EC6 | Laverty2022 | 0 | | | | | | not assayed |
| LIN28A | 5UDZ | Jolma2020 | 0 | | | | | | not assayed |
| LIN28A | 5UDZ | Laverty2022 | 20 | 0.197 (0.398 / 0.333) | 1.0 | 0.998 | no | 0 / 0 | no corresponding candidate |
| RBMY1A1 | 2FY1 | Jolma2020 | 0 | | | | | | not assayed |
| RBMY1A1 | 2FY1 | Laverty2022 | 20 | 0.625 (0.798 / 0.625) | 0.023 | 0.6713 | no | 0 / 0 | corresponding, no stem predicted |
| SNU13 | 2OZB | Jolma2020 | 0 | | | | | | not assayed |
| SNU13 | 2OZB | Laverty2022 | 20 | 0.263 (0.447 / 0.333) | 0.3477 | 0.6174 | yes | 0 / 0 | no corresponding candidate |
| NOVA1 | 2ANR | Jolma2020 | 0 | | | | | | not assayed |
| NOVA1 | 2ANR | Laverty2022 | 20 | 0.2 (0.43 / 0.2) | 1.0 | 0.6573 | no | 0 / 0 | no corresponding candidate |
| SLBP | 4L8R | Jolma2020 | 0 | | | | | | not assayed |
| SLBP | 4L8R | Laverty2022 | 20 | 0.5 (0.759 / 0.5) | 0.7193 | 0.6663 | yes | 0 / 0 | no corresponding candidate |
| Vts1 | 2ESE | Jolma2020 | 0 | | | | | | not assayed |
| Vts1 | 2ESE | Laverty2022 | 20 | 0.517 (0.638 / 0.8) | 0.5894 | 0.9281 | no | 3 / 3 | no corresponding candidate |
