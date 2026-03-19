# RaptCouple
Code of RaptCouple, an unsupervised machine learning framework for SELEX data. RaptCouple learns structure and fitness information from SELEX data.

# Environment
```
mamba env create -f environment.yaml
```
## Install plmc
After installing [plmc](https://github.com/debbiemarkslab/plmc) as the instruction, please edit `PLMC_TO_PATH` variable in `src/plmc.py`.


# Description
`example/Ishida2020/Ishida2020.ipynb` contains the whole workflow described below.

## Examples
Pre-computed outputs (MSA, model parameters, coupling scores) are available under the `example/` directory:
```
example/
├── Ishida2020/       # Ishida et al. 2020
├── Adachi2024/       # Adachi et al. 2024
├── Jolma2020/        # Jolma et al. 2020
├── Laverty2023/      # Laverty et al. 2023
├── Methylation/
├── PRJDB19138/
├── PRJDB19139/
├── PRJDB19140/
└── PRJDB19141/
```
Each example contains a config YAML file and an `outputs/` directory with pre-computed results. Raw SELEX FASTA files are not included; place them under the path specified by `data_dir` in the config file.

## Data preparation
Raw FASTQ files can be downloaded from [SRA](https://www.ncbi.nlm.nih.gov/sra). After downloading, preprocess them according to the config YAML file to obtain FASTA files for each SELEX round.
```
data_dir/
├── 2nd_round.fa
├── 3rd_round.fa
├── 4th_round.fa
├── 5th_round.fa
└── 6th_round.fa
```
Config file example:
```yaml
N_random: 40
adapter_3: TATGTGCGCATACATGGATCCTC
adapter_5: TAATACGACTCACTATAGGGAGAACTTCGACCAGAAG
data_dir: ./where_the_data_is
fasta_annotation:
  2nd_round.fa: 2R
  3rd_round.fa: 3R
  4th_round.fa: 4R
  5th_round.fa: 5R
  6th_round.fa: 6R
# remove_lowcount:  # optional: remove sequences with count smaller than threshold
#   2nd_round.fa: 1
#   3rd_round.fa: 1
```

## Preprocessing
`python scripts/merge_and_cutadapt_all_rounds.py` performs preprocessing based on the config file.
This script performs:
1. cutadapt & fastaptamer_count
2. sequence merging
3. remove seqs of small count (optional)

```
python scripts/merge_and_cutadapt_all_rounds.py --config ./example/Ishida2020/config_6R_rank1.yaml
```
## MSA construction
Set the MSA parameters (jackhmmer) in config.yaml as follows:
```yaml
MSA_parameters:
  all_fasta: ./example/Ishida2020/data/Ishida2020.count.ann.all_selex.unique.fa
  target_id: Ishida2020-6R-1-2626-55264.43-0
  save_dir: ./example/Ishida2020/outputs
  prefix: ""
  iters: 10
  F1: 0.02
  F2: 0.001 # 1e-3
  F3: 0.0001 # 1e-4
  T: 5
  domT: 5
  incT: 5
  incdomT: 5
  print_result: true
```
Generate a multiple sequence alignment (MSA) using jackhmmer:
```
python ./scripts/run_jackhmmer.py --config ./example/Ishida2020/config_6R_rank1.yaml
```
Note: We found these parameters work for most SELEX data. But, if the MSA depth is insufficient, consider relaxing the jackhmmer parameters (iters, F1, F2, F3, T, domT, incT, incdomT). For further details, please refer to the HMMER3 user guide.

## Potts model training
Set Potts model parameters in config.yaml:
```yaml
Potts_parameters:
  input_fasta: ./example/Ishida2020/outputs/Ishida2020-6R-1-2626-55264.43-0.msa
  sim_threshold: 0.05 # theta
  vocab: AUGC.
  iters: 200
  suffix: ""
  print_result: true
```
sim_threshold is re-weighting parameters of each sequence in MSA. If the sequences are highly similar, smaller sim_threshold may be suitable.
Train the Potts model by running:
```
python scripts/train_potts.py --config ./example/Ishida2020/config_6R_rank1.yaml
```

## Model parameters file format
The `.model_params` file produced by plmc is a binary format based on the [MATLAB reader in the plmc repository](https://github.com/debbiemarkslab/plmc/blob/master/scripts/read_params.m). A Python reader is provided in `src/plmc.py` as `read_params()`:
```python
from src.plmc import read_params
params = read_params("example/Ishida2020/outputs/Ishida2020-6R-1-2626-55264.43-0.model_params")
# params contains: target_seq, alphabet, hi (fields), Jij (couplings), FN_apc (Frobenius norms), etc.
```

## Folding with coupling scores
Once you have obtained coupling scores from the Potts model training, predict the 2D structure by using the coupling information. For example:
```
python scripts/fold_by_coupling.py --coupling ./example/Ishida2020/outputs/Ishida2020-6R-1-2626-55264.43-0.model_params --min_loop_len 3 --z_threshold 2 --output ./example/Ishida2020/outputs/fold.yaml
```

## Prediction of mutation effects
Evaluate the impact of mutations on sequence fitness and structure with:
```
python scripts/predict_mutation_effects.py --param_file example/Ishida2020/outputs/Ishida2020-6R-1-2626-55264.43-0.model_params --mutations_file ./example/Ishida2020/variants/mutations.txt > ./example/Ishida2020/variants/mutations_effect_prediction.txt
```
or
```
python scripts/predict_mutation_effects.py --param_file example/Ishida2020/outputs/Ishida2020-6R-1-2626-55264.43-0.model_params --mutations G1A,A21.
```

`mutations.txt` should list mutations in a standard format (e.g., A15G).
The script outputs predicted effects for each mutation, facilitating the analysis of mutation impact.

## Sampling and Annealing Scripts

### Gibbs Sampling

Generate sequences via Gibbs sampling and output them in FASTA format along with energy values. For example, run the following command:
```
python scripts/gibbs_sampling.py --param_file example/Ishida2020/outputs/Ishida2020-6R-1-2626-55264.43-0.model_params > ./example/Ishida2020/outputs/gibbs_sampling_output.fa
```

### Simulated Annealing

Generate sequences via simulated annealing and output them in FASTA format along with energy values. For example, run the following command:
```
python scripts/simulated_annealing.py --param_file ./example/Ishida2020/outputs/Ishida2020-6R-1-2626-55264.43-0.model_params > ./example/Ishida2020/outputs/simulated_annealing_output.fa
```

# Citation
If you use this code, please cite the following paper:

```bibtex
@article{sumi2025raptcouple,
  title={Discovering structural and functional landscapes of nucleic acids through in vitro evolution},
  author={Sumi, Shunsuke and Kawahara, Daiki and Hada, Yuki and Yoshii, Tatsuyuki and Adachi, Tatsuo and Saito, Hirohide and Hamada, Michiaki},
  journal={submitted},
  volume={XX},
  number={YY},
  pages={ZZ-ZZ},
  year={2025},
  note={Correspondence should be addressed to: mhamada@waseda.jp, hirosaito@iqb.u-tokyo.ac.jp}
}
