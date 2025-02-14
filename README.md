# raptcouple_test
Code of RaptCouple, a unsupervised machine learning of SELEX data.  


# Install
```
mamba env create -f environment.yaml
```
# Description
## Data preparation
```
.
├── LHa4-1R_S5_L001_R1_001.fa
├── LHa4-2R_S6_L001_R1_001.fa
├── LHa4-3R_S7_L001_R1_001.fa
├── LHa4-4R_S8_L001_R1_001.fa
└── config.yaml

```
## Preprocessing
`python scripts/merge_and_cutadapt_all_rounds.py` performs preprocessing based on `config.yaml`.  
Preprocessing is done in the following order: 1.cutadapt, 2.sequence merging, 3.remove seqs of small count.  
`config.yaml` should be in this format. 
```
N_random: 20  # length of random region
adapter_3: TGCTCGAGAAAAACAAGACGTAA  # 3' adapter sequence
adapter_5: AATACGACTCACTATAGGGCA # 5' adapter sequence
data_dir: ./SELEXdata/Ribomic/RAPT27
fasta_annotation: # fasta files to be analyzed
  LHa4-1R_S5_L001_R1_001.fa: RAPT27-1R
  LHa4-2R_S6_L001_R1_001.fa: RAPT27-2R
  LHa4-3R_S7_L001_R1_001.fa: RAPT27-3R
  LHa4-4R_S8_L001_R1_001.fa: RAPT27-4R
remove_lowcount: # remove sequences which count is smaller than mincount.
  LHa4-1R_S5_L001_R1_001.fa: 1
  LHa4-2R_S6_L001_R1_001.fa: 1
```

This is an example of preprocessing.
```
python scripts/merge_and_cutadapt_all_rounds.py  \
    --config ./SELEXdata/Ribomic/RAPT27/config.yaml \
    --merged_fasta ./SELEXdata/Ribomic/RAPT27/RAPT27.count.ann.all_selex.fa
```
## MSA constraction by jackhmmer
```
python ./scripts/run_jackhmmer.py --help
usage: run_jackhmmer.py [-h] --fasta FASTA --target TARGET --save_dir SAVE_DIR [--prefix PREFIX] [--iters ITERS] [--F1 F1] [--F2 F2] [--F3 F3] [--T T] [--domT DOMT] [--incT INCT] [--incdomT INCDOMT] [--print_result]

Run jackhmmer and plmc on a fasta file

options:
  -h, --help           show this help message and exit
  --fasta FASTA        fasta file
  --target TARGET      target id of fasta file
  --save_dir SAVE_DIR  directory to save results
  --prefix PREFIX      prefix of output files
  --iters ITERS        number of iterations for jackhmmer (default: 10)
  --F1 F1              threshold F1 for jackhmmer (default: 0.02)
  --F2 F2              threshold F2 for jackhmmer (default: 1e-3)
  --F3 F3              threshold F3 for jackhmmer (default: 1e-4)
  --T T                threshold T for jackhmmer (default: 5)
  --domT DOMT          threshold domT for jackhmmer (default: 5)
  --incT INCT          threshold incT for jackhmmer (default: 5)
  --incdomT INCDOMT    threshold incdomT for jackhmmer (default: 5)
  --print_result
```
## Potts model training

```
python scripts/train_potts.py --help
usage: train_potts.py [-h] --input_fasta INPUT_FASTA [--target TARGET] [--threshold THRESHOLD] [--vocab VOCAB] [--iters ITERS] [--suffix SUFFIX] [--print_result]

options:
  -h, --help            show this help message and exit
  --input_fasta INPUT_FASTA
                        input fasta file (=alignment)
  --target TARGET       specified target
  --threshold THRESHOLD
                        similarity threshold for plmc (default: 0.05)
  --vocab VOCAB         vocabulary for plmc (default: AUGC.)
  --iters ITERS         number of iterations for plmc (default: 200)
  --suffix SUFFIX       suffix for output files (default: '')
  --print_result        print result of plmc
```
## Folding with coupling scores

## Prediction of mutation effects

