# raptcouple_test
Code of RaptCouple, a unsupervised machine learning of SELEX data. Raptcouple learns structure and fitness information from SELEX data.

# Environment
```
mamba env create -f environment.yaml
```
## install plmc
After installing [plmc](https://github.com/debbiemarkslab/plmc) as the instruction, please edit `PLMC_TO_PATH` variable in `src/coupling.py`.  


# Description
`example/Ishida2020/Ishida2020.ipynb` contains the whole workflow described below.  
## Data preparation
```
.
├── DRR201870.fa
├── DRR201871.fa
├── DRR201872.fa
├── DRR201873.fa
└── config.yaml

```
## Preprocessing
`python scripts/merge_and_cutadapt_all_rounds.py` performs preprocessing based on `config.yaml`.  
This script performs:  
1. cutadapt  
2. sequence merging  
3. remove seqs of small count  
`config.yaml` should follow this format:   
```
N_random: 40
adapter_3: TATGTGCGCATACATGGATCCTC
adapter_5: TAATACGACTCACTATAGGGAGAACTTCGACCAGAAG
data_dir: ./example/Ishida2020/data
fasta_annotation:
  DRR201870.fa: Ishida2020-3R
  DRR201871.fa: Ishida2020-4R
  DRR201872.fa: Ishida2020-5R
  DRR201873.fa: Ishida2020-6R
# remove_lowcount: # remove sequences which count is smaller than mincount.
#   DRR201870.fa: 1
#   DRR201871.fa: 1
```

This is an example of preprocessing.
```
python scripts/merge_and_cutadapt_all_rounds.py  \
    --config ./example/data/Ishida2020/config.yaml \
    --merged_fasta ./example/data/Ishida2020/Ishida2020.count.ann.all_selex.fa
```
## MSA constraction by jackhmmer
Generate a multiple sequence alignment (MSA) using jackhmmer:
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
Once you have obtained coupling scores from the Potts model training, predict the 2D structure by using the coupling information. For example: 
```

```


## Prediction of mutation effects
Evaluate the impact of mutations on sequence fitness and structure with:  
```

```

`mutations.txt` should list mutations in a standard format (e.g., A15G).
The script outputs predicted effects for each mutation, facilitating the analysis of mutation impact.

# Citation
If you use this code, please cite the following paper:

```bibtex
@article{
  title={Learning structure and fitness of RNA discovered by SELEX},
  author={Sumi, Shunsuke and Kawahara, Daiki and Hada, Yuki and Yoshii, Tatsuyuki and Adachi, Tatsuo and Saito, Hirohide and Hamada, Michiaki},
  journal={Journal Name},        % 論文掲載ジャーナル名に置き換えてください
  volume={XX},                    % 巻番号に置き換えてください
  number={YY},                    % 号番号に置き換えてください
  pages={ZZ-ZZ},                  % ページ番号に置き換えてください
  year={2025},                    % 発行年に置き換えてください
  note={Correspondence should be addressed to: mhamada@waseda.jp, hirosaito@iqb.u-tokyo.ac.jp}
}