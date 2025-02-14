import subprocess 
import yaml
import os 
import pandas as pd 
from Bio import SeqIO
import re

class FASTAProcessor:
    def __init__(self, config_path):
        self.config_path = config_path
        self.config = yaml.safe_load(open(config_path, "r"))
        self.fasta_count = False
        self.fasta_count_ann = False
        self.fasta_count_ann_merge = False
        self.fasta_count_ann_merge_reunique = False

    def cutadapt_of_all_fasta(self):
        for fasta_file in self.config["fasta_annotation"].keys():
            adapter_5 = self.config["adapter_5"]
            adapter_3 = self.config["adapter_3"]
            N_random = self.config["N_random"]

            fasta_trim_tmp_file = os.path.join(self.config["data_dir"], fasta_file.replace(".fa", ".trim.tmp.fa"))
            fasta_trim_file = os.path.join(self.config["data_dir"], fasta_file.replace(".fa", ".trim.fa"))
            
            cmd_5 = f"""cutadapt --discard-untrimmed -g {adapter_5} {os.path.join(self.config["data_dir"], fasta_file)} -o {fasta_trim_tmp_file}"""
            res = subprocess.run(cmd_5, shell=True, capture_output=True)
            print(res.stdout.decode("utf-8"))
            print(res.stderr.decode("utf-8"))

            cmd_3 = f"""cutadapt --minimum-length {N_random} --discard-untrimmed -a {adapter_3} {fasta_trim_tmp_file} -o {fasta_trim_file}"""
            print(cmd_3)
            res = subprocess.run(cmd_3, shell=True, capture_output=True)
            print(res.stdout.decode("utf-8"))
            print(res.stderr.decode("utf-8"))

        return 

    def _remove_duplicate_of_trimfasta_by_fastaptamer(self, fasta_file):
        fasta_trim_file = fasta_file.replace(".fa", ".trim.fa")
        fasta_trim_count_file = fasta_file.replace(".fa", ".trim.count.fa")

        cmd = f"""fastaptamer_count -i {fasta_trim_file} -o {fasta_trim_count_file}"""
        res = subprocess.run(cmd, shell=True, capture_output=True)

        print(res.stdout.decode("utf-8"))
        print(res.stderr.decode("utf-8"))

        return fasta_trim_count_file
    
    def remove_duplicate_of_all_trimfasta(self):
        print("Removing duplicates in a fasta")
        for fasta_file in self.config["fasta_annotation"].keys():
            print(f"Removing duplicates in {fasta_file}")
            self._remove_duplicate_of_trimfasta_by_fastaptamer(os.path.join(self.config["data_dir"], fasta_file))
        self.fasta_count = True
        return

    def _add_info_to_trimcount_fasta_by_config(self, ann, trimcountfasta_file):
        cmd = f"""sed 's/>/>{ann}-/g' {trimcountfasta_file} > {trimcountfasta_file.replace(".fa", ".ann.fa")}"""
        res = subprocess.run(cmd, shell=True, capture_output=True)
        # print(res.stdout.decode("utf-8"))
        # print(res.stderr.decode("utf-8"))
        return res
    
    def add_info_to_all_trimcount_fasta(self):
        assert self.fasta_count, "Run remove_duplicate_of_all_fasta"
        for fasta_file in self.config["fasta_annotation"].keys():
            self._add_info_to_trimcount_fasta_by_config(
                ann = self.config["fasta_annotation"][fasta_file],
                trimcountfasta_file = os.path.join(self.config["data_dir"], fasta_file.replace(".fa", ".trim.count.fa"))
                )
        self.fasta_count_ann = True
        return 

    def merge_all_fasta(self, fasta_merged_file):
        assert self.fasta_count_ann, "Run add_info_to_all_counted_fasta"
        self.fasta_merged_file = fasta_merged_file
        print(fasta_merged_file)
        with open(fasta_merged_file, 'w') as fasta_trim_count_ann_all:
            for fasta_file in self.config["fasta_annotation"].keys():
                fasta_trim_count_ann_file = fasta_file.replace(".fa", ".trim.count.ann.fa")
                with open(os.path.join(self.config["data_dir"], fasta_trim_count_ann_file), "r") as fasta_trim_count_ann:
                    fasta_trim_count_ann_all.write(fasta_trim_count_ann.read())
        self.fasta_count_ann_merge = True
        return fasta_merged_file

    # optional
    def remove_count1_from_trimcountann_fasta(self):
        if not "remove_lowcount" in self.config:
            return 
    
        else: 
            records = SeqIO.parse(self.fasta_merged_file, "fasta")
            
            for fasta_file, annot in self.config["fasta_annotation"].items():
                if fasta_file in self.config["remove_lowcount"].keys():
                    count = self.config["remove_lowcount"][fasta_file]
                    rm_count_range = "|".join([str(c) for c in range(1, count+1)])
                    records = [record for record in records if re.search(f"{annot}-.*-[{rm_count_range}]-.*", record.id) == None]
                    print(len(records))

            SeqIO.write(
                records,
                os.path.join(self.fasta_merged_file.replace(".fa", ".rmlow.fa")),
                "fasta"
                )
            self.fasta_merged_file = self.fasta_merged_file.replace(".fa", ".rmlow.fa")
            return

    # def reuniquenize(self):
    #     assert self.fasta_count_ann_merge, "Run merge_all_fasta"

    #     records = SeqIO.parse(self.fasta_merged_file, "fasta")    
    #     records_dicts = {}
    #     for record in records:
    #         records_dicts[str(record.seq)] = record.id

    #     df_records = pd.DataFrame.from_dict(records_dicts, columns=["id"], orient="index").reset_index().rename(columns={"index": "seq"})
    #     df_records = df_records.sort_values(by="id").reset_index(drop=True)

    #     # add int to id if duplicated
    #     id_cumc = ["-".join([idx, str(cumc)]) for idx, cumc in zip(df_records["id"], df_records.groupby("id").cumcount())]
    #     df_records["id"] = id_cumc

    #     self.fasta_merge_reunique_file = self.fasta_merged_file.replace(".fa", ".unique.fa")
    #     with open(self.fasta_merge_reunique_file, "w") as f:
    #         for idx, row in df_records.iterrows():
    #             f.write(f">{row['id']}\n{row['seq']}\n")
    #     self.fasta_count_ann_merge_reunique = True

    #     return 
    
    def reuniquenize(self):
        """
        retain duplicated ids.
        v2
        """
        assert self.fasta_count_ann_merge, "Run merge_all_fasta"

        records = SeqIO.parse(self.fasta_merged_file, "fasta")    
        records_dicts = {}
        for record in records:
            if str(record.seq) in records_dicts:
                records_dicts[str(record.seq)].append(record.id)
            else:
                records_dicts[str(record.seq)] = [record.id]

        df_records = pd.DataFrame([(seq, ids) for seq, ids in records_dicts.items()], columns=["seq", "ids"])
        df_records = df_records.sort_values(by="ids").reset_index(drop=True)

        # Merge IDs for identical sequences
        df_records["merged_id"] = df_records["ids"].apply(lambda x: "_".join(x))

        # add int to id if duplicated
        id_cumc = ["-".join([idx, str(cumc)]) for idx, cumc in zip(df_records["merged_id"], df_records.groupby("merged_id").cumcount())]
        df_records["merged_id_"] = id_cumc

        self.fasta_merge_reunique_file = self.fasta_merged_file.replace(".fa", ".unique.fa")
        with open(self.fasta_merge_reunique_file, "w") as f:
            for idx, row in df_records.iterrows():
                f.write(f">{row['merged_id_']}\n{row['seq']}\n")
        self.fasta_count_ann_merge_reunique = True

        return 

    
    def print_proc_info(self):
        print(f"Config: {self.config_path}")
        print(f"Final trim/count/ann/merged/reunique fasta: {self.fasta_count_ann_merge_reunique}")
        return


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=str, help="selex config file path")
    parser.add_argument("--merged_fasta", required=True, type=str, help="file name of merged fasta")

    args = parser.parse_args()

    proc = FASTAProcessor(args.config)
    proc.cutadapt_of_all_fasta()
    proc.remove_duplicate_of_all_trimfasta()
    proc.add_info_to_all_trimcount_fasta()
    proc.merge_all_fasta(args.merged_fasta)
    proc.remove_count1_from_trimcountann_fasta()
    proc.reuniquenize()
    proc.print_proc_info()
    