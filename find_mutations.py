#!/usr/bin/env python3
from pysam import VariantFile, tabix_index

import sys
import os
from collections import defaultdict
from datetime import date
import shutil

def read_gene_coordinates(gene_file):
    gene_coordinates = dict()
    with open(gene_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.strip().split(",")
            if len(line) < 4:
                print(f"Invalid line: {line}. Skipping...")
                continue
            gene_id = line[0]
            chromosome = line[1]
            start = int(line[2])
            end = int(line[3])
            gene_coordinates[gene_id] = (chromosome, start, end)
    return gene_coordinates

def agregate_gene_mutations(vcf_file, sample_id, gene_coordinates, samp_mut_index, mut_samp_index):
    # rec.samples: dict[str, pysam.libcbcf.VariantRecordSample]
    # VariantRecordSample: dict[str, Any]
        # for example [('GT', (0, 1)), ('GQ', 100), ('AD', (14671, 462)), ('VF', 0.030500000342726707), ('NL', 20), ('SB', -100.0), ('GQX', 100)]
    # where:
        # GQX "Minimum of {Genotype quality assuming variant position,Genotype quality assuming non-variant position}"
        # GT "Genotype"
        # GQ "Genotype Quality"
        # AD "Allele Depth"
        # VF "Variant Frequency" -- corresponds to the "coverage allele-fraction threshold"
        # NL "Applied BaseCall Noise Level"
        # SB "StrandBias Score"
    # rec.info contains "DP", which should be the depth of the reads at the position
    with VariantFile(vcf_file) as vcf:
        for gene_id, pos in gene_coordinates.items():
            chromosome, start, end = pos
            for rec in vcf.fetch(chromosome, start, end):
                mut_position = f"{chromosome}:{rec.pos}"
                mut_samp_index[mut_position].add(sample_id)
                samp_mut_index[sample_id].append({"sample_id": sample_id,
                                                    "gene": gene_id,
                                                    "position": f"{chromosome}:{rec.pos}",
                                                    "ref": rec.ref,
                                                    "alt": rec.alts,
                                                    "var_freqs": [round(var_rec["VF"], 3) for var_rec in rec.samples.values()]})
    return samp_mut_index, mut_samp_index

def asses_mutations(samp_mut_index, mut_samp_index):
    num_samples = len(samp_mut_index)
    mut_counts = {mut_position: len(mut_samp_index[mut_position]) for mut_position in mut_samp_index}
    for samp_mut_infos in samp_mut_index.values():
        for samp_mut_info in samp_mut_infos:
            mut_position = samp_mut_info["position"]
            if mut_position in mut_counts:
                samp_mut_info["mut_count"] = mut_counts[mut_position]
                samp_mut_info["mut_freq"] = round(mut_counts[mut_position]/num_samples, 3)
    return samp_mut_index

def write_mutations(samp_mut_index, out_file):
    with open(out_file, "w+") as f:
        for sample_id, mut_infos in sorted(samp_mut_index.items(), key=lambda x: x[0]):
            for mut_info in sorted(mut_infos, key=lambda x: (x["sample_id"], x["mut_count"])):
                chromosome, position = mut_info["position"].split(":")
                f.write(f"{mut_info['sample_id']},{mut_info['gene']},{chromosome},{position},{mut_info['ref']},{mut_info['alt']},{mut_info['var_freqs']},{mut_info['mut_count']},{mut_info['mut_freq']}\n")

if __name__ == "__main__":
    try:
        vcf_folder = sys.argv[1]
        gene_list_file = sys.argv[2]
        indexed = bool(int(sys.argv[3]))
    except:
        print("Invalid arguments.")
        print("Usage: python3 find_mutations.py <vcf_files_folder> <gene_list_file> <already_indexed>")
        print("Example: python3 find_mutations.py /path/to/vcf/files/folder /path/to/gene_list_file.csv 0")
        print("Make sure you are providing correct arguments and try again.")
        exit()
    
    gene_list_file_path = os.path.normpath(os.path.join(os.getcwd(), gene_list_file))
    vcf_folder_path = os.path.normpath(os.path.join(os.getcwd(), vcf_folder))

    gene_coordinates = read_gene_coordinates(gene_list_file_path)

    if not indexed:
        for file in os.listdir(vcf_folder_path):
            if not file.endswith(".vcf"):
                continue
            print(f"Indexing {file}")
            file_path = os.path.join(vcf_folder_path, file)
            # Make a copy of original file, tabix_index would destroy it
            tmp_file_path = os.path.join(vcf_folder_path, file+".tmp")
            shutil.copy(file_path, tmp_file_path)
            tabix_index(file_path, preset="vcf", force=True)
            # Rename the original file back
            os.rename(tmp_file_path, file_path)

    sample_mutation_index, mutation_sample_index = defaultdict(list), defaultdict(set)
    for file in os.listdir(vcf_folder_path):
        if not file.endswith(".vcf.gz"):
            continue
        print(f"Processing {file}")
        file_path = os.path.join(vcf_folder_path, file)
        try:
            sample_id = int(file.split("_")[0])
        except:
            print(f"Can't process {file}, infered sample_id is not an integer. Skipping...")
            continue
        if sample_id in sample_mutation_index and sample_mutation_index[sample_id]:
            print(f"Duplicate sample_id found: {sample_id}. Skipping...")
            continue
        sample_mutation_index, mutation_sample_index = agregate_gene_mutations(file_path, sample_id, gene_coordinates, sample_mutation_index, mutation_sample_index)
        
    print(sample_mutation_index)
    print(mutation_sample_index)
    
    sample_mutation_index = asses_mutations(sample_mutation_index, mutation_sample_index)
    print(sample_mutation_index)

    out_name = f"{str(date.today())}-{os.path.basename(vcf_folder)}-{os.path.basename(gene_list_file)}-mutations.csv"
    write_mutations(sample_mutation_index, out_name)
    print(f"Output written to {out_name}")
    
    # vcf_file = "1_S1.vcf"
    # chromosome = "chr14"
    # pos = (81609282,81610974)
    # pos = (pos[0]-1, pos[1])




