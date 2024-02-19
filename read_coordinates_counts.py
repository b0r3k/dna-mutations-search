#!/usr/bin/env python3
import pysam
import pandas as pd

import sys
import os
from datetime import date

def get_read_count(bam_file, chromosome, position):
    """
    Get read count at given position in given bam file.
    
    :param bam_file: str, path to the bam file
    :param chromosome: str, chromosome name
    :param position: int, position in the chromosome

    :return: int, number of reads at the position, -1 if error
    """
    try:
        with pysam.AlignmentFile(bam_file, "rb", index_filename=bam_file+".bai") as bam:
            return bam.count(chromosome, position-1, position)
    except:
        print(f"Error while processing {bam_file}, either it can't be opened or `.bai` index file is missing.")
        return -1

def get_base_counts(bam_file, chromosome, position, min_base_quality=13):
    """
    Efficiently get base counts at given position in given bam file. Beware of the `min_base_quality` parameter, it can significantly affect the results! Position is expected 1-based as the graphical editors use it.
    
    :param bam_file: str, path to the bam file
    :param chromosome: str, chromosome name
    :param position: int, position in the chromosome
    :param min_base_quality: int, minimum base quality to consider the read into results, default is 13
    
    :return: dict[str, int], base as key and count as value
    """
    with pysam.AlignmentFile(bam_file, "rb", index_filename=bam_file+".bai") as bam:
        # Use position-1 to get 0-based position
        for column in bam.pileup(chromosome, position-1, position, truncate=True, stepper="nofilter", min_base_quality=min_base_quality):
            counts = {base : 0 for base in "ATCGND"}
            for read in column.pileups:
                if read.is_del:
                    counts["D"] += 1
                base = read.alignment.seq[read.query_position]
                counts[base] += 1
    return counts

if __name__ == "__main__":
    try:
        bam_folder = sys.argv[1]
        chromosome = sys.argv[2]
        position = int(sys.argv[3])
    except:
        print("Invalid arguments.")
        print("Usage: python3 read_coordinates_counts.py <bam_files_folder> <chromosome> <position>")
        print("Example: python3 read_coordinates_counts.py /path/to/bam/files/folder chr17 40857116")
        print("Make sure you are providing correct arguments and try again.")
        exit()
    
    bam_folder_path = os.path.normpath(os.path.join(os.getcwd(), bam_folder))
    number_of_measurments = dict()

    for file in os.listdir(bam_folder_path):
        if not file.endswith(".bam"):
            continue
        print(f"Processing {file}")
        file_path = os.path.join(bam_folder_path, file)
        try:
            sample_id = int(file.split("_")[0])
        except:
            print(f"Can't process {file}, infered sample_id is not an integer. Skipping...")
            continue
        n = get_read_count(file_path, chromosome, position)
        if sample_id in number_of_measurments and number_of_measurments[sample_id] != -1:
            print(f"Duplicate sample_id found: {sample_id}. Skipping...")
            continue
        number_of_measurments[sample_id] = n    
    res = pd.Series(number_of_measurments).sort_index()
    out_name = f"{str(date.today())}-{chromosome}-{position}-counts.csv"
    res.to_csv(out_name, header=False, sep=";")
    print(f"Output written to {out_name}")
    
    # bam_file = "1_S1.bam"
    # chromosome = "chr17"
    # position = 40857116
    # "chr17:40,857,075-40,857,156"
    # "40 857 116"
    # "expected out: total 1381, G=3, T=1378"
