# dna-mutations-search

## Intro
I have created this repo when I found out that my friend did a tedious job of looking through files and writing some numbers out of them into a table, when it could be easily automated. **Be sure to reach out (for example in Issues) if you have any problem or request! I will happily keep working on this repo.**

If you work with results of DNA sequencing in `.bam`, `.bam.bai` and `.vcf` files, you could find it useful. Currently it contains two scripts. Also, I found `pysam` documentation quite incomplete in parts, so if you struggle with it as well, you might find inspiration in my code.

## Installation
```bash
git clone https://github.com/b0r3k/dna-mutations-search.git
cd dna-mutations-search
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Usage

### `read_coordinates_counts.py`
This script reads `.bam` and `.bam.bai` files and counts the number of reads that cover given position in the genome. It expects that you have a folder with all your `.bam` and `.bam.bai` files named as `sampleid_*.bam` and `sampleid_*.bam.bai` respectively. The `sampleid` is used to identify the sample in the output table.

```bash
python3 read_coordinates_counts.py <bam_files_folder> <chromosome> <position>
```

Where `bam_files_folder` is a path to your folder (as `/data/bams`), `chromosome` is a string (as `chr17`) and `position` is an integer of the position as you would see it in some graphical program (as `40857116`).

The output is a `.csv` file, where there is a line for each sample from your folder and the number of reads that cover the position specified, i.e. the structure of the line is `sampleid;reads_count`. Semicolon (`;`) is used as a separator.

### `find_mutations.py`
This script reads `.vcf` files and finds mutations on given genes (regions) in the genome. For each mutation, it counts how many samples have it. Furthermore it uses the function from `read_coordinates_counts.py` to count the number of reads that cover the position of the mutation. From this you should be able to easily see, which mutations are really interesting mutations (only a few people from your sample have it) and which mutations are just common mutations (many people from your sample have it). Also you are able to see how many reads found the mutation for that sample and how many reads covered the position of the mutation; which should allow you to decide if it is a real mutation or just some measurement error.

```bash
python3 find_mutations.py <vcf_files_folder> <gene_list_file> <already_indexed>
```

Where `vcf_files_folder` is a path to your folder (as `/data/vcfs`), `gene_list_file` is a path to a file with a list of genes (regions) you are interested in (as `gene_list.csv`), and `already_indexed` specifies if you have also indexed versions of your files (`.vcf.gz` and `.vcf.gz.tbi`) –– if you set it to `0`, the files are also indexed, if to `1`, it is expected that the indexed files are already present. The `gene_list_file` should contain a list of genes (regions) in the genome you are interested in. The format of each line is `gene,chromosome,start_position,end_position`, where `gene` is any identifier you choose for the gene, `chromosome` is a string (as `chr17`), `start_position` and `end_position` are integers of the start and end position of the gene in the genome as you would see it in some graphical program (as `40857116` and `40857987`).

The output is again semicolon (`;`) separated `.csv` file, where there is a line for every mutation found in each sample in given regions. The format of each line is following:
```
sample_id;gene;chromosome;position;reference;alternatives;variant_frequencies;read_count;mutation_freqency;mutation_count
```
Where `reference` and `alternatives` are strings of the reference and alternative alleles, `variant_frequencies` is a string of frequencies of the alternative alleles (in what ratio of reads was the mutation confirmed), `read_count` is the number of reads that cover the position of the mutation (how many reads were on that position overall), `mutation_freqency` is the frequency of the mutation across the samples, and `mutation_count` is the number of samples that have given mutation (how many samples exist such that at least 1 read found the mutation for that sample).