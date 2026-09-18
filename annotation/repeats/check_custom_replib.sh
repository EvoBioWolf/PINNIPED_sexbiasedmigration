#!/bin/bash

#SBATCH --export=none
#SBATCH --get-user-env
#SBATCH --clusters=***
#SBATCH --partition=***
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1Gb
#SBATCH --time=48:00:00
#SBATCH --mail-user=***
#SBATCH --job-name=check_repeats
#SBATCH --output=slurm_logs/%x-%j.out

echo "///////////////////////////////////////////"; 
echo "GET FILTERED VARIANTS FOR Y CHROMOSOME";
echo "$(date '+%Y-%m-%d %H:%M:%S')";
echo "///////////////////////////////////////////";
echo "";
source ~/.bashrc;
conda activate rm;

basepath=""

ZalCalRMDB_fasta=${basepath}/Assembly/Repeats/RM/ZalCalRMDB/ZalCalRMDB-families.fa
dedup_fasta=${ZalCalRMDB_fasta}.cd98
uniprot_fasta=${basepath}/Assembly/Repeats/Uniprot/uniprot_sprot.fasta
prot_hits=${basepath}/Assembly/Repeats/RM/ZalCalRMDB/ZalCalRMDB-families.fa.cd98.blastx_sprot_hits.tsv
prot_hits_85=${basepath}/Assembly/Repeats/RM/ZalCalRMDB/ZalCalRMDB-families.fa.cd98.blastx_sprot_hits_pid85.tsv
gene_ids=${basepath}/Assembly/Repeats/RM/ZalCalRMDB/ZalCalRMDB-families.fa.cd98.blastx_sprot_hits_pid85.likely_gene_ids.txt
final_fasta=${ZalCalRMDB_fasta}.truereps

# Deduplicate sequences
cd-hit-est -i $ZalCalRMDB_fasta -o $dedup_fasta -c 0.98 -n 10 -aS 0.9 -T 8 -M 0

# Make Swiss-Prot BLAST database
conda deactivate;
conda activate tools;
cd ...../tools/blast/ncbi-blast-2.15.0+/bin/;

./makeblastdb -in $uniprot_fasta -dbtype prot

# BLASTX: search for protein-coding contaminants
./blastx -query $dedup_fasta -db $uniprot_fasta -evalue 1e-5 -seg yes \
       -outfmt '6 qseqid sseqid pident length evalue bitscore' \
       -max_target_seqs 10 -num_threads 8 > $prot_hits
awk '$3 > 85' $prot_hits > $prot_hits_85
cut -f1 $prot_hits_85 | sort | uniq > $gene_ids


# Remove genes-like sequences from repeats database
conda activate rm;
# seqtk subseq -l 80 $dedup_fasta $gene_ids | seqtk seq -A > $genelike_fasta
seqtk seq -A $dedup_fasta | seqkit grep -v -f <(cut -f1 $gene_ids) > $final_fasta

# RepeatClassifier -consensi $ZalCalRMDB_fasta
