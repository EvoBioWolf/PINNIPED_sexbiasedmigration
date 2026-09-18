#!/bin/bash

#SBATCH --export=none
#SBATCH --get-user-env
#SBATCH --clusters=***
#SBATCH --partition=***
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=2Gb
#SBATCH --time=12:00:00
#SBATCH --mail-user=***
#SBATCH --job-name=blast
#SBATCH --output=_slurm_logs/%x-%j.out

echo "///////////////////////////////////////////"; 
hostname; 
echo "`date "+%Y-%m-%d %H:%M:%S"`";
source ~/.bashrc;
conda activate tools;
echo "///////////////////////////////////////////";


mkdir -p $outputdir;
basepath=""
ref="${basepath}/Annotation/CM019820.2_xtranssingle.fasta"
query="${basepath}/Assembly/split_by_chr/GCA_009762305.2_mZalCal1.pri.v2_genomic.CM019819.2.fasta.masked"
ref_name="Y_XtransSingle"
query_name="X"
outputdir="${basepath}/Annotation"
echo " ";
echo "-----------------------START-----------------------";
echo "`date "+%H:%M:%S"`";

echo "";
echo "REF="$ref_name;
echo "QUERY="$query_name;
echo "";
echo "REF FASTA="$ref;
echo "QUERY FASTA="$query;
echo "";
echo "OUTPUT="${outputdir};
echo "";


cd ...../tools/blast/ncbi-blast-2.15.0+/bin/;

if ! test -f ${ref}.db.ndb; then
echo " ";
echo "--------------------makeblastdb---------------------";
echo "`date "+%H:%M:%S"`";
echo "makeblastdb -in $ref -dbtype nucl -out ${ref}.db";
./makeblastdb -in $ref -dbtype nucl -out ${ref}.db;
fi


for threshhold in 0 90 95; do
output_name=${outputdir}/homologs.r${ref_name}.q${query_name}.IDT${threshhold}.W20
if ! test -f ${output_name}.csv; then
    echo " ";
    echo "----------------------blastn----------------------";
    echo "`date "+%H:%M:%S"`";
    echo "blastn -query $query -db ${ref}.db -out ${output_name}.csv \
    -outfmt 6 -num_threads $SLURM_CPUS_PER_TASK \
    -perc_identity $threshhold -word_size 20 -evalue 1e-3";
    ./blastn -query $query -db ${ref}.db -out ${output_name}.csv \
    -outfmt 6 -num_threads $SLURM_CPUS_PER_TASK \
    -perc_identity $threshhold -word_size 20 -evalue 1e-3 ;
fi
done;

echo "-----------------------END-----------------------";
echo "`date "+%Y-%m-%d %H:%M:%S"`";
