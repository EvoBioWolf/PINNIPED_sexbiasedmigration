#!/bin/bash

#SBATCH --export=none
#SBATCH --get-user-env
#SBATCH --clusters=***
#SBATCH --partition=***
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2Gb
#SBATCH --time=02:00:00
#SBATCH --mail-user=***
#SBATCH --job-name=cov
#SBATCH --output=_slurm_logs/%x-%j.out

echo "`date "+%Y-%m-%d %H:%M:%S"`";
source ~/.bashrc;
echo "conda activate coverage;"
conda activate coverage;

basepath=""
filelist="${basepath}/bams_list_ssl_HC_test.txt"
outdir="${basepath}/5.realign_gatk_md_coverage"

while IFS= read -r bam; do

    [[ -z "$bam" ]] && continue

    bam_filename="${bam##*/}";
    prefix=$outdir/${bam_filename%.*};

    if ! test -f ${prefix}.mosdepth.summary.txt; then
        echo "mosdepth --threads 8 --mapq 20 --min-frag-len 20 $prefix $bam"
        mosdepth --threads 8 --mapq 20 --min-frag-len 20 $prefix $bam
    fi

done < "$filelist"

echo "END";
echo "`date "+%Y-%m-%d %H:%M:%S"`";