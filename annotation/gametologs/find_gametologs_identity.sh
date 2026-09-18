#!/bin/bash

#SBATCH --export=none
#SBATCH --get-user-env
#SBATCH --clusters=***
#SBATCH --partition=***
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2Gb
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=***
#SBATCH --job-name=fdgam
#SBATCH --output=_slurm_logs/%x-%j.out

hostname; 
echo "`date "+%Y-%m-%d %H:%M:%S"`";
source ~/.bashrc;
conda activate emboss

basepath=""
regions="${basepath}/Annotation/X_Y_regions.tsv"
xref="${basepath}/Assembly/split_by_chr/GCA_009762305.2_mZalCal1.pri.v2_genomic.CM019819.2.fasta.masked"
yref="${basepath}/Assembly/split_by_chr/GCA_009762305.2_mZalCal1.pri.v2_genomic.CM019820.2.fasta.masked"
outdir="${basepath}/Annotation/gametologs"

mkdir -p "$outdir"

# samtools faidx "$xref"
# samtools faidx "$yref"

while IFS=$'\t' read -r -u 3 \
    locidy gene_biotype ys ye strand yregion ygene xgene xs xe; do

    echo "`date "+%H:%M:%S"`";
    printf 'Ygene=%q Xgene=%q ys=%q ye=%q xs=%q xe=%q\n' \
        "$ygene" "$xgene" "$ys" "$ye" "$xs" "$xe"

    prefix="${outdir}/${locidy}_${ygene}_${xgene}"
    xopts=()
    yopts=()

    # if start < end reverse
    if (( xs > xe )); then
        tmp=$xs; xs=$xe; xe=$tmp
        xopts=(-i)
    fi

    if (( ys > ye )); then
        tmp=$ys; ys=$ye; ye=$tmp
        yopts=(-i)
    fi

    samtools faidx "${xopts[@]}" "$xref" "CM019819.2:${xs}-${xe}" > "${prefix}_X.fasta"
    samtools faidx "${yopts[@]}" "$yref" "CM019820.2:${ys}-${ye}" > "${prefix}_Y.fasta"

    stretcher \
        -asequence "${prefix}_X.fasta" \
        -bsequence "${prefix}_Y.fasta" \
        -gapopen 16 \
        -gapextend 4 \
        -outfile "${prefix}_alignment.txt" \
        -auto </dev/null

    # needle \
    #     -asequence "${prefix}_X.fasta" \
    #     -bsequence "${prefix}_Y.fasta" \
    #     -gapopen 10 \
    #     -gapextend 0.5 \
    #     -outfile "${prefix}_alignment.txt" \
    #     -auto

done 3< "$regions"
