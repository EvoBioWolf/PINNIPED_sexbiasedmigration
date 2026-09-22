#!/bin/bash

#SBATCH --export=none
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=1Gb
#SBATCH --time=00:10:00
#SBATCH --mail-type END
#SBATCH --mail-type FAIL
#SBATCH --output=logs/01_relatedness_%j.out

echo "///////////////////////////////////////////"
echo "relatedness analysis"
echo "$(date '+%Y-%m-%d %H:%M:%S')"
echo "///////////////////////////////////////////"
echo ""
source ~/.bashrc
conda activate myEnv
######################################
##     INPUT & OUTPUT PARAMETERS    ##
######################################
# $1 = vcf file
sp=$(echo $1 | awk -F/ '$0=$NF' | sed 's/_.*$//')
chr=$(echo $1 | sed -e 's/^.*CM019819.2.*$/Xchr/' -e 's/^.*CM019810.2.*$/chr9/')
outpath=resutls/01_relatedness
out=$outpath/01_relatedness_${chr}__$sp

######################################
##           OUTPUT HEADER          ##
######################################
echo "NODE: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}"
echo ""

echo "INPUTS:";
echo "  - VCF:                  $1"
echo ""
echo "OUTPUTS:"
echo "  - plink files:        ${out}.bed/.bim/.fam"
echo "  - relatedness matrix: ${out}.kin0"
echo ""

echo "-----------------------START-----------------------"
echo "$(date '+%Y-%m-%d %H:%M:%S')"
echo " "
######################################
##             WORKFLOW             ##
######################################
mkdir -p $outpath

# convert vcf to plink
echo "plink --vcf $1 --double-id --allow-extra-chr --set-missing-var-ids @:# --keep-allele-order --make-bed --out $out"
plink --vcf $1 --double-id --allow-extra-chr --set-missing-var-ids @:# --keep-allele-order --make-bed --out $out
echo " "

# relatedness matrix
echo "plink --bfile $out --allow-extra-chr --make-king square --out $out.mat"
plink2 --bfile $out --allow-extra-chr --make-king square --out $out.mat
echo " "

# paste sample IDs
echo "awk 'NR>1 {print $1}' $out.mat.king.id | paste -d'\t' - $out.mat.king > $out.mat.king.nm"
awk 'NR>1 {print $1}' $out.mat.king.id | paste -d'\t' - $out.mat.king > $out.mat.king.nm



echo " "
echo "-----------------------END-----------------------"
echo "$(date '+%Y-%m-%d %H:%M:%S')"

