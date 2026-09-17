#!/bin/bash

#SBATCH --export=none
#SBATCH --get-user-env
#SBATCH --clusters=***
#SBATCH --partition=***
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1Gb
#SBATCH --time=24:00:00
#SBATCH --mail-user=***
#SBATCH --job-name=pop
#SBATCH --output=slurm_logs/%x-%j.out

echo "///////////////////////////////////////////"; 
echo "GET FILTERED VARIANTS FOR Y CHROMOSOME";
echo "$(date '+%Y-%m-%d %H:%M:%S')";
echo "///////////////////////////////////////////";
echo "";
source ~/.bashrc;
conda activate tools;

######################################
##     INPUT & OUTPUT PARAMETERS    ##
######################################

DIR="$(dirname "${VCF}")"
FILE="$(basename "${VCF}")"
FILE="${FILE%.vcf.gz}"
mkdir -p ${DIR}/RENAMED/;
mkdir -p ${DIR}/QC/;
mkdir -p ${DIR}/PCA/;
mkdir -p ${DIR}/ADMIX/;


RENAMED_FILE=${DIR}/RENAMED/${FILE}.rn.vcf.gz
bcftools annotate --rename-chrs $CHRMAP -Oz -o $RENAMED_FILE $VCF

conda activate pop;
tabix -p vcf $RENAMED_FILE
PLINK_FILE=${DIR}/RENAMED/${FILE}.rn
plink --vcf $RENAMED_FILE --double-id --keep-allele-order --set-missing-var-ids @:# --make-bed --out $PLINK_FILE

PLINK_FILE_QC=${DIR}/QC/${FILE}.rn.qc
plink --bfile $PLINK_FILE --mind 0.1 --geno 0.05 --maf 0.01 --make-bed --out $PLINK_FILE_QC

# PCA
PLINK_FILE_PCAld=${DIR}/PCA/${FILE}.rn.pcald
PLINK_FILE_PCAin=${DIR}/PCA/${FILE}.rn.pcain
PLINK_FILE_PCA=${DIR}/PCA/${FILE}.rn.PCA
plink --bfile $PLINK_FILE_QC --indep-pairwise 200 50 0.1 --out $PLINK_FILE_PCAld
plink --bfile $PLINK_FILE_QC --extract $PLINK_FILE_PCAld.prune.in --make-bed --out $PLINK_FILE_PCAin
plink2 --bfile $PLINK_FILE_PCAin --freq --out $PLINK_FILE_PCAin
plink --bfile $PLINK_FILE_PCAin --recode A --out $PLINK_FILE_PCAin.raw

# ADMIXTURE
PLINK_FILE_ADMIXld=${DIR}/ADMIX/${FILE}.rn.admixld
PLINK_FILE_ADMIXin=${DIR}/ADMIX/${FILE}.rn.admixin
plink --bfile $PLINK_FILE_QC --indep-pairwise 50 5 0.2 --out $PLINK_FILE_ADMIXld
plink --bfile $PLINK_FILE_QC --extract ${PLINK_FILE_ADMIXld}.prune.in --make-bed --out $PLINK_FILE_ADMIXin
cd ${DIR}/ADMIX/
for i in $(seq 2 5);
do
    admixture --bootstrap=1000 --cv=10 ${PLINK_FILE_ADMIXin}.bed $i
done
