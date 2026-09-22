#!/bin/bash

#SBATCH --export=none
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=12Gb
#SBATCH --time=04:00:00
#SBATCH --mail-type END
#SBATCH --mail-type FAIL
#SBATCH --job-name=02_angsd_pseudo
#SBATCH --output=logs/%x-%j.out

echo "///////////////////////////////////////////"
echo "pseudohaploidisation in angsd"
echo "$(date '+%Y-%m-%d %H:%M:%S')"
echo "///////////////////////////////////////////"
echo ""
source ~/.bashrc
conda activate myEnv

######################################
##              USAGE               ##
######################################
# $1 = SAMPLING_SUFFIX (ALL, CSL, GSL)
# $2 = sites file
# sbatch angsd_varCall_Xchr.sh ALL

######################################
##  SET INPUT & OUTPUT PARAMETERS  ##
######################################
## inputs
PROJECT_PATH="/dss/dsslegfs01/pr53da/pr53da-dss-0019/projects/2024__SexBias_GeneFlow/Autosomes"
FASTA="/dss/dsslegfs01/pr53da/pr53da-dss-0019/assemblies/Zalophus.californianus/genome/VGP.v3/GCA_009762305.2_mZalCal1.pri.v2_genomic.fna"
BAMLIST="${PROJECT_PATH}/data/bamlist_$1.txt"
chr=$(echo $2 | sed -e 's/^.*Xchr.*$/Xchr/' -e 's/^.*chr9.*$/chr9/')
region=$(awk '{print $1}' $2 | uniq)
### Outputs
OUTPUT_PREF="${PROJECT_PATH}/data/angsd/02_angsd_${chr}__$1"


######################################
##           OUTPUT HEADER          ##
######################################


echo "NODE: $(hostname)";
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}";
echo "";

echo "INPUTS:";
echo "  - BAMLIST:                   ${BAMLIST}";
echo "  - sitesfile:               $2"
echo "";
echo "OUTPUTS:";
echo "  - ANGSD outputs              $OUTPUT_PREF";
echo "";

echo " ";
echo "---------------------START-----------------------";
echo "$(date '+%Y-%m-%d %H:%M:%S')";
echo " ";

######################################
##             WORKFLOW             ##
######################################

echo "angsd sites index $2"
angsd sites index $2

echo "angsd -bam $BAMLIST -doMajorMinor 1 -out ${OUTPUT_PREF} \
    -nThreads 4 -doGeno 4 -doMaf 1 -doPost 1 -gl 2 -dohaplocall 2 \
    -doCounts 1 -doGlf 2 -r $region: -sites $2"
    
echo " "
angsd -bam $BAMLIST -doMajorMinor 1 -out ${OUTPUT_PREF} \
    -nThreads 4 -doMaf 1 -doPost 1 -gl 2 -dohaplocall 2 \
    -doCounts 1 -r $region: -sites $2


## convert to plink format
haplo2plink=/dss/dsshome1/08/ra74duv/angsd/misc/haploToPlink

# 1st, haplo to plink input
echo ""
echo "$haplo2plink $OUTPUT_PREF.haplo.gz $OUTPUT_PREF"
$haplo2plink $OUTPUT_PREF.haplo.gz $OUTPUT_PREF


# 2nd, plink to vcf

echo ""
echo "plink --missing-genotype N --tped $OUTPUT_PREF.tped --tfam $OUTPUT_PREF.tfam --threads 2 --allow-extra-chr --recode vcf --out $OUTPUT_PREF"
plink --missing-genotype N --tped $OUTPUT_PREF.tped --tfam $OUTPUT_PREF.tfam --threads 2 --allow-extra-chr --recode vcf --out $OUTPUT_PREF


echo " ";
echo "-----------------------END-----------------------";
echo "$(date '+%Y-%m-%d %H:%M:%S')";
