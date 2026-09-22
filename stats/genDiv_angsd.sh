#!/bin/bash

#SBATCH --export=none
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=4Gb
#SBATCH --time=04:00:00
#SBATCH --mail-type END
#SBATCH --mail-type FAIL
#SBATCH --job-name=test_angsd_pi
#SBATCH --output=logs/%x-%j.out

echo "///////////////////////////////////////////"
echo "pi and thetaW in angsd"
echo "$(date '+%Y-%m-%d %H:%M:%S')"
echo "///////////////////////////////////////////"
echo ""
source ~/.bashrc
conda activate myEnv

######################################
##              USAGE               ##
######################################
# $1 = SAMPLING_SUFFIX (ALL, CSL, GSL)
# $2 = bamfiles
# sbatch angsd_varCall_Xchr.sh ALL

######################################
##  SET INPUT & OUTPUT PARAMETERS  ##
######################################
## inputs
# PROJECT_PATH=
FASTA="GCA_009762305.2_mZalCal1.pri.v2_genomic.fna"
angsd="$HOME/angsd/misc"
region_file="$PROJECT_PATH/data/region_file_angsdPi"
bamlist="$PROJECT_PATH/data/bamlist_${1}_females.txt"
### Outputs
OUTPUT_PREF="$PROJECT_PATH/resutls/25_pi/angsd_test/angsd_pi_FEM_ONLY__$1"


######################################
##           OUTPUT HEADER          ##
######################################

echo " angsd pi for $1"

echo "NODE: $(hostname)";
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}";
echo "";

echo "INPUTS:";
echo "  - BAMLIST:                   $bamlist"
echo "  - region file:               $region_file"
echo "";
echo "OUTPUTS:";
echo "  - ANGSD outputs              $OUTPUT_PREF"
echo "";

echo " ";
echo "---------------------START-----------------------";
echo "$(date '+%Y-%m-%d %H:%M:%S')";
echo " ";

######################################
##             WORKFLOW             ##
######################################


# substract repeats from complete chromosome
bedtools subtract -a data/chr_size_CM019810.2.txt -b /dss/dsslegfs01/pr53da/pr53da-dss-0019/projects/2024__SexBias_GeneFlow/Ychromosome/VARCALL/hom_rep_par_bedfiles/CM019810.2_hom_repOtZC.bed > data/angsd_pi_sites_ansgdFilter_CM019810.2.txt

bedtools subtract -a data/chr_size_CM019819.2.txt -b /dss/dsslegfs01/pr53da/pr53da-dss-0019/projects/2024__SexBias_GeneFlow/Ychromosome/VARCALL/hom_rep_par_bedfiles/CM019819.2_hom_repOtZC_par.bed > data/angsd_pi_sites_ansgdFilter_CM019819.2.txt

# get minimum coverage and nind
nind=$(wc -l $bamlist | awk '{print int($1*0.9)}')
minD=$(wc -l $bamlist | awk '{print int($1*2)}')
maxD=$(wc -l $bamlist | awk '{print int($1*15)}')

for chr in CM019810.2 CM019819.2
do
echo "---------------------------------------"
echo "angsd pi for $1 - chromosome $chr"

region_file="data/angsd_pi_sites_ansgdFilter_$chr.txt"

### with min depth filter
OUTPUT_PREF="$PROJECT_PATH/resutls/25_pi/angsd_test/angsd_pi_sites_FEM_ONLY_angsdFilter_MQ40_MD1000_minD2_NA01__${1}_$chr"
echo " "
echo "angsd -bam $bamlist -doSaf 1 -anc $FASTA -minMapQ 40 -GL 2 -nInd $nind -doCounts 1 -setMaxDepth 1000 -setMinDepth $minD -P 8 -out $OUTPUT_PREF -r $chr: -sites $region_file"
angsd -bam $bamlist -doSaf 1 -anc $FASTA -minMapQ 40 -GL 2 -nInd $nind -doCounts 1 -setMaxDepth 1000 -setMinDepth $minD -P 8 -out $OUTPUT_PREF -r $chr: -sites $region_file

echo " "
echo "$angsd/realSFS $OUTPUT_PREF.saf.idx -P 8 -fold 1 > $OUTPUT_PREF.sfs"
$angsd/realSFS $OUTPUT_PREF.saf.idx -P 8 -fold 1 > $OUTPUT_PREF.sfs

echo " "
echo "$angsd/realSFS saf2theta $OUTPUT_PREF.saf.idx -P 8 -outname $OUTPUT_PREF -sfs $OUTPUT_PREF.sfs -fold 1"
$angsd/realSFS saf2theta $OUTPUT_PREF.saf.idx -P 8 -outname $OUTPUT_PREF -sfs $OUTPUT_PREF.sfs -fold 1

echo " "
echo "angsd/thetaStat do_stat $OUTPUT_PREF.thetas.idx -win 50000 -step 10000  -outnames $OUTPUT_PREF"
$angsd/thetaStat do_stat $OUTPUT_PREF.thetas.idx -win 50000 -step 10000  -outnames $OUTPUT_PREF



done




echo " ";
echo "-----------------------END-----------------------";
echo "$(date '+%Y-%m-%d %H:%M:%S')";



