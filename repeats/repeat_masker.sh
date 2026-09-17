#!/bin/bash

#SBATCH --export=none
#SBATCH --get-user-env
#SBATCH --clusters=***
#SBATCH --partition=***
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1Gb
#SBATCH --time=12:00:00
#SBATCH --mail-user=***
#SBATCH --job-name=repmask
#SBATCH --output=slurm_logs/%x-%j.out

echo "///////////////////////////////////////////"; 
echo "GET REPEAT GFFs";
echo "`date "+%Y-%m-%d %H:%M:%S"`";
source ~/.bashrc;
echo "conda activate rm;";
conda activate rm;
echo "///////////////////////////////////////////";



######################################
## INPUT & OUTPUT PATHS & FILENAMES ##
######################################
PROJECT_PATH="";
RMPATH=${PROJECT_PATH}/Assembly/Repeats/RM;
REPEATS_OUTPUT=${RMPATH}/Repeats_gff;

FASTA=${PROJECT_PATH}/Assembly/Y/Y.fa;
FASTA_PREFIX=$(basename "${FASTA%.*}");



######################################
##           OUTPUT HEADER          ##
######################################
echo "$(date '+%Y-%m-%d %H:%M:%S')";
echo "";

echo "NODE: $(hostname)";
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}";
echo "";

echo "FASTA: ${FASTA}";

echo "";



######################################
##             WORKFLOW             ##
######################################

declare -A LIBS OUTPUTS PREFIXES
LIBS=(
    ["CANIFORMIA"]=$RMPATH/Dfam/Caniformia/dfam_caniformia.fa
    ["OTARIIDAE"]=$RMPATH/Dfam/Otariidae/dfam_otariidae.fa
    ["PHOCIDAE"]=$RMPATH/Dfam/Phocidae/dfam_phocidae.fa
    ["ZALCAL"]=$RMPATH/ZalCalRMDB/ZalCalRMDB-families.fa
)
OUTPUTS=(
    ["CANIFORMIA"]=$REPEATS_OUTPUT/Caniformia
    ["OTARIIDAE"]=$REPEATS_OUTPUT/Otariidae
    ["PHOCIDAE"]=$REPEATS_OUTPUT/Phocidae
    ["ZALCAL"]=$REPEATS_OUTPUT/ZalCal
)
PREFIXES=(
    ["CANIFORMIA"]=${REPEATS_OUTPUT}/Caniformia/${FASTA_PREFIX}.RM_Caniformia_sm
    ["OTARIIDAE"]=${REPEATS_OUTPUT}/Otariidae/${FASTA_PREFIX}.RM_Otariidae_sm
    ["PHOCIDAE"]=${REPEATS_OUTPUT}/Phocidae/${FASTA_PREFIX}.RM_Phocidae_sm
    ["ZALCAL"]=${REPEATS_OUTPUT}/ZalCal/${FASTA_PREFIX}.RM_ZalCal_sm
)

mkdir -p $REPEATS_OUTPUT ${OUTPUTS["CANIFORMIA"]} ${OUTPUTS["OTARIIDAE"]} ${OUTPUTS["PHOCIDAE"]} ${OUTPUTS["ZALCAL"]}

for SPECIES in CANIFORMIA OTARIIDAE PHOCIDAE ZALCAL; do
    LIB=${LIBS[$SPECIES]}
    OUTPUT=${OUTPUTS[$SPECIES]}
    PREFIX=${PREFIXES[$SPECIES]}

    if ! test -f ${OUTPUT}/${FASTA_PREFIX}.fa.masked; then
        echo "";
        echo "----------------------REPEATMASKER: $SPECIES----------------------";
        echo "";
        echo "`date "+%H:%M:%S"`";
        
        RepeatMasker -lib $LIB $FASTA -dir $OUTPUT -gff -xsmall -parallel $((${SLURM_CPUS_PER_TASK}/4));
        mv ${OUTPUT}/${FASTA_PREFIX}.fa.out.gff ${PREFIX}.gff;
        awk 'BEGIN{OFS="\t"} !/^#/ {if ($4 > $5) print $1, $5-1, $4; else print $1, $4-1, $5}' ${PREFIX}.gff > ${PREFIX}.bed;
        rm ${OUTPUT}/${FASTA_PREFIX}.fa.*;
    fi
done



echo " "
echo "-----------------------END-----------------------";
echo "`date "+%Y-%m-%d %H:%M:%S"`";
