#!/bin/bash

#SBATCH --export=none
#SBATCH --get-user-env
#SBATCH --clusters=***
#SBATCH --partition=***
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2Gb
#SBATCH --time=48:00:00
#SBATCH --mail-user=***
#SBATCH --job-name=repmod
#SBATCH --output=slurm_logs/%x-%j.out

echo "///////////////////////////////////////////"; 
echo "MASK REPEATS: REPEATMODELER + REPEATMASKER";
echo "`date "+%Y-%m-%d %H:%M:%S"`";
source ~/.bashrc;
conda activate rm;
echo "///////////////////////////////////////////";



######################################
## INPUT & OUTPUT PATHS & FILENAMES ##
######################################
PROJECT_PATH="";
OUTPUTPATH=${PROJECT_PATH}/Assembly/Repeats/RM;

FASTA=${PROJECT_PATH}/Assembly/wgs/GCA_009762305.2_mZalCal1.pri.v2_genomic.fa; # Input raw fasta
MASKED_FASTA=${FASTA%.*}.softmasked.fa; # Output soft masked fasta
HARDMASKED_FASTA=${FASTA%.*}.hardmasked.fa; # Output hard masked fasta

DFAM_PATH=${OUTPUTPATH}/Dfam; # Dfam library of repeats
DFAM_DB="${DFAM_PATH}/Caniformia/dfam_caniformia.fa ${DFAM_PATH}/Otariidae/dfam_otariidae.fa ${DFAM_PATH}/Phocidae/dfam_phocidae.fa";

RMDB_NAME=ZalCalRMDB; #Zalophus californianus RepeateModeler database name
RMDB_PATH=$OUTPUTPATH/$RMDB_NAME;
RMDB=${RMDB_PATH}/${RMDB_NAME}-families.fa; # RepeatModeler library of repeats

LIB=${OUTPUTPATH}/rm_library.fa; # Full library of repeats Dfam+RepeatModeler

logpath=${PROJECT_PATH}/Assembly/Repeats/script/_logs/; # Path to the log folder
log_file_suffix="script.log";
log=$logpath/$log_file_suffix;
mkdir -p $logpath;



######################################
##           OUTPUT HEADER          ##
######################################
echo "$(date '+%Y-%m-%d %H:%M:%S')" > $log;
echo "" >> $log;

echo "NODE: $(hostname)" >> $log;
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}" >> $log;
echo "" >> $log;

echo "FASTA: ${FASTA}" >> $log;
echo "SOFT MASKED FASTA: ${MASKED_FASTA}" >> $log;
echo "HARD MASKED FASTA: ${HARDMASKED_FASTA}" >> $log;
echo "" >> $log;



######################################
##             WORKFLOW             ##
######################################
if ! test -f $RMDB_PATH/$RMDB_NAME.nsq; then # .nhr/.nin/.njs/.nnd/.nni/.nog/.nsq/.translation
    echo "" | tee -a $log;
    echo "----------------------REPEATMODELER: BUILD DATABASE----------------------" | tee -a $log;
    echo "";
    echo "`date "+%H:%M:%S"`";
    echo "mkdir -p $RMDB_PATH; cd $RMDB_PATH; BuildDatabase -name $RMDB_NAME $FASTA;" | tee -a $log;
    mkdir -p $RMDB_PATH; cd $RMDB_PATH; BuildDatabase -name $RMDB_NAME $FASTA;
fi

if ! test -f $RMDB; then
    echo "" | tee -a $log;
    echo "----------------------REPEATMODELER: MODEL DATABASE----------------------" | tee -a $log;
    echo "";
    echo "`date "+%H:%M:%S"`";
    echo "cd $RMDB_PATH; RepeatModeler -database $RMDB_NAME -threads $SLURM_CPUS_PER_TASK;" | tee -a $log;
    cd $RMDB_PATH; RepeatModeler -database $RMDB_NAME -threads $SLURM_CPUS_PER_TASK;
fi

if ! test -f $LIB; then
    echo "";
    echo "----------------------MERGE ALL DATABASES----------------------";
    echo "";
    echo "`date "+%H:%M:%S"`";
    echo "cat $DFAM_DB $RMDB > $LIB;" | tee -a $log;
    cat $DFAM_DB $RMDB > $LIB;
fi

if ! test -f $MASKED_FASTA; then
    echo "";
    echo "----------------------REPEATMASKER: SOFT MASK----------------------";
    echo "";
    echo "`date "+%H:%M:%S"`";
    echo "RepeatMasker -lib $LIB $FASTA -dir $OUTPUTPATH -gff -xsmall -parallel $((${SLURM_CPUS_PER_TASK}/4)); mv $OUTPUTPATH/$(basename $FASTA).masked $MASKED_FASTA;" | tee -a $log;
    RepeatMasker -lib $LIB $FASTA -dir $OUTPUTPATH -gff -xsmall -parallel $((${SLURM_CPUS_PER_TASK}/4)); mv $OUTPUTPATH/$(basename $FASTA).masked $MASKED_FASTA;
fi

if ! test -f $HARDMASKED_FASTA; then
    echo "";
    echo "----------------------REPEATMASKER: HARD MASK----------------------";
    echo "";
    echo "`date "+%H:%M:%S"`";
    echo "awk '{ if ($0 ~ /^>/) { print $0 } else { gsub(/[acgt]/, "N", $0); print $0 }}' $MASKED_FASTA > $HARDMASKED_FASTA;" | tee -a $log;
    awk '{ if ($0 ~ /^>/) { print $0 } else { gsub(/[acgt]/, "N", $0); print $0 }}' $MASKED_FASTA > $HARDMASKED_FASTA;
fi


if ! test -f $HARDMASKED_FASTA.fai; then
    echo "";
    echo "----------------------INDEX HARD MASKED FASTA----------------------";
    echo "";
    echo "`date "+%H:%M:%S"`";
    echo "samtools faidx $HARDMASKED_FASTA;" | tee -a $log;
    samtools faidx $HARDMASKED_FASTA;
fi


# Rename log file with unique name and move to the folder with experiment
new_log="${logpath}/repmod.$(date '+%Y-%m-%d.%H:%M').job-${SLURM_JOB_ID}.${log_file_suffix}"
mv $log $new_log

echo " "
echo "-----------------------END-----------------------";
echo "LOGS: $new_log"
echo "`date "+%Y-%m-%d %H:%M:%S"`" | tee -a $new_log;
