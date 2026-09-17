#!/bin/bash

#SBATCH --export=none
#SBATCH --get-user-env
#SBATCH --clusters=***
#SBATCH --partition=***
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1Gb
#SBATCH --time=01:00:00
#SBATCH --mail-user=***
#SBATCH --job-name=filter_bams
#SBATCH --output=slurm_logs/%x-%j.out

echo "///////////////////////////////////////////"; 
echo "FILTER BAMS";
echo "$(date '+%Y-%m-%d %H:%M:%S')";
echo "///////////////////////////////////////////";
echo "";
source ~/.bashrc;
conda activate tools;

basepath=""
CHROMOSOME=$1; # CM019810.2=9chr CM019819.2=X CM019820.2=Y CM019821.1=MT
SPECIES=$2; # CSL GSL SSL
SEX=$3; # MF M
BAMLIST="${basepath}/lists/bamlists/bamlist_${SPECIES}_${SEX}.txt";
FBAMDIR="${basepath}/VARCALL/bam/bam_${CHROMOSOME}_${SPECIES}";


######################################
##           OUTPUT HEADER          ##
######################################

echo "INPUTS:";
echo "    CHR=${CHROMOSOME}";
echo "    BAMLIST:                     ${BAMLIST}";
echo "";
echo "OUTPUTS:";
echo "  - DIR FOR FILTERED BAM FILES:  ${FBAMDIR}"
echo "";

######################################
##             WORKFLOW             ##
######################################


mkdir -p $FBAMDIR;

while IFS=$'\t' read -r SAMPLE BAM; do
    FBAM="${SAMPLE}.REA.MD.${CHROMOSOME}.Q20.F3844.bam";
    if ! test -f $FBAMDIR/$FBAM; then
        echo "samtools view --threads 4 -h -b -F 3844 -q 20 $BAM $CHROMOSOME > $FBAMDIR/$FBAM";
        samtools view --threads 4 -h -b -F 3844 -q 20 $BAM $CHROMOSOME > $FBAMDIR/$FBAM;

        echo "samtools index $FBAMDIR/$FBAM";
        samtools index $FBAMDIR/$FBAM;
    fi
done < $BAMLIST;


echo " ";
echo "-----------------------END-----------------------";
echo "$(date '+%Y-%m-%d %H:%M:%S')";


