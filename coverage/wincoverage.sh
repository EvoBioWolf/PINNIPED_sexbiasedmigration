#!/bin/bash

#SBATCH --export=none
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=500Mb
#SBATCH --time=00:10:00
#SBATCH --mail-user=yakupova@bio.lmu.de
#SBATCH --job-name=wincov
#SBATCH --output=slurm_logs/%x-%j.out

echo "///////////////////////////////////////////"; 
hostname; 
echo "`date "+%Y-%m-%d %H:%M:%S"`";
source ~/.bashrc;
echo "conda activate tools";
conda activate tools;
echo "///////////////////////////////////////////";

threads=8;
min_mapping_quality=20;

echo " ";
echo "-----------------------START-----------------------";
echo "`date "+%H:%M:%S"`";
echo "SAMPLE="$SAMPLE;
echo "CHROMOSOME="$CHR;
echo "BAM="$BAM;
echo "OUTPUT="${OUTPUT_PREFIX};

if ! test -f ${OUTPUT_PREFIX}.per-base.bed.gz; then
echo " ";
echo "--------------------mosdepth---------------------";
echo "`date "+%H:%M:%S"`";
echo "mosdepth -t ${threads} --mapq ${min_mapping_quality} ${OUTPUT_PREFIX} ${BAM} ${USE_MEDIAN} ${BY} ${CHR};";
mosdepth -t ${threads} --mapq ${min_mapping_quality} ${OUTPUT_PREFIX} ${BAM} $USE_MEDIAN $BY $CHR;

fi

echo " ";
echo "-----------------------END-----------------------";
echo "`date "+%Y-%m-%d %H:%M:%S"`";