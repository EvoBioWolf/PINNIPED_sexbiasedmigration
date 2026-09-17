#!/bin/bash

#SBATCH --export=none
#SBATCH --get-user-env
#SBATCH --clusters=***
#SBATCH --partition=***
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8Gb
#SBATCH --time=48:00:00
#SBATCH --mail-user=***
#SBATCH --job-name=trm
#SBATCH --output=_slurm_logs/%x-%j.out


echo "///////////////////////////////////////////"; 
echo "TRIMMING"
echo "`date "+%Y-%m-%d %H:%M:%S"`";
source ~/.bashrc;
echo "conda activate tools;"
conda activate tools;
echo "///////////////////////////////////////////";

logpath=_logs/;
log=${logpath}/${log_file_suffix};
mkdir -p $logpath;

threads=$SLURM_CPUS_PER_TASK

echo "$(date '+%Y-%m-%d %H:%M:%S')" > $log
echo "" >> $log

echo "NODE: $(hostname)" >> $log
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}" >> $log
echo "" >> $log

mkdir -p ${trimmedpath}/tmp/${sample};
mkdir -p ${trimmedpath}/${sample};

echo "SAMPLE: $sample" | tee -a $log;

runs=$(ls "${rawpath}/${sample}"  \
| grep -E '_[12]\.fq\.gz$'  \
| sed -E 's/(_1\.fq\.gz|_2\.fq\.gz)$//'  \
| sort -u);
# runs=$(ls -1 "${rawpath}/${sample}" \
# | grep -E '_(R)?[12](_[^.]*)?\.fq\.gz$' \
# | sed -E 's/_(R)?[12](_[^.]*)?\.fq\.gz$//' \
# | sort -u)

for run in $runs; do

    echo "" | tee -a $log;
    echo "RUN: $run" | tee -a $log;
    
    if ! test -f ${trimmedpath}/${sample}/${run}*.fq.gz; then
        echo "trim_galore -j $threads --paired -q 0 --gzip --length 40 --stringency 3 \
                    --fastqc --fastqc_args "-t $threads" \
                    -o ${trimmedpath}/tmp/${sample} \
                    ${rawpath}/${sample}/${run}*.fq.gz" | tee -a $log;
        trim_galore -j $threads --paired -q 0 --gzip --length 40 --stringency 3 \
                    --fastqc --fastqc_args "-t $threads" \
                    -o ${trimmedpath}/tmp/${sample} \
                    ${rawpath}/${sample}/${run}*.fq.gz

        
        echo "trim_galore -j $threads --paired -q 0 --gzip --length 40 --stringency 3 \
                    --fastqc --fastqc_args "-t $threads" \
                    -a A{10} -a2 A{10} \
                    -a G{10} -a2 G{10} \
                    -o ${trimmedpath}/${sample} \
                    ${trimmedpath}/tmp/${sample}/${run}*.fq.gz" | tee -a $log;
        trim_galore -j $threads --paired -q 0 --gzip --length 40 --stringency 3 \
                    --fastqc --fastqc_args "-t $threads" \
                    -a A{10} -a2 A{10} \
                    -a G{10} -a2 G{10} \
                    -o ${trimmedpath}/${sample} \
                    ${trimmedpath}/tmp/${sample}/${run}*.fq.gz

        echo "rm -rf ${trimmedpath}/tmp/${sample};" | tee -a $log;
        rm -rf ${trimmedpath}/tmp/${sample};
    fi
done
