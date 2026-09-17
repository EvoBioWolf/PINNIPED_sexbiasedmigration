#!/bin/bash

#SBATCH --export=none
#SBATCH --get-user-env
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_highmem
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16Gb
#SBATCH --time=5-00:00:00
#SBATCH --mail-user=yakupova@bio.lmu.de
#SBATCH --job-name=alg
#SBATCH --output=_slurm_logs/%x-%j.out


echo "///////////////////////////////////////////"; 
echo "ALIGNMENT"
echo "`date "+%Y-%m-%d %H:%M:%S"`";
source ~/.bashrc;
echo "conda activate tools;"
conda activate tools;
echo "///////////////////////////////////////////";

logpath=_logs/;
log=${logpath}/${log_file_suffix};
mkdir -p $logpath;

threads=$SLURM_CPUS_PER_TASK
SAMPLE=$sample

P_TRM=$trimmedpath
P_ALG=$outpath/2.align
P_MG=$outpath/3.merge
P_MG_REAG=$outpath/4.realign_gatk
P_MG_REAG_MD=$outpath/5.realign_gatk_md
P_MG_REAG_MD_F3844=$outpath/6.realign_gatk_md_F3844
P_ASSEMBLY=$(dirname $fasta)

fasta_file="${fasta##*/}"

mkdir -p ${P_ALG}/${SAMPLE};
mkdir -p ${P_MG}/${SAMPLE};
mkdir -p ${P_MG_REAG}/${SAMPLE};
mkdir -p ${P_MG_REAG_MD}/${SAMPLE};
mkdir -p ${P_MG_REAG_MD_F3844}/${SAMPLE};

echo "$(date '+%Y-%m-%d %H:%M:%S')" > $log
echo "" >> $log

echo "NODE: $(hostname)" >> $log
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}" >> $log
echo "" >> $log

echo "INPUT:" >> $log
echo "  FASTA:               ${fasta}" >> $log # specified in INPUT
echo "  READS:               ${P_TRM}" >> $log # specified in INPUT

echo "" >> $log

echo "OUTPUT:" >> $log
echo "  ALIGNMENTS:          ${P_ALG}" >> $log
echo "  MERGED ALIGNMENTS:   ${P_MG}" >> $log
echo "   + REALIGNMENT GATK: ${P_MG_REAG}" >> $log
echo "    + MARKDUP:         ${P_MG_REAG_MD}" >> $log

echo "" >> $log

echo "TR=trimmed" >> $log
echo "MG=merged" >> $log
echo "REAG=realigned near indels with GATK" >> $log
echo "MD=markdup" >> $log
echo "SRT=sorted" >> $log

echo "" >> $log

echo "SAMPLE: $SAMPLE" | tee -a $log;
runs=$(ls "${P_TRM}/${SAMPLE}" | grep -E '_(R)?[12](_[^.]*)?_val_' | sed -E 's/_(R)?[12](_[^.]*)?_val_.*$//' | sort -u);
for run in $runs; do
    echo "" | tee -a $log;
    echo "RUN: $run" | tee -a $log;

    FQF=${P_TRM}/${SAMPLE}/${run}*1.fq.gz;
    FQR=${P_TRM}/${SAMPLE}/${run}*2.fq.gz;

    # Align
    BAM=${P_ALG}/${SAMPLE}/${run}.TR.SRT.bam;
    if ! test -f $BAM; then
        echo " ";
        echo "----------------------bwa mem----------------------"
        echo "`date "+%H:%M:%S"`";
        echo "bwa mem  $fasta $FQF $FQR -M -R "@RG\tID:${SAMPLE}\tPU:x\tSM:${SAMPLE}\tPL:NovaSeqXplus\tLB:x" | samtools fixmate -@ $threads -m - - | samtools sort -@ $threads -m 8G -o $BAM;" | tee -a $log;
        bwa mem  $fasta $FQF $FQR -M -R "@RG\tID:${SAMPLE}\tPU:x\tSM:${SAMPLE}\tPL:NovaSeqXplus\tLB:x" | \
        samtools fixmate -@ $threads -m - - | \
        samtools sort -@ $threads -m 8G -o $BAM;
    fi
    
    # Index bam
    if ! test -f ${BAM}.bai; then
        echo " ";
        echo "----------------------samtools index----------------------"
        echo "`date "+%H:%M:%S"`";
        echo "samtools index -@ $threads $BAM" | tee -a $log;
        samtools index -@ $threads $BAM;
    fi

done

# Merge all bams of one sample
BAM_MG=${P_MG}/${SAMPLE}/${SAMPLE}.TR.MG.SRT.bam;     
if ! test -f $BAM_MG; then
    echo "`date "+%H:%M:%S"`" | tee -a $log;
    echo "picard MergeSamFiles \
        $(printf 'I=%s ' "${P_ALG}/${SAMPLE}"/*.bam) \
        O=$BAM_MG;" | tee -a $log;
    picard MergeSamFiles \
        $(printf 'I=%s ' "${P_ALG}/${SAMPLE}"/*.bam) \
        O=$BAM_MG;
fi;

# Index bam
if ! test -f ${BAM_MG}.bai; then
    echo " ";
    echo "----------------------samtools index----------------------"
    echo "`date "+%H:%M:%S"`";
    echo "samtools index -@ $threads $BAM_MG" | tee -a $log;
    samtools index -@ $threads $BAM_MG;
fi

# Realign indels
BAM_MG_REAG=${P_MG_REAG}/${SAMPLE}/${SAMPLE}.TR.MG.REAG.bam;
if ! test -f $BAM_MG_REAG; then
    echo "`date "+%H:%M:%S"`" | tee -a $log;
    module load /lrz/sys/share/modules/extfiles/singularity/3.8.3;
    # build gatk3 singularity image and determine path to .jar:
    # singularity build gatk3.sif docker://broadinstitute/gatk3:3.8-1;
    # singularity shell gatk3.sif;
    # find / -name "GenomeAnalysisTK.jar" 2>/dev/null --> /usr/GenomeAnalysisTK.jar
    echo "singularity exec \
        -B ${P_ASSEMBLY}:/assembly \
        -B ${P_MG}:/bam \
        -B ${P_MG_REAG}:/reabam \
        /dss/dsshome1/08/ra63rof/tools/gatk3.8.1.0/gatk3.sif \
        java -jar /usr/GenomeAnalysisTK.jar \
            -T RealignerTargetCreator \
            -nt $threads \
            -R /assembly/${fasta_file} \
            -I /bam/${SAMPLE}/${SAMPLE}.TR.MG.SRT.bam \
            -o /reabam/${SAMPLE}/${SAMPLE}.TR.MG.targets.list;" | tee -a $log;
    
    singularity exec \
        -B ${P_ASSEMBLY}:/assembly \
        -B ${P_MG}:/bam \
        -B ${P_MG_REAG}:/reabam \
        /dss/dsshome1/08/ra63rof/tools/gatk3.8.1.0/gatk3.sif \
        java -jar /usr/GenomeAnalysisTK.jar \
            -T RealignerTargetCreator \
            -nt $threads \
            -R /assembly/${fasta_file} \
            -I /bam/${SAMPLE}/${SAMPLE}.TR.MG.SRT.bam \
            -o /reabam/${SAMPLE}/${SAMPLE}.TR.MG.targets.list;

    echo "`date "+%H:%M:%S"`" | tee -a $log;
    echo "singularity exec \
        -B ${P_ASSEMBLY}:/assembly \
        -B ${P_MG}:/bam \
        -B ${P_MG_REAG}:/reabam \
        /dss/dsshome1/08/ra63rof/tools/gatk3.8.1.0/gatk3.sif \
        java -jar /usr/GenomeAnalysisTK.jar \
            -T IndelRealigner \
            -R /assembly/${fasta_file} \
            -I /bam/${SAMPLE}/${SAMPLE}.TR.MG.SRT.bam \
            -targetIntervals /reabam/${SAMPLE}/${SAMPLE}.TR.MG.targets.list \
            -o /reabam/${SAMPLE}/${SAMPLE}.TR.MG.REAG.bam;" | tee -a $log;
    
    singularity exec \
        -B ${P_ASSEMBLY}:/assembly \
        -B ${P_MG}:/bam \
        -B ${P_MG_REAG}:/reabam \
        /dss/dsshome1/08/ra63rof/tools/gatk3.8.1.0/gatk3.sif \
        java -jar /usr/GenomeAnalysisTK.jar \
            -T IndelRealigner \
            -R /assembly/${fasta_file} \
            -I /bam/${SAMPLE}/${SAMPLE}.TR.MG.SRT.bam \
            -targetIntervals /reabam/${SAMPLE}/${SAMPLE}.TR.MG.targets.list \
            -o /reabam/${SAMPLE}/${SAMPLE}.TR.MG.REAG.bam;

    samtools index $BAM_MG_REAG;
fi

# Mark duplicates
BAM_MG_REAG_MD=${P_MG_REAG_MD}/${SAMPLE}/${SAMPLE}.TR.MG.REAG.MD.bam;
metrix=${P_MG_REAG_MD}/${SAMPLE}/${SAMPLE}.TR.MG.REAG.MD.metrics.txt;
if ! test -f $BAM_MG_REAG_MD; then
    echo "`date "+%H:%M:%S"`" | tee -a $log;
    echo "picard MarkDuplicates \
        I=$BAM_MG_REAG \
        O=$BAM_MG_REAG_MD \
        M=$metrix \
        CREATE_INDEX=true \
        VALIDATION_STRINGENCY=SILENT;" | tee -a $log;
    picard MarkDuplicates \
        I=$BAM_MG_REAG \
        O=$BAM_MG_REAG_MD \
        M=$metrix \
        CREATE_INDEX=true \
        VALIDATION_STRINGENCY=SILENT;
fi


BAM_MG_REAG_MD_F3844=${P_MG_REAG_MD_F3844}/${SAMPLE}/${SAMPLE}.TR.MG.REAG.MD.F3844.bam;
if ! test -f ${BAM_MG_REAG_MD_F3844}.bai; then
    echo "`date "+%H:%M:%S"`" | tee -a $log;
    echo "samtools view --threads $threads -h -b -F 3844 -q 20 $BAM_MG_REAG_MD > $BAM_MG_REAG_MD_F3844; samtools index $BAM_MG_REAG_MD_F3844;" | tee -a $log;
    samtools view --threads $threads -h -b -F 3844 -q 20 $BAM_MG_REAG_MD > $BAM_MG_REAG_MD_F3844;
    samtools index $BAM_MG_REAG_MD_F3844;
fi
    
echo "`date "+%H:%M:%S"`" | tee -a $log;
echo "END" | tee -a $log;