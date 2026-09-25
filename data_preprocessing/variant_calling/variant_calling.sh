#!/bin/bash

#SBATCH --export=none
#SBATCH --get-user-env
#SBATCH --clusters=***
#SBATCH --partition=***
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1Gb
#SBATCH --time=24:00:00
#SBATCH --mail-user=***
#SBATCH --job-name=vc
#SBATCH --output=slurm_logs/%x-%j.out

echo "///////////////////////////////////////////"; 
if [[ $VARIANTS == "v" ]]; then
    echo "SNP CALLER";
else
    echo "SEGREGATING AND INVARIANT SITES CALLER";
fi
echo "$(date '+%Y-%m-%d %H:%M:%S')";
echo "///////////////////////////////////////////";
echo "";
source ~/.bashrc;
conda activate tools;

######################################
##     INPUT & OUTPUT PARAMETERS    ##
######################################
PROJECT_PATH=""

REP_FEMHOM_BED="${PROJECT_PATH}/lists/filtered_sites/Ychromosome__merged.rep_Otariidae_ZalCal.fem90.bed";

MPILEUP_ADDITIONAL_PARAMS="-q $MQ -Q $BQ --annotate AD,ADF,ADR,DP,SP,INFO/AD,INFO/ADF,INFO/ADR"
MPILEUP_ADDITIONAL_SUFFIX="MQ${MQ}.BQ${BQ}"

VCF_OUTPUT_DIR="${PROJECT_PATH}/VARCALL/${VCF_DIR}${VARIANTS}/${SPECIES}_${CHR}";
RVCF=${VCF_OUTPUT_DIR}/${SPECIES}_${CHR}.${MPILEUP_ADDITIONAL_SUFFIX}.vcf.gz;
SVCF=${VCF_OUTPUT_DIR}/${SPECIES}_${CHR}.${MPILEUP_ADDITIONAL_SUFFIX}.snps.vcf.gz;

######################################
##           OUTPUT HEADER          ##
######################################


echo "NODE: $(hostname)";
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}";
echo "";

echo "INPUTS:";
echo "  - FASTA:                  ${FASTA}";
echo "  - BAMLIST:                ${BAMLIST}";
echo "  - REPEAT BED:             ${REP_FEMHOM_BED}";
echo "    - CHROMOSOME:           ${CHR}";
echo "    - SPECIES:              ${SPECIES}";
echo "";
echo "OUTPUTS:";
echo "  - VCF WITH FILTERED SNPS: ${FVCF}";
echo "";

######################################
##             WORKFLOW             ##
######################################
mkdir -p $VCF_OUTPUT_DIR;

# CALL VARIANTS
if ! test -f $RVCF; then
    echo "bcftools mpileup --threads $SLURM_CPUS_PER_TASK -Ou -r $CHR -f $FASTA -b $BAMLIST $MPILEUP_ADDITIONAL_PARAMS | \
    bcftools call --threads $SLURM_CPUS_PER_TASK --ploidy-file $PLOIDY_FILE --samples-file $SAMPLE_FILE -m${VARIANTS} -Oz -o $RVCF";

    bcftools mpileup --threads $SLURM_CPUS_PER_TASK -Ou -r $CHR -f $FASTA -b $BAMLIST $MPILEUP_ADDITIONAL_PARAMS | \
    bcftools call --threads $SLURM_CPUS_PER_TASK --ploidy-file $PLOIDY_FILE --samples-file $SAMPLE_FILE -m${VARIANTS} -Oz -o $RVCF;
fi

# EXTRACT SNPs
if [[ $VARIANTS == "v" ]]; then
    if ! test -f $SVCF; then
        echo "bcftools view --types snps $RVCF -o $SVCF";
        bcftools view --types snps $RVCF -o $SVCF;
    fi
fi

# SUBTRACT REPEATS AND HOMOLOGOUS REGIONS
# GO TO exclude_regions_from_vcf.sh INSTEAD
if ! test -f $FVCF; then
    echo "(grep '^#' $SVCF; bedtools subtract -a $SVCF -b $REP_FEMHOM_BED | grep -v '^#') > $FVCF";
    (grep '^#' $SVCF; bedtools subtract -a $SVCF -b $REP_FEMHOM_BED | grep -v '^#') > $FVCF;
fi

echo " ";
echo "-----------------------END-----------------------";
echo "$(date '+%Y-%m-%d %H:%M:%S')";


