#!/bin/bash


#### NOTE: This file contains a collection of scripts used for testing pixy and angsd on a single species and chromosome. It is not meant to be run as a batch job and was not run as is, rather each step was run separately. The purpose of this file is to document the commands used for testing genetic diversity using different methods.

source ~/.bashrc
conda activate pixyEnv
module load bcftools
######################################
##     INPUT & OUTPUT PARAMETERS    ##
######################################
# $1 = invariant vcf
sp=$(echo $1 | awk -F/ '$0=$NF' | sed -e 's|^.*_inv/||' -e 's/_.*$//' )
chr=$(echo $1 | sed -e 's/^.*19820.*$/chrY/' -e 's/^.*19819.*$/chrX/' -e 's/^.*19810.*$/chr9/' -e 's/^.*19821.*$/chrM/')
bedfile="data/GCA_009762305.2_mZalCal1.pri.v2_genomic_${chr}.fa.fai.bed"
K=K3
echo "species: $sp, chr: $chr, K = $K"


######################################
##       pixy tests workflow        ##  
######################################

# Index VCF if not already indexed
if [[ ! -f $1.gz.tbi && $chr != "chrX" ]]
then
    echo "Indexing invariants VCF file with bgzip and tabix"
    echo "bgzip --threads 4 --keep $1"
    echo "tabix --threads 4 -p vcf $1.gz"
    bgzip --threads 4 --keep $1
    tabix --threads 4 -p vcf $1.gz
    echo " "
fi


# have to split chrX by sex
if [[ $chr == "chrX" ]]
then
    females=$(awk '{print $1}' data/rooksList_${sp}_females_$K | tr '\n' ',' | sed 's/,$//')
    out=$(echo $1 | sed 's/\.vcf//')
    if [[ ! -f "${out}_females.vcf.gz.tbi" ]]
    then
    echo "chrX: splitting by sex on invariant file"
    echo "bcftools view --threads 4 -s $females $1 > ${out}_females.vcf"
    echo "bgzip --threads 4 --keep ${out}_females.vcf"
    echo "tabix --threads 4 -f -p vcf ${out}_females.vcf.gz"   
    bcftools view --threads 4 -s $females $1 > ${out}_females.vcf
    bgzip --threads 4 -f --keep ${out}_females.vcf
    tabix --threads 4 -f -p vcf ${out}_females.vcf.gz
    echo " "

    males=$(awk '{print $1}' data/rooksList_${sp}_males_$K | tr '\n' ',' | sed 's/,$//')
    echo "bcftools view --threads 4 -s $males $1 > ${out}_males.vcf"
    echo "bgzip --threads 4 --keep ${out}_males.vcf"
    echo "tabix --threads 4 -f -p vcf ${out}_males.vcf.gz"
    bcftools view --threads 4 -s $males $1 > ${out}_males.vcf
    bgzip --threads 4 -f --keep ${out}_males.vcf
    tabix --threads 4 -f -p vcf ${out}_males.vcf.gz
    echo " "

    fi

fi

# start analyses
pop="stock"
## set up popfiles
popfile=data/${pop}List_${sp}_$K

if [[ $chr == "chrX" ]]
then
    # output prefix
    out=$(echo $1 | sed 's/\.vcf//')
    # reset popfile
    popfile="data/${pop}List_${sp}_males_$K"
    # males invariants
    echo "running pixy on males invariant sites"
    echo "pixy --stats pi fst dxy watterson_theta tajima_d --n_cores 4 --vcf ${out}_males.vcf.gz --populations $popfile --output_folder resutls/25_pi/25_pixy_$K --output_prefix 25_pixy_invar_${pop}_${chr}__${sp}_males --bed_file $bedfile --fst_type hudson"
    pixy --stats pi fst dxy watterson_theta tajima_d --n_cores 4 --vcf ${out}_males.vcf.gz --populations $popfile --output_folder resutls/25_pi/25_pixy_$K --output_prefix 25_pixy_invar_${pop}_${chr}__${sp}_males --bed_file $bedfile --fst_type hudson
    echo " "

    # reset popfile
    popfile="data/${pop}List_${sp}_females_$K"
    # females invariants
    echo "running pixy on females invariant sites"
    echo "pixy --stats pi fst dxy watterson_theta tajima_d --n_cores 4 --vcf ${out}_females.vcf.gz --populations $popfile --output_folder resutls/25_pi/25_pixy_$K --output_prefix 25_pixy_invar_${pop}_${chr}__${sp}_females --bed_file $bedfile --fst_type hudson"
    pixy --stats pi fst dxy watterson_theta tajima_d --n_cores 4 --vcf ${out}_females.vcf.gz --populations $popfile --output_folder resutls/25_pi/25_pixy_$K --output_prefix 25_pixy_invar_${pop}_${chr}__${sp}_females --bed_file $bedfile --fst_type hudson
    echo " "


else
    # all indivs invariants
    echo "running pixy on all individuals invariants"
    echo "pixy --stats pi fst dxy watterson_theta tajima_d --n_cores 4 --vcf $1.gz --populations $popfile --output_folder resutls/25_pi/25_pixy_$K --output_prefix 25_pixy_invar_${pop}_${chr}__$sp --bed_file $bedfile --fst_type hudson"
    pixy --stats pi fst dxy watterson_theta tajima_d --n_cores 4 --vcf $1.gz --populations $popfile --output_folder resutls/25_pi/25_pixy_$K --output_prefix 25_pixy_invar_${pop}_${chr}__$sp --bed_file $bedfile --fst_type hudson
    echo " "

fi


######################################
##       ANGSD tests workflow       ##  
######################################
## inputs
bamlist="$PROJECT_PATH/data/bamlist_${2}.txt" # $2 is species abbreviation (GSL ot SSL)
### Outputs

## ANGSD no filtring
OUTPUT_PREF="$PROJECT_PATH/resutls/25_pi/angsd_test/angsd_pi_noFilt__$2"
echo "angsd -bam $bamlist -doSaf 1 -anc $FASTA -GL 2 -P 8 -out $OUTPUT_PREF -rf $region_file"
angsd -bam $bamlist -doSaf 1 -anc $FASTA -GL 2 -P 8 -out $OUTPUT_PREF -rf $region_file

echo " "
echo "$angsd/realSFS $OUTPUT_PREF.saf.idx -P 8 -fold 1 > $OUTPUT_PREF.sfs"
$angsd/realSFS $OUTPUT_PREF.saf.idx -P 8 -fold 1 > $OUTPUT_PREF.sfs

echo " "
echo "$angsd/realSFS saf2theta $OUTPUT_PREF.saf.idx -P 8 -outname $OUTPUT_PREF -sfs $OUTPUT_PREF.sfs -fold 1"
$angsd/realSFS saf2theta $OUTPUT_PREF.saf.idx -P 8 -outname $OUTPUT_PREF -sfs $OUTPUT_PREF.sfs -fold 1

echo " "
echo "angsd/thetaStat do_stat $OUTPUT_PREF.thetas.idx -win 50000 -step 10000  -outnames $OUTPUT_PREF.thetasWindow.gz"
$angsd/thetaStat do_stat $OUTPUT_PREF.thetas.idx -win 50000 -step 10000  -outnames $OUTPUT_PREF.thetasWindow.gz



for chr in CM019810.2 CM019819.2
do

## ANGSD sites from variant calling only
echo "angsd pi for $2 - chromosome $chr"

region_file="$PROJECT_PATH/data/sites_file_angsdPi_${2}_$chr.txt"
OUTPUT_PREF="$PROJECT_PATH/resutls/25_pi/angsd_test/angsd_pi_sites_FEM_ONLY__${2}_$chr"

echo "angsd sites index $region_file"
angsd sites index $region_file

echo "angsd -bam $bamlist -doSaf 1 -anc $FASTA -GL 2 -P 8 -out $OUTPUT_PREF -r $chr: -sites $region_file"
angsd -bam $bamlist -doSaf 1 -anc $FASTA -GL 2 -P 8 -out $OUTPUT_PREF -r $chr: -sites $region_file

echo " "
echo "$angsd/realSFS $OUTPUT_PREF.saf.idx -P 8 -fold 1 > $OUTPUT_PREF.sfs"
$angsd/realSFS $OUTPUT_PREF.saf.idx -P 8 -fold 1 > $OUTPUT_PREF.sfs

echo " "
echo "$angsd/realSFS saf2theta $OUTPUT_PREF.saf.idx -P 8 -outname $OUTPUT_PREF -sfs $OUTPUT_PREF.sfs -fold 1"
$angsd/realSFS saf2theta $OUTPUT_PREF.saf.idx -P 8 -outname $OUTPUT_PREF -sfs $OUTPUT_PREF.sfs -fold 1

echo " "
echo "angsd/thetaStat do_stat $OUTPUT_PREF.thetas.idx -win 50000 -step 10000  -outnames $OUTPUT_PREF"
$angsd/thetaStat do_stat $OUTPUT_PREF.thetas.idx -win 50000 -step 10000  -outnames $OUTPUT_PREF


## ANGSD sites from variant calling only + MQ40
OUTPUT_PREF="$PROJECT_PATH/resutls/25_pi/angsd_test/angsd_pi_sites_FEM_ONLY_angsdFilter_MQ40_MD1000_NA01__${2}_$chr"

echo "angsd -bam $bamlist -doSaf 1 -anc $FASTA -minMapQ 40 -GL 2 -minInd $nind -doCounts 1 -setMaxDepth 1000 -P 8 -out $OUTPUT_PREF -r $chr: -sites $region_file"
angsd -bam $bamlist -doSaf 1 -anc $FASTA -minMapQ 40 -GL 2 -minInd $nind -doCounts 1 -setMaxDepth 1000 -P 8 -out $OUTPUT_PREF -r $chr: -sites $region_file

echo " "
echo "$angsd/realSFS $OUTPUT_PREF.saf.idx -P 8 -fold 1 > $OUTPUT_PREF.sfs"
$angsd/realSFS $OUTPUT_PREF.saf.idx -P 8 -fold 1 > $OUTPUT_PREF.sfs

echo " "
echo "$angsd/realSFS saf2theta $OUTPUT_PREF.saf.idx -P 8 -outname $OUTPUT_PREF -sfs $OUTPUT_PREF.sfs -fold 1"
$angsd/realSFS saf2theta $OUTPUT_PREF.saf.idx -P 8 -outname $OUTPUT_PREF -sfs $OUTPUT_PREF.sfs -fold 1

echo " "
echo "angsd/thetaStat do_stat $OUTPUT_PREF.thetas.idx -win 50000 -step 10000  -outnames $OUTPUT_PREF"
$angsd/thetaStat do_stat $OUTPUT_PREF.thetas.idx -win 50000 -step 10000  -outnames $OUTPUT_PREF

## ANGSD sites from variant calling only + MQ40 + max depth
OUTPUT_PREF="$PROJECT_PATH/resutls/25_pi/angsd_test/angsd_pi_sites_FEM_ONLY_angsdFilter_MQ40_MDindiv15_NA01__${2}_$chr"

echo "angsd -bam $bamlist -doSaf 1 -anc $FASTA -minMapQ 40 -GL 2 -minInd $nind -doCounts 1 -setMaxDepth $maxD -P 8 -out $OUTPUT_PREF -r $chr: -sites $region_file"
angsd -bam $bamlist -doSaf 1 -anc $FASTA -minMapQ 40 -GL 2 -minInd $nind -doCounts 1 -setMaxDepth $maxD -P 8 -out $OUTPUT_PREF -r $chr: -sites $region_file

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


