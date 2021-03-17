#!/usr/bash
###################################################################
#0. Resources			                                  #
###################################################################
#ANNOTATIONS=($PROJECT_DIR'data/meta/annotation/Human_annotation.bed' $PROJECT_DIR'data/meta/annotation/Mouse_annotation.bed')
PROJECT_DIR='/crex/proj/snic2019-30-25/private/UserDirectories/SMS_4882_19_Prostate_Bulk_RNA_Seq/'
CODE=$PROJECT_DIR'code/GSNAP_Human/'
SAMPLES=$CODE'Samples_read.txt'
INTERMEDIATE_DIR=$PROJECT_DIR'intermediate/'
RESULTS_DIR=$PROJECT_DIR'results'
RAW_READS=''
TRIMMED_READS=$PROJECT_DIR'data/Trimmed-reads/'
SAMPLE_FILE=$PROJECT_DIR'Samples.txt' 
GSNAP=$INTERMEDIATE_DIR'GSNAP_Human/'
GSNAP_QC=$INTERMEDIATE_DIR'GSNAP_Human_QC/'
GENOME_DIR=$PROJECT_DIR'data/meta/reference/Human/'
GENOME='GSNAPIndex_Human'
GTF_FILE=$PROJECT_DIR'data/meta/annotation/Human_annotation.gff2'
ANNOTATION_GFF2=$PROJECT_DIR'data/meta/annotation/Human_annotation.gtf'
GENOME_SIZE=$PROJECT_DIR'data/meta/reference/Human_genome.fai'

XENOFILTER_QC=$INTERMEDIATE_DIR'XenofilteR_GSNAP_Human_QC/'
XENOFILTER=$INTERMEDIATE_DIR'XenofilteR_GSNAP/'

THREADS=20
READ_LENGTH=150

mkdir $GSNAP_QC
SBATCH="#!/bin/bash -l
#SBATCH -A snic2020-15-15
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nimarafati@gmail.com"

###################################################################
#1. Align the reads to genome					  #
###################################################################
###################################################################
#2. Index bam files		                                  #       
###################################################################
rm -f run_align.sh
while read -r sample R1 R2
do
	echo "$SBATCH
#SBATCH -p node
#SBATCH -t 72:00:00
#SBATCH -J ${sample}_gsnap_Human

module load bioinfo-tools gmap-gsnap/2017-09-11 samtools zlib

#Alignign the reads
cd \$SNIC_TMP
cp $TRIMMED_READS$sample/$R1 \$SNIC_TMP &
cp $TRIMMED_READS$sample/$R2 \$SNIC_TMP &
wait

mkdir -p $GSNAP/$sample
gsnap --gunzip -m 0.1 -D $GENOME_DIR -d $GENOME --orientation FR -B 5 -N 1 -n 30 -E 4 --nthreads=$THREADS --gmap-mode=all -A sam -J 33 -O --quiet-if-excessive \
--read-group-id=$sample --read-group-name=$sample --read-group-library=200PE --read-group-platform=Illumina \
$R1 $R2 | samtools view -bSt $GENOME_SIZE  - >${sample}.bam

sort by coordinate
samtools sort -@ $THREADS -O bam -o ${sample}.sort.bam ${sample}.bam
Indexing bam files
samtools index -@ $THREADS ${sample}.sort.bam
ln -s  ${sample}.sort.bam.bai ${sample}.sort.bai
mv  ${sample}.sort* $GSNAP/$sample
cd $CODE
sbatch Sbatch_QoRTs_${sample}.script" >Sbatch_align_${sample}.script
echo "sbatch Sbatch_align_${sample}.script" >>run_align.sh
done<$SAMPLES

###################################################################
#3. QC by QoRTs		                                          #       
###################################################################
rm -f run_QoRTs.sh
THREADS=10
while read -r sample R1 R2
do
        echo "$SBATCH
#SBATCH -p core -n 10
#SBATCH -t 24:00:00
#SBATCH -J ${sample}_QC_GSNAP_Human
module load bioinfo-tools samtools

cd $GSNAP/$sample
samtools stats -@ 20 ${sample}.sort.bam >${sample}.sort.stats
samtools flagstat -@ 20 ${sample}.sort.bam >${sample}.sort.flagstat

java -Xmx50G -jar ~/git/QoRTs/QoRTs.jar QC --prefilterImproperPairs --generatePlots --stranded --generatePlots --numThreads $THREADS --maxReadLength 141 --minMAPQ 20 --title ${sample} \
--addFunctions JunctionCalcs,writeJunctionSeqCounts,writeSpliceExon,writeNovelSplices,annotatedSpliceExonCounts,makeAllBrowserTracks,FPKM,writeGeneCounts,writeDESeq,writeDEXSeq,writeGeneBody,writeGenewiseGeneBody,calcDetailedGeneCounts,makeJunctionBed,writeDocs,writeGeneBodyIv --chromSizes $GENOME_SIZE $GSNAP$sample/${sample}.sort.bam \
$ANNOTATION_GFF2 $GSNAP_QC${sample}/ " >Sbatch_QoRTs_${sample}.script
echo "sbatch Sbatch_QoRTs_${sample}.script" >>run_QoRTs.sh
done<$SAMPLES





rm -f run_qualimap.sh
THREADS=10
while read -r sample R1 R2
do
        echo "$SBATCH
#SBATCH -p node
#SBATCH -t 24:00:00
#SBATCH -J ${sample}_QC_GSNAP_Human
module load bioinfo-tools QualiMap

qualimap bamqc -nt 20 -c \
--gff $ANNOTATION_GFF2 \
--bam $GSNAP$sample/${sample}.sort.bam \
--outdir  $GSNAP_QC${sample}/ " >Sbatch_QualiMap_${sample}.script
echo "sbatch Sbatch_QualiMap_${sample}.script" >>run_qualimaps.sh
done<$SAMPLES



rm -f run_featurecounts.sh
SAMPLES=$CODE'Samples_read_XenofilteR_available.txt'
THREADS=10
while read -r sample R1 R2
do
        echo "$SBATCH
#SBATCH -p node
#SBATCH -t 24:00:00
#SBATCH -J ${sample}_FC_GSNAP_H

cd \$SNIC_TMP
cp /crex/proj/snic2019-30-25/private/UserDirectories/SMS_4882_19_Prostate_Bulk_RNA_Seq/intermediate/XenofilteR_GSNAP/${sample}/Filtered_bams_Human/${sample}_Filtered.bam \$SNIC_TMP
featureCounts_path='/crex/proj/snic2019-30-25/private/UserDirectories/SMS_4882_19_Prostate_Bulk_RNA_Seq/intermediate/featurecounts/XenofilteR/GSNAP_Human/$sample'
output='count-s-2'
annotation='/crex/proj/snic2019-30-25/private/UserDirectories/SMS_4882_19_Prostate_Bulk_RNA_Seq/data/meta/annotation/Human_annotation.gff3'
mkdir -p \$featureCounts_path
~/glob/Software/subread-2.0.0-source/bin/featureCounts -s 2 -t exon -g gene_id -F GTF -C -T $THREADS -p -B --primary -Q 20 -o \$output -a \$annotation ${sample}_Filtered.bam
mv \$output* \$featureCounts_path " >Sbatch_featurecounts_${sample}_XenofilteR_Human.script
echo "sbatch Sbatch_${sample}_featurecounts_XenofilteR_Human.script" >>run_featurecounts.sh
done<$SAMPLES



###################################################################
#. QoRTs: XenofilteR                                              #       
###################################################################

SAMPLES=$CODE'Samples_read.txt'

rm -f run_QoRTs.sh
THREADS=10
while read -r sample R1 R2
do
        echo "$SBATCH
#SBATCH -p core -n 10
#SBATCH -t 24:00:00
#SBATCH -J ${sample}_QC_XE_Human
module load bioinfo-tools samtools

#cd $GSNAP/$sample
#samtools stats -@ 20 ${sample}.sort.bam >${sample}.sort.stats
#samtools flagstat -@ 20 ${sample}.sort.bam >${sample}.sort.flagstat

java -Xmx50G -jar ~/git/QoRTs/QoRTs.jar QC --prefilterImproperPairs --generatePlots --stranded --generatePlots --numThreads $THREADS --maxReadLength 141 --minMAPQ 20 --title ${sample} \
--addFunctions JunctionCalcs,writeJunctionSeqCounts,writeSpliceExon,writeNovelSplices,annotatedSpliceExonCounts,makeAllBrowserTracks,FPKM,writeGeneCounts,writeDESeq,writeDEXSeq,writeGeneBody,writeGenewiseGeneBody,calcDetailedGeneCounts,makeJunctionBed,writeDocs,writeGeneBodyIv --chromSizes $GENOME_SIZE $XENOFILTER$sample/Filtered_bams_Human/${sample}_Filtered.bam \
$ANNOTATION_GFF2 $XENOFILTER_QC${sample}/ " >Sbatch_QoRTs_XenofilteR_${sample}.script
echo "sbatch Sbatch_QoRTs_XenofilteR_${sample}.script" >>run_QoRTs.sh
done<$SAMPLES
