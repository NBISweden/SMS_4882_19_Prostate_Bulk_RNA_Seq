#!/usr/bash
###################################################################
#0. Resources			                                  #
###################################################################
#ANNOTATIONS=($PROJECT_DIR'data/meta/annotation/BALB_CJ_annotation.bed' $PROJECT_DIR'data/meta/annotation/Mouse_annotation.bed')
PROJECT_DIR='/crex/proj/snic2019-30-25/nobackup/private/UserDirectories/project_template/'
#SAMPLES=$PROJ_DIR'Samples_read_Clone1N.txt'
SAMPLES=$PROJ_DIR'Samples_read.txt'
INTERMEDIATE_DIR=$PROJECT_DIR'intermediate/'
RESULTS_DIR=$PROJECT_DIR'results'
CODE='/domus/h1/nimar/SMS_4882_19/nobackup/private/UserDirectories/project_template/code/'
RAW_READS='/proj/snic2019-30-25/private/data/Demultiplexing_2020/2019_61_R1/2019_61_R1/'
TRIMMED_READS=$PROJECT_DIR'data/Trimmed-reads/'
SAMPLE_FILE=$PROJECT_DIR'Samples.txt' 
STAR=$INTERMEDIATE_DIR'STAR_BALB_CJ/'
STAR_QC=$INTERMEDIATE_DIR'STAR_BALB_CJ_QC/'
FEATURECOUNTS_DIR=$INTERMEDIATE_DIR'featureCounts/STAR_BALB_CJ'
#GENOME_DIR=$PROJECT_DIR'data/meta/reference/STARIndex_BALB_CJ/'
GENOME_DIR=$PROJECT_DIR'data/meta/reference/STARIndex_BALB_CJ/'
#GTF_FILE=$PROJECT_DIR'data/meta/annotation/Merged_annotation.gff3'
GTF_FILE=$PROJECT_DIR'data/meta/annotation/BALB_CJ/Mus_musculus_balbcj.BALB_cJ_v1.99.chr.gtf' #This is for BALB_CJ genome
GFF_FILE=$PROJECT_DIR'data/meta/annotation/BALB_CJ/Mus_musculus_balbcj.BALB_cJ_v1.99.chr.gff3'

HUMAN_REF=$PROJECT_DIR'data/meta/reference/BALB_CJ_genome.bed'
MOUSE_REF=$PROJECT_DIR'data/meta/reference/Mouse_genome.bed'

HUMAN_SIZE=$PROJECT_DIR'data/meta/reference/BALB_CJ_genome.fai'
MOUSE_SIZE=$PROJECT_DIR'data/meta/reference/Mouse_genome.fai'

HUMAN_ANNOTATION_GFF2=$PROJECT_DIR'data/meta/annotation/BALB_CJ_annotation.gff2'
MOUSE_ANNOTATION_GFF2=$PROJECT_DIR'data/meta/annotation/Mouse_annotation.gff2'
THREADS=20

SBATCH="#!/bin/bash -l
#SBATCH -A snic2019-8-295
#SBATCH --mail-type=all
#SBATCH --mail-user=nimarafati@gmail.com"



rm -f run_trim.sh
while read -r sample R1 R2
do
THREADS=2
        echo "$SBATCH
#SBATCH -J ${sample}_trim
#SBATCH -p core -n $THREADS
#SBATCH -t 24:00:00
module load bioinfo-tools trimmomatic
mkdir $TRIMMED_READS/${sample}
cp $RAW_READS/$sample/*fastq.gz \$SNIC_TMP/
cd \$SNIC_TMP
cp /home/nimar/glob_old/Contamination/adapters.fa adapters.fa 
trimmomatic PE -threads $THREADS -phred33 -trimlog ${sample}.log \
$R1 $R2 \
${sample}_R1_paired.fq ${sample}_R1_unpaired.fq \
${sample}_R2_paired.fq ${sample}_R2_unpaired.fq \
ILLUMINACLIP:adapters.fa:2:40:15:8:true SLIDINGWINDOW:4:15 LEADING:10 TRAILING:10 MINLEN:36 HEADCROP:10
#Compress paired reads
gzip ${sample}_R1_paired.fq &
gzip ${sample}_R2_paired.fq &
wait 
#Transfer only paired reads
mv ${sample}_R1_paired.fq.gz $TRIMMED_READS/${sample} &
mv ${sample}_R2_paired.fq.gz $TRIMMED_READS/${sample} &
wait
cd $CODE/
sbatch Sbatch_align_${sample}.script" >Sbatch_trim_${sample}.script
        echo "sbatch Sbatch_trim_${sample}.script" >>run_trim.log
done<$SAMPLES



###################################################################
#1. Align the reads to genome					  #
###################################################################
###################################################################
#2. Index bam files		                                  #       
###################################################################
rm -f run_align.sh
THREADS=20
while read -r sample R1 R2
do
	echo "$SBATCH
#SBATCH -p node
#SBATCH -t 5:00:00
#SBATCH -J ${sample}_align

module load bioinfo-tools star samtools

#Alignign the reads
mkdir $INTERMEDIATE_DIR/STAR_BALB_CJ/${sample}_demultiplexed_200825
cd $INTERMEDIATE_DIR/STAR_BALB_CJ/${sample}_demultiplexed_200825
STAR --genomeDir $GENOME_DIR \
--sjdbGTFfile $GTF_FILE \
--readFilesIn $TRIMMED_READS$sample/$R1 $TRIMMED_READS$sample/$R2 \
--runThreadN  $THREADS \
--twopassMode Basic \
--outWigType bedGraph \
--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 64324509440 \
--readFilesCommand zcat \
--runDirPerm All_RWX  \
--quantMode GeneCounts \
--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 20 \
--outFileNamePrefix $sample --outSAMattrRGline ID:$sample 'SM:${sample}' 

#Simplify name of the bam and count files 
mv ${sample}Aligned.sortedByCoord.out.bam ${sample}.sort.bam
mv ${sample}.ReadsPerGene.out.tab ${sample}.counts

#Indexing bam files
samtools index -@ $THREADS ${sample}.sort.bam
ln -s  ${sample}.sort.bam.bai ${sample}.sort.bai" >Sbatch_align_${sample}.script
echo "sbatch Sbatch_align_${sample}.script" >>run_align.sh
done<$SAMPLES


###################################################################
# QC by QoRTs, qualimap                                           #
###################################################################
#Runnign QoRTs
rm -f run_QoRTs.sh
THREADS=20
while read -r sample R1 R2
do
        echo "$SBATCH
#SBATCH -p node
#SBATCH -t 4:00:00
#SBATCH -J QoRTS_${sample}
#module load bioinfo-tools QualiMap

java -Xmx50G -jar ~/git/QoRTs/QoRTs.jar QC --prefilterImproperPairs --generatePlots --stranded --generatePlots --numThreads $THREADS --maxReadLength $READ_LENGTH  --minMAPQ 20 --title ${sample} \
--addFunctions JunctionCalcs,writeJunctionSeqCounts,writeSpliceExon,writeNovelSplices,annotatedSpliceExonCounts,makeAllBrowserTracks,FPKM,writeGeneCounts,writeDESeq,writeDEXSeq,writeGeneBody,writeGenewiseGeneBody,calcDetailedGeneCounts,makeJunctionBed,writeDocs,writeGeneBodyIv --chromSizes $GENOME_SIZE $STAR/$sample/${sample}.sort.bam \
$GTF_FILE $STAR_QC${sample}_QoRTs/ 
#qualimap bamqc -nt 20 -bam $STAR/$sample/${sample}.sort.bam -c -gff $GFF_FILE -outdir $STAR_QC${sample}_qualimap -p strand-specific-reverse  -sd " >Sbatch_QoRTs_${sample}.script
echo "sbatch Sbatch_QoRTs_${sample}.script" >>run_QoRTs.sh
done<$SAMPLES


###################################################################
# featurecounts                                                   #
###################################################################
rm -f run_featureCounts.sh
while read -r sample R1 R2
do
        echo "$SBATCH
#SBATCH -p core -n 16
#SBATCH -t 1:30:00
#SBATCH -J ${sample}_FC
cd \$SNIC_TMP
cp $STAR/$sample/${sample}.sort.bam \$SNIC_TMP
featureCounts_path=\"$FEATURECOUNTS_DIR/${sample}\"
output=\"count-s-2\"
annotation=\"$GTF_FILE\"
mkdir -p \$featureCounts_path
~/glob/Software/subread-2.0.0-source/bin/featureCounts -s 2 -t exon -g gene_id -F GFF -C -T $THREADS -p -B --primary -Q 20 \
-o \$output \
-a \$annotation \
${sample}.sort.bam 
mv \$output* \$featureCounts_path" >Sbatch_featurecounts_${sample}.script
echo "sbatch Sbatch_featurecounts_${sample}.script " >>run_featureCounts.sh
done<$SAMPLES

