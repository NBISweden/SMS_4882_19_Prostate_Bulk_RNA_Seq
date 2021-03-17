#!/usr/bash
###################################################################
#0. Resources			                                  #
###################################################################
#ANNOTATIONS=($PROJECT_DIR'data/meta/annotation/Human_annotation.bed' $PROJECT_DIR'data/meta/annotation/Mouse_annotation.bed')
PROJECT_DIR='/crex/proj/snic2020-16-70/nobackup/private/tmp_4882/'
CODE=$PROJECT_DIR'code/GSNAP_BALB_CJ/'
SAMPLES=$CODE'Samples_read_R2_R3.txt' #Samples_read_XenofilteR.txt' #Samples_read.txt'
echo $SAMPLES
INTERMEDIATE_DIR=$PROJECT_DIR'intermediate/'
RESULTS_DIR=$PROJECT_DIR'results'
RAW_READS=''
TRIMMED_READS=$PROJECT_DIR'data/Trimmed-reads/'
SAMPLE_FILE=$PROJECT_DIR'Samples.txt' 
GSNAP=$INTERMEDIATE_DIR'GSNAP_BALB_CJ/'
GSNAP_HUMAN=$INTERMEDIATE_DIR'GSNAP_Human/'
GSNAP_QC=$INTERMEDIATE_DIR'GSNAP_BALB_CJ_QC/'
GENOME_DIR=$PROJECT_DIR'data/meta/reference/BALB_CJ/'
GENOME='GSNAPIndex_BALB_CJ'
GTF_FILE=$PROJECT_DIR'data/meta/annotation/Human_annotation.gff2'
ANNOTATION_GFF2=$PROJECT_DIR'data/meta/annotation/BALB_CJ/Mus_musculus_balbcj.BALB_cJ_v1.99.chr.gtf'
GENOME_SIZE=$PROJECT_DIR'data/meta/reference/BALB_CJ/Mus_musculus_balbcj.BALB_cJ_v1.fa.fai'

XENOFILTER=$INTERMEDIATE_DIR'XenofilteR_GSNAP/'
DISAMBIGUATE=$INTERMEDIATE_DIR'Disambiguate_GSNAP'

XENOFILTER_QC=$INTERMEDIATE_DIR'XenofilteR_GSNAP_BALB_CJ_QC/'

THREADS=20
READ_LENGTH=150

mkdir $GSNAP_QC $DISAMBIGUATE $XENOFILTER
SBATCH="#!/bin/bash -l
#SBATCH -A snic2019-8-295
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
#SBATCH -J ${sample}_gsnap_BALB

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

#sort by coordinate
samtools sort -@ $THREADS -O bam -o ${sample}.sort.bam ${sample}.bam
#Indexing bam files
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
#Runnign QoRTs
rm -f run_QoRTs.sh
THREADS=10
while read -r sample R1 R2
do
        echo "$SBATCH
#SBATCH -p core -n 10
#SBATCH -t 24:00:00
#SBATCH -J ${sample}_QC_GSNAP_BALB_CJ
module load bioinfo-tools samtools

cd $GSNAP/$sample
samtools stats -@ 20 ${sample}.sort.bam >${sample}.sort.stats
samtools flagstat -@ 20 ${sample}.sort.bam >${sample}.sort.flagstat

java -Xmx50G -jar ~/git/QoRTs/QoRTs.jar QC --prefilterImproperPairs --generatePlots --stranded --generatePlots --numThreads $THREADS --maxReadLength 141 --minMAPQ 20 --title ${sample} \
--addFunctions JunctionCalcs,writeJunctionSeqCounts,writeSpliceExon,writeNovelSplices,annotatedSpliceExonCounts,makeAllBrowserTracks,FPKM,writeGeneCounts,writeDESeq,writeDEXSeq,writeGeneBody,writeGenewiseGeneBody,calcDetailedGeneCounts,makeJunctionBed,writeDocs,writeGeneBodyIv --chromSizes $GENOME_SIZE $GSNAP$sample/${sample}.sort.bam \
$ANNOTATION_GFF2 $GSNAP_QC${sample}/ " >Sbatch_QoRTs_${sample}.script
echo "sbatch Sbatch_QoRTs_${sample}.script" >>run_QoRTs.sh
done<$SAMPLES

###################################################################
# Disambiguate                                                    #
###################################################################
rm -f run_disambiguate.sh
while read -r sample R1 R2
do
        echo "$SBATCH
#SBATCH -p core -n 1
#SBATCH -t 12:00:00
#SBATCH -J Disambg_$sample
module load bioinfo-tools samtools

cp $GSNAP/$sample/${sample}.sort.bam \$SNIC_TMP/${sample}_BALB_CJ_original.sort.bam &
cp $GSNAP_HUMAN/$sample/${sample}.sort.bam \$SNIC_TMP/${sample}_Human_original.sort.bam &
wait

cd \$SNIC_TMP
~/git/disambiguate_conda/bin/ngs_disambiguate -s $sample -a star \
-o ${sample} \
${sample}_BALB_CJ_original.sort.bam ${sample}_Human_original.sort.bam 
mkdir -p $DISAMBIGUATE/$sample/

mv $sample/${sample}.disambiguatedSpeciesA.bam $sample/${sample}_BALB_CJ.sort.bam
mv $sample/${sample}.disambiguatedSpeciesB.bam $sample/${sample}_Human.sort.bam

mv $sample $DISAMBIGUATE/

cd $CODE/  ">Sbatch_Disambiguate_${sample}.script
echo "sbatch Sbatch_Disambiguate_${sample}.script" >>run_disambiguate.sh
done<$SAMPLES


###################################################################
# XenofiltR                                                       #
###################################################################
rm -f run_xenofiltR.sh
while read -r sample R1 R2
do
        echo "$SBATCH
#SBATCH -p core -n 10
#SBATCH -t 6:00:00
#SBATCH -J XenofiltR_${sample}_GSNAP
module load bioinfo-tools samtools R_packages

mkdir -p $XENOFILTER/$sample/
cd $XENOFILTER/$sample/

ln -s  $GSNAP/$sample/${sample}.sort.bam ${sample}_BALB_CJ.sort.bam 
ln -s $GSNAP_HUMAN/$sample/${sample}.sort.bam ${sample}_Human.sort.bam 


#Human
echo \"library(XenofilteR)
library(Rsamtools)
library(GenomicAlignments)
library(BiocParallel)
library(futile.logger)

bp.param <- SnowParam(workers = 10, type = 'SOCK')

samples.list <-data.frame(Human = '${sample}_Human.sort.bam', 
                         BALB_CJ = '${sample}_BALB_CJ.sort.bam')
output <- c('$sample')

XenofilteR(sample.list = samples.list, destination.folder = '$XENOFILTER/$sample/',  bp.param = bp.param, output.names = output)
\" >cmd_Human.Rscript

R --vanilla -q <cmd_Human.Rscript

mv $XENOFILTER$sample/Filtered_bams/ \
$XENOFILTER$sample/Filtered_bams_Human/

#BALB_CJ
echo \"library(XenofilteR)
library(Rsamtools)
library(GenomicAlignments)
library(BiocParallel)
library(futile.logger)

bp.param <- SnowParam(workers = 10, type = 'SOCK')

samples.list <-data.frame(BALB_CJ = '${sample}_BALB_CJ.sort.bam',
                         Human = '${sample}_Human.sort.bam')
output <- c('$sample')

XenofilteR(sample.list = samples.list, destination.folder = '$XENOFILTER/$sample/',  bp.param = bp.param, output.names = output)
\" >cmd_BALB_CJ.Rscript

R --vanilla -q < cmd_BALB_CJ.Rscript 

mv $XENOFILTER$sample/Filtered_bams/ \
$XENOFILTER$sample/Filtered_bams_BALB_CJ/

samtools view -f8 $XENOFILTER$sample/Filtered_bams_BALB_CJ/${sample}_Filtered.bam \
| sed 's/.*NM:i://' | cut -f1 \
> $XENOFILTER$sample/Filtered_bams_BALB_CJ/${sample}_Filtered_nM.txt &

samtools view -f8 $XENOFILTER$sample/Filtered_bams_Human/${sample}_Filtered.bam \
| sed 's/.*NM:i://' | cut -f1 \
> $XENOFILTER$sample/Filtered_bams_Human/${sample}_Filtered_nM.txt &

#samtools view -f8 $GSNAP/$sample/${sample}.sort.bam \
| sed 's/.*NM:i://' | cut -f1 | grep -v \"NB501986\" \
> $GSNAP/$sample/${sample}_nM.txt &

#samtools view -f8 $GSNAP_HUMAN/$sample/${sample}.sort.bam \
| sed 's/.*NM:i://' | cut -f1 | grep -v \"NB501986\" \
> $GSNAP_HUMAN/$sample/${sample}_nM.txt &
wait
cd $CODE
sbatch Sbatch_featurecounts_${sample}_Xenofilter_Human.script
sbatch Sbatch_featurecounts_${sample}_Xenofilter_BALB_CJ.script " >Sbatch_XenofiltR_${sample}.script
echo "sbatch Sbatch_XenofiltR_${sample}.script" >>run_xenofiltR.sh
done<$SAMPLES


###################################################################
#. QoRTs: XenofilteR                                              #       
###################################################################


#Runnign QoRTs
rm -f run_QoRTs.sh
THREADS=10
while read -r sample R1 R2
do
        echo "$SBATCH
#SBATCH -p core -n 10
#SBATCH -t 24:00:00
#SBATCH -J ${sample}_QC_GSNAP_XE_BALB_CJ
module load bioinfo-tools samtools

cd $GSNAP/$sample
samtools stats -@ 20 ${sample}.sort.bam >${sample}.sort.stats
samtools flagstat -@ 20 ${sample}.sort.bam >${sample}.sort.flagstat

java -Xmx50G -jar ~/git/QoRTs/QoRTs.jar QC --prefilterImproperPairs --generatePlots --stranded --generatePlots --numThreads $THREADS --maxReadLength 141 --minMAPQ 20 --title ${sample} \
--addFunctions JunctionCalcs,writeJunctionSeqCounts,writeSpliceExon,writeNovelSplices,annotatedSpliceExonCounts,makeAllBrowserTracks,FPKM,writeGeneCounts,writeDESeq,writeDEXSeq,writeGeneBody,writeGenewiseGeneBody,calcDetailedGeneCounts,makeJunctionBed,writeDocs,writeGeneBodyIv --chromSizes $GENOME_SIZE $XENOFILTER$sample/Filtered_bams_BALB_CJ/${sample}_Filtered.bam \
$ANNOTATION_GFF2 $XENOFILTER_QC${sample}/ " >Sbatch_QoRTs_XenofilteR_${sample}.script
echo "sbatch Sbatch_QoRTs_XenofilteR_${sample}.script" >>run_QoRTs.sh
done<$SAMPLES



###################################################################
#. Featurecount                                                  #       
###################################################################
rm -f run_featureCounts.sh
#BALB_CJ
FEATURECOUNTS_DIR=$INTERMEDIATE_DIR'featurecounts/XenofilteR/GSNAP_BALB_CJ/'
GFF_FILE=$PROJECT_DIR'data/meta/annotation/BALB_CJ/Mus_musculus_balbcj.BALB_cJ_v1.99.chr.gff3'

while read -r sample R1 R2
do
        echo "$SBATCH
#SBATCH -p core -n 10
#SBATCH -t 3:00:00
#SBATCH -J ${sample}_FC_BALB_CJ
cd \$SNIC_TMP
cp $XENOFILTER$sample/Filtered_bams_BALB_CJ/${sample}_Filtered.bam \$SNIC_TMP
featureCounts_path=\"$FEATURECOUNTS_DIR/${sample}\"
output=\"count-s-2\"
annotation=\"$ANNOTATION_GFF2\"
mkdir -p \$featureCounts_path
~/glob/Software/subread-2.0.0-source/bin/featureCounts -s 2 -t exon -g gene_id -F GFF -C -T $THREADS -p -B --primary -Q 20 \
-o \$output \
-a \$annotation \
${sample}_Filtered.bam
mv \$output* \$featureCounts_path" >Sbatch_featurecounts_${sample}_Xenofilter_BALB_CJ.script
echo "sbatch Sbatch_featurecounts_${sample}_Xenofilter_BALB_CJ.script " >>run_featureCounts.sh
done<$SAMPLES

#Human
FEATURECOUNTS_DIR=$INTERMEDIATE_DIR'featurecounts/XenofilteR/GSNAP_Human/'
GFF_FILE='/crex/proj/snic2019-30-25/nobackup/private/UserDirectories/project_template/data/meta/annotation/Human_annotation.gff3'

while read -r sample R1 R2
do
        echo "$SBATCH
#SBATCH -p core -n 10
#SBATCH -t 3:30:00
#SBATCH -J ${sample}_FC_Human
cd \$SNIC_TMP
cp $XENOFILTER$sample/Filtered_bams_Human/${sample}_Filtered.bam \$SNIC_TMP
featureCounts_path=\"$FEATURECOUNTS_DIR/${sample}\"
output=\"count-s-2\"
annotation=\"$GFF_FILE\"
mkdir -p \$featureCounts_path
~/glob/Software/subread-2.0.0-source/bin/featureCounts -s 2 -t exon -g gene_id -F GFF -C -T $THREADS -p -B --primary -Q 20 \
-o \$output \
-a \$annotation \
${sample}_Filtered.bam
mv \$output* \$featureCounts_path" >Sbatch_featurecounts_${sample}_Xenofilter_Human.script
echo "sbatch Sbatch_featurecounts_${sample}_Xenofilter_Human.script " >>run_featureCounts.sh
done<$SAMPLES


#overrepresentation
echo "Sample_Bone10	Sample_Bone10_13_S5_R1_001-paired_fastqc
Sample_Bone11	Sample_Bone11_14_S9_R1_001-paired_fastqc
Sample_Bone12	Sample_Bone12_15_S3_R1_001-paired_fastqc
Sample_Bone2	Sample_Bone2_6_S3_R1_001-paired_fastqc
Sample_Bone3	Sample_Bone3_7_S3_R1_001-paired_fastqc
Sample_Bone4	Sample_Bone4_8_S6_R1_001-paired_fastqc
Sample_Bone5	Sample_Bone5_9_S8_R1_001-paired_fastqc
Sample_Bone6	Sample_Bone6_10_S4_R1_001-paired_fastqc
Sample_Bone7	Sample_Bone7_11_S7_R1_001-paired_fastqc
Sample_Bone9	Sample_Bone9_12_S2_R1_001-paired_fastqc
Sample_Clone1N	Sample_Clone1N_16_S6_R1_001-paired_fastqc
Sample_Clone1o	Sample_Clone1o_19_S7_R1_001-paired_fastqc
Sample_Clone2N	Sample_Clone2N_17_S1_R1_001-paired_fastqc
Sample_Clone2o	Sample_Clone2o_20_S6_R1_001-paired_fastqc
Sample_Clone3N	Sample_Clone3N_18_S4_R1_001-paired_fastqc
Sample_Clone3o	Sample_Clone3o_21_S9_R1_001-paired_fastqc
Sample_liver2	Sample_liver2_1_S1_R1_001-paired_fastqc
Sample_liver4	Sample_liver4_2_S2_R1_001-paired_fastqc
Sample_liver5	Sample_liver5_3_S5_R1_001-paired_fastqc
Sample_liver6	Sample_liver6_4_S2_R1_001-paired_fastqc
Sample_liver9	Sample_liver9_5_S5_R1_001-paired_fastqc
Sample_Scr1N	Sample_Scr1N_22_S8_R1_001-paired_fastqc
Sample_Scr1o	Sample_Scr1o_25_S9_R1_001-paired_fastqc
Sample_Scr2N	Sample_Scr2N_23_S4_R1_001-paired_fastqc
Sample_Scr2o	Sample_Scr2o_26_S8_R1_001-paired_fastqc
Sample_Scr3N	Sample_Scr3N_24_S7_R1_001-paired_fastqc
Sample_Scr3o	Sample_Scr3o_27_S1_R1_001-paired_fastqc
Sample_Bone10	Sample_Bone10_13_S5_R2_001-paired_fastqc
Sample_Bone11	Sample_Bone11_14_S9_R2_001-paired_fastqc
Sample_Bone12	Sample_Bone12_15_S3_R2_001-paired_fastqc
Sample_Bone2	Sample_Bone2_6_S3_R2_001-paired_fastqc
Sample_Bone3	Sample_Bone3_7_S3_R2_001-paired_fastqc
Sample_Bone4	Sample_Bone4_8_S6_R2_001-paired_fastqc
Sample_Bone5	Sample_Bone5_9_S8_R2_001-paired_fastqc
Sample_Bone6	Sample_Bone6_10_S4_R2_001-paired_fastqc
Sample_Bone7	Sample_Bone7_11_S7_R2_001-paired_fastqc
Sample_Bone9	Sample_Bone9_12_S2_R2_001-paired_fastqc
Sample_Clone1N	Sample_Clone1N_16_S6_R2_001-paired_fastqc
Sample_Clone1o	Sample_Clone1o_19_S7_R2_001-paired_fastqc
Sample_Clone2N	Sample_Clone2N_17_S1_R2_001-paired_fastqc
Sample_Clone2o	Sample_Clone2o_20_S6_R2_001-paired_fastqc
Sample_Clone3N	Sample_Clone3N_18_S4_R2_001-paired_fastqc
Sample_Clone3o	Sample_Clone3o_21_S9_R2_001-paired_fastqc
Sample_liver2	Sample_liver2_1_S1_R2_001-paired_fastqc
Sample_liver4	Sample_liver4_2_S2_R2_001-paired_fastqc
Sample_liver5	Sample_liver5_3_S5_R2_001-paired_fastqc
Sample_liver6	Sample_liver6_4_S2_R2_001-paired_fastqc
Sample_liver9	Sample_liver9_5_S5_R2_001-paired_fastqc
Sample_Scr1N	Sample_ScR2N_22_S8_R2_001-paired_fastqc
Sample_Scr1o	Sample_ScR2o_25_S9_R2_001-paired_fastqc
Sample_Scr2N	Sample_Scr2N_23_S4_R2_001-paired_fastqc
Sample_Scr2o	Sample_Scr2o_26_S8_R2_001-paired_fastqc
Sample_Scr3N	Sample_Scr3N_24_S7_R2_001-paired_fastqc
Sample_Scr3o	Sample_Scr3o_27_S1_R2_001-paired_fastqc" >fastqc.list 
 echo "library(dplyr)
library(magrittr)" >Extract_Overrepresented_sequences_R1.Rscript
while read -r sample fastqc
do
	echo "
setwd('~/SMS_4882_19/nobackup/private/UserDirectories/project_template/intermediate/FastQC_Trimmed/${sample}_QC/')
fastqc <- '$fastqc/fastqc_data.txt' 
fastqcr::qc_read(fastqc)\$overrepresented_sequences %>%
     dplyr::mutate(name=paste('>',1:n(),'-',Count,'-','$sample',sep=''),fa=paste(name,Sequence,sep='\n')) %>%
     dplyr::pull(fa) %>% 
     readr::write_lines('overrepresented_${fastqc}.fa')" >>Extract_Overrepresented_sequences_R1.Rscript
done<fastqc.list

echo "$SBATCH
#SBATCH -p node -n 1
#SBATCH -t 48:00:00
#SBATCH -J blast_overrepresented
module load bioinfo-tools blast
mkdir $INTERMEDIATE_DIR/blast_overrepresented/
blastn -query ~/SMS_4882_19/nobackup/private/UserDirectories/project_template/intermediate/FastQC_Trimmed/All_overrepresented.fa -db /sw/data/blast_databases/nt -out $INTERMEDIATE_DIR/blast_overrepresented/all_overrepresented.blast -evalue 0.00001 -outfmt '6 qseqid sseqid stitle pident evalue scomname staxids ssciname length mismatch qlen slen' -num_alignments 1 -html -sorthits 0 -num_threads 16
rm -f fastqc.list" >Sbatch_blast_overrepresented.script
