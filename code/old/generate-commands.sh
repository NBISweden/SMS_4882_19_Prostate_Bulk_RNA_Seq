#!/usr/bash
###################################################################
#0. Resources			                                  #
###################################################################
#ANNOTATIONS=($PROJECT_DIR'data/meta/annotation/Human_annotation.bed' $PROJECT_DIR'data/meta/annotation/Mouse_annotation.bed')
PROJECT_DIR='/crex/proj/snic2019-30-25/nobackup/private/UserDirectories/project_template/'
SAMPLES=$PROJ_DIR'Samples_read.txt'
INTERMEDIATE_DIR=$PROJECT_DIR'intermediate/'
RESULTS_DIR=$PROJECT_DIR'results'
CODE=$PROJECT_DIR'code/'
RAW_READS=$PROJECT_DIR'data/raw_internal/'
TRIMMED_READS=$PROJECT_DIR'data/Trimmed-reads/'
SAMPLE_FILE=$PROJECT_DIR'Samples.txt' 
STAR=$INTERMEDIATE_DIR'STAR/'
STAR_QC=$INTERMEDIATE_DIR'STAR_QC/'
GENOME_DIR=$PROJECT_DIR'data/meta/reference/STARIndex/'
#GTF_FILE=$PROJECT_DIR'data/meta/annotation/Merged_annotation.gff3'
GTF_FILE=$PROJECT_DIR'data/meta/annotation/BALB_CJ/Mus_musculus_balbcj.BALB_cJ_v1.99.chr.gtf' #This is for BALB_CJ genome

HUMAN_REF=$PROJECT_DIR'data/meta/reference/Human_genome.bed'
MOUSE_REF=$PROJECT_DIR'data/meta/reference/Mouse_genome.bed'

HUMAN_SIZE=$PROJECT_DIR'data/meta/reference/Human_genome.fai'
MOUSE_SIZE=$PROJECT_DIR'data/meta/reference/Mouse_genome.fai'
BALB_CJ_SIZE=$PROJECT_DIR'data/meta/reference/BALB_CJ/Mus_musculus_balbcj.BALB_cJ_v1.fa.fai'

HUMAN_ANNOTATION_GFF2=$PROJECT_DIR'data/meta/annotation/Human_annotation.gff2'
MOUSE_ANNOTATION_GFF2=$PROJECT_DIR'data/meta/annotation/Mouse_annotation.gff2'
BALB_CJ_ANNOTATION_GFF2=$PROJECT_DIR'data/meta/annotation/BALB_CJ/Mus_musculus_balbcj.BALB_cJ_v1.99.chr.gtf'

THREADS=20

SBATCH="#!/bin/bash -l
#SBATCH -A snic2019-8-295
#SBATCH --mail-type=all
#SBATCH --mail-user=nimarafati@gmail.com"

###################################################################
# FastQC                                                          #
###################################################################
rm -rf run_fastqc.txt
while read -r sample R1 R2
do
        echo "$SBATCH
#SBATCH -J ${sample}_FastQC
#SBATCH -p core -n 4
#SBATCH -t 2:00:00
module load bioinfo-tools FastQC

cd \$SNIC_TMP
cp $RAW_READS/$sample/*fastq.gz  \$SNIC_TMP/ 
mkdir  ${sample}_QC
mkdir -p $INTERMEDIATE_DIR/FastQC/
fastqc -t 4 -o ${sample}_QC -f fastq -t 4 *fastq.gz
mv ${sample}_QC $INTERMEDIATE_DIR/FastQC/" >Sbatch_fastqc_${sample}.script
echo "sbatch Sbatch_fastqc_${sample}.script" >>run_fastqc.txt
done<$SAMPLES
#Samples_reads.txt

###################################################################
# FastQC after trimming                                           #
###################################################################
rm -rf run_fastqc_trimmed.txt
while read -r sample R1 R2
do
        echo "$SBATCH
#SBATCH -J ${sample}_FastQC
#SBATCH -p core -n 4
#SBATCH -t 2:00:00
module load bioinfo-tools FastQC

cd \$SNIC_TMP
cp $TRIMMED_READS/$sample/*fq.gz  \$SNIC_TMP/ 
mkdir  ${sample}_QC
mkdir -p $INTERMEDIATE_DIR/FastQC_Trimmed/
fastqc -t 4 -o ${sample}_QC -f fastq *fq.gz
mv ${sample}_QC $INTERMEDIATE_DIR/FastQC_Trimmed/" >Sbatch_fastqc_trimmed_${sample}.script
echo "sbatch Sbatch_fastqc_trimmed_${sample}.script" >>run_fastqc_trimmed.txt
done<$SAMPLES


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
#SBATCH -t 5:00:00
#SBATCH -J ${sample}_align

module load bioinfo-tools star samtools

#Alignign the reads
cd $INTERMEDIATE_DIR/STAR/$sample
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
#3. Separate alignments by genome                                 #       
###################################################################
rm -f run_split.sh
THREADS=10
while read -r sample R1 R2
do
        echo "$SBATCH
#SBATCH -p node
#SBATCH -t 5:00:00
#SBATCH -J ${sample}_split

module load bioinfo-tools samtools
cd $INTERMEDIATE_DIR/STAR/$sample
samtools view -@ $THREADS -h -L $MOUSE_REF ${sample}.sort.bam | samtools view -bS - >${sample}_mouse.sort.bam & 
samtools view -@ $THREADS -h -L $HUMAN_REF ${sample}.sort.bam | samtools view -bS - >${sample}_human.sort.bam &
wait 
samtools index -@ $THREADS ${sample}_mouse.sort.bam &
samtools index -@ $THREADS ${sample}_human.sort.bam &
wait
ln -s ${sample}_mouse.sort.bam.bai ${sample}_mouse.sort.bai
ln -s ${sample}_human.sort.bam.bai ${sample}_human.sort.bai " >Sbatch_split_genomes_${sample}.script
echo "sbatch Sbatch_split_genomes_${sample}.script" >>run_split.sh
done<$SAMPLES

###################################################################
#4. QC of the alignments                                          #       
###################################################################
#To generate decoder file for QC we need to have number of all aligned reads and multi_mapped reads:
#For this we first generate stats by samtools and parse the output to extract these numbers 

rm -f run_read_stats.sh
THREADS=10
while read -r sample R1 R2
do
        echo "$SBATCH
#SBATCH -p node
#SBATCH -t 00:30:00
#SBATCH -J ${sample}_stats
module load bioinfo-tools samtools
python samtools_stats_parser_V1.py --threads $THREADS --bam $STAR$sample/${sample}_human.sort.bam & 
python samtools_stats_parser_V1.py --threads $THREADS --bam $STAR$sample/${sample}_mouse.sort.bam &
wait" >Sbatch_read_stats_${sample}.script
echo "sbatch Sbatch_read_stats_${sample}.script" >>run_read_stats.sh
done<$SAMPLES

#Runnign QoRTs
rm -f run_QoRTs.sh
THREADS=10
while read -r sample R1 R2
do
        echo "$SBATCH
#SBATCH -p node
#SBATCH -t 12:00:00
#SBATCH -J ${sample}_QC_post_mapping
module load bioinfo-tools QoRTs
cd $INTERMEDIATE_DIR/STAR_BALB_CJ/$sample
#java -Xmx50G -jar ~/git/QoRTs/
QoRTs QC --generatePlots --stranded --generatePlots --numThreads $THREADS --maxReadLength 141 --minMAPQ 20 --title ${sample}_BALB_CJ \
--addFunctions JunctionCalcs,writeJunctionSeqCounts,writeSpliceExon,writeNovelSplices,annotatedSpliceExonCounts,makeAllBrowserTracks,FPKM,writeGeneCounts,writeDESeq,writeDEXSeq,writeGeneBody,writeGenewiseGeneBody,calcDetailedGeneCounts,makeJunctionBed,writeDocs,writeGeneBodyIv --chromSizes $BALB_CJ_SIZE $INTERMEDIATE_DIR/STAR_BALB_CJ/$sample/${sample}.sort.bam \
$BALB_CJ_ANNOTATION_GFF2 $INTERMEDIATE_DIR/STAR_QC_BALB_CJ/$sample/ &

cd $INTERMEDIATE_DIR/STAR_Human/$sample
#java -Xmx50G -jar ~/git/QoRTs/
QoRTs QC --generatePlots --stranded --generatePlots --numThreads $THREADS --maxReadLength 141 --minMAPQ 20 --title ${sample}_Human \
--addFunctions JunctionCalcs,writeJunctionSeqCounts,writeSpliceExon,writeNovelSplices,annotatedSpliceExonCounts,makeAllBrowserTracks,FPKM,writeGeneCounts,writeDESeq,writeDEXSeq,writeGeneBody,writeGenewiseGeneBody,calcDetailedGeneCounts,makeJunctionBed,writeDocs,writeGeneBodyIv --chromSizes $HUMAN_SIZE $INTERMEDIATE_DIR/STAR_Human/$sample/${sample}.sort.bam \
$HUMAN_ANNOTATION_GFF2 $INTERMEDIATE_DIR/STAR_QC_Human/$sample/ &

wait " >Sbatch_QoRTs_${sample}.script
echo "sbatch Sbatch_QoRTs_${sample}.script" >>run_QoRTs.sh
done<$SAMPLES
###################################################################
#5. Featurecount                                                  #       
###################################################################
rm -f run_featureCounts.sh
#GTF_FILE=$PROJECT_DIR'data/meta/annotation/Merged_annotation_Human_BALB_cJ.gtf'
THREADS=20
#S1
echo "$SBATCH
#SBATCH -p node
#SBATCH -t 24:00:00
#SBATCH -J featurecounts_S1

featureCounts_path=\"$INTERMEDIATE_DIR/featureCounts/Disambiguate_STAR_BALB_CJ_s1\"
output=\"count-s-1\"
annotation=\"$GTF_FILE\"
mkdir -p \$featureCounts_path
~/glob/Software/subread-2.0.0-source/bin/featureCounts -s 1 -t exon -g gene_id -F GTF \
-C -T $THREADS -p -B --primary -Q 20 \
-o \$output \
-a \$annotation \\" >Sbatch_featurecounts_Disambiguate_STAR_s1.script
while read -r sample R1 R2
do
	echo -n "$INTERMEDIATE_DIR/Disambiguate/${sample}_STAR_BALB_CJ/${sample}.disambiguatedSpeciesA.bam $INTERMEDIATE_DIR/Disambiguate/${sample}_STAR_BALB_CJ/${sample}.disambiguatedSpeciesB.bam " >>Sbatch_featurecounts_Disambiguate_STAR_s1.script
done<$SAMPLES

echo "sbatch Sbatch_featurecounts_Disambiguate_STAR_s1.script ">>run_featureCounts.sh
#S1
echo "$SBATCH
#SBATCH -p node
#SBATCH -t 24:00:00
#SBATCH -J featurecounts_S2

featureCounts_path=\"$INTERMEDIATE_DIR/featureCounts/Disambiguate_STAR_BALBA_CJ_s2\"
output=\"count-s-2\"
annotation=\"$GTF_FILE\"
mkdir -p \$featureCounts_path
cd \$featureCounts_path
~/glob/Software/subread-2.0.0-source/bin/featureCounts -s 2 -t exon -g gene_id -F GTF -C -T $THREADS -p -B --primary -Q 20 \
-o \$output \
-a \$annotation \\" >Sbatch_featurecounts_Disambiguate_STAR_s2.script
while read -r sample R1 R2
do
        echo -n "$INTERMEDIATE_DIR/Disambiguate/${sample}_STAR_BALB_CJ/${sample}.disambiguatedSpeciesA.bam $INTERMEDIATE_DIR/Disambiguate/${sample}_STAR_BALB_CJ/${sample}.disambiguatedSpeciesB.bam " >>Sbatch_featurecounts_Disambiguate_STAR_s2.script
done<$SAMPLES

echo "sbatch Sbatch_featurecounts_Disambiguate_STAR_s2.script ">>run_featureCounts.sh

#Individually S2 BALB_CJ
while read -r sample R1 R2
do
	echo "$SBATCH
#SBATCH -p node
#SBATCH -t 2:00:00
#SBATCH -J ${sample}_S2_BALB_CJ

featureCounts_path=\"$INTERMEDIATE_DIR/featureCounts/Disambiguate_STAR_BALBA_CJ_s2/${sample}\"
output=\"count-s-2_BALB_CJ\"
annotation=\"$GTF_FILE\"
mkdir -p \$featureCounts_path
cd \$featureCounts_path
~/glob/Software/subread-2.0.0-source/bin/featureCounts -s 2 -t exon -g gene_id -F GTF -C -T $THREADS -p -B --primary -Q 20 \
-o \$output \
-a \$annotation \
$INTERMEDIATE_DIR/Disambiguate/${sample}_STAR_BALB_CJ/${sample}.disambiguatedSpeciesA.bam " >Sbatch_featurecounts_Disambiguate_STAR_${sample}_s2_BALB_CJ.script
echo "sbatch Sbatch_featurecounts_Disambiguate_STAR_${sample}_s2_BALB_CJ.script " >>run_featureCounts.sh
done<$SAMPLES

#Individually S2 Human
GTF_FILE='/crex/proj/snic2019-30-25/nobackup/private/UserDirectories/project_template/data/meta/annotation/Human_v32_no_chr.gff3'
while read -r sample R1 R2
do
        echo "$SBATCH
#SBATCH -p node
#SBATCH -t 24:00:00
#SBATCH -J ${sample}_S2_Human
#SBATCH -M Snowy 

featureCounts_path=\"$INTERMEDIATE_DIR/featureCounts/Disambiguate_STAR_BALBA_CJ_s2/${sample}\"
output=\"count-s-2_Human\"
annotation=\"$GTF_FILE\"
mkdir -p \$featureCounts_path
cd \$featureCounts_path
~/glob/Software/subread-2.0.0-source/bin/featureCounts -s 2 -t exon -g gene_id -F GTF -C -T $THREADS -p -B --primary -Q 20 \
-o \$output \
-a \$annotation \
$INTERMEDIATE_DIR/Disambiguate/${sample}_STAR_BALB_CJ/${sample}.disambiguatedSpeciesB.bam " >Sbatch_featurecounts_Disambiguate_STAR_${sample}_s2_Human.script
echo "sbatch Sbatch_featurecounts_Disambiguate_STAR_${sample}_s2_Human.script " >>run_featureCounts.sh
done<$SAMPLES


#Extract multimapped reads
while read -r sample R1 R2
do
        echo "$SBATCH
#SBATCH -p core -n 1
#SBATCH -t 4:00:00
#SBATCH -J ${sample}_Human_multim
module load bioinfo-tools samtools
samtools view -h -f258 $INTERMEDIATE_DIR/STAR_Human/$sample/${sample}.sort.bam | samtools view -bS - >$INTERMEDIATE_DIR/STAR_Human/$sample/${sample}_multimapped.sort.bam
samtools index $INTERMEDIATE_DIR/STAR_Human/$sample/${sample}_multimapped.sort.bam
sbatch Sbatch_featurecount_multimapped_${sample}.script" >Sbatch_multimapped_bam_${sample}.script

THREADS=20
echo "$SBATCH
#SBATCH -p core -n 4
#SBATCH -t 2:00:00
#SBATCH -J ${sample}_FC_Human_multim

featureCounts_path=\"$INTERMEDIATE_DIR/featureCounts/Multi_mapped_Human/$sample/\"
output=\"count-s-2\"
annotation=\"$GTF_FILE\"
mkdir -p \$featureCounts_path
cd \$featureCounts_path
~/glob/Software/subread-2.0.0-source/bin/featureCounts -s 2 -t exon -g gene_id -F GTF -C -T $THREADS -p -B -M \
-o \$output \
-a \$annotation \
$INTERMEDIATE_DIR/STAR_Human/$sample/${sample}_multimapped.sort.bam " >Sbatch_featurecount_multimapped_${sample}.script

done<$SAMPLES

while read -r sample R1 R2
do

echo "$SBATCH
#SBATCH -p core -n 1
#SBATCH -t 12:00:00
#SBATCH -J ${sample}_depth

mkdir -p $INTERMEDIATE_DIR/Depth/Multi_mapped_Human/${sample}
module load bioinfo-tools samtools
samtools depth -g 258 \
-b /crex/proj/snic2019-30-25/nobackup/private/UserDirectories/project_template/data/meta/reference/Human_GRCh38_1k_window_non_coding.bed \
$INTERMEDIATE_DIR/STAR_Human/$sample/${sample}_multimapped.sort.bam >$INTERMEDIATE_DIR/Depth/Multi_mapped_Human/$sample/${sample}_multimapped_non_coding_1kb.depth " >Sbatch_depth_multimapped_${sample}.script

done<$SAMPLES




##bbduk_Human for rRNA

echo "$SBATCH
#SBATCH -p node
#SBATCH -t 16:00:00
#SBATCH -J bbduk_Human 
module load bioinfo-tools bbmap" >Sbatch_bbduk_Human.script

cntr=0
#ribo_reference=$PROJECT_DIR'data/meta/annotation/BALB_CJ/Mus_musculus_balbcj.BALB_cJ_v1.99.chr_rRNA.fa'
#ribo_reference=$PROJECT_DIR'data/meta/annotation/Human_annotation_rRNA.fa'
ribo_reference=$PROJECT_DIR'data/meta/annotation/Human_Ensembl_rRNA.fa'
while read -r sample R1 R2
do
        if [[ $cntr -lt 8 ]]
        then
                echo "mkdir -p $INTERMEDIATE_DIR/bbduk_Human/$sample/
cd $INTERMEDIATE_DIR/bbduk_Human/$sample/
bbduk.sh \
in=$TRIMMED_READS/$sample/$R1 \
in2=$TRIMMED_READS/$sample/$R2 \
outm=ribo.fa ref=$ribo_reference \
qhist=qhist \
aqhist=aqhist \
lhist=lhist \
gchist=gchist \
rpkm=rpkm \
stats=stats \
1>err 2>${sample}_bbduk_Human.txt &" >>Sbatch_bbduk_Human.script
                 cntr=$((cntr+1))
        else
                echo "wait
mkdir -p $INTERMEDIATE_DIR/bbduk_Human/$sample/
cd $INTERMEDIATE_DIR/bbduk_Human/$sample/
bbduk.sh \
in=$TRIMMED_READS/$sample/$R1 \
in2=$TRIMMED_READS/$sample/$R2 \
outm=ribo.fa ref=$ribo_reference \
qhist=qhist \
aqhist=aqhist \
lhist=lhist \
gchist=gchist \
rpkm=rpkm \
stats=stats \
1>err 2>${sample}_bbduk_Human.txt &" >>Sbatch_bbduk_Human.script
                cntr=0
        fi
done<$SAMPLES
echo "wait" >>Sbatch_bbduk_Human.script


#Report bbduk_Human results in a table
rm -f Parse_bbduk_Humant.sh
while read sample  R1 R2
do
        echo "cd $INTERMEDIATE_DIR/bbduk_Human/$sample/
python $CODE/parse_bbduk.py -file ${sample}_bbduk_Human.txt" >>Parse_bbduk_Humant.sh
done<$SAMPLES

#while read sample
#do
#  echo "bam_stat.py -i $sample &
#clipping_profile.py -i $sample -o $sample -s \"PE\" & "
#  for ann in $ANNOTATIONS
#  do
#    echo "
#geneBody_coverage.py -r $ann -i $sample -o $sample &
#infer_experiment.py -r $ann -i $sample &
#inner_distance.py -r $ann -i $sample -o $sample &
#junction_annotation.py -r $ann -i $sample -o $sample &
#junction_saturation.py -r $ann -i $sample -o $sample &
#read_distribution.py -r $ann -i $sample &
#read_duplication.py -i $sample -o $sample &
#read_GC.py -i $sample -o $sample &
#read_NVC.py -i $sample -o $sample &
#read_quality.py -i $sample -o $sample &
#wait"
#  done
#done<samples.txt
