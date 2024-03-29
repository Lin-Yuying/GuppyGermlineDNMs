# This script align reads and generate genotype files for guppy DNM project.
# Author: Y.Lin 

######################################################
#1. index the genome
# GATL, bwa, samtools
######################################################
ref=$1 #QuH-F25.ref.fa
# gatk index
java1.8 -Xmx64G -jar gatk-package-4.1.2.0-local.jar CreateSequenceDictionary \
                -R ${ref}
# bwa index
bwa index ${ref}
# samtools 
samtools faidx ${ref}


######################################################
#2. bwa, align reads to reference genome
######################################################
prefix=$2
bwa mem -M -R "@RG\tPL:Illumina\tID:${prefix}\tSM:${prefix}" \
		   ${ref} \
		   ${prefix}.r1.fq.gz \
		   ${prefix}.r2.fq.gz 1> ${prefix}.sam 2> ${prefix}.samlog

#sam2bam
samtools view -bS ${prefix}.sam > ${prefix}.bam


######################################################
#3. sort, mark duplications
######################################################
java1.8 -jar "-Djava.io.tmpdir=/tmp" picard.jar SortSam \
	I=${prefix}.bam \
	O=${prefix}.sorted.bam \
	SORT_ORDER=coordinate \
	CREATE_INDEX=true

java1.8 -jar "-Djava.io.tmpdir=/tmp/" picard.jar MarkDuplicates \
	I=${prefix}.sorted.bam \
	O=${prefix}.mark_duplicates.bam \
	M=${prefix}.marked_dup_metrics.txt

java1.8 -jar picard.jar BuildBamIndex \
      I=${prefix}.mark_duplicates.bam

# average depth -> outDep.txt
samtools depth ${prefix}.mark_duplicates.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' > ${prefix}.dep.txt 



