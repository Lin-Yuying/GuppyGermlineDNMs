# This script reconstructs reference genome for guppy DNM project.
# Author: Y.Lin


######################################################
#1. longranger
#https://support.10xgenomics.com/genome-exome/software/pipelines/latest/using/wgs 
######################################################
longranger wgs \
--reference=refdata-Guppy_superscaffold \
--id=QuH-F25 \
--fastqs=/path/QuH-F25_Chromium/SI-GA-E12_1,/path/QuH-F25_Chromium/SI-GA-E12_2,/path/QuH-F25_Chromium/SI-GA-E12_3,/path/QuH-F25_Chromium/SI-GA-E12_4 \
--sex=f \
--vcmode=gatk:gatk-package-4.0.3.0-local.jar 


######################################################
#2.running last, remove redundant sequence. 
#https://gitlab.com/mcfrith/last/-/tree/main 
######################################################
SEQ=QuH-F25
lastdb -uNEAR -c ${SEQ} ${SEQ}.fa

lastal -C2 ${SEQ} ${SEQ}.fa|last-split -fMAF+ > ${SEQ}_${SEQ}.maf

last-split ${SEQ}_${SEQ}.maf | last-postmask > ${SEQ}_${SEQ}.2.maf 

python3 02.dedupRef.py -i ${SEQ}_${SEQ}.2.maf -g ${SEQ}.dedup.fa


######################################################
#3. ARCS, k-mer based approach, choose the best k = 40
#https://github.com/bcgsc/arcs
######################################################
for i in {40,60,80,100}
do
arcs-make arcs draft=QuH-F25.dedup reads=QuH-F25 k=${i} t=10 -o QuH-F25.dedup.k${i}
done


######################################################
#4. RagTag -> QuH-F25_ref.fa
#https://github.com/malonge/RagTag
######################################################
SCORE=0.3
GAP=100
CPUS=6

ragtag.py scaffold -s ${SCORE} -t ${CPUS} -g ${GAP} EnsembleGenome.fa QuH-F25.dedup.k40.scaffolds.fa -o ./


######################################################
#5. BUSCO
# https://busco.ezlab.org/
######################################################
#cmd, auto lineage selection
mode=genome
lineage=Cyprinodontiformes_odb10 #eukaryota_odb10
busco -i QuH-F25_ref.fa --auto-lineage -m ${mode} -l ${lineage}


######################################################
#6. QUAST
# https://busco.ezlab.org/
######################################################
quast.py QuH-F25_ref.fa 

