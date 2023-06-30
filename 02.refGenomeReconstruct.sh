# This script reconstructs reference genome for guppy DNM project.
# Author: Y.Lin


######################################################
#1. longranger
#https://support.10xgenomics.com/genome-exome/software/pipelines/latest/using/wgs 
######################################################
seq=$1 #QuH-F25
longranger wgs \
--reference=refdata-Guppy_superscaffold \
--id=${seq} \
--fastqs=${seq}_Chromium/SI-GA-E12_1,${seq}_Chromium/SI-GA-E12_2,${seq}_Chromium/SI-GA-E12_3,${seq}_Chromium/SI-GA-E12_4 \
--sex=f \
--vcmode=gatk:gatk-package-4.0.3.0-local.jar 


######################################################
#2.running last, remove redundant sequence. 
#https://gitlab.com/mcfrith/last/-/tree/main 
######################################################

lastdb -uNEAR -c ${seq} ${seq}.fa

lastal -C2 ${seq} ${seq}.fa|last-split -fMAF+ > ${seq}_${seq}.maf

last-split ${seq}_${seq}.maf | last-postmask > ${seq}_${seq}.2.maf 

python dedupRef.py -i ${seq}_${seq}.2.maf -g ${seq}.dedup.fa


######################################################
#3. ARCS, k-mer based approach, choose the best k = 40
#https://github.com/bcgsc/arcs
######################################################
for i in {40,60,80,100}
do
arcs-make arcs draft=${seq}.dedup reads=${seq} k=${i} t=10 -o ${seq}.dedup.k${i}
done


######################################################
#4. RagTag -> QuH-F25.ref.fa
#https://github.com/malonge/RagTag
######################################################
score=0.3
gap=100
cpus=6
ref=$2 #EnsembleGenome.fa
ragtag.py scaffold -s ${score} -t ${cpus} -g ${gap} ${ref}.fa ${seq}.dedup.k40.scaffolds.fa -o ./


######################################################
#5. BUSCO
# https://busco.ezlab.org/
######################################################
#cmd, auto lineage selection
mode=genome
lineage=Cyprinodontiformes_odb10 #eukaryota_odb10
busco -i ${seq}.ref.fa --auto-lineage -m ${mode} -l ${lineage}


######################################################
#6. QUAST
# https://busco.ezlab.org/
######################################################
quast.py ${seq}.ref.fa 

