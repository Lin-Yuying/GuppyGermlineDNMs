# This script outputs repeat gff file. 
prefix=$1
output=$2
ref=$3
genomesize=$4

#1. RepeatModeler
BuildDatabase -name ${ref} ${ref}.ref.fa

RepeatModeler -database ${ref} -threads 20 -LTRStruct > run.out &

#2. RepeatMasker
RepeatMasker -lib consensi.classified.fa -dir ${output} -s -a -nolow -html -gff ${ref}
# align2divsum 
# calcDivergenceFromAlign.pl -s ${prefix}.divsum ${prefix}.align
# Repeat Landscape
# createRepeatLandscape.pl -div ${prefix}.divsum -g ${genomesize} > ${prefix}.html