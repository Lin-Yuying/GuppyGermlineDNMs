# This script will insert mutations into bam file
# next step run the whole genotypeing and filtering pipeline

python3 bamsurgeon/bin/addsnv.py \
        -v snp.list \
        -f ${prefix}.sorted.bam \
        -r QuH-F25.ref.fa \
        -o ${prefix}.sim.bam \
        --ignorepileup \
        --maxopen 100000 \
        --picardjar picard.jar \
        --aligner mem \
        --seed 1234 2> sim.log

samtools sort --write-index -o ${prefix}.sim.sorted.bam ${prefix}.sim.bam 
samtools index ${prefix}.sim.sorted.bam

