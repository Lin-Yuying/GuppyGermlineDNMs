ref=$1
prefix=$2
# run GATK4 in a BP resolution 
gatk HaplotypeCaller \
-R ${ref} \
-ERC BP_RESOLUTION \
--min-base-quality-score 30 \
--read-filter MappingQualityReadFilter \
--minimum-mapping-quality 30 \
-I ${prefix}.bam \
-XL LGs.repeats.bed \
-O ${prefix}.BP.vcf.gz

# calculate callable genome size 
bcftools query -f "%CHROM\t%POS[\t%SAMPLE\t%AD]\n" ${prefix}.BP.vcf.gz> ${prefix}.stat.txt

python3 callableGenomeSize.py ${prefix}.stat.txt outDep.txt > ${prefix}.call.txt


