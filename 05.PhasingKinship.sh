#1. kinship coefficent
fam=$1
#family-based genotyping
ref=$2 #QuH-F25.ref.fa


gatk CombineGVCFs \
        -R ${ref} \
        --variant indv1.g.vcf.gz \
        --variant indv2.g.vcf.gz \
        ... \
        #other indvs \ 
        --variant indv12.g.vcf.gz \
        -O ${fam}.all.g.vcf.gz

gatk GenotypeGVCFs \
        -R $ref \
        -V ${fam}.all.g.vcf.gz \
        -O ${fam}.all.out.vcf.gz

# hard filtering 
gatk VariantFiltration \
    -V  ${fam}.all.out.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O ${fam}.all.out.filtered.vcf.gz

# kinship 
vcftools --gzvcf ${fam}.all.out.filtered.vcf.gz --plink --out ${fam}
plink --file ${fam} --make-bed --out ${fam}
king -b ${fam}.bed --kinship --prefix ${fam}

# phasing
whatshap phase -o ${fam}.phased.vcf \
				--reference=${ref} \
				${fam}.all.out.filtered.vcf.gz \
				indv1.bam indv1.bam ... indv12.bam 2> FR26.whap.log 