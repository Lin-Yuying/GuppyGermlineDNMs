# This script generate trio-based genotype files for guppy DNM project.
# Author: Y.Lin 

######################################################
#1. trio-based genotyping, GATK
######################################################
#(1) GATK
father=$1
mother=$2
prefix=$3
ref=$4
gatk CombineGVCFs \
	-R ${ref} \
	--variant ${father}.g.vcf.gz \
	--variant ${mother}.g.vcf.gz \
	--variant ${prefix}.20k.g.vcf.gz \
        --annotation PossibleDeNovo \
        --annotation StrandBiasBySample \
	-O ${prefix}.cohort.g.vcf.gz

gatk GenotypeGVCFs \
        -R $ref \
        -V ${prefix}.cohort.g.vcf.gz \
        --annotation PossibleDeNovo \
        --annotation StrandBiasBySample \
        -O ${prefix}.out.vcf.gz

#SNPs and indels
gatk SelectVariants \
        -V ${prefix}.out.vcf.gz \
        -select-type SNP \
        -O ${prefix}.gatk.snps.vcf.gz

gatk SelectVariants \
        -V ${prefix}.out.vcf.gz \
        -select-type INDEL \
        -O ${prefix}.gatk.indels.vcf.gz    

gatk VariantsToTable \
        -V ${prefix}.snps.vcf.gz \
        -F CHROM -F POS -F TYPE -F QD -F FS -F SOR -F MQ -F MQRankSum -F ReadPosRankSum -GF AD \
        -O ${prefix}.summary.vcf.table

gatk VariantFiltration \
    -R ${ref} \
    -V  ${prefix}.gatk.snps.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -filter "DP < 10" --filter-name "LowDepth" \
    --genotype-filter-expression "GQ < 30" --genotype-filter-name "GQ30" \
    -O ${prefix}.gatk.hardfilter.vcf.gz



######################################################
#2. trio-based genotyping, BCFtools, hard filtering 
######################################################
#BCFTools
bcftools mpileup -a FORMAT/SP,FORMAT/ADF,FORMAT/ADR,INFO/SCR,FORMAT/QS,FORMAT/AD,FORMAT/DP,FORMAT/SCR,FORMAT/NMBZ \
                 --min-MQ 30 --min-BQ 30 \
                 -f ${ref} \
                 -Ou ${prefix}.bam ${father}.bam ${mother}.bam | bcftools call -mv -Oz -o ${prefix}.bcftools.vcf.gz

#bcftools +trio-dnm2 -p ${prefix},${father},${mother} \
#         --dnm-tag DNM:prob -Oz -o ${prefix}.bcftools.vcf.gz \
#         ${prefix}.bcftools.dnm.vcf.gz

bcftools view -e 'QUAL <= 100 || MQBZ < -3 || RPBZ < -3 || RPBZ > 3 || FORMAT/SP > 32 || SCBZ > 3 || TYPE="INDEL" ' ${prefix}.bcftools.vcf.gz |bgzip > ${prefix}.bcftools.hardfilter.vcf.gz
######################################################
#3.individual-level filtering to both dataset
######################################################
# AB allele, Depth, etc.
python SNPfiltering.py ${prefix}.bcftools.hardfilter.vcf.gz outDep.txt pedigree.ped
python SNPfiltering.py ${prefix}.gatk.hardfilter.vcf.gz outDep.txt pedigree.ped
# remove snps with reads with gaps or without properly mapped mate pair
python BAMfilter.py ${prefix}.AB_DP.bcftools.vcf.gz pedigree.ped
python BAMfilter.py ${prefix}.AB_DP.gatk.vcf.gz pedigree.ped


######################################################
#4.intersect two datasets
######################################################
bcftools isec ${prefix}.AB_DP_BAM.bcftools.vcf.gz ${prefix}.AB_DP_BAM.gatk.vcf.gz

