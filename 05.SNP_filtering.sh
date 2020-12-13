'''
####  for non human samples, we are unbale to use VQSR (recalibration) to filter variants ###### 
####### solution: Hard filtering ########
########## filtering stardard ###########
1.	QualByDepth (QD)    > 2
2.	FisherStrand (FS)       > 60
3.	StrandOddsRatio (SOR)    > 3
4.	RMSMappingQuality (MQ)    > 59 & < 61, around 60
5.	MappingQualityRankSumTest (MQRankSum)   +: alt > ref;  - :  ref > alt; 0 is the best 

see: https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
'''
#################### 00. setup #######################
gatk=/Linux/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar

############################ 01. (1)VariantToTable plot ########################
java1.8 -Xmx64G -jar $gatk VariantsToTable \
-V input.vcf \
-F CHROM -F POS -F TYPE -F QD -F FS -F SOR -F MQ -F MQRankSum -F ReadPosRankSum -GF AD \
-O summary.vcf.table

# then visualize the each annotation variants using R and set filtering threshold

##############  01.(2)alternatively, using bcftools for info extraction ############## 
# 1.	QualByDepth (QD)    > 2
#bcftools query -f "%QD\n" out.vcf.gz > my_snp.QD &
# 2.	FisherStrand (FS)       > 60
#bcftools query -f "%FS\n" out.vcf.gz > my_snp.FS &

# 3.	StrandOddsRatio (SOR)    > 3
#bcftools query -f "%SOR\n" out.vcf.gz > my_snp.SOR
# 4.	RMSMappingQuality (MQ)    > 59 & < 61, around 60
#bcftools query -f "%MQ\n" out.vcf.gz > my_snp.MQ
# 5.	MappingQualityRankSumTest (MQRankSum)   +: alt > ref;  - :  ref > alt; 0 is the best 
#bcftools query -f "%MQRankSum\n" out.vcf.gz > my_snp.MQRankSum
# 6.	ReadPosRankSumTest (ReadPosRankSum) , position 
#bcftools query -f "%ReadPosRankSum\n" out.vcf.gz > my_snp.ReadPosRankSum

############## 02. site-level annotation_SelectVariants ###########
# using `INFO` field annotation

# 2.1 SNPs-only callset
gatk SelectVariants \
    -V cohort.vcf.gz \
    -select-type SNP \
    -O snps.vcf.gz
# 2.2 indels-only callset 
gatk SelectVariants \
    -V cohort.vcf.gz \
    -select-type INDEL \
    -O indels.vcf.gz
    
# 2.3 hard filter SNPs on multi-expressions using VariantFiltration
# provide each expression separately with the `-filter` parameter followed by the `–-filter-name`
gatk VariantFiltration \
    -V snps.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O snps_filtered.vcf.gz
    
# site-level evaluation ???
#gatk CollectVariantCallingMetrics \
#    -I filtered.vcf.gz \
#    --DBSNP Homo_sapiens_assembly38.dbsnp138.vcf \
#    -SD Homo_sapiens_assembly38.dict \
#    -O metrics

############## 3. genotype filtering   ###########################
#e.g. heterozygous genotype call, single sample or multi-sample callsets, using `FORMAT` field
###### 3.1 annotate genotype using VariantFiltaration ###
#gatk VariantFiltration \
#-V trio.vcf \
#-O trio_VF.vcf \
#--genotype-filter-expression "isHet == 1" \
#--genotype-filter-name "isHetFilter"
#output
#GT:AD:DP:FT:GQ:PL 0/1:17,15:32:isHetFilter:99:399,0,439 0/1:11,12:23:isHetFilter:99:291,0,292 1/1:0,30:30:PASS:90:948,90,0

#### 3.2 tranform filtered genotype to no call 
#gatk SelectVariants \
#-V trio_VF.vcf \
#--set-filtered-gt-to-nocall \
#-O trioGGVCF_VF_SV.vcf
# output
# GT:AD:DP:FT:GQ:PL ./.:17,15:32:isHetFilter:99:399,0,439   ./.:11,12:23:isHetFilter:99:291,0,292   1/1:0,30:30:PASS:90:948,90,0

########### 4. sample-level filtering, using `GenotypeConcordance` ################
#To be finished



