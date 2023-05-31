#!/usr/bin/env python3 
# -*- coding: UTF-8 -*-

'''
This script will filter SNPs based on following criteria:

Done! **(1) DNMs: loci in which genotypes are heterozygous in parents but is homozygous in child.**
Done! **(2) biallelic loci.**
Done! **(3) No missing data.**
Done! **(4) Hard filtering following GATK best practice pipeline (done with GATK).**
Done! **(5) Remove SNPs in low complexity genomic regions.** 
            RepeatModeler & Repeatmasker identified low complexity genomic regions and filter SNPs in those regions. 

**(6) Individual site depth filtering.**
remove SNPs that are with read depth < 0.5x or > 2x of average read depth in the indv.
#Note: not using the following one. 
#Two-sided Poisson test, assuming lambda is equal to mean read depth of the individual, 
#assesses if the read depth of a given locus is deviated from the mean read depth of that individual,  
#The cut-off p-value was set as 2 x 10-4.
#For each trio, all three members need to pass this criteria.

**(7) AB allelic balance.**
False discovery rate could arise when heterozygous loci of parents are miscalled as homozygous. 
Thus, loci in parents were required to have AD = 0. 
We performed 2-sided binomial test on heterozygous site in offspring, 
under null hypothesis, freq(A/B) = 0.5.
SNPs with p-value > 0.05 will be marked as 'AB0.05'.

# ped file needed

Usage: python3 SNPfiltering.py [in.vcf.gz] [outDep.txt] [pedigree.ped]
Author: Y.Lin
Date: May 16, 2022
Update: Jun 30, 2022
Update V2: Aug 9, 2022
Update V3: Sep 1, 2022

'''

from scipy import stats
import sys, os
import gzip
import vcfpy

def main(file=sys.argv[1], DepFile=sys.argv[2], PedFile=sys.argv[3]):
	try:
		reader = vcfpy.Reader.from_path(file)
		#print(reader.header)
	except:
		sys.exit('Cannot find VCF file')

	#name of paternal, maternal, child
	ped_list = getPED(PedFile)
	for i in reader.header.samples.names:
		try:
			if ped_list[i]:
				c = i
				p, m = ped_list[i]
		except:
			pass

	# sample column
	c_col = reader.header.samples.names.index(c)
	p_col = reader.header.samples.names.index(p)
	m_col = reader.header.samples.names.index(m)
	

	writer = vcfpy.Writer.from_path('./{}.AB_DP.vcf'.format(c), reader.header)
	writer_PASS_only = vcfpy.Writer.from_path('./{}_pass.AB_DP.vcf'.format(c), reader.header)

	for record in reader:
		# SNP quality
		if record.QUAL < 100:
			record.add_filter("LowQual")

		# Biallele
		if len(record.ALT) != 1:
			record.add_filter("BiAlleles")
		
		# SNV
		if not record.is_snv():
			continue

		# genotype
		pGT = [call.data.get('GT') for call in record.calls][p_col]
		mGT = [call.data.get('GT') for call in record.calls][m_col]
		cGT = [call.data.get('GT') for call in record.calls][c_col]
		# allele depth
		pAD = [call.data.get('AD') for call in record.calls][p_col]
		mAD = [call.data.get('AD') for call in record.calls][m_col]
		cAD = [call.data.get('AD') for call in record.calls][c_col]
		# Alt and ref allele depth
		pDP = [call.data.get('DP') for call in record.calls][p_col]
		mDP = [call.data.get('DP') for call in record.calls][m_col]
		cDP = [call.data.get('DP') for call in record.calls][c_col]

		# Read Depth
		depDict = getAvgDep(DepFile)
		if ReadDep(int(pDP), depDict[p]) and ReadDep(int(mDP), depDict[m]) and ReadDep(int(cDP), depDict[c]):
			record = record
		else:
			record.add_filter("PoissonDep")

		# AB005 
		c_is_het = [call.is_het for call in record.calls][c_col]
		p_is_het = [call.is_het for call in record.calls][p_col]
		m_is_het = [call.is_het for call in record.calls][m_col]

		if ABTest(c_is_het, cAD):
			record = record
		else:
			record.add_filter("cAB005")

		if ABTest(p_is_het, pAD):
			record = record
		else:
			record.add_filter("pAB005")

		if ABTest(m_is_het, mAD):
			record = record
		else:
			record.add_filter("mAB005")


		# output all SNPs		
		writer.write_record(record)

		# output PASS only 
		if record.FILTER == ["PASS"]:
			#print(record)
			writer_PASS_only.write_record(record)
		else:
			continue


def getPED(PedFile):
	'''
	read DNG file and return a dictionary containing info of parents for each child
	{child:[father, mother]}
	'''
	ped_list = {}
	with open(PedFile,'r') as inf:
		for line in inf:
			if line.startswith("#"):
				pass
			elif line.split():
				child, father, mother, *_ = line.split()
				if father != '.' and mother != '.':
					ped_list[child] = [father, mother]
	return ped_list


def getAvgDep(DepFile):
	'''
	Read dep.txt file and return a dictionary = {'indv' : Dep}.
	'''
	avgDepDict = {}
	with open(DepFile, 'r') as inf:
		for line in inf:
			line = line.rstrip()
			indv, dep = line.split('\t')
			avgDepDict[indv] = float(dep)
	return avgDepDict

def ReadDep(dep, avgDep):
	'''
	#Two-sided Poisson distribution.
	0.5x < Read Depth < 2x
	'''
	pval = stats.poisson.pmf(int(dep), int(avgDep))
	#check pmf, prob mass func
	#if pval >= 2 * 1e-4:
	if 2 * int(avgDep) >= int(dep) >= 0.5 * int(avgDep):
		return True
	else:
		return False

def ABTest(is_het,AD):
	'''
	*** if genotype in child is het, then do the AB test, 
	*** else if homo, then the alt allele coverage need to be < 2.
	ABTest(True, AD)
	'''
	if len(AD) == 2:
		refAD, altAD = AD
		if is_het:
			# cumulative distribution function
			pval = stats.binom.pmf(int(refAD), int(refAD) + int(altAD), 0.5)
			if pval >= 0.05:
				return True
			else:
				return False
		else:
			if min(int(refAD), int(altAD)) <= 2:
				return True
			else:
				return False
	else:
		return False


if '__main__':
	main()
