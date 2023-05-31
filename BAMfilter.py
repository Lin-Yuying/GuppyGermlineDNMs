'''
1. check directions of reads containing DNM 
2. if there are indels in reads, if any, could be misaligned reads
Author: Y.Lin
'''

import pysam, vcfpy
import sys, os


def main():
	VCF = sys.argv[1]
	PEDfile = sys.argv[2]
	#child, father, mother = getPED(PED)
	#CHROM, POS = getVCFInfo (VCF)

	try:
		reader = vcfpy.Reader.from_path(VCF)
		#print(reader.header)
	except:
		sys.exit('Cannot find VCF file')

	ped_list = getPED(PedFile)
	for i in reader.header.samples.names:
		try:
			if ped_list[i]:
				c = i
				p, m = ped_list[i]
		except:
			pass

	pBAM = p + '.dedup.20k.bam'
	mBAM = m + '.dedup.20k.bam'
	cBAM = c + '.dedup.20k.bam'

	writer = vcfpy.Writer.from_path('./{}.AB_DP_BAM.filter.vcf'.format(c), reader.header)
	writer_PASS_only = vcfpy.Writer.from_path('./{}_pass.AB_DP_BAM.filter.vcf'.format(c), reader.header)


	for record in reader:
		CHROM = record.CHROM
		POS = record.POS
		for read in AlignedRead(cBAM, CHROM, POS):
			if checkGaps(read):
				record = record
			else:
				record.add_filter("cGaps")
			if secMappedReads(read):
				record = record
			else:
				record.add_filter("csecMapped")

		for read in AlignedRead(pBAM, CHROM, POS):
			if checkGaps(read):
				record = record
			else:
				record.add_filter("pGaps")
			if secMappedReads(read):
				record = record
			else:
				record.add_filter("psecMapped")

		for read in AlignedRead(mBAM, CHROM, POS):
			if checkGaps(read):
				record = record
			else:
				record.add_filter("mGaps")
			if secMappedReads(read):
				record = record
			else:
				record.add_filter("msecMapped")

		 
		# output all sites
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


def getVCFInfo(VCF):
	try:
		reader = vcfpy.Reader.from_path(VCF)
		#print(reader.header)
	except:
		sys.exit('Cannot find VCF file')

	for record in reader:
		CHROM = record.CHROM
		POS = record.POS
		
	return CHROM, POS


def AlignedRead(BAM, CHROM, POS):
	'''
	remove duplicate reads and mapQ < 30 reads
	'''
	samfile = pysam.AlignmentFile(BAM,'rb')
	for read in samfile.fetch(CHROM, POS - 1, POS + 1):
		if not read.is_duplicate and read.mapping_quality >= 30 and read.is_proper_pair:
			yield read


def checkGaps(read):
	# check cigar, CIGAR: 3M1I3M1D5M -> 3 matches, 1 insertion (not exist in ref seq), 3 matches, 1 deletion (not exist in quary seq), 5 matches
	# read.cigarstring
	read_count = 0
	gap_count = 0
	cigars = [i for i in read.cigarstring]
	if 'I' in cigars or 'D' in cigars:
		gap_count += 1
		read_count += 1
	else:
		read_count += 1
	return False if gap_count / read_count > 0.5 else True


def secMappedReads(read):
	read_count = 0
	secMapped_count = 0
	if read.is_secondary and read.is_proper_pair:
		secMapped_count += 1
		read_count += 1
	else:
		read_count += 1
	return False if secMapped_count / read_count > 0.5 else True


if "__main__":
	main()










