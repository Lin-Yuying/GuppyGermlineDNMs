'''
This script takes output of 
`bcftools query -f "%CHROM\t%POS\t%AD[\t%SAMPLE\t%AD]\n" P492_2.bcftools.vcf > P492_2.stat.txt`
which has been removed repetitive sequences with `bcftools view -T ^repetitive.txt`.

Author: Y.Lin 
Date: Oct 25, 2022 Version 1
Update: Feb 17, 2023, Version 2
Usage: python callable_single.py [outDep.txt] [avgDep]
'''

import sys,os
import pandas as pd

#def ReadDepCutoff(avgDep): #-> return maxDep and minDep 
def ReadDepCutoff(avgDep): #-> return maxDep and minDep 
    '''
    #Two-sided Poisson distribution.
    0.5x < Read Depth < 2x
    '''
    #pval = stats.poisson.pmf(int(dep), int(avgDep))
    #check pmf, prob mass func
    #if pval >= 2 * 1e-4:
    minDep = 0.5 * avgDep
    maxDep = 2 * avgDep
    return maxDep, minDep

if __name__ == '__main__':
	# input 1: individual depth file output by bcftools
	DepFile = sys.argv[1]
	name = sys.argv[1].split(".")[0]
	#average depth 
	avgDep = 30 #FR26
	# min depth and min depth
	maxDep, minDep = ReadDepCutoff(avgDep)
	# define the chunksize to process
	chunksize=10000
	with pd.read_csv(DepFile, dtype={'DEPTH': 'float'}, chunksize=chunksize, sep = "\t", names=["CHROM","POS","t", "INDV","DEPTH", "_"], na_values='.') as reader:
		'''
        LG1     678     P484    12      
        LG1     679     P484    12      
        LG1     680     P484    12      
        LG1     681     P484    12      
        LG1     1106    P484    51      
        LG1     1107    P484    51      
        LG1     1895    P484    108
		'''
		for chunk in reader:
			df = chunk
			df_na_1 = df.where(df!='.',0)
			
			df_na = df_na_1.drop("_", axis=1)
			#print(df_na)
			#df_na = df.dropna().copy()
			df_dep = df_na[(df_na['DEPTH'] >= int(minDep)) & (df_na['DEPTH'] <= int(maxDep))]
			df_dep.to_csv(name + ".callable.csv", index=False, header=False, mode="a")



