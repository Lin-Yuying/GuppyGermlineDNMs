#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
This script will filter short reads with 95% overlaps in length mapping to longer reads. 
And if the 
Usage: python dedup.py -i [QuH-F25_QuH-F25.maf] -g [genome.fa] 
Author: Yuying Lin 
Date: Mar 1st, 2022
'''
import argparse
import sys
from collections import defaultdict
from Bio import SeqIO
from time import sleep

def parseMAF(maf_file):
	'''
    Alignment Block Lines (lines starting with "a" -- parameters for a new alignment block)
    "i" -- before and after this block in the aligning species.
    "s" -- sequence alignment
    "q" -- quality
    check http://genome.ucsc.edu/FAQ/FAQformat.html#format5
    '''
	linestarts = ['a','s','q','i','e'] 
	with open(maf_file, 'r') as inf:
		block = []
		for line in inf:
			line = line.strip() 
			if line.startswith('#'):
				pass
			elif line.startswith('a'):
				if block:
					yield block
				block = [line] # reset block if
			elif line.startswith('s'):
				block.append(line)
			else: # if line is empty, continue
				continue
		yield block # last block

def parseBlock(block) -> 'list':
	'''
	read each block from MAF file；
	# a score=13458 mismap=1e-09,
	# s 1 0 2243 + 2243 ATTTTAattttactt...
	# s 1 0 2243 + 2243 ATTTTAattttactt...
	return 
	[[chrR, startR, endR, sizeR],[chrQ, startQ, endQ, sizeQ]]
	# [['1', '0', 2243, '2243'], ['1', '0', 2243, '2243']]
	'''
	header, Ref, Qry = block
	# filter header later
	_, chrR, startR, mapSizeR, strandR, sizeR, *_ = Ref.split()
	endR = str(int(startR) + int(mapSizeR))
	_, chrQ, startQ, mapSizeQ, strandQ, sizeQ, *_ = Qry.split()
	endQ = str(int(startQ) + int(mapSizeQ))
	# we assume the longer seqs are the reference seqs
	if int(sizeR) > int(sizeQ):
		return [[chrR, startR, endR, mapSizeR, sizeR],[chrQ, startQ, endQ, mapSizeQ, sizeQ]]
	elif int(sizeR) == int(sizeQ):
		if int(chrR) < int(chrQ):
			return [[chrR, startR, endR, mapSizeR, sizeR],[chrQ, startQ, endQ, mapSizeQ, sizeQ]]
		else:
			return [[chrQ, startQ, endQ, mapSizeQ, sizeQ],[chrR, startR, endR, mapSizeR, sizeR]]
	else:
		return [[chrQ, startQ, endQ, mapSizeQ, sizeQ],[chrR, startR, endR, mapSizeR, sizeR]]

def filterSeqs(blockList):
	'''
	[[chrR, startR, endR, sizeR],[chrQ, startQ, endQ, sizeQ]]
	'''
	splitAln = []
	ref, qry = blockList
	if ref[0] == qry[0]:
		return None
	else:
		if int(ref[3]) >= 0.95 * int(ref[4]):
			return qry[0]
		else:
			if  int(ref[-1]) == int(qry[-1]):
				return qry[0]
			else:
				return None
			#if qry[0] in splitAln

def rmSeqs(genome,rmList):
	for seq_record in SeqIO.parse(genome,'fasta'):
		ids = str(seq_record.id)
		seqs = str(seq_record.seq)
		if ids not in rmList:
			if checkLen(seqs):
				if checkNs(seqs):
					yield '>' + ids + '\n' + seqs
				else:
					continue
			else:
				continue
		else:
			continue


def checkNs(seq):
    if seq.count('N') >= 0.5 * len(seq):
    	return False
    else:
    	return True
    

def checkLen(seq):
	if len(seq) <= 200:
		return True
	else:
		return True


def getArgs():
    '''
    set up parameters
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str,
                        help="Mapping file from Last v")
    parser.add_argument("-o", "--output", type=str,
                        help="Output file ")
    # This checks if the user supplied any arguments. If not, help is printed.
    #if len(sys.argv) != 3:
    #    parser.print_help()
    #    sys.exit(1)
    # Shorten the access to command line arguments.
    parser.add_argument("-g", "--genome", type=str,
                        help="Ref genome")
    args = parser.parse_args()
    return args

# exe the script
if '__main__':
	args = getArgs()
	rmList = [] # remove sequences chr ID
	for block in parseMAF(args.input):
		dedupSeqs = filterSeqs(parseBlock(block))
		if dedupSeqs:
			rmList.append(dedupSeqs)
	#print(rmList, len(rmList))
	for i in rmSeqs(args.genome, rmList):
		print(i)
	


