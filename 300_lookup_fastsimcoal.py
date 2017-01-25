# PJM is merging 6 individs w/ highest coverage (hal, cro, aust. lyr) in 1 vcf file.
# This script will loop through the merged vcf, parse through the lines, 
# and ask if there's a genotype present in at least 4 of the individuals.
# It will ultimately create a lookup key from the merged vcf.

'''
Suggestions from Filip:

input SNP matrix:
- sites present in all indivs of the selected pops., min X-times covered 
- re-polarized
- spits DSFS / (optionally?? pairwise joint_SFS for 2pops, 3pops, 4pops ...)

- non-parametric bootstrap with replacement (??? LD-block bootstrap, as suggested on p. 49 in manual)

input params
- N of bootstraps of the original input matrix (optional)
- min X- depth per site/sample
- re-polarization reference table (optional)
- input VCF
- output DSFS.obs (? optionally multiple "pairwise" joint_SFS)

(into in the script - names of pops, names of indivs, ploidy)


output
- re-polarized SFS 
!! the order of pops is reversed in the SFS (see p. 46 and 47 of the manual)

- with a two-line header:

1 observations. No. of demes and sample sizes are on next line
3	52	32	16 
(n pops, n chromosomes pop1, pop2, pop3)
'''

import os, sys, subprocess, argparse, random, numpy, csv, gzip

parser = argparse.ArgumentParser(description = 'This program takes a merged vcf with highest-coverage individuals from specific pops as input, and outputs a lookup key from the merged vcf')
parser.add_argument('-v', type = str, metavar = 'vcf_path', required = True, help = 'path to directory with vcf(s)')
parser.add_argument('-gz', type = str, metavar = 'gzipped?', required = True, help = 'are vcfs gzipped (true) or not (false)')
parser.add_argument('-o', type = str, metavar = 'Output_Prefix', required = True, help = 'Vcfs retain original name but the concatenated lookup key will be a text file specified by output')
parser.add_argument('-b', type = str, metavar = 'bootstrap?', required = True, help = 'do you want to create an LD-block, non-parametric bootstrap with replacement')

args = parser.parse_args()

if os.path.exists(args.v + '/lookup_table/') == False: #Create folder for output if it doesn't already exist
    os.mkdir(args.v + '/lookup_table/')
vcf_list = []

for file in os.listdir(args.v): #get names of vcf files in args.v directory
    if args.gz == 'true':
        if file[-3:] == '.gz':
            vcf_list.append(file)

    elif args.gz == 'false':   
        if file[-3:] == 'vcf':
            vcf_list.append(file)

    else:
        print 'error'

count = 0

for file in os.listdir(args.v):
	for line in file:

	lookuptempfile = open(args.v+ '/lookup_table/'+args.o+".lookup_table.txt", 'w')
    lookup_table= open(args.v+ '/lookup_table/'+args.o+".lookup_table.txt",'w'). #.par since fastsimcoal uses .par as input file format

    for i, vcf in enumerate(vcf_list):
    	if args.gz=='true':
    		newVCF = open(args.v + '/lookup_table/'+ args.o+ )













