# PJM is merging 6 individs w/ highest coverage (hal, cro, aust. lyr) in 1 vcf file.
# This script will loop through the merged vcf, parse through the lines, 
# and ask if there's a genotype present in at least 4 of the individuals. Then take sites from these.
# It will ultimately create a lookup key from the merged vcf.

#loop over /Repolarization/ vcf 

import os, sys, subprocess, argparse, random, numpy, csv, gzip

parser = argparse.ArgumentParser(description = 'This program takes a merged vcf with highest-coverage individuals from specific pops as input, and outputs a lookup key from the merged vcf')
parser.add_argument('-v', type = str, metavar = 'vcf_path', required = True, help = 'path to vcf')
parser.add_argument('-gz', type = str, metavar = 'gzipped?', required = True, help = 'are vcfs gzipped (true) or not (false)')
parser.add_argument('-o', type = str, metavar = 'Output_Prefix', required = True, help = 'Vcfs retain original name but the concatenated lookup key will be a text file specified by output')

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

with open(args.v) as vcf:
    for line_idx, line in enumerate(vcf): #Cycle over lines in the VCF file
        cols = line.replace('\n', '').split('\t')  #Split each line of vcf
        if len(cols) < 2:               ## This should be info just before header
            pass
        elif cols[0] == "#CHROM": #This should be header
            
            # for i in range(len(cols)):
            #     print(i,cols[i])         # printing the header name & its index position. Won't keep in final code; just to help build it
            
            newVCF.write(line)
            for j in cols[9:]: #get names of individuals in vcf
                names.append(j)

        else: #
            scaff = cols[0]               #parse important info from each line     
            position = int(cols[1])
            ref_base = cols[3]
            alt_base = cols[4]      # parsing all these things
            info = cols[7].split(";")
            AN = float(info[2].split("=")[1])
            AC = float(info[0].split("=")[1])
            newsite=[]if line_num < 10:             # trouble shooting; look at first 10 lines of code. remove later. make it larger to start.
            num_ind = 0
            alt_ind = 0
            for ind in cols[9:]:
                ind = ind.split(":")
                dp = ind[2]
                gt = ind[0]
                gt = gt.split("/")
                if gt[0] != ".":
                    num_ind += 1
                    if sum([int(x) for x in gt]) == 2:
                        alt_ind += 1

            if num_ind >= min_ind and float(alt_ind)/float(num_ind) >= min_prop_alt:
                out.write(scaff + "\t" + pos + "\n")














