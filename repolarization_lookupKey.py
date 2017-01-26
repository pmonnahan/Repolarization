# PJM is merging 6 individs w/ highest coverage (hal, cro, aust. lyr) in 1 vcf file.
# This script will loop through the merged vcf, parse through the lines, 
# and ask if there's a genotype present in at least 4 of the individuals. Then take sites from these.
# It will ultimately create a lookup key from the merged vcf.

#loop over /Repolarization/ vcf 

import os, sys, subprocess, argparse, random, numpy, csv, gzip

parser = argparse.ArgumentParser(description = 'This program takes a merged vcf with highest-coverage individuals from specific pops as input, and outputs a lookup key from the merged vcf')
parser.add_argument('-v', type = str, metavar = 'vcf_path', required = True, help = 'path to vcf')
parser.add_argument('-gz', type = str, metavar = 'gzipped?', required = False, default = 'false', help = 'are vcfs gzipped (true) or not (false)')
parser.add_argument('-o', type = str, metavar = 'Output_Prefix', required = True, help = 'Vcfs retain original name but the concatenated lookup key will be a text file specified by output')
parser.add_argument('-mi', type = int, metavar = 'minimum amount individuals', required = True, help = 'minimum amount of individuals required to support alternative alleles')
parser.add_argument('-mp', type = int, metavar = 'min. proportion alt alleles', required = True, help = 'minimum proportion of alternative alleles to allow')
parser.add_argument('-nm', type = int, metavar = 'number_missing_alleles', required = True, help = 'amount of missing data (integer corresponding to # of alleles) to allow')
parser.add_argument('-ly', type = str, metavar = 'lyrata_only?', required = False, default = 'false', help = 'do you want to include lyrata only (true) or not (false)?')

args = parser.parse_args()

#if os.path.exists('/repolarization_lookup/') == False: #Create folder for output if it doesn't already exist
 #   os.mkdir('/repolarization_lookup/')

if args.gz == 'true' and args.v[-3:] == '.gz':
    gzip.gunzip(args.v)

# I need help below
   # elif args.gz == 'false' and args.v[-3:] == 'vcf':   
lookup_table_file = open(args.v+args.o+"repolarized.lookupKey.minAlleles_"+str(args.mi)+".txt", 'w')
# is the above line necessary, or is that redundant due to line 38 opening args.v as 'vcf' ?

 #   else:
  #      print 'error'

count = 0

with open(args.v) as vcf:
    for line_idx, line in enumerate(vcf): #Cycle over lines in the VCF file
        # Do I need the below line, or is it redundant?
        # newVCF = open(args.v + '/repolarization_lookup/'+vcf[:-3]+"repolarized.lookupKey.minAlleles_"+str(args.mi)+".vcf", 'w')  
        cols = line.replace('\n', '').split('\t')  #Split each line of vcf
        if len(cols) < 2:               ## This should be info just before header
            pass
        elif cols[0] == "#CHROM": #This should be header
            
            # for i in range(len(cols)):
            #     print(i,cols[i])         # printing the header name & its index position. Won't keep in final code; just to help build it
            
            # Help- I'm overthinking this; should I do lookup_table_file.write or one of the next two options?
            lookup_table_file.write(line)
            # vcf.write(line)
            # newVCF.write(line)
            
            names = [] #list to append pop names to (may not be necessary..)
            for j in cols[9:]: #get names of individuals in vcf
                names.append(j)

        else: 
            #if line_num < 10:             # trouble shooting; look at first 10 lines of code. remove later. make it larger to start.
            scaff = cols[0]               #parse important info from each line     
            position = int(cols[1])
            ref_base = cols[3]
            alt_base = cols[4]      # parsing all these things
            info = cols[7].split(";")
            AN = float(info[2].split("=")[1])
            AC = float(info[0].split("=")[1])
            newsite=[]
            num_ind = 0
            alt_ind = 0
            min_ind = args.mi
            min_prop_alt = args.mp

            for ind in cols[9:]:
                ind = ind.split(":")
                dp = ind[2]  # Help- this now prints out dp (correct; I compared to vcf) but breaks after a while and says 'listindex out of range'
                print dp
                print ind
                gt = ind[0]
                gt = gt.split("/")

# below is the 'Lyrata Only' section.

                if args.ly == 'true':
                    if gt[0] != ".":
                        num_ind += 1
                        if sum([int(x) for x in gt]) == 2:
                            alt_ind += 1
                            if gt == '1/1':
                                print 'blah2'
                                #add code to save the info/write to file
# above section needs work; need to think more.

                elif args.ly != 'true':
                    if gt[0] != ".":
                        num_ind += 1
                        if sum([int(x) for x in gt]) == 2:
                            alt_ind += 1

            if num_ind >= min_ind and float(alt_ind)/float(num_ind) >= min_prop_alt:
                #out.write(scaff + "\t" + str(position) + "\n")
                lookup_table_file.write(scaff + "\t" + str(position) + "\n")
                # Help- is above line correct, or should it be vcf.write(.....)
