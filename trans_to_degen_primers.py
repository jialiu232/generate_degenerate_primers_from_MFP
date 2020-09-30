#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 27 17:28:55 2020

@author: liujia
This tool takes the output file `prepped.....out` from MetaFunPrimer,
and generates a fasta file for degenerate primers from those single primers
and a summarize file with each degenerate primer and the inexclusive genes that they can target. 

Input: 
    - `prepped.....out`: MetaFunPrimer output

Output:
    - `output.fa`: a fasta file contains the degenrate primers
    - `output.fa.sum`: a file that summarize each degenerate primer length, number of 
        unique genes it targets, and the unique genes it target
Note:
    - The output file `output.fa.sum` summarizes each degenerate primer, the product length it can target,
        the total number of inexlusive genes that it can target, and all the genes that it can target are
        listed below the primer.
    - If you want to choose the top several degenerate primer pairs with the largest gene coverage,
        it's not wise to refer on the total number of inexlusive genes that it can target, 
        and the listed genes below the primer. Since those genes are inexclusive which means
        that top several primers chosen from this file may target overlapped genes
    - RECOMMENDATION for choosing the top several genes:
        Use `output.fa.sum` as an INTERMEDIATE file, and run:
        python choose_top_degen_primers.py -i output.fa.sum -o final_output.fa.sum -n 10
  
Usage:
    python trans_to_degen_primers.py -i prepped.mfpOutput.out -o degen_prepped.mfpOutput.fa

Example: 
    python trans_to_degen_primers.py -i prepped.fungene_9.6_amoA_AOB_1205_unaligned_nucleotide_seqs.fa.primers.out -o degen_prepped.fungene_9.6_amoA_AOB_1205_unaligned_nucleotide_seqs.fa.primers.fa
"""

# Generate command line tool
#import dict
import argparse

parser = argparse.ArgumentParser(description = 'Summarize MetaFunPrimer single primers output to degenerate primers')
parser.add_argument('-i', '--input', type = str, metavar = '', required = True, help = 'Input file: MetaFunPrimer output file with each primer pairs and their targeting genes information')
parser.add_argument('-o', '--output', type=str, metavar='', required=True, help='define the output file name')
#parser.add_argument('-g', '--geneNum', type=int, metavar='', required=False, help='the total number of genes these primers were designed based on')
#TODO: Add an argument of the total number of genes these primers were built based on (e.g., lamp abundant amoA-AOB genes)
args = parser.parse_args()



# Build 4 info dictionaries:
#   f_dic and r_dic contains information of `clusterID: primersIDwithin cluster`;
#   fseq_dic and rseq_dic contains `primerID: sequence`
f_dic = {}
r_dic = {}
fseq_dic = {}
rseq_dic = {}

# some variables to help summarize each degenerate primers target situation
tarDic = {} # key: degenerate primer pair (e.g., C001); value: the total number of genes each primer can target
prim = ""
numb = 0
pg = {} # key: primer pair; value: a string combining all the genes each degenerate primer pair can target

with open(args.input, 'r') as fin:
    for line in fin:
        line = line.strip()
        if line[0] == 'F':
            lli = line.split()
            #print(lli)
            fn = lli[0]     # fn is forward primer in the line
            rn = lli[2]     # rn is reverse primer in the line
            #print(fn)
            #print(rn)
            # Get the cluster ID (key for f_dic and r_dic)
            if len(rn) <= len(fn):
                cn = rn.split('.')[1] + "." + str(lli[4])
            else:
                cn = fn.split('.')[1] + "." + str(lli[4])
            #print(cn)
            # Add the checked fn into values of f_dic under certain clusterID
            if (cn in f_dic) and (fn not in f_dic.get(cn)):
                f_dic[cn] = f_dic.get(cn) + ";" + fn
            elif cn not in f_dic:
                f_dic[cn] = fn
            # Add the checked rn into values of r_dic under certain clusterID
            if (cn in r_dic) and (rn not in r_dic.get(cn)):
                r_dic[cn] = r_dic.get(cn) + ";" + rn
            elif cn not in r_dic:
                r_dic[cn] = rn
            
            fs = lli[1]
            fs = fs.upper()
            rs = lli[3]
            rs = rs.upper()
            if fn not in fseq_dic:
                fseq_dic[fn] = fs
            if rn not in rseq_dic:
                rseq_dic[rn] = rs
            
            # update targeting info for each degenerate primer pair
            plen = lli[4]
            prim = cn + "\tpro_length=" + str(plen) # cn for example: C001
            numb = 0
        elif line[0] == ">":
            if prim not in pg:
                pg[prim] = line
                numb += 1
            else:
                if line not in pg.get(prim):
                    pg[prim] = pg.get(prim) + "\n" + line
                    numb += 1
        else:
            if prim not in tarDic:
                tarDic[prim] = numb
            else:
                tarDic[prim] = int(tarDic.get(prim)) + numb
                
#print(f_dic)
#print(r_dic)
#print(len(fseq_dic))
#print(len(rseq_dic))
#print(tarDic)
#print(pg)

# Sort the dictionary based on values (number of primer pairs a primer can target) in descending order
# The sort_tarDic is not dictionary but tuple 
sort_tarDic = sorted(tarDic.items(), key=lambda x: x[1], reverse = True)
#print(sort_tarDic)

outF1 = open(args.output + ".sum", "w")

for u in range(len(sort_tarDic)):
    outF1.write(sort_tarDic[u][0] + "\t" + str(sort_tarDic[u][1]) + "\n")
    outF1.write(pg[sort_tarDic[u][0]] + "\n")
    
outF1.close()


# helper method to generate a degenerate base based on a list of given different nucleotide bases
# Input: a list of different nucleotide bases (character list)
# Output: a new generated base
def getBase(basel):
    if len(basel) == 2:
        if ("A" in basel) and ("G" in basel):
            newb = "R"
        elif ("C" in basel) and ("T" in basel):
            newb = "Y"
        elif ("G" in basel) and ("C" in basel):
            newb = "S"
        elif ("A" in basel) and ("T" in basel):
            newb = "W"
        elif ("G" in basel) and ("T" in basel):
            newb = "K"
        elif ("A" in basel) and ("C" in basel):
            newb = "M"
    elif len(basel) == 3:
        if ("C" in basel) and ("G" in basel) and ("T" in basel):
            newb = "B"
        elif ("A" in basel) and ("G" in basel) and ("T" in basel):
            newb = "D"
        elif ("A" in basel) and ("C" in basel) and ("T" in basel):
            newb = "H"
        elif ("A" in basel) and ("C" in basel) and ("G" in basel):
            newb = "V"        
    elif len(basel) == 4:
        newb = "N"
    
    return newb
        
    
    





# Generate degenerate primers
# iterate through f_dic
final_fdic = {}
final_rdic = {}

for fkey in f_dic:  # fkey: C001; 
    fli = f_dic[fkey].split(";")    # a list of forward primer IDs under each C001
    #print(fli)
    fsli = []   # generate a list for all the forward degenerate primer sequences under each cluster 
    for fi in fli:
        fp = fseq_dic[fi]
        fsli.append(fp)
#        print(fp)
    #print(fsli)
    fsli.sort(key = len)    # sort the degenerate sequences based on length
    #print(fsli)
    
    # Adjust all the forward degenerate sequence lengths to be the same as the shortest one
    slen = len(fsli[0])
    for m in range(len(fsli)):
        fsli[m] = fsli[m][0:slen]
    #print(fsli)
    
    # Start to generate the each base for the new primers 
    # based on all the sequences in fsli
    newPrim = ""
    for i in range(slen):   # loop through all the positions of forward sequences
        nucl = []
        for j in range(len(fsli)):  # loop through all the forward degenerate primers 
            if fsli[j][i] not in nucl:
                nucl.append(fsli[j][i])
        if len(nucl) == 1:  # if all the base at postion i are the same
            newbase = nucl[0]
        else:   # if there are different base at the same position, use a helper method
            newbase = getBase(nucl)
        newPrim += newbase
    final_fdic[fkey] = newPrim
#print(final_fdic)
        
for rkey in r_dic:  # rkey: C001; 
    rli = r_dic[rkey].split(";")    # a list of reverse primer IDs under each C001
    #print(rli)
    rsli = []   # generate a list for all the reverse degenerate primer sequences under each cluster
    for ri in rli:
        rp = rseq_dic[ri]
        rsli.append(rp)
    #print(rsli)
    rsli.sort(key = len)    # sort the degenerate sequences based on length
    
    # Adjust all the reverse degenerate sequence lengths to be the same as the shortest one
    rslen = len(rsli[0])
    for n in range(len(rsli)):
        rsli[n] = rsli[n][0:rslen]
    #print(rsli)
    
    # Start to generate the each base for the new primers 
    # based on all the sequences in rsli
    newRPrim = ""
    for a in range(rslen):  # loop through all the positions of the reverse degenerate primer
        rnucl = []  # a list contains all the different bases from degenerate primers appear at position a
        for b in range(len(rsli)):  # loop through all the reverse degenerate primers
            if rsli[b][a] not in rnucl:
                rnucl.append(rsli[b][a])
        if len(rnucl) == 1:
            newbaseR = rnucl[0]
        else:
            newbaseR = getBase(rnucl)
        newRPrim += newbaseR
    final_rdic[rkey] = newRPrim
print("Total number of degenerate primers is: " + str(len(final_rdic)))

outF2 = open(args.output, "w")
for key in final_fdic:  # key is C002 for example
    fID = ">F:" + key + ".F"
    fSeq = final_fdic.get(key)
    rID = ">R:" + key + ".R"
    rSeq = final_rdic.get(key)
    #outF2.write(fID + "\n" + fSeq + "\n" + rID + "\n" + rSeq + "\n")
    outF2.write('{}\n{}\n{}\n{}\n'.format(fID, fSeq, rID, rSeq))
    



outF2.close()

            
