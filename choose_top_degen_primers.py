#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 14:08:11 2020

@author: liujia
This script takes the output file `output.fa.sum` from `trans_to_degen_primers.py`, 
and choose the top n (user defined) best degenerate primers with biggest exclusive
targeting rate against unique genes that can be targeted

Input:
    - `degen_output.fa.sum`: output file from `trans_to_degen_primers.py`
    
Output:
    - `final_degen_output.fa.sum`: summarizes from the best primer, it's product length,
        total number of exclusive genes it can target, and eclusive genes it can target
        are listed below

Usage:
    python choose_top_degen_primers.py -i degen_output.fa.sum -o final_degen_output.fa.sum -n 10


Example:
	python choose_top_degen_primers.py -i degen_4_prepped.fungene_9.6_amoA_AOB_1205_unaligned_nucleotide_seqs.fa.primers.fa.sum -o final_degen_4_prepped.fungene_9.6_amoA_AOB_1205_unaligned_nucleotide_seqs.fa.primers.fa.sum -n 32
"""

import argparse
import sys

## Read in file and build a list with all unique genes that can be targeted by degenerate primers
#   and a dictionary with key (degerate primer ID) and value (unique genes the key primer can target)
def read_file(file):
    genes = []
    prim_dic = {}
    curr_prim = ""
    with open(file, 'r') as f:
        for line in f:          
            line = line.strip()
            if line[0] == "C":
                prim_dic[line.split()[0] + "\t" + line.split()[1]] = []
                curr_prim = line.split()[0] + "\t" + line.split()[1]
            elif line[0] == ">":
                if line not in prim_dic.get(curr_prim):
                    prim_dic[curr_prim].append(line)
                if line not in genes:
                    genes.append(line)
    return genes, prim_dic
                
            
# This is a helper method to get the best primer with the largest coverage
def get_best_primer(prim_dic):
    max = 0
    prim = ""
    for key in prim_dic:
        gn = prim_dic.get(key)  # get the total number of genes primer "key" can target
        gnl = len(gn)
        if max < gnl:
            max = gnl
            prim = key  # define the best primer as key since it targets most of the genes
    return prim
        

def main():
    parser = argparse.ArgumentParser(description = 'Choose the top n primer pairs')
    parser.add_argument('-i', '--input', type = str, metavar = '', required = True, help = 'trans_to_degen_primers.py output file with each primer pairs and their targeting genes information')
    parser.add_argument('-n', '--num', type=int, metavar='', required=True, help='the top n primer pairs that you want')
    parser.add_argument('-o', '--output', type=str, metavar='', required=True, help='define the output file name')
    args = parser.parse_args()
    
    # Read file and generate a list of unique genes and a dictionary with key (primer ID) and value (unique genes the certain primer can target)
    genes, prim_dic = read_file(args.input)
    num_genes = len(genes)
    #print(len(genes))
    #print(prim_dic)
    #for key in prim_dic:
    #    print(len(prim_dic.get(key)))
    
    if args.num > len(prim_dic):
        print("Given number of top primers ("  + str(args.num) + ") is bigger than the total number of primers (" +
               str(len(prim_dic)) + ").")
        sys.exit()
        
    # Get the n (user defined) best primer pairs with the largest coverage
    final_prim = {}
    total_tar = {}
    outF = open(args.output, 'w')
    
    for i in range(args.num):
        # Get the best primer with the largest coverage
        best_prim = get_best_primer(prim_dic)
        #print("Number " + str(i+1) + " best primer is: " + best_prim)
        outF.write("Number " + str(i+1) + " best primer is: " + best_prim + "\n")
        
        bestp_genes = prim_dic.get(best_prim)   # a list of genes that can be targeted by best_prim; already excluded genes that can be targeted by previous best primers 
        #print("Number of exclusive genes that can be targeted by <" + best_prim + "> is: " + str(len(bestp_genes)))
        outF.write("Total number of exclusive genes that can be targeted by <" + best_prim + "> is: " + str(len(bestp_genes)) + "\n")
        # get the number of genes from 'genes' list been targeted by max_prim 
        # and add to final primer list if the number is bigger than 0;
        # update the genes list (by removing the genes that targeted by max_prim from genes list)
        total_tar[best_prim] = []
        num_tar = 0 
        for g in prim_dic.get(best_prim):
            if g in genes:
                genes.remove(g)
                num_tar += 1
                if g not in total_tar[best_prim]:
                    total_tar[best_prim].append(g)

        if num_tar > 0:
            final_prim[best_prim] = num_tar
        
        #print("\n".join(total_tar[best_prim]))
        #print("total_tar[best_prim] number is: " + str(len(total_tar[best_prim])))
        outF.write("\n".join(total_tar[best_prim]) + "\n")
        
        #del prim_dic[best_prim]
        # update the target genes of each primer in prim_dic 
        # (by removing the genes that being targeted by max_prim within each primer targeting genes list)
        for j in prim_dic:
            gene_li = prim_dic.get(j)
            #print("primer " + j + " relatively original coverage: " + str(len(gene_li)))
            #gene_li = [gs for gs in gene_li if gs not in prim_dic.get(best_prim)]
            newli = []
            for gg in gene_li:
                if gg not in bestp_genes:
                    newli.append(gg)
            prim_dic[j] = newli
            #print("eliminate the targeted primers for " + j  + " " + str(len(prim_dic.get(j))) + "\n")
        #del prim_dic[best_prim]
            
    #print(final_prim)
    #print(len(final_prim))
    
    sum_tar = 0
    for t in total_tar:
        sum_tar += len(total_tar[t])
    #print(sum_tar)
    print("The selected top " + str(args.num) + " degenerated primer pairs can target " 
               + str(sum_tar) + " out of " + str(num_genes) + " number of genes.\nTargeting rate: "
               + str(round(sum_tar/num_genes, 5)))
    outF.write("The selected top " + str(args.num) + " degenerated primer pairs can target " 
               + str(sum_tar) + " out of " + str(num_genes) + " number of genes.\nTargeting rate: "
               + str(round(sum_tar/num_genes, 5)))
    
    outF.close()
    
if __name__ == '__main__':
    main()
