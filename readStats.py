# Import libraries
import pysam
import numpy as np
import pandas as pd
import math
import sys
import os
import csv
from collections import defaultdict
from math import log
from pathlib import Path

bamName=sys.argv[1]
bamBase = os.path.basename(bamName)
bamBase = bamBase.replace(".bam","_")
print(bamBase)
path = Path(bamName)
svPath = path.parent.absolute()
path = svPath.parent.absolute()

#sys.stdout = open(str(path) + '/quickAlStats/' + bamBase + 'quickAlStats_read_arrays.txt', 'w')

# Define variables
readCount=0
suppleCount=0
secondCount=0
suppleBps=0
primaryReads=0
primaryBps=0
bps=0
unmappedReads=0
unmappedbps=0
alignedbps=0
supplealignedbps=0
readLenList=[]
readNames=[]
readSeqs=[]
qualityList=[]
alignQuals=[]
alignQualsSNV=[]
alignQualsDel=[]
alignQualsIns=[]
alignQualsNn=[]
percentIdent=0.0
qualityListPrim=[]
d = {}
dqual = {}
temp_quality_vals=[]

def errs_tab(n):
    return [10**(q / -10) for q in range(n+1)]

tab=errs_tab(128)

# Import BAM file
bamfile = pysam.AlignmentFile(bamName, "rb")
    
# Look for unmapped reads (ONLY WORKS IN THIS LOOP W/ PYSAM UNTIL_EOF
for r in bamfile.fetch(until_eof=True):
    if r.is_unmapped:
        unmappedReads += 1
        unmappedbps += r.query_length                        

# Second read in bamfile to return all relevant variables
for read in bamfile.fetch():

    if read.is_unmapped:
        pass
    else:
        readCount += 1
        bps += read.query_length
        key = read.query_name
        if d.__contains__(key):
            print("Dict already contains "+key)
            tempReadLen=read.query_length + d[key][0]
            temp_cig_stats_counts = read.get_cigar_stats()[0]
            temp_nn=d[key][6]
            temp_nm=d[key][7] + temp_cig_stats_counts[10]
            temp_ins=d[key][8] + temp_cig_stats_counts[1]
            temp_dels=d[key][9] + temp_cig_stats_counts[2]
            temp_matches=d[key][10] + temp_cig_stats_counts[0]
            try:
                temp_nn += read.get_tag('nn')
            except KeyError:
                temp_nn = 0
            temp_mismatches = temp_nm - temp_dels - temp_ins - temp_nn
            temp_matches_calc = temp_matches - temp_mismatches
            temp_accuracy = temp_matches_calc / (temp_matches_calc + temp_nm) * 100
            # Accuracy only SNVs/mismatches and no deletions and insertions
            temp_accuracy_SNV = float(temp_mismatches) / (temp_matches_calc + temp_nm) * 100
            temp_accuracy_Del = float(temp_dels) / (temp_matches_calc + temp_nm) * 100
            temp_accuracy_Ins = float(temp_ins) / (temp_matches_calc + temp_nm) * 100
            temp_accuracy_Nn = float(temp_nn) / (temp_matches_calc + temp_nm) * 100
            # Calculate quality
            dqual[key].extend(read.query_qualities)
            for l in range(len(dqual[key])):
                temp_quality_vals.append(dqual[key][l])
            sum_prob = 0.0
            if temp_quality_vals:
                temp_mq = -10 * log(sum([tab[q] for q in temp_quality_vals]) / len(temp_quality_vals), 10)
            d[key] = [tempReadLen, temp_accuracy,temp_accuracy_SNV,temp_accuracy_Del,temp_accuracy_Ins,temp_accuracy_Nn,temp_nn,temp_nm,temp_ins,temp_dels,temp_matches,temp_mq]
        else:
            d.setdefault(key, [])
            dqual.setdefault(key, [])
            
            # Calculate accuracy
            cig_stats_counts = read.get_cigar_stats()[0]
            try:
                nn = read.get_tag('nn')
            except KeyError:
                nn = 0
            nm = cig_stats_counts[10]
            ins = cig_stats_counts[1]
            dels = cig_stats_counts[2]
            matches = cig_stats_counts[0]
            mismatches = nm - dels - ins - nn
            matches_calc = matches - mismatches
            accuracy = matches_calc / (matches_calc + nm) * 100
            # Accuracy only SNVs/mismatches and no deletions and insertions
            accuracy_SNV = float(mismatches) / (matches_calc + nm) * 100
            accuracy_Del = float(dels) / (matches_calc + nm) * 100
            accuracy_Ins = float(ins) / (matches_calc + nm) * 100
            accuracy_Nn = float(nn) / (matches_calc + nm) * 100
            # Calculate quality
            quality = read.query_qualities
            dqual[key]=read.query_qualities
            sum_prob = 0.0
            if quality:
                mq = -10 * log(sum([tab[q] for q in quality]) / len(quality), 10)
                d[key].append(mq)
            d[key]=[read.query_length,accuracy,accuracy_SNV,accuracy_Del,accuracy_Nn,accuracy_Ins,nn,nm,ins,dels,matches,mq]
        
    if read.is_unmapped:
        pass
        
    elif read.is_supplementary:
        suppleCount += 1
        suppleBps += read.query_length
        supplealignedbps += read.query_alignment_length
    elif read.is_secondary:
        secondCount += 1
    else:
        primaryReads += 1
        primaryBps += read.query_length
        alignedbps += read.query_alignment_length
        qualityPrim = read.query_qualities
        sum_prob = 0.0
        if qualityPrim:
            mqPrim = -10 * log(sum([tab[q] for q in qualityPrim]) / len(qualityPrim), 10)
            qualityListPrim.append(mqPrim)


with open (str(path) + '/quickAlStats/' + bamBase + 'quickAlStats_read_arrays.csv',"w") as csv_file:
    csv.writer(csv_file).writerow(["Read_Names", "Read_Length","Alignment_Identity","Error_Rate_SNV","Error_Rate_Del","Error_Rate_Nn","Error_Rate_INS","NN","NM","INS","DEL","Matches","Quals"])
    csv.writer(csv_file,delimiter = ",").writerows([k, *v] for k,v in d.items())
