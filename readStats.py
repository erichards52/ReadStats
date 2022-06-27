#!/resources/tools/apps/software/lang/Python/3.7.4-GCCcore-8.3.0/bin/python
# Import libraries
import pysam
import numpy as np
import pandas as pd
import math
import sys
import os
import csv
from collections import defaultdict
from sqlitedict import SqliteDict
from math import log
from pathlib import Path
import pickle
from sqlitedict import decode

bamName=sys.argv[1]
bamBase = os.path.basename(bamName)
bamBase = bamBase.replace(".bam","_")
print(bamBase)
path = Path(bamName)
svPath = path.parent.absolute()
path = svPath.parent.absolute()

#sys.stdout = open(str(path) + '/quickAlStats/' + bamBase + 'quickAlStats_read_arrays.txt', 'w')

# Define variables
readLenList=[]
readNames=[]
readSeqs=[]
percentIdent=0.0
qualityListPrim=[]
keyNames=[]
d = SqliteDict("read_stats_upgrade.sqlite", tablename="reads")
dqual = SqliteDict("read_quals_upgrade.sqlite")
dpqual = SqliteDict("prim_read_quals_upgrade.sqlite")
temp_quality_vals=[]
ptemp_quality_vals=[]
accuracy=0
accuracy_SNV=0
accuracy_Del=0
accuracy_Nn=0
accuracy_Ins=0
nn=0
nm=0
ins=0
dels=0
matches=0
mq=0
suppleBps=0
supplealignedbps=0
primaryBps=0
primaryalignedbps=0
mqPrim=0
temp_suppleBps=0
temp_supplealignedbps=0

def errs_tab(n):
    return [10**(q / -10) for q in range(n+1)]

tab=errs_tab(128)

# Import BAM file
bamfile = pysam.AlignmentFile(bamName, "rb")
    
# Look for unmapped reads (ONLY WORKS IN THIS LOOP W/ PYSAM UNTIL_EOF
# for r in bamfile.fetch(until_eof=True):
 #    if r.is_unmapped:
  #       unmappedReads += 1
   #      unmappedbps += r.query_length                        

# Second read in bamfile to return all relevant variables
for read in bamfile.fetch():
    accuracy=0
    accuracy_SNV=0
    accuracy_Del=0
    accuracy_Nn=0
    accuracy_Ins=0
    nn=0
    nm=0
    ins=0
    dels=0
    matches=0
    mq=0
    suppleBps=0
    supplealignedbps=0
    primaryBps=0
    primaryalignedbps=0
    mqPrim=0
    temp_suppleBps=0
    temp_supplealignedbps=0
    temp_mq=0
    if read.is_unmapped:
        pass
    else:
        key = read.query_name
#        if d.__contains__(key):
#        if key in d.keys():
        if key in d:
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
            sum_prob = 0.0
            temp_mq = -10 * log(sum([tab[q] for q in dqual[key]]) / len(dqual[key]), 10)
            if read.is_unmapped:
                pass

            elif read.is_supplementary:
                temp_suppleBps = read.query_length + d[key][12]
                temp_supplealignedbps = read.query_alignment_length + d[key][13]
                d[key]= tempReadLen, temp_accuracy,temp_accuracy_SNV,temp_accuracy_Del,temp_accuracy_Ins,temp_accuracy_Nn,temp_nn,temp_nm,temp_ins,temp_dels,temp_matches,temp_mq,temp_suppleBps,temp_supplealignedbps,d[key][14],d[key][15],d[key][16],d[key][17]
                
            elif read.is_secondary:
                pass
            else:
                primaryBps = read.query_length + d[key][14]
                primaryalignedbps = read.query_alignment_length + d[key][15]
                qualityPrim = read.query_qualities
                dpqual[key]=read.query_qualities
                sum_prob = 0.0
                mqPrim = -10 * log(sum([tab[q] for q in dpqual[key]]) / len(dpqual[key]), 10)
                d[key]=tempReadLen, temp_accuracy,temp_accuracy_SNV,temp_accuracy_Del,temp_accuracy_Ins,temp_accuracy_Nn,temp_nn,temp_nm,temp_ins,temp_dels,temp_matches,temp_mq,d[key][12],d[key][13],primaryBps,primaryalignedbps,mqPrim,read.mapping_quality
                print(d[key])
        else:
#            dqual.setdefault(key, [])
#            dpqual.setdefault(key, [])   
                     
            readLen=read.query_length
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
            mq = -10 * log(sum([tab[q] for q in dqual[key]]) / len(dqual[key]), 10)
            if read.is_unmapped:
                pass
        
            elif read.is_supplementary:
                suppleBps = read.query_length
                supplealignedbps = read.query_alignment_length
                d[key]=readLen,accuracy,accuracy_SNV,accuracy_Del,accuracy_Nn,accuracy_Ins,nn,nm,ins,dels,matches,mq,suppleBps,supplealignedbps,primaryBps,primaryalignedbps,mqPrim,read.mapping_quality
            elif read.is_secondary:
                pass
            else:
                primaryBps = read.query_length
                primaryalignedbps = read.query_alignment_length
                qualityPrim = read.query_qualities
                dpqual[key]=qualityPrim
                sum_prob = 0.0
                mqPrim = -10 * log(sum([tab[q] for q in dpqual[key]]) / len(dpqual[key]), 10)
                d[key]=readLen,accuracy,accuracy_SNV,accuracy_Del,accuracy_Nn,accuracy_Ins,nn,nm,ins,dels,matches,mq,suppleBps,supplealignedbps,primaryBps,primaryalignedbps,mqPrim,read.mapping_quality
                print(d[key])

with open (str(path) + bamBase + '_per_read_stats.csv','w', newline='') as csv_file:
    csv.writer(csv_file).writerows([k, *v] for k,v in d.items())

d.close()
