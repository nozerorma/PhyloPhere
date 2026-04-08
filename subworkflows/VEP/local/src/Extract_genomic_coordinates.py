#!/usr/bin/python3

######################################################################
#ALEJANDRO VALENZUELA - RETRIEVAL GENOMIC COORDINATES TEMPLATE SCRIPT#
######################################################################


import pandas as pd
import numpy as np
import sys
from itertools import combinations
#import seaborn as sns
#from scipy import stats
#import pickle
from collections import Counter
import copy
#from scipy.stats import sem, t, chi2_contingency, norm
#from scipy import mean
import re
import os
import gzip
import time
from Bio import AlignIO


"""
Here we create the function for checking the input parameters and saving
them in different variables, error if the usage is not good
"""

if len(sys.argv) == 9:
    alignment_path = sys.argv[1]
    fasta_path = sys.argv[2]
    tracking_file = sys.argv[3]
    gff_file = sys.argv[4]
    caas_file = sys.argv[5]
    output_coordinates = sys.argv[6]
    gene_equivalences = sys.argv[7]
    gene_name = sys.argv[8]
else:
    sys.exit("The usage shoud be: length 8 files")


#align = AlignIO.read(alignment_path, "fasta")
selected_positions = []
removed_positions = []


#FUNCTION A) --> GET INDEX OF FILTERED POSITIONS IN ALIGNMENTS (html files)
with open(tracking_file, "r") as in_fh:
    for line in in_fh:
        line=line.rstrip()
        if "selected:" in line:
            # print(line)
            intervals = line.split(" ")
            for i in range(4,len(intervals)):
                if "-" in intervals[i]:
                    init = int(intervals[i].split("-")[0])
                    end = int(intervals[i].split("-")[1])
                    nums = list(range(init, end+1))
                    for num in nums:
                        selected_positions.append(num)
                else:
                    selected_positions.append(int(intervals[i]))
        elif "removed:" in line:
            # print(line)
            intervals = line.split(" ")
            for i in range(5,len(intervals)):
                if "-" in intervals[i]:
                    init = int(intervals[i].split("-")[0])
                    end = int(intervals[i].split("-")[1])
                    nums = list(range(init, end+1))
                    for num in nums:
                        removed_positions.append(num)
                else:
                    removed_positions.append(int(intervals[i]))



#FUNCTION B.1) --> GET INFORMATION FROM EQUIVALENCE BETWEEN GENE NAMES AND PROTEINS
protein_id = None
with open(gene_equivalences, "r") as in_fh2a:
    print(f"gene_name '{gene_name}' found in {gene_equivalences}")
    target_sequence = False
    aligned_CDS = ""
    for line in in_fh2a:
        line = line.rstrip()
        parts = line.split()
        if parts and parts[0] == gene_name:
            protein_id = parts[1]
            break

if protein_id is None:
    # Gracefully skip when the gene is absent from the equivalence table
    print(f"gene_name '{gene_name}' not found in {gene_equivalences}; skipping", file=sys.stderr)
    sys.exit(0)

#print(transcript_name)

#FUNCTION B) --> GET CORRESPONDING HUMAN CDS POSITIONS ALIGNED

with open(alignment_path, "r") as in_fh2:
    target_sequence = False
    aligned_CDS = ""
    for line in in_fh2:
        line = line.rstrip()
        if line.startswith(">Homo"):
            target_sequence = True
        elif line.startswith(">") and not line.startswith(">Homo"):
                target_sequence = False
        elif not line.startswith(">") and target_sequence:
            aligned_CDS += line

#FUNCTION C) read fasta sequence

#print(fasta_path)
with gzip.open(fasta_path, "rt") as in_fh3:
    fasta_header = protein_id
    #print(fasta_header)
    real_CDS = ""
    for line in in_fh3:
        line = line.rstrip()
        if line.startswith(">"+fasta_header):
            target_sequence = True
        elif line.startswith(">") and not line.startswith(">"+fasta_header):
            target_sequence = False
        elif not line.startswith(">") and target_sequence:
            real_CDS += line

#print(aligned_CDS)
#print(real_CDS)



#FUNCTION D) Get indexes of human sequence based on alignment

counter_nt = 3
counter_codons = 0

prev_ambiguity=False
n_ambiguities=0
ref_aligned_position = 0
indexes = {}
warned_N = False
warned_mismatch = False
for real_position in range(0,len(real_CDS),3):
    codon = real_CDS[real_position:real_position+3]
    # print(codon, "real")
    for aligned_position in range(ref_aligned_position,len(aligned_CDS),3):
        aligned_codon = aligned_CDS[aligned_position:aligned_position +3]
        # print(aligned_codon, "aligned")
        if aligned_CDS[aligned_position] == "-":
            counter_codons += 1
            indexes[counter_codons] = 0
        else:
            if "N" in aligned_codon.upper() and not warned_N:
                warned_N = True
                print(f"WARN {gene_name}: masked codon (N) in aligned human CDS; position {aligned_position//3}", file=sys.stderr)
            if aligned_codon != codon and not warned_mismatch:
                warned_mismatch = True
                print(f"WARN {gene_name}: mismatch aligned {aligned_codon} vs CDS {codon} at codon {aligned_position//3}", file=sys.stderr)
            # Accept mismatches/unknowns: advance counters even if codons differ
            counter_nt += 3
            counter_codons += 1
            break
    ref_aligned_position = 3*counter_codons
    indexes[counter_codons] = counter_nt-3


selected_positions_filt = selected_positions
# Write selected_positions_filt to a file
with open("selected_positions_filt.txt", "w") as in_fh7:
    for i in selected_positions_filt:
        print(i, file=in_fh7)

selected_indexes = indexes

#FUNCTION F) GET coordinate info from gff gff_file

starts = []
ends = []
CDS_positions = []
with open(gff_file, "r") as in_fh4:
    fasta_header = protein_id
    real_CDS = ""
    for line in in_fh4:
        line = line.rstrip()
        if fasta_header in line and line.split("\t")[2] == "CDS":
            #print(line)
            fields = line.split("\t")
            chrom = fields[0]
            start = int(fields[3])
            starts.append(start)
            end = int(fields[4])
            ends.append(end)
            strand = fields[6]
            # print(strand)
    if strand == "+":
        init_coord = starts[0]
        for i in range(0, len(starts)):
            for j in range(starts[i], ends[i]+1):
                #print(j)
                CDS_positions.append(j)
    elif strand == "-":
        init_coord = ends[-1]
        for i in range(0, len(starts)):
            curr_end = ends[-(i+1)]
            curr_start = starts[-(i+1)]
            for j in range(curr_end, curr_start-1, -1):
                #print(j)
                CDS_positions.append(j)

#####################################################
#Translate selected positions to output_coordinates##
#####################################################

#print(CDS_positions)
#print(strand)

with open(caas_file, "r") as in_fh5:
    for line in in_fh5:
        fields=line.rstrip().split("\t")
        if fields[0] == gene_name and int(fields[1]) < len(selected_positions_filt):
            caas_position = int(fields[1])
            print(caas_position)
            prefiltered_position = selected_positions_filt[caas_position]
            if prefiltered_position <= len(selected_indexes):
                prefiltered_nt = selected_indexes[prefiltered_position]
#LAST ELEMENT ADDING
                if prefiltered_nt == len(CDS_positions):
                    if strand == "+":
                        current_coord = CDS_positions[-1]+1
                        real_coord = CDS_positions[-1]
                    elif strand == "-":
                            current_coord = CDS_positions[-1]-1
                            real_coord = CDS_positions[-1]
                    with open(output_coordinates, "a") as in_fh6:
                            print("{}\t{}\t{}\t{}\t{}\t{}".format(fields[0], fields[1], fields[2], chrom, current_coord, real_coord), file=in_fh6)
                            # Change fields from number to column name, note that fields[0] is named as gene_name
#NON-GAPPY filtered positions in reference
                elif prefiltered_nt != 0:
                    #print(prefiltered_nt, type(prefiltered_nt), "NON-O")
                    current_coord = CDS_positions[prefiltered_nt]
                    real_coord = CDS_positions[prefiltered_nt-1]
                    with open(output_coordinates, "a") as in_fh6:
                            print("{}\t{}\t{}\t{}\t{}\t{}".format(fields[0], fields[1], fields[2], chrom, current_coord, real_coord), file=in_fh6)
