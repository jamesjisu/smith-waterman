#!/usr/bin/python
__author__ = "James Han"
__email__ = "james.han@yale.edu"
__copyright__ = "Copyright 2023"
__license__ = "GPL"
__version__ = "1.0.0"

### Usage: python hw1.py -i <input file> -s <score file>
### Example: python hw1.py -i input.txt -s blosum62.txt
### Note: Smith-Waterman Algorithm

import argparse
import numpy as np
import pandas as pd

### This is one way to read in arguments in Python. 
parser = argparse.ArgumentParser(description='Smith-Waterman Algorithm')
parser.add_argument('-i', '--input', help='input file', required=True)
parser.add_argument('-s', '--score', help='score file', required=True)
parser.add_argument('-o', '--opengap', help='open gap', required=False, default=-2)
parser.add_argument('-e', '--extgap', help='extension gap', required=False, default=-1)
parser.add_argument('-r', '--result', help='result file', required=False, default='output.txt')
args = parser.parse_args()



### Implement your Smith-Waterman Algorithm

def traceback(score_matrix, direction_matrix): ### Traceback algorithm using final scoring and direction matrix
    alignment = []
    traceback_vecs = np.array([[-1,-1], [0,-1], [-1,0]]) # using direction matrix, index of traceback_vecs corresponds to case in direction matrix
    i, j = np.unravel_index(score_matrix.argmax(), score_matrix.shape)
    while i > 0 and j > 0:
        if direction_matrix[i, j] == 3:
            break
        alignment.append([i,j])
        i, j = np.array([i,j]) + traceback_vecs[direction_matrix[i, j]]
    alignment.reverse()
    return alignment ### Returns information necessary for alignment
        
    

def runSW(input_file, score_file, open_gap, extension_gap, output_file):
    with open(input_file, 'r') as f: 
        seq_1, seq_2 = [x.strip() for x in f.readlines()] ### Read input sequences

    with open(output_file, 'w') as f: ### Write sequences to output
        f.write("-----------\n|Sequences|\n-----------\nsequence1\n")
        f.write(seq_1)
        f.write("\nsequence2\n")
        f.write(seq_2)
        f.write("\n--------------\n|Score Matrix|\n--------------\n")
    

        

    weight_matrix = pd.read_csv(score_file, delim_whitespace=True) ### Read in the weight matrix

    score_matrix = np.zeros([len(seq_2) + 1, len(seq_1) + 1]) ### Initialize the scoring matrix
    score_matrix[0, ] = 0
    score_matrix[:, 0] = 0

    direction_matrix = np.zeros([len(seq_2) + 1, len(seq_1) + 1]) ### 0 if match, 1 if gap in seq_1, 2 if gap in seq_2, 3 if 0 is max
    
    for j in range(len(seq_1)): ### j is column number, corresponding to sequence 1
        for i in range(len(seq_2)): ### i is row number, corresponding to sequence 2
            no_gap, gap_seq1, gap_seq2 = [score_matrix[i,j] + weight_matrix.loc[seq_1[j], seq_2[i]], # score added if no gap is introduced
                                            np.max(score_matrix[i+1,:j+1] - np.arange(-j * extension_gap - open_gap, -open_gap - 1, extension_gap)), # add a gap in seq_1 (j-indexed), adding penalties for extension and open gaps
                                            np.max(score_matrix[:i+1,j+1] - np.arange(-i * extension_gap - open_gap, -open_gap - 1, extension_gap))] # add a gap in seq_2 (i-indexed), adding penalties for extension and open gaps
            
            score_matrix[i+1,j+1] = np.max([no_gap, gap_seq1, gap_seq2, 0])
            direction_matrix[i+1, j+1] = np.argmax([no_gap, gap_seq1, gap_seq2, 0])

            # vals = [np.max(score_matrix[:i+1,j+1] - np.arange(-i * extension_gap - open_gap, -open_gap - 1, extension_gap)), # gap is introduced/extended in seq_2
            #                                 np.max(score_matrix[i+1,:j+1] - np.arange(-j * extension_gap - open_gap, -open_gap - 1, extension_gap)), # gap is introduced/extended in seq_1
            #                                 score_matrix[i,j] + weight_matrix.loc[seq_1[j], seq_2[i]], # match or mismatch, no gap introduced or extended
            #                                 0]
            
            # max_val = np.argwhere(vals == np.amax(vals)).flatten().tolist()
            # if len(max_val) > 1:
            #     print('SGHWI')
            
                                            
                                            
    score_matrix = score_matrix.astype(int)
    direction_matrix = direction_matrix.astype(int)

    alignment = traceback(score_matrix, direction_matrix) ### Use custom function traceback()

    aligned_seq_1 = "" ### Initialize string for aligned seq_1
    aligned_seq_2 = "" ### Initialize string for aligned seq_2
    aligned_markers = "" ### Intiialize string to show matches in alignment

    if alignment[0][1] > alignment[0][0]: ### Case where seq_2 begins before seq_1 in the local alignment
        aligned_seq_2 += " "*(alignment[0][1] - alignment[0][0])
    elif alignment[0][0] > alignment[0][1]: ### Case where seq_1 begins before seq_2 in the local alignment
        aligned_seq_1 += " "*(alignment[0][0] - alignment[0][1])

    aligned_markers += " "*max(alignment[0])

    aligned_seq_1 += seq_1[0:alignment[0][1] - 1]
    aligned_seq_2 += seq_2[0:alignment[0][0] - 1]

    aligned_seq_1 += '('
    aligned_seq_2 += '('

    prev_i, prev_j = 0, 0 ### Code to fill in the alignment
    for i, j in alignment:
        if prev_j == j:
            aligned_seq_1 += '-'
        else:
            aligned_seq_1 += seq_1[j-1]

        if prev_i == i:
            aligned_seq_2 += '-'
        else:
            aligned_seq_2 += seq_2[i-1]

        if seq_1[j-1] == seq_2[i-1]:
            aligned_markers += '|'
        else:
            aligned_markers += ' '
        
        prev_i, prev_j = i, j

    aligned_seq_1 += ')'
    aligned_seq_2 += ')'

    aligned_seq_1 += seq_1[alignment[-1][1]:len(seq_1)]
    aligned_seq_2 += seq_2[alignment[-1][0]:len(seq_2)]

    ### Write to output file
    pd.DataFrame(score_matrix, columns = [' '] + [x for x in seq_1], index = [' '] + [x for x in seq_2]).to_csv(output_file, sep = '\t', mode = 'a') ### Write score matrix
    with open(output_file, 'a') as f:
        f.write("----------------------\n|Best Local Alignment|\n----------------------\n")
        f.write("Alignment Score:" + str(score_matrix.max()))
        f.write("\nAlignment Results:\n")
        f.write('\n'.join([aligned_seq_1, aligned_markers, aligned_seq_2]))
        f.write('\n')

### Run your Smith-Waterman Algorithm
# runSW(args.input, args.score, args.opengap, args.extgap, args.result)