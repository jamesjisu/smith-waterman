# smith-waterman
An implementation of the Smith-Waterman local alignment algorithm in Python

## Introduction
The Smith-Waterman algorithm is a dynamic programming approach to the sequence/string alignment problem: succintly, given two strings, how can we align matching regions of the two strings to have the greatest similarity between the two? The [Wikipedia article](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm) about the Smith-Waterman algorithm provides an excellent explanation of the algorithm

## Installation
### Requirements
The only packages required to run this script are NumPy and Pandas, which come installed with [Anaconda](https://www.anaconda.com/).
### Installation
To install the package, run `pip install git+https://github.com/jamesjisu/smith-waterman.git`. This will install the package locally

## Usage
The package has two main functions with different use-cases: `SW` and `runSW`.
### `SW`
This function simply performs the Smith-Waterman algorithm on two strings given a weight matrix and returns the score matrix, the direction matrix, and the aligned sequences
```
SW(seq_1, seq_2, weight_matrix, open_gap = -2, extension_gap = -1)
```
returns
```
(score_matrix, direction_matrix, aligned_seq_1, aligned_seq_2)
```
### `runSW`
This function performs the function associated with assignment 1 for S&DS 352: Biomedical Data Mining. Given an input file, score file, gap penalties, and an output file, it will provide the input sequences, score matrix, and the aligned sequences. The input file must follow the following format
```
seq_1
seq_2
```
and will output the following
```
-----------
|Sequences|
-----------
sequence1
seq_1
sequence2
seq_2
--------------
|Score Matrix|
--------------
SCORE_MATRIX
----------------------
|Best Local Alignment|
----------------------
Alignment Score: BEST ALIGNMENT SCORE
Alignment Results: ALIGNMENT
```
The function within Python is called as such
```
runSW(input_file_path, score_file_path, output_file_path, open_gap, extension_gap)
```
and will output a file to `output_file_path` following the above format.

## Examples
Attached in the repository are two samples: `sample-input1.txt` and `sample-input2.txt`, as well as their corresponding outputs. You can run `runSW` to compare your installation with the repository.
