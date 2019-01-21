# Encode alignment file using a weighting scheme derived from Krogh et all (1995).
# Duplicates are filtered and external gaps are skipped.
# No gap is considered to calculate the weight of the sequence.
# Clipping of the AA in the query sequence is implemented as in https://doi.org/10.1101/289033.

# Usage: python3 process-alignment.py alignment_file file_extension
# Expected input: 1 line header followed by primary sequences (and "." gaps) only.

import os
import sys
import math
from numpy import ravel as np

pattern = sys.argv[2] # alignment extension
pattern1 = ".ann" # encoded alignment extension
alignment = sys.argv[1] # alignment path

# create dictionary to represent the primary structure
aa = {"A": 0, "C": 1, "D": 2, "E": 3, "F": 4, "G": 5, "H": 6, "I": 7, "K": 8, "L": 9, "M": 10, "N": 11, "P": 12,
        "Q": 13, "R": 14, "S": 15, "T": 16, "V": 17, "W": 18, "Y": 19,
        "B": 20, "J": 20, "O": 20, "U": 20, "X": 20, "Z": 20,
        ".": 21}


# open and read file to parse
f = open(alignment, "r")
lines_raw = f.readlines()

# filter and delete duplicates
lines = []
lines_set = set()
for l in lines_raw:
    if l not in lines_set:
        lines.append(l)
    lines_set.add(l)
# update number of sequences
sequences = lines[0] = len(lines) - 1

# delete pattern from the filename
pid = alignment.replace("."+pattern, "")
# create file to fill up
f = open(alignment + pattern1, "w")

# write header and filename on the first 3 rows
f.write("1\n22 3\n"+ pid + "\n")

# check length and write it on the second row
length = len(lines[1].strip())
f.write(str(length) + "\n")

# create a list of lists to archive the frequencies,
frequencies = [[0] * 22 for i in range(length)]
# a list of lists to archive the extremities,
extremities = [[0] * 2 for i in range(sequences)]
# a list of list to archive the profiles,
result = [[0] * 22 for i in range(length)]
# and a list for the weighting scheme
weights = [0] * int(lines[0])


if sequences > 1:
    # check each column/AA to find the right and left extremity per sequence
    for i in range(sequences):

        for j in range(length):
            if lines[i+1][j] != ".":
                extremities[i][0] = j
                break

        for j in range(length):
            if lines[i+1][length-j-1] != ".":
                extremities[i][1] = length-j-1
                break


    # check each column/AA to calculate its percentages
    for j in range(length):
        gaps = 0

        # how many gaps and invalid sequences in the column?
        for i in range(sequences):
            if lines[i+1][j] == ".":
                gaps += 1

        # calculate the relative frequencies per column
        add = 1 / (sequences - gaps)
        for i in range(sequences):
            if lines[i+1][j] != ".":  # No gap is considered to calculate the weight of the sequence
                frequencies[j][aa[lines[i+1][j]]] += add

        # calculate weights
        for i in range(sequences):
            if lines[i+1][j] != ".": # No gap is considered to calculate the weight of the sequence
                weights[i] = weights[i] - math.log(frequencies[j][aa[lines[i+1][j]]])

    if weights[0] > 0:
        # use the weights to create the profile
        for j in range(length):
            weightsSum = 0

            # calculate weights sum and entropy
            for i in range(sequences):
                if not j < extremities[i][0] and not j > extremities[i][1]:
                    result[j][aa[lines[i+1][j]]] += weights[i]
            for i in range(21): # Gaps are normalised separately
                weightsSum += result[j][i]

            # normalize result
            for i in range(21):
                if result[j][i] != 0:
                    result[j][i] = result[j][i] / weightsSum
            if result[j][21] != 0:
                result[j][21] = result[j][21] / (weightsSum + result[j][21])
    elif sequences > 1:
        # 1 aligned sequence diverge only with gaps so no weights/entropy but gaps found
        result = frequencies

# clip to 1 the frequence of the AA presents in the profiled protein
for j in range(length):
    result[j][aa[lines[1][j]]] = 1.0

# convert result (a list of lists) in a string and write out the result
f.write(" ".join(map(str, np(result))) + "\n")
f.close()
