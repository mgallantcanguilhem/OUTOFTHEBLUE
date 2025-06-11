#!/usr/bin/env python2

import numpy as np

# Read the first line to count the number of columns
with open("meff_matrix.raw", "r") as f:
    first_line = f.readline().strip().split()
    num_cols = len(first_line)
    
# Load PLINK raw file (skip first 6 columns)
data = np.loadtxt("meff_matrix.raw", dtype=str, delimiter=" ", usecols=range(6, num_cols))

# Transpose the matrix so SNPs are rows
data_transposed = data.T

# Remove the first row (SNP names)
data_transposed_clean = data_transposed[:, 1:]

# Find rows where all values are the same
rows_to_keep = ~np.all(data_transposed_clean == data_transposed_clean[:, 0][:, None], axis=1)

# Keep only polymorphic SNPs
filtered_data = data_transposed_clean[rows_to_keep, :]

# Save in space-separated format
np.savetxt("genotype_matrix.txt", filtered_data, fmt="%s", delimiter=" ")
