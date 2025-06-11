#!/usr/bin/env python3

import sys
import argparse
import allel
import numpy as np

# Function to get command-line options
def get_option_value(option):
    if option in sys.argv:
        option_pos = sys.argv.index(option)
        if option_pos + 1 < len(sys.argv):
            return sys.argv[option_pos + 1]
        else:
            print(f"\nWarning, option {option} not specified.\n", file=sys.stderr)
            sys.exit(1)
    else:
        print(f"\nplease specify {option} option\n", file=sys.stderr)
        sys.exit(1)

# Parse input and output file arguments
file_name = get_option_value("-i")
output_prefix = get_option_value("-o")

vcf_data = allel.read_vcf(file_name)

# Extract genotype data
gt = allel.GenotypeArray(vcf_data['calldata/GT'])

# Compute heterozygosity per site
het = gt.count_het(axis=1) / gt.count_called(axis=1)

# Get variant positions
positions = vcf_data['variants/POS']

# Save results in a two-column format: Position | Heterozygosity
output_file = f"{output_prefix}_heterozygosity.txt"
np.savetxt(output_file, np.column_stack((positions, het)), fmt="%d %.6f", header="POS HET", comments='')

print(f"Heterozygosity values saved to {output_file}")
