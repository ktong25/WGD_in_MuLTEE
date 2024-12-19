#!/usr/bin/env python3
"""
Author: g. ozan bozdag

This script calculates allele frequencies from VCF-derived TAB files.
It processes TAB files containing variant information and computes the frequency
of alternative alleles using the AD (Allele Depth) and DP (Total Depth) fields.

Input: TAB files in the current directory
Output: New TAB files with '_freq.tab' suffix containing original data plus allele frequencies

The script:
1. Reads TAB files containing variant information
2. Extracts AD (Allele Depth, for the reference allele) and DP (Total Depth: total number of reads covering the position) values
3. Calculates alternative allele frequency as a percentage
4. Writes results to new files with '_freq.tab' suffix

IMPORTANT NOTE: Before running this script, please verify that:
1. Input files must contain FORMAT and sample columns with AD and DP fields
2. Check that columns 8 and 9 in your input file contain 
   the correct AD and DP information. Column positions are zero-based, 
   so these are the 9th and 10th columns in the file.
"""

import glob

def compute_allele_frequency(input_file, output_file):
    """
    Calculate allele frequencies from variant data in TAB format.
    
    Args:
        input_file (str): Path to input TAB file containing variant data
        output_file (str): Path to output file for frequency results
        
    The function:
    - Extracts AD (Allele Depth: comma-separated read counts for ref,alt alleles) 
      and DP (Total Depth: total number of reads at position) values
    - Calculates alternative allele frequency as percentage
    - Handles potential errors in data format or missing values
    - Writes original data plus calculated frequency to output file
    """
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            fields = line.strip().split('\t')
            # Extracting AD and DP values from FORMAT and sample columns
            info_fields = fields[8].split(":")  # FORMAT column (9th column)
            info_values = fields[9].split(":")  # Sample column (10th column)
            try:
                # Find indices for AD and DP in FORMAT field
                ad_index = info_fields.index("AD")
                dp_index = info_fields.index("DP")
                
                # Parse reference and alternative allele counts
                ref_count, alt_count = map(int, info_values[ad_index].split(","))
                total_depth = int(info_values[dp_index])
                
                # Calculate alternative allele frequency as percentage
                freq = (100 * alt_count) / total_depth
                
                # Write original data plus frequency to output
                f_out.write(line.strip() + f"\t{freq:.2f}%\n")
            
            except (ValueError, IndexError):
                # Handle parsing errors or missing fields
                print(f"Error processing line: {line.strip()} in {input_file}")
                continue

def main():
    """
    Main function to process all TAB files in current directory.
    """
    # Get list of all TAB files in current directory
    tab_files = glob.glob('*.tab')
    
    # Process each TAB file
    for tab_file in tab_files:
        output_file = tab_file.rsplit('.', 1)[0] + '_freq.tab'
        compute_allele_frequency(tab_file, output_file)
    
    print(f"Processed {len(tab_files)} TAB files and calculated allele frequencies.")

if __name__ == "__main__":
    main()