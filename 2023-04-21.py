#!/usr/bin/env python3


""" Usage: primerMinder.py --reference_genome= --vcf= --marker_file= --primer_range= > output_file.txt
Dependencies: pandas, Primer3-py (conda: Primer3-py), pyfaidx, yaml (conda: pyyaml).
For a given target marker, looks for adjacent snps that may interfere with primer design.
 Provided a list of markers, it will sort the markers by number of adjacent snps (lowest to highest).
Reference_genome must be the same assembly used to generate the vcf file.
Marker_file must be formatted Chromosome:base_pair; chromosome names must be the same across marker_file, reference_genome and vcf.
 One target marker per line.
Need help? kim.vertacnik@uky.edu"""


import sys
import argparse
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
"filtered_SNPs = vcf_df.loc[(vcf_df['CHROM'] == target_snp_dict[key]['target_chromosome']) & (vcf_df['POS']).between(target_snp_dict[key]['target_lower_range'], target_snp_dict[key]['target_upper_range'], inclusive=True)]"
import pandas as pd
import primer3
import yaml
from pyfaidx import Fasta
from collections import defaultdict


parser = argparse.ArgumentParser()
parser.add_argument('--reference_genome',
                    type=str,
                    required=True,
                    metavar='<Same assembly used to create the vcf file>'
                    )
parser.add_argument('--vcf',
                    type=str,
                    required=True
                    )
parser.add_argument('--marker_file',
                    type=str,
                    required=True,
                    metavar='chromosome:position format'
                    )
parser.add_argument('--range',
                    type=int,
                    required=True,
                    metavar='Value is applied up and downstream from the target marker'
                    )
parser.add_argument('--primer3_settings',
                    type=str,
                    required=True,
                    metavar='.yaml file (see primer3_default.yaml)'
                    )
args = parser.parse_args()


reference_genome = args.reference_genome
assembly = Fasta(reference_genome)
vcf_file = args.vcf
marker_file = args.marker_file
target_range = args.range
primer3_settings_file = args.primer3_settings
target_snp_dict = defaultdict(dict)
adjacent_snp_dict = defaultdict(dict)
primer3_output_dict = defaultdict(dict)


def insert_letter(snp_string,index,alt_allele):
    insert_size = -abs(len(alt_allele))
    snp_string = snp_string[:index] + alt_allele + snp_string[index:]
    snp_string = snp_string[:insert_size]
    return snp_string


def degenerate_nucleotide_conversion(multiallelic):
    if '*' in multiallelic:
        return '*'
    elif 'A' and 'T' in multiallelic: 
        return 'W'
    elif 'A' and 'G' in multiallelic:
        return 'R'
    elif 'A' and 'C' in multiallelic:
        return 'M'
    elif 'T' and 'G' in multiallelic:
        return 'K'
    elif 'T' and 'C' in multiallelic:
        return 'Y'
    elif 'G' and 'C' in multiallelic:
        return 'S'


def rev_comp(seq):
    complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    return "".join(complement[base] for base in reversed(seq))


# Read in vcf file as a pandas database. Add column names. Convert 'POS' column from string to integer.
vcf_df = pd.read_csv(vcf_file, sep="\t", comment='#', header=None, usecols=[0,1,3,4])
vcf_df.columns = ['CHROM', 'POS', 'REF', 'ALT']
vcf_df[["POS"]] = vcf_df[["POS"]].apply(pd.to_numeric)


# Read in target snp file
with open(marker_file) as f:
    for line in f:
        if not line.strip(): # 'not' reverses conditional; empty line is considered False
            continue
        else:
            target_snp_key = line.strip()
            target_snp_coords = line.strip().split(':')
            target_chromosome = target_snp_coords[0]
            target_snp_dict[target_snp_key]['target_chromosome'] = target_chromosome
            target_snp_dict[target_snp_key]['target_position'] = int(target_snp_coords[1])
            target_snp_dict[target_snp_key]['target_index'] = target_range + 1
            target_snp_dict[target_snp_key]['target_lower_range'] = int(target_snp_coords[1]) - target_range
            target_snp_dict[target_snp_key]['target_upper_range'] = int(target_snp_coords[1]) + target_range
            target_snp_dict[target_snp_key]['target_ref_seq']= assembly[target_chromosome][int(target_snp_coords[1])-target_range-1:int(target_snp_coords[1])+target_range].seq
            target_snp_dict[target_snp_key]['target_primer3'] = {}


# For each target snp, find adjacent snps within the user-specified range & append to the target snp's dictionary as a nested dictionary
#  Pandas df is filtered based on chromosome and +- primer_range from the target's position); the filtered rows are converted to a dictionary and appended
for target in target_snp_dict:
    filtered_SNPs = vcf_df.loc[(vcf_df['CHROM'] == target_snp_dict[target]['target_chromosome']) & (vcf_df['POS']).between(target_snp_dict[target]['target_lower_range'], target_snp_dict[target]['target_upper_range'], inclusive='both')]
    adjacent_snp_dict = filtered_SNPs.to_dict(orient='index')
    for adjacent in adjacent_snp_dict:
        # Get the adjacent snp's positon
        if target_snp_dict[target]['target_position'] <= adjacent_snp_dict[adjacent]['POS']:
            index = abs(target_snp_dict[target]['target_position'] - adjacent_snp_dict[adjacent]['POS']) + target_range
        else:
            index = target_range - abs(target_snp_dict[target]['target_position'] - adjacent_snp_dict[adjacent]['POS'])
        adjacent_snp_dict[adjacent]['index'] = index
        # Classify polymorphism: target_snp=T, other_snp=s, multiallelic=m, insertion=i, deletion=d
        if len(adjacent_snp_dict[adjacent]['ALT']) == 1:
            adjacent_snp_dict[adjacent]['type'] = 's'
            if target_snp_dict[target]['target_position'] == adjacent_snp_dict[adjacent]['POS']:
                adjacent_snp_dict[adjacent]['type'] = 'T'
        elif ',' in adjacent_snp_dict[adjacent]['ALT']:
            adjacent_snp_dict[adjacent]['type'] = 'm'
            adjacent_snp_dict[adjacent]['ALT'] = degenerate_nucleotide_conversion(adjacent_snp_dict[adjacent]['ALT'])
        elif len(adjacent_snp_dict[adjacent]['ALT']) > 1:
            adjacent_snp_dict[adjacent]['type'] = 'i'
        elif '*' in adjacent_snp_dict[adjacent]['ALT']:
            adjacent_snp_dict[adjacent]['type'] = 'd'
    target_snp_dict[target]['adjacent_snps'] = adjacent_snp_dict


# For each target snp, generate the alignment string
for target in target_snp_dict:
    target_ref_seq_length = len(target_snp_dict[target]['target_ref_seq'])
    snp_string = '-' * (target_ref_seq_length)
    for adjacent in target_snp_dict[target]['adjacent_snps']:
        index = target_snp_dict[target]['adjacent_snps'][adjacent]['index']
        alt_allele = target_snp_dict[target]['adjacent_snps'][adjacent]['ALT']
        snp_string = insert_letter(snp_string,index,alt_allele)
    target_snp_dict[target]['snp_string'] = snp_string


# For each target snp, generate the type string
for target in target_snp_dict:
    target_ref_seq_length = len(target_snp_dict[target]['target_ref_seq'])
    snp_string = ' ' * (target_ref_seq_length)
    for adjacent in target_snp_dict[target]['adjacent_snps']:
        index = target_snp_dict[target]['adjacent_snps'][adjacent]['index']
        polymorph_type = target_snp_dict[target]['adjacent_snps'][adjacent]['type']
        snp_string = insert_letter(snp_string,index,polymorph_type)
    target_snp_dict[target]['type_string'] = snp_string


# Read in Primer3 settings
with open(primer3_settings_file, "r") as stream:
    try:
        settings = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)


# Run Primer3. Export all output to file (file name = markerName_primer3.out). For predicted primers, capture sequence and index position info
for target in target_snp_dict:
    primer3_output_dict = primer3.bindings.design_primers(
       seq_args = {'SEQUENCE_ID': target_snp_dict[target],
            'SEQUENCE_TEMPLATE': target_snp_dict[target]['target_ref_seq'],
            'SEQUENCE_TARGET': [target_snp_dict[target]['target_index'],1]
        },
        global_args = settings['global_args']
    )

    with open(target.replace(":", "_") + "_primer3.out", 'w') as f:
        yaml.dump(primer3_output_dict, f)
    
    how_many_primers = primer3_output_dict['PRIMER_PAIR_NUM_RETURNED']
    if how_many_primers == 0:
        target_snp_dict[target]['target_primer3'] = 'No predicted primers'
        continue
    else:
        primer_list = list(range(how_many_primers))
        for i in primer_list:
            target_snp_dict[target]['target_primer3'][f'primer_{i}'] = {}
            selected_keys = [f'PRIMER_LEFT_{i}_SEQUENCE', f'PRIMER_LEFT_{i}',f'PRIMER_RIGHT_{i}_SEQUENCE', f'PRIMER_RIGHT_{i}']
            for key in primer3_output_dict.keys():
                if key in selected_keys:
                    target_snp_dict[target]['target_primer3'][f'primer_{i}'][key] = primer3_output_dict[key]


# For each primer pair, generate the alignment string
for target in target_snp_dict:
        for primer in target_snp_dict[target]['target_primer3']:
            if target_snp_dict[target]['target_primer3'] == 'No predicted primers':
                continue
            else:
                primer_num = primer.split('_')[1]
                target_ref_seq_length = len(target_snp_dict[target]['target_ref_seq'])
                snp_string = ' ' * (target_ref_seq_length)
            
                left_index = target_snp_dict[target]['target_primer3'][f'primer_{primer_num}'][f'PRIMER_LEFT_{primer_num}'][0]
                left_seq =  target_snp_dict[target]['target_primer3'][f'primer_{primer_num}'][f'PRIMER_LEFT_{primer_num}_SEQUENCE']
                snp_string = insert_letter(snp_string,left_index,left_seq)
            
                right_seq = rev_comp(target_snp_dict[target]['target_primer3'][f'primer_{primer_num}'][f'PRIMER_RIGHT_{primer_num}_SEQUENCE'])
                right_index = (target_snp_dict[target]['target_primer3'][f'primer_{primer_num}'][f'PRIMER_RIGHT_{primer_num}'][0]) - len(right_seq) + 1
                snp_string = insert_letter(snp_string,right_index,right_seq)
            
            
                target_snp_dict[target]['target_primer3'][primer]['primer_string'] = snp_string


# Get a list of target_snps sorted by number of adjacent snps (low to high)
target_snps_sorted_by_num_adj_snps = sorted(target_snp_dict, key=lambda k: len(target_snp_dict[k]['adjacent_snps']))


# Print all alignment strings to stdout
print("(0-index numbering)\n")
for target in target_snps_sorted_by_num_adj_snps:
    snp_pos_list = []
    for adjacent in target_snp_dict[target]['adjacent_snps']:
        snp_pos_list.append(target_snp_dict[target]['adjacent_snps'][adjacent]['POS'])
    print(f"Marker: {target}")
    print(f"Sequence: {target_snp_dict[target]['target_lower_range']}:{target_snp_dict[target]['target_upper_range']} bp")
    print(f"Variable site(s) at: {snp_pos_list}")
    print("Polymorphism labels: target_snp=T, other_snp=s, multiallelic=m, insertion=i, deletion=d")
    print(target_snp_dict[target]['type_string'])
    print(target_snp_dict[target]['snp_string'])
    print(target_snp_dict[target]['target_ref_seq'])
    if target_snp_dict[target]['target_primer3'] == 'No predicted primers':
        print("No predicted primers")
    else:
        for primer_pair in primer_list:
            print(target_snp_dict[target]['target_primer3'][f'primer_{primer_pair}']['primer_string'])
    print('\n')


