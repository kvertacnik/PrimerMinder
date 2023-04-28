# PrimerMinder

PrimerMinder (as in mind your surrounding polymorphisms) facilitates primer design by visualizing sequence variation that may interfere with primer annealing. 


Inputs are a VCF file, the genome assembly the VCF is derived from, a list of target loci, Primer3 settings, and a user-specified flanking region. PrimerMinder searches the VCF for variants in the region, uses Primer3 to generate primers, and produces a visual alignment of the flanking region, target marker, adjacent snps, and predicted primers. Provided a list of markers, it will sort the markers by number of adjacent snps (lowest to highest).


![Example output](/repository/image.png?raw=true "Example output")


## Installation

Dependencies: pandas, Primer3-py, pyfaidx, yaml

Import environment:
'conda env create -f PrimerMinder_env.yml'


## Usage
'./PrimerMinder.py --reference_genome= --vcf= --marker_file= --range= --primer3_settings= > output.txt'


- Reference_genome must be the same assembly used to generate the vcf file
- Marker_file must be formatted Chromosome:base_pair; one target marker per line
- Chromosome names must be consistent across the assembly, VCF, and marker_file
- Modify Primer3 settings using the Primer3_settings_default.yaml file
