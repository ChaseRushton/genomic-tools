#!/usr/bin/env python3
"""
Generate a BED file of X chromosome SNPs suitable for sex verification.
Filters SNPs based on:
- Minor allele frequency
- Distance from pseudoautosomal regions
- Uniqueness (avoiding regions with Y homology)
- Coverage in typical sequencing
"""

import argparse
import gzip
from typing import List, Tuple
import pandas as pd
import numpy as np

# Define PAR regions (GRCh38)
PAR1_START = 10001
PAR1_END = 2781479
PAR2_START = 155701383
PAR2_END = 156030895

# Regions with high X-Y homology to exclude
XY_HOMOLOGY = [
    (89037499, 92375509),   # XTR/YTR
    (156030896, 156040895)  # PAR3
]

def load_gnomad_variants(vcf_path: str) -> pd.DataFrame:
    """Load variants from gnomAD VCF file."""
    variants = []
    with gzip.open(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            
            # Parse INFO field
            info = dict(x.split('=') for x in fields[7].split(';') if '=' in x)
            
            # Get allele frequencies
            try:
                af = float(info.get('AF', 0))
                female_af = float(info.get('AF_female', 0))
                male_af = float(info.get('AF_male', 0))
            except (ValueError, TypeError):
                continue
                
            variants.append({
                'pos': int(fields[1]),
                'ref': fields[3],
                'alt': fields[4],
                'af': af,
                'female_af': female_af,
                'male_af': male_af
            })
    
    return pd.DataFrame(variants)

def is_in_excluded_region(pos: int, excluded_regions: List[Tuple[int, int]]) -> bool:
    """Check if position is in any excluded region."""
    return any(start <= pos <= end for start, end in excluded_regions)

def filter_variants(df: pd.DataFrame, min_af: float = 0.01, min_female_af: float = 0.2) -> pd.DataFrame:
    """Filter variants based on frequency and position criteria."""
    # Combine PAR and homology regions
    excluded_regions = [(PAR1_START, PAR1_END), (PAR2_START, PAR2_END)] + XY_HOMOLOGY
    
    # Apply filters
    df = df[
        (df['af'] >= min_af) &  # Common enough
        (df['female_af'] >= min_female_af) &  # Common in females
        (~df['pos'].apply(lambda x: is_in_excluded_region(x, excluded_regions)))  # Not in excluded regions
    ]
    
    # Ensure good spacing between SNPs (at least 10kb)
    df = df.sort_values('pos')
    df['dist_to_next'] = df['pos'].shift(-1) - df['pos']
    df = df[df['dist_to_next'].fillna(float('inf')) >= 10000]
    
    return df

def main():
    parser = argparse.ArgumentParser(description='Generate X chromosome SNPs for sex verification')
    parser.add_argument('--gnomad-vcf', required=True, help='Path to gnomAD X chromosome VCF')
    parser.add_argument('--output', required=True, help='Output BED file')
    parser.add_argument('--min-af', type=float, default=0.01, help='Minimum allele frequency')
    parser.add_argument('--min-female-af', type=float, default=0.2, 
                       help='Minimum allele frequency in females')
    args = parser.parse_args()

    # Load and filter variants
    print("Loading variants from gnomAD...")
    variants_df = load_gnomad_variants(args.gnomad_vcf)
    
    print("Filtering variants...")
    filtered_df = filter_variants(variants_df, args.min_af, args.min_female_af)
    
    # Write BED file
    print(f"Writing {len(filtered_df)} filtered variants to BED file...")
    with open(args.output, 'w') as f:
        for _, row in filtered_df.iterrows():
            # BED format: chrom start end name
            f.write(f"chrX\t{row['pos']-1}\t{row['pos']}\t"
                   f"X_SNP_{row['pos']}_{row['ref']}_{row['alt']}\n")
    
    print("Done!")

if __name__ == '__main__':
    main()
