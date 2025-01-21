#!/usr/bin/env python3

import pysam
import argparse
import sys
from collections import defaultdict

def check_variant_support(bam_file, vcf_file, min_depth=10, min_mapping_quality=20, min_base_quality=20, min_allele_freq=0.2):
    """
    Validate variants in VCF file against evidence in BAM file.
    
    Args:
        bam_file (str): Path to BAM file
        vcf_file (str): Path to VCF file
        min_depth (int): Minimum read depth required
        min_mapping_quality (int): Minimum mapping quality for reads
        min_base_quality (int): Minimum base quality for variant position
        min_allele_freq (float): Minimum allele frequency required
    """
    try:
        # Open BAM and VCF files
        bam = pysam.AlignmentFile(bam_file, "rb")
        vcf = pysam.VariantFile(vcf_file)
        
        print("Starting validation...")
        
        # Counter for validation results
        results = defaultdict(int)
        
        # Iterate through variants in VCF
        for variant in vcf:
            chrom = variant.chrom
            pos = variant.pos - 1  # Convert to 0-based coordinate
            ref = variant.ref
            alt = variant.alts[0] if variant.alts else None
            
            if not alt:
                continue
                
            # Get pileup at variant position
            pileup = bam.pileup(chrom, pos, pos + 1, truncate=True, min_mapping_quality=min_mapping_quality, min_base_quality=min_base_quality)
            
            for pileupcolumn in pileup:
                if pileupcolumn.pos == pos:
                    # Count bases at this position
                    bases = defaultdict(int)
                    total_depth = 0
                    alt_depth = 0
                    
                    for read in pileupcolumn.pileups:
                        if not read.is_del and not read.is_refskip:
                            base = read.alignment.query_sequence[read.query_position]
                            bases[base] += 1
                            total_depth += 1
                            if base == alt[0]:  # For SNVs
                                alt_depth += 1
                    
                    # Calculate allele frequency
                    allele_freq = alt_depth / total_depth if total_depth > 0 else 0
                    
                    # Check validation criteria
                    support_found = (
                        total_depth >= min_depth and  # Check depth
                        allele_freq >= min_allele_freq and  # Check allele frequency
                        (
                            (len(alt) == len(ref) and bases[alt] > 0) or  # SNV
                            (len(alt) != len(ref) and bases[alt[0]] > 0)  # INDEL
                        )
                    )
                    
                    if support_found:
                        results['supported'] += 1
                    else:
                        results['unsupported'] += 1
                        print(f"Warning: Variant at {chrom}:{pos+1} {ref}>{alt}")
                        print(f"  Depth: {total_depth} (min: {min_depth})")
                        print(f"  Allele Frequency: {allele_freq:.3f} (min: {min_allele_freq})")
                        print(f"  Base counts: {dict(bases)}")
            
        # Print summary
        print("\nValidation Summary:")
        print(f"Total variants processed: {sum(results.values())}")
        print(f"Variants with support: {results['supported']}")
        print(f"Variants without support: {results['unsupported']}")
        
    except Exception as e:
        print(f"Error during validation: {str(e)}", file=sys.stderr)
        sys.exit(1)
    finally:
        if 'bam' in locals():
            bam.close()
        if 'vcf' in locals():
            vcf.close()

def main():
    parser = argparse.ArgumentParser(description='Validate VCF variants against BAM file evidence')
    parser.add_argument('--bam', required=True, help='Input BAM file')
    parser.add_argument('--vcf', required=True, help='Input VCF file')
    parser.add_argument('--min-depth', type=int, default=10, help='Minimum read depth required (default: 10)')
    parser.add_argument('--min-mapping-quality', type=int, default=20, help='Minimum mapping quality (default: 20)')
    parser.add_argument('--min-base-quality', type=int, default=20, help='Minimum base quality (default: 20)')
    parser.add_argument('--min-allele-freq', type=float, default=0.2, help='Minimum allele frequency required (default: 0.2)')
    
    args = parser.parse_args()
    
    check_variant_support(
        args.bam, 
        args.vcf,
        args.min_depth,
        args.min_mapping_quality,
        args.min_base_quality,
        args.min_allele_freq
    )

if __name__ == '__main__':
    main()
