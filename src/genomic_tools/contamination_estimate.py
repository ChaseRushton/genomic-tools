#!/usr/bin/env python3

import pysam
import argparse
import numpy as np
from collections import defaultdict
import pandas as pd
from typing import Dict, List, Tuple
import json
from scipy import stats

class ContaminationEstimator:
    """Estimate sample contamination using allele fractions at common SNP sites."""
    
    def __init__(self, bam_path: str, reference_path: str):
        """Initialize with paths to BAM and reference files."""
        self.bam = pysam.AlignmentFile(bam_path, "rb")
        self.reference = pysam.FastaFile(reference_path)
        
    def calculate_allele_fractions(self, snp_sites: str, min_base_qual: int = 20,
                                 min_mapping_qual: int = 20) -> Dict:
        """
        Calculate allele fractions at common SNP sites.
        
        Args:
            snp_sites: Path to BED file containing common SNP positions
            min_base_qual: Minimum base quality to consider
            min_mapping_qual: Minimum mapping quality to consider
            
        Returns:
            Dictionary containing allele fractions and related metrics
        """
        allele_counts = defaultdict(lambda: defaultdict(int))
        total_sites = 0
        usable_sites = 0
        
        # Process each SNP site
        with open(snp_sites) as f:
            for line in f:
                chrom, pos, *_ = line.strip().split('\t')
                pos = int(pos)
                total_sites += 1
                
                # Get reference base
                ref_base = self.reference.fetch(chrom, pos-1, pos).upper()
                
                # Process pileup at this position
                for pileup in self.bam.pileup(chrom, pos-1, pos, truncate=True,
                                            min_base_quality=min_base_qual,
                                            min_mapping_quality=min_mapping_qual):
                    if pileup.pos == pos-1:  # Found our position
                        bases = defaultdict(int)
                        total_depth = 0
                        
                        # Count bases
                        for read in pileup.pileups:
                            if not read.is_del and not read.is_refskip:
                                base = read.alignment.query_sequence[read.query_position].upper()
                                qual = read.alignment.query_qualities[read.query_position]
                                if qual >= min_base_qual and read.alignment.mapping_quality >= min_mapping_qual:
                                    bases[base] += 1
                                    total_depth += 1
                        
                        # Calculate allele fractions if we have sufficient depth
                        if total_depth >= 20:  # Minimum depth threshold
                            usable_sites += 1
                            sorted_bases = sorted(bases.items(), key=lambda x: x[1], reverse=True)
                            if len(sorted_bases) >= 2:  # At least two alleles present
                                major_allele, major_count = sorted_bases[0]
                                minor_allele, minor_count = sorted_bases[1]
                                af = minor_count / total_depth
                                
                                # Store allele fraction
                                key = f"{chrom}:{pos}"
                                allele_counts[key] = {
                                    'major_allele': major_allele,
                                    'minor_allele': minor_allele,
                                    'major_count': major_count,
                                    'minor_count': minor_count,
                                    'depth': total_depth,
                                    'allele_fraction': af
                                }
        
        return {
            'allele_counts': dict(allele_counts),
            'total_sites': total_sites,
            'usable_sites': usable_sites
        }
    
    def estimate_contamination(self, allele_data: Dict) -> Dict:
        """
        Estimate contamination level based on allele fraction distribution.
        
        Args:
            allele_data: Dictionary containing allele fraction data
            
        Returns:
            Dictionary containing contamination estimates and confidence metrics
        """
        allele_fractions = []
        for site_data in allele_data['allele_counts'].values():
            af = site_data['allele_fraction']
            allele_fractions.append(af)
            
        if not allele_fractions:
            return {
                'contamination_estimate': 'NA',
                'confidence': 'NA',
                'warning': 'No usable sites found'
            }
            
        afs = np.array(allele_fractions)
        
        # Expected peaks for homozygous ref (0.0), het (0.5), and homozygous alt (1.0)
        # Calculate deviation from expected
        hom_ref_dev = np.mean(afs[afs < 0.15]) if len(afs[afs < 0.15]) > 0 else 0
        het_dev = np.abs(np.mean(afs[(afs >= 0.4) & (afs <= 0.6)]) - 0.5) if len(afs[(afs >= 0.4) & (afs <= 0.6)]) > 0 else 0
        hom_alt_dev = np.abs(1 - np.mean(afs[afs > 0.85])) if len(afs[afs > 0.85]) > 0 else 0
        
        # Calculate contamination estimate
        # For heterozygous sites, deviation from 0.5 indicates contamination
        het_sites = afs[(afs >= 0.25) & (afs <= 0.75)]
        if len(het_sites) < 10:
            return {
                'contamination_estimate': 'NA',
                'confidence': 'Low',
                'warning': 'Too few heterozygous sites'
            }
            
        # Calculate contamination estimate based on deviation from expected ratios
        contamination_estimate = 2 * np.mean([hom_ref_dev, het_dev, hom_alt_dev])
        
        # Calculate confidence based on number of sites and consistency
        n_sites = len(allele_fractions)
        consistency = 1 - np.std([hom_ref_dev, het_dev, hom_alt_dev])
        
        confidence = 'High' if n_sites >= 100 and consistency >= 0.9 else \
                    'Medium' if n_sites >= 50 and consistency >= 0.7 else 'Low'
                    
        # Additional metrics
        metrics = {
            'contamination_estimate': round(contamination_estimate * 100, 2),  # Convert to percentage
            'confidence': confidence,
            'number_of_sites': n_sites,
            'consistency_score': round(consistency, 3),
            'homozygous_deviation': round(np.mean([hom_ref_dev, hom_alt_dev]), 3),
            'heterozygous_deviation': round(het_dev, 3),
            'warning': None
        }
        
        # Add warnings if necessary
        if n_sites < 50:
            metrics['warning'] = 'Low number of informative sites'
        elif consistency < 0.7:
            metrics['warning'] = 'Low consistency across sites'
            
        return metrics
    
    def run_analysis(self, snp_sites: str, min_base_qual: int = 20,
                    min_mapping_qual: int = 20) -> Dict:
        """Run complete contamination analysis."""
        # Calculate allele fractions
        allele_data = self.calculate_allele_fractions(snp_sites, min_base_qual, min_mapping_qual)
        
        # Estimate contamination
        contamination_metrics = self.estimate_contamination(allele_data)
        
        # Combine results
        results = {
            'sample_metrics': {
                'total_snp_sites': allele_data['total_sites'],
                'usable_sites': allele_data['usable_sites']
            },
            'contamination_metrics': contamination_metrics
        }
        
        return results
    
    def close(self):
        """Close open file handles."""
        self.bam.close()
        self.reference.close()

def main():
    parser = argparse.ArgumentParser(
        description="Estimate sample contamination using allele fractions at common SNP sites"
    )
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--reference", required=True, help="Reference genome FASTA")
    parser.add_argument("--snps", required=True, help="BED file with common SNP positions")
    parser.add_argument("--min-base-qual", type=int, default=20,
                       help="Minimum base quality (default: 20)")
    parser.add_argument("--min-mapping-qual", type=int, default=20,
                       help="Minimum mapping quality (default: 20)")
    parser.add_argument("--output", help="Output JSON file (optional, default: stdout)")
    
    args = parser.parse_args()
    
    # Run contamination estimation
    estimator = ContaminationEstimator(args.bam, args.reference)
    results = estimator.run_analysis(args.snps, args.min_base_qual, args.min_mapping_qual)
    estimator.close()
    
    # Format results for output
    output_data = {
        'sample_name': args.bam,
        'analysis_parameters': {
            'min_base_quality': args.min_base_qual,
            'min_mapping_quality': args.min_mapping_qual
        },
        'results': results
    }
    
    # Output results
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(output_data, f, indent=2)
    else:
        print(json.dumps(output_data, indent=2))

if __name__ == "__main__":
    main()
