#!/usr/bin/env python3

import os
import glob
import argparse
import pandas as pd
import numpy as np
from typing import Dict, List, Tuple
import json
from contamination_estimate import ContaminationEstimator
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations
import pysam

class BatchContaminationAnalyzer:
    """Analyze potential cross-sample contamination in a batch of BAM files."""
    
    def __init__(self, bam_dir: str, reference_path: str):
        """Initialize with directory containing BAM files and reference genome path."""
        self.bam_dir = bam_dir
        self.reference_path = reference_path
        self.bam_files = []
        for ext in ['*.bam', '*.cram']:
            self.bam_files.extend(glob.glob(os.path.join(bam_dir, ext)))
    
    def analyze_single_sample(self, bam_path: str, snp_sites: str,
                            min_base_qual: int = 20,
                            min_mapping_qual: int = 20) -> Dict:
        """Analyze contamination for a single sample."""
        try:
            estimator = ContaminationEstimator(bam_path, self.reference_path)
            results = estimator.run_analysis(snp_sites, min_base_qual, min_mapping_qual)
            estimator.close()
            return {
                'sample': os.path.basename(bam_path),
                'results': results
            }
        except Exception as e:
            return {
                'sample': os.path.basename(bam_path),
                'error': str(e)
            }
    
    def get_sample_alleles(self, bam_path: str, snp_sites: str,
                          min_base_qual: int = 20,
                          min_mapping_qual: int = 20) -> Dict[str, float]:
        """Get allele fractions at SNP sites for a sample."""
        allele_fractions = {}
        bam = pysam.AlignmentFile(bam_path, "rb")
        
        with open(snp_sites) as f:
            for line in f:
                chrom, pos, *_ = line.strip().split('\t')
                pos = int(pos)
                
                # Process pileup at this position
                for pileup in bam.pileup(chrom, pos-1, pos, truncate=True,
                                       min_base_quality=min_base_qual,
                                       min_mapping_quality=min_mapping_qual):
                    if pileup.pos == pos-1:
                        bases = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
                        total_depth = 0
                        
                        for read in pileup.pileups:
                            if not read.is_del and not read.is_refskip:
                                base = read.alignment.query_sequence[read.query_position].upper()
                                qual = read.alignment.query_qualities[read.query_position]
                                if qual >= min_base_qual and read.alignment.mapping_quality >= min_mapping_qual:
                                    if base in bases:
                                        bases[base] += 1
                                        total_depth += 1
                        
                        if total_depth >= 20:
                            # Store allele fractions for major and minor alleles
                            sorted_bases = sorted(bases.items(), key=lambda x: x[1], reverse=True)
                            if len(sorted_bases) >= 2 and sorted_bases[1][1] > 0:
                                site_key = f"{chrom}:{pos}"
                                allele_fractions[site_key] = {
                                    'major': sorted_bases[0][0],
                                    'minor': sorted_bases[1][0],
                                    'major_frac': sorted_bases[0][1] / total_depth,
                                    'minor_frac': sorted_bases[1][1] / total_depth
                                }
        
        bam.close()
        return allele_fractions
    
    def compare_samples(self, sample1_alleles: Dict, sample2_alleles: Dict) -> Dict:
        """Compare allele patterns between two samples to detect potential contamination."""
        shared_sites = set(sample1_alleles.keys()) & set(sample2_alleles.keys())
        if len(shared_sites) < 100:
            return {
                'correlation': 0,
                'shared_sites': len(shared_sites),
                'potential_contamination': False,
                'confidence': 'Low',
                'warning': 'Insufficient shared sites'
            }
        
        correlations = []
        contamination_signals = []
        
        for site in shared_sites:
            s1 = sample1_alleles[site]
            s2 = sample2_alleles[site]
            
            # Check if the samples share the same major/minor alleles
            if (s1['major'] == s2['major'] and s1['minor'] == s2['minor']):
                correlations.append(abs(s1['minor_frac'] - s2['minor_frac']))
                
                # Look for intermediate allele fractions suggesting contamination
                if (0.1 <= s1['minor_frac'] <= 0.4 and 0.1 <= s2['minor_frac'] <= 0.4) or \
                   (0.6 <= s1['minor_frac'] <= 0.9 and 0.6 <= s2['minor_frac'] <= 0.9):
                    contamination_signals.append(1)
                else:
                    contamination_signals.append(0)
        
        if not correlations:
            return {
                'correlation': 0,
                'shared_sites': len(shared_sites),
                'potential_contamination': False,
                'confidence': 'Low',
                'warning': 'No comparable sites'
            }
        
        mean_correlation = 1 - np.mean(correlations)
        contamination_score = np.mean(contamination_signals) if contamination_signals else 0
        
        return {
            'correlation': round(mean_correlation, 3),
            'shared_sites': len(shared_sites),
            'potential_contamination': contamination_score > 0.1,
            'contamination_score': round(contamination_score * 100, 2),
            'confidence': 'High' if len(shared_sites) >= 500 else 'Medium' if len(shared_sites) >= 200 else 'Low'
        }
    
    def run_batch_analysis(self, snp_sites: str, min_base_qual: int = 20,
                          min_mapping_qual: int = 20, threads: int = 1) -> Dict:
        """Run contamination analysis on all samples and compare them."""
        print(f"Analyzing {len(self.bam_files)} samples...")
        
        # First get individual sample metrics
        individual_results = {}
        with ProcessPoolExecutor(max_workers=threads) as executor:
            future_to_bam = {
                executor.submit(self.analyze_single_sample, bam, snp_sites,
                              min_base_qual, min_mapping_qual): bam
                for bam in self.bam_files
            }
            
            for future in as_completed(future_to_bam):
                result = future.result()
                individual_results[result['sample']] = result.get('results', {'error': 'Analysis failed'})
        
        # Then get allele fractions for cross-sample comparison
        print("Collecting allele fractions for cross-sample comparison...")
        sample_alleles = {}
        with ProcessPoolExecutor(max_workers=threads) as executor:
            future_to_bam = {
                executor.submit(self.get_sample_alleles, bam, snp_sites,
                              min_base_qual, min_mapping_qual): bam
                for bam in self.bam_files
            }
            
            for future in as_completed(future_to_bam):
                bam = future_to_bam[future]
                sample_alleles[os.path.basename(bam)] = future.result()
        
        # Compare all sample pairs
        print("Comparing sample pairs...")
        cross_sample_results = []
        for sample1, sample2 in combinations(sample_alleles.keys(), 2):
            comparison = self.compare_samples(sample_alleles[sample1], sample_alleles[sample2])
            if comparison['potential_contamination']:
                cross_sample_results.append({
                    'sample1': sample1,
                    'sample2': sample2,
                    'comparison_results': comparison
                })
        
        return {
            'individual_results': individual_results,
            'cross_sample_contamination': cross_sample_results
        }

def main():
    parser = argparse.ArgumentParser(
        description="Analyze potential cross-sample contamination in a batch of BAM files"
    )
    parser.add_argument("--bam-dir", required=True, help="Directory containing BAM files")
    parser.add_argument("--reference", required=True, help="Reference genome FASTA")
    parser.add_argument("--snps", required=True, help="BED file with common SNP positions")
    parser.add_argument("--min-base-qual", type=int, default=20,
                       help="Minimum base quality (default: 20)")
    parser.add_argument("--min-mapping-qual", type=int, default=20,
                       help="Minimum mapping quality (default: 20)")
    parser.add_argument("--threads", type=int, default=1,
                       help="Number of threads to use (default: 1)")
    parser.add_argument("--output", help="Output JSON file (optional, default: stdout)")
    
    args = parser.parse_args()
    
    # Run batch analysis
    analyzer = BatchContaminationAnalyzer(args.bam_dir, args.reference)
    results = analyzer.run_batch_analysis(
        args.snps,
        args.min_base_qual,
        args.min_mapping_qual,
        args.threads
    )
    
    # Format results for output
    output_data = {
        'analysis_parameters': {
            'min_base_quality': args.min_base_qual,
            'min_mapping_quality': args.min_mapping_qual,
            'number_of_samples': len(analyzer.bam_files)
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
