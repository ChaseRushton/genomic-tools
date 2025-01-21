#!/usr/bin/env python3

import pysam
import argparse
import numpy as np
from collections import defaultdict
import json
from typing import Dict, List, Tuple
import pandas as pd

class SexVerifier:
    """Verify sample sex using multiple genomic indicators."""
    
    def __init__(self, bam_path: str, reference_path: str):
        """Initialize with paths to BAM and reference files."""
        self.bam = pysam.AlignmentFile(bam_path, "rb")
        self.reference = pysam.FastaFile(reference_path)
        
        # Default SRY region (can be overridden)
        self.sry_region = {
            "chrY": [(2786854, 2787699)]  # GRCh38 coordinates
        }
        
        # PAR regions (pseudoautosomal regions) - GRCh38
        self.par_regions = {
            "chrX": [(10001, 2781479), (155701383, 156030895)],
            "chrY": [(10001, 2781479), (56887903, 57217415)]
        }
        
    def calculate_chromosome_coverage(self, chrom: str, exclude_regions: List[Tuple[int, int]] = None) -> Dict:
        """Calculate average coverage for a chromosome, excluding specified regions."""
        total_depth = 0
        total_bases = 0
        chunk_size = 1000000  # Process in 1Mb chunks
        
        # Get chromosome length
        chrom_length = self.reference.get_reference_length(chrom)
        
        for start in range(0, chrom_length, chunk_size):
            end = min(start + chunk_size, chrom_length)
            
            # Skip excluded regions
            if exclude_regions:
                skip = False
                for ex_start, ex_end in exclude_regions:
                    if start < ex_end and end > ex_start:
                        skip = True
                        break
                if skip:
                    continue
            
            # Calculate coverage
            for pileup in self.bam.pileup(chrom, start, end, truncate=True):
                total_depth += pileup.n
                total_bases += 1
        
        if total_bases == 0:
            return {
                "mean_coverage": 0,
                "bases_covered": 0
            }
            
        return {
            "mean_coverage": total_depth / total_bases,
            "bases_covered": total_bases
        }
    
    def calculate_autosomal_coverage(self) -> float:
        """Calculate mean coverage across autosomes (chr1-22)."""
        autosomal_coverage = []
        
        for chrom in [f"chr{i}" for i in range(1, 23)]:
            try:
                cov = self.calculate_chromosome_coverage(chrom)
                if cov["mean_coverage"] > 0:
                    autosomal_coverage.append(cov["mean_coverage"])
            except Exception:
                continue
        
        return np.mean(autosomal_coverage) if autosomal_coverage else 0
    
    def check_sry_presence(self) -> Dict:
        """Check for presence of SRY gene."""
        total_depth = 0
        covered_bases = 0
        
        for chrom, regions in self.sry_region.items():
            for start, end in regions:
                try:
                    for pileup in self.bam.pileup(chrom, start, end, truncate=True):
                        if pileup.n >= 5:  # Minimum coverage threshold
                            total_depth += pileup.n
                            covered_bases += 1
                except Exception:
                    continue
        
        coverage = total_depth / (end - start) if end > start else 0
        
        return {
            "sry_present": covered_bases > 0.7 * (end - start),  # At least 70% coverage
            "mean_coverage": coverage,
            "covered_bases": covered_bases
        }
    
    def calculate_x_heterozygosity(self, snp_sites: str,
                                 min_base_qual: int = 20,
                                 min_mapping_qual: int = 20) -> Dict:
        """Calculate X chromosome heterozygosity using SNP sites."""
        het_sites = 0
        total_sites = 0
        
        with open(snp_sites) as f:
            for line in f:
                chrom, pos, *_ = line.strip().split('\t')
                if not chrom == "chrX":
                    continue
                    
                pos = int(pos)
                
                # Skip PAR regions
                skip = False
                for start, end in self.par_regions["chrX"]:
                    if start <= pos <= end:
                        skip = True
                        break
                if skip:
                    continue
                
                total_sites += 1
                
                # Analyze site
                for pileup in self.bam.pileup(chrom, pos-1, pos, truncate=True,
                                            min_base_quality=min_base_qual,
                                            min_mapping_quality=min_mapping_qual):
                    if pileup.pos == pos-1:
                        bases = defaultdict(int)
                        total_depth = 0
                        
                        for read in pileup.pileups:
                            if not read.is_del and not read.is_refskip:
                                base = read.alignment.query_sequence[read.query_position].upper()
                                qual = read.alignment.query_qualities[read.query_position]
                                if qual >= min_base_qual and read.alignment.mapping_quality >= min_mapping_qual:
                                    bases[base] += 1
                                    total_depth += 1
                        
                        if total_depth >= 20:  # Minimum depth threshold
                            sorted_bases = sorted(bases.items(), key=lambda x: x[1], reverse=True)
                            if len(sorted_bases) >= 2:
                                major_frac = sorted_bases[0][1] / total_depth
                                if 0.2 <= major_frac <= 0.8:  # Heterozygous site
                                    het_sites += 1
        
        return {
            "heterozygous_sites": het_sites,
            "total_sites": total_sites,
            "heterozygosity_rate": het_sites / total_sites if total_sites > 0 else 0
        }
    
    def determine_sex(self, snp_sites: str = None) -> Dict:
        """Determine sample sex using multiple methods."""
        # Calculate coverage ratios
        auto_cov = self.calculate_autosomal_coverage()
        x_cov = self.calculate_chromosome_coverage("chrX", self.par_regions["chrX"])["mean_coverage"]
        y_cov = self.calculate_chromosome_coverage("chrY", self.par_regions["chrY"])["mean_coverage"]
        
        x_ratio = x_cov / auto_cov if auto_cov > 0 else 0
        y_ratio = y_cov / auto_cov if auto_cov > 0 else 0
        
        # Check SRY
        sry_results = self.check_sry_presence()
        
        # Calculate X heterozygosity if SNP sites provided
        x_het = None
        if snp_sites:
            x_het = self.calculate_x_heterozygosity(snp_sites)
        
        # Determine sex based on multiple indicators
        indicators = {
            "xx_female": x_ratio >= 0.8 and y_ratio < 0.1 and not sry_results["sry_present"],
            "xy_male": 0.4 <= x_ratio <= 0.6 and y_ratio >= 0.1 and sry_results["sry_present"],
            "x0_turner": x_ratio >= 0.8 and y_ratio < 0.1 and not sry_results["sry_present"],
            "xxy_klinefelter": x_ratio >= 0.8 and y_ratio >= 0.1 and sry_results["sry_present"]
        }
        
        # Determine confidence
        confidence_metrics = []
        if x_ratio > 0:
            confidence_metrics.append(1)
        if y_ratio > 0 or y_ratio < 0.05:
            confidence_metrics.append(1)
        if sry_results["mean_coverage"] > 0:
            confidence_metrics.append(1)
        if x_het and x_het["total_sites"] > 100:
            confidence_metrics.append(1)
            
        confidence = "High" if len(confidence_metrics) >= 3 else \
                    "Medium" if len(confidence_metrics) >= 2 else "Low"
        
        # Determine most likely sex
        predicted_sex = None
        for sex, condition in indicators.items():
            if condition:
                predicted_sex = sex
                break
        
        if not predicted_sex:
            predicted_sex = "Inconclusive"
            confidence = "Low"
        
        return {
            "predicted_sex": predicted_sex,
            "confidence": confidence,
            "metrics": {
                "x_autosome_ratio": round(x_ratio, 3),
                "y_autosome_ratio": round(y_ratio, 3),
                "sry_present": sry_results["sry_present"],
                "sry_coverage": round(sry_results["mean_coverage"], 2),
                "x_heterozygosity": round(x_het["heterozygosity_rate"], 3) if x_het else None,
                "autosomal_coverage": round(auto_cov, 2)
            }
        }
    
    def close(self):
        """Close open file handles."""
        self.bam.close()
        self.reference.close()

def main():
    parser = argparse.ArgumentParser(
        description="Verify sample sex using multiple genomic indicators"
    )
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--reference", required=True, help="Reference genome FASTA")
    parser.add_argument("--snps", help="BED file with SNP positions (optional)")
    parser.add_argument("--output", help="Output JSON file (optional, default: stdout)")
    
    args = parser.parse_args()
    
    # Run sex verification
    verifier = SexVerifier(args.bam, args.reference)
    results = verifier.determine_sex(args.snps)
    verifier.close()
    
    # Format results for output
    output_data = {
        'sample_name': args.bam,
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
