#!/usr/bin/env python3

import pysam
import argparse
from collections import defaultdict
import numpy as np
import json
from typing import Dict, List, Tuple
import pandas as pd

class BamMetricsCalculator:
    """Calculate comprehensive quality metrics from BAM files."""
    
    def __init__(self, bam_path: str):
        """Initialize with path to BAM file."""
        self.bam_path = bam_path
        self.bam = pysam.AlignmentFile(bam_path, "rb")
        
    def calculate_basic_stats(self) -> Dict:
        """Calculate basic alignment statistics including sequencing and mapping metrics."""
        total_reads = 0
        mapped_reads = 0
        seq_reads = 0  # Reads that pass basic sequencing QC
        
        for read in self.bam.fetch(until_eof=True):  # until_eof=True to include unmapped reads
            total_reads += 1
            
            # Count reads that pass basic sequencing QC
            if not read.is_qcfail and not read.is_secondary:
                seq_reads += 1
                
            if not read.is_unmapped:
                mapped_reads += 1
        
        percentage_lost = ((total_reads - seq_reads) / total_reads * 100) if total_reads > 0 else 0
        percentage_mapping = (mapped_reads / seq_reads * 100) if seq_reads > 0 else 0
        
        return {
            "total_reads": total_reads,
            "seq_reads": seq_reads,
            "percentage_lost": percentage_lost,
            "reads_mapped": mapped_reads,
            "percentage_read_mapping": percentage_mapping
        }
    
    def calculate_target_metrics(self, bed_file: str) -> Dict:
        """Calculate on-target metrics using a BED file of target regions."""
        if not bed_file:
            return {}
            
        # Read target regions
        targets = defaultdict(list)
        with open(bed_file) as f:
            for line in f:
                chrom, start, end = line.strip().split("\t")[:3]
                targets[chrom].append((int(start), int(end)))
        
        # Initialize counters
        total_on_target = 0
        total_on_target_filtered = 0  # Reads passing additional quality filters
        coverage_by_position = defaultdict(int)
        
        # Process reads
        for chrom in targets:
            for start, end in targets[chrom]:
                for pileup in self.bam.pileup(chrom, start, end, truncate=True):
                    pos = pileup.pos
                    coverage_by_position[f"{chrom}:{pos}"] = pileup.n
                    
                    for read in pileup.pileups:
                        if not read.is_del and not read.is_refskip:
                            total_on_target += 1
                            # Count high-quality reads
                            if (not read.alignment.is_duplicate and 
                                read.alignment.mapping_quality >= 20):
                                total_on_target_filtered += 1
        
        # Calculate coverage metrics
        coverages = list(coverage_by_position.values())
        mean_coverage = np.mean(coverages) if coverages else 0
        
        # Calculate percentage of bases above thresholds
        total_positions = len(coverage_by_position)
        bases_above_250 = sum(1 for cov in coverages if cov >= 250)
        bases_above_1000 = sum(1 for cov in coverages if cov >= 1000)
        
        # Calculate exons below thresholds
        exons_below_250 = []
        exons_below_150 = []
        
        for chrom in targets:
            for start, end in targets[chrom]:
                exon_coverages = []
                for pos in range(start, end):
                    key = f"{chrom}:{pos}"
                    if key in coverage_by_position:
                        exon_coverages.append(coverage_by_position[key])
                
                if exon_coverages:
                    mean_exon_cov = np.mean(exon_coverages)
                    if mean_exon_cov < 250:
                        exons_below_250.append(f"{chrom}:{start}-{end}")
                    if mean_exon_cov < 150:
                        exons_below_150.append(f"{chrom}:{start}-{end}")
        
        # Calculate percentages
        basic_stats = self.calculate_basic_stats()
        seq_reads = basic_stats["seq_reads"]
        
        return {
            "total_on_target_reads": total_on_target,
            "percentage_on_target": (total_on_target / seq_reads * 100) if seq_reads > 0 else 0,
            "total_on_target_filter_reads": total_on_target_filtered,
            "percentage_usable": (total_on_target_filtered / seq_reads * 100) if seq_reads > 0 else 0,
            "mean_coverage": mean_coverage,
            "percent_bases_above_250": (bases_above_250 / total_positions * 100) if total_positions > 0 else 0,
            "percent_bases_above_1000": (bases_above_1000 / total_positions * 100) if total_positions > 0 else 0,
            "exons_below_250": exons_below_250,
            "exons_below_150": exons_below_150
        }
    
    def calculate_all_metrics(self, bed_file: str = None) -> Dict:
        """Calculate all available metrics."""
        basic_stats = self.calculate_basic_stats()
        target_metrics = self.calculate_target_metrics(bed_file)
        
        metrics = {
            "basic_stats": basic_stats,
            "target_metrics": target_metrics
        }
        return metrics
    
    def close(self):
        """Close the BAM file."""
        self.bam.close()

def main():
    parser = argparse.ArgumentParser(description="Calculate quality metrics from BAM files")
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--bed", help="BED file with target regions")
    parser.add_argument("--output", help="Output JSON file (optional, default: stdout)")
    parser.add_argument("--tsv", help="Output TSV file with key metrics (optional)")
    
    args = parser.parse_args()
    
    # Calculate metrics
    calculator = BamMetricsCalculator(args.bam)
    metrics = calculator.calculate_all_metrics(args.bed)
    calculator.close()
    
    # Output JSON results
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(metrics, f, indent=2)
    else:
        print(json.dumps(metrics, indent=2))
    
    # Output TSV format if requested
    if args.tsv:
        basic = metrics["basic_stats"]
        target = metrics["target_metrics"]
        
        # Write main metrics
        with open(args.tsv, 'w') as f:
            # Write key metrics
            df = pd.DataFrame({
                "Metric": [
                    "TotalReads", "SeqReads", "PercentageLost", "ReadsMapped", 
                    "PercentageReadMapping", "TotalOnTargetReads", "PercentageOnTarget",
                    "TotalOnTargetFilterReads", "PercentageUsable", "MeanCoverage",
                    "PercentBasesAbove250", "PercentBasesAbove1000",
                    "ExonsBelow250x", "ExonsBelow150x"
                ],
                "Value": [
                    basic["total_reads"],
                    basic["seq_reads"],
                    basic["percentage_lost"],
                    basic["reads_mapped"],
                    basic["percentage_read_mapping"],
                    target.get("total_on_target_reads", "NA"),
                    target.get("percentage_on_target", "NA"),
                    target.get("total_on_target_filter_reads", "NA"),
                    target.get("percentage_usable", "NA"),
                    target.get("mean_coverage", "NA"),
                    target.get("percent_bases_above_250", "NA"),
                    target.get("percent_bases_above_1000", "NA"),
                    len(target.get("exons_below_250", [])),
                    len(target.get("exons_below_150", []))
                ]
            })
            df.to_csv(f, sep='\t', index=False)
        
        # Write exons below thresholds to separate files if they exist
        if target.get("exons_below_250"):
            with open(args.tsv + ".exons_below_250.txt", 'w') as f:
                f.write("\n".join(target["exons_below_250"]))
                
        if target.get("exons_below_150"):
            with open(args.tsv + ".exons_below_150.txt", 'w') as f:
                f.write("\n".join(target["exons_below_150"]))

if __name__ == "__main__":
    main()
