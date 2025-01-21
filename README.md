# Genomic Tools

A comprehensive suite of Python tools for genomic data quality control and validation.

## Overview

This toolkit provides essential utilities for analyzing and validating genomic data, with a focus on quality control and sample verification. It includes tools for BAM file analysis, variant validation, contamination detection, and sex verification.

## Table of Contents
- [Features](#features)
- [Installation](#installation)
- [Requirements](#requirements)
- [Usage](#usage)
  - [VCF Validator](#vcf-validator)
  - [BAM Quality Metrics](#bam-quality-metrics)
  - [Contamination Estimator](#contamination-estimator)
  - [Batch Contamination Analyzer](#batch-contamination-analyzer)
  - [Sex Verification](#sex-verification)
- [Common Workflows](#common-workflows)
- [Contributing](#contributing)
- [License](#license)

## Features

### Core Tools

- **VCF Validator**
  - Validates variants against BAM evidence
  - Checks read depth and mapping quality
  - Verifies base quality
  - Confirms allele frequency
  - Ideal for variant call verification

- **BAM Quality Metrics**
  - Comprehensive BAM file analysis
  - Coverage statistics and uniformity
  - Mapping quality distribution
  - Read depth analysis
  - Exon coverage verification
  - Outputs clear, tabulated results

- **Contamination Estimator**
  - Sample-level contamination detection
  - Uses common SNP sites
  - Calculates allele fraction patterns
  - Provides confidence scores
  - Detailed contamination metrics

- **Batch Contamination Analyzer**
  - Multi-sample contamination detection
  - Cross-sample comparison
  - Parallel processing support
  - Detailed pair-wise analysis
  - Identifies potential sample swaps

- **Sex Verification**
  - Multi-method sex determination
  - Coverage-based analysis
  - SRY gene detection
  - X chromosome heterozygosity
  - Aneuploidy detection (X0, XXY)

## Installation

1. Clone the repository:
```bash
git clone https://github.com/ChaseRushton/genomic-tools.git
cd genomic-tools
```

2. Create a virtual environment (recommended):
```bash
python -m venv venv
source venv/bin/activate  # Linux/Mac
venv\\Scripts\\activate   # Windows
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

## Requirements

### Python Dependencies
- Python ≥ 3.8
- pysam ≥ 0.22.0
- numpy ≥ 1.21.0
- pandas ≥ 1.3.0
- scipy ≥ 1.7.0

### External Requirements
- Samtools (for BAM/CRAM support)
- Reference genome in FASTA format
- Index files (.fai, .bai) for reference and BAM files

## Usage

### VCF Validator

Validate variant calls against supporting evidence in BAM files:

```bash
python src/genomic_tools/vcf_validator.py --bam input.bam --vcf variants.vcf
```

Parameters:
- Required:
  - `--bam`: Input BAM file
  - `--vcf`: VCF file to validate
- Optional:
  - `--min-depth`: Minimum read depth (default: 10)
  - `--min-mapping-quality`: Minimum mapping quality (default: 20)
  - `--min-base-quality`: Minimum base quality (default: 20)
  - `--min-allele-freq`: Minimum allele frequency (default: 0.2)

### BAM Quality Metrics

Generate comprehensive quality metrics for BAM files:

```bash
python src/genomic_tools/bam_metrics.py --bam input.bam --bed targets.bed --tsv output.tsv
```

Key Metrics:
- Basic Statistics
  - Total/Mapped reads
  - Mapping rates
  - Coverage statistics
- Target Analysis
  - On-target rate
  - Coverage uniformity
  - Problematic regions
- Output Formats
  - TSV summary file
  - Detailed JSON report
  - Coverage histograms

### Contamination Analysis

#### Single Sample
```bash
python src/genomic_tools/contamination_estimate.py --bam input.bam --reference ref.fa --snps snps.bed
```

#### Multiple Samples
```bash
python src/genomic_tools/batch_contamination.py --bam-dir /path/to/bams --reference ref.fa --snps snps.bed --threads 4
```

### Sex Verification

Determine biological sex using multiple genomic indicators:

```bash
python src/genomic_tools/sex_verify.py --bam input.bam --reference ref.fa --snps snps.bed
```

Interpretation Guidelines:
| Pattern | XX Female | XY Male | X0 (Turner) | XXY (Klinefelter) |
|---------|-----------|---------|-------------|-------------------|
| X/Auto  | ~1.0     | ~0.5    | ~1.0        | ~1.0             |
| Y/Auto  | ~0       | >0.1    | ~0          | >0.1             |
| SRY     | No       | Yes     | No          | Yes              |
| X-het   | High     | Low     | Low         | High             |

## Common Workflows

### Basic QC Pipeline
```bash
# 1. Run BAM metrics
python src/genomic_tools/bam_metrics.py --bam sample.bam --bed targets.bed --tsv qc.tsv

# 2. Check contamination
python src/genomic_tools/contamination_estimate.py --bam sample.bam --reference ref.fa --snps snps.bed

# 3. Verify sex
python src/genomic_tools/sex_verify.py --bam sample.bam --reference ref.fa --snps snps.bed
```

### Batch Analysis
```bash
# Analyze multiple samples for cross-contamination
python src/genomic_tools/batch_contamination.py --bam-dir /path/to/bams --reference ref.fa --snps snps.bed --threads 4
```

## Contributing

We welcome contributions! Please follow these steps:

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Run tests
5. Submit a pull request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For questions or support, please open an issue on GitHub.
