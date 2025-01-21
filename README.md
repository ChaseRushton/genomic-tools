# Genomic Tools

A collection of tools for genomic data analysis and validation.

## VCF Validator

The VCF Validator is a Python tool that validates VCF (Variant Call Format) files against their corresponding BAM files. It ensures that variants called in the VCF have sufficient supporting evidence in the alignment data.

### Features

- Validates variants against BAM file evidence
- Checks read depth at variant positions
- Filters by mapping quality and base quality
- Calculates and validates allele frequencies
- Provides detailed validation reports

### Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/genomic-tools.git
cd genomic-tools

# Install dependencies
pip install -r requirements.txt
```

### Usage

```bash
python -m genomic_tools.vcf_validator --bam input.bam --vcf input.vcf
```

Optional parameters:
- `--min-depth`: Minimum read depth required (default: 10)
- `--min-mapping-quality`: Minimum mapping quality (default: 20)
- `--min-base-quality`: Minimum base quality (default: 20)
- `--min-allele-freq`: Minimum allele frequency required (default: 0.2)

### Example

```bash
python -m genomic_tools.vcf_validator \
    --bam sample.bam \
    --vcf variants.vcf \
    --min-depth 15 \
    --min-mapping-quality 30 \
    --min-base-quality 25 \
    --min-allele-freq 0.1
```

## License

MIT License
