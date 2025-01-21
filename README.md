# Genomic Tools

A collection of Python tools for genomic data analysis and validation.

## Features

- **VCF Validator**: Validates variants in VCF files against supporting evidence in BAM files
  - Checks read depth
  - Validates mapping quality
  - Verifies base quality
  - Confirms allele frequency

## Installation

```bash
git clone https://github.com/ChaseRushton/genomic-tools.git
cd genomic-tools
pip install -r requirements.txt
```

## Requirements

- Python 3.x
- pysam >= 0.22.0
- argparse >= 1.4.0

## Usage

### VCF Validator

Validate variants in a VCF file against a BAM file:

```bash
python src/genomic_tools/vcf_validator.py --bam input.bam --vcf variants.vcf
```

Optional parameters:
- `--min-depth`: Minimum read depth required (default: 10)
- `--min-mapping-quality`: Minimum mapping quality for reads (default: 20)
- `--min-base-quality`: Minimum base quality for variant position (default: 20)
- `--min-allele-freq`: Minimum allele frequency required (default: 0.2)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is open source and available under the MIT License.

## Contact

- GitHub: [@ChaseRushton](https://github.com/ChaseRushton)
