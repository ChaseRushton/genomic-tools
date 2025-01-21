from setuptools import setup, find_packages

setup(
    name="genomic_tools",
    version="0.1.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        "pysam>=0.22.0",
        "argparse>=1.4.0",
    ],
    author="Your Name",
    author_email="your.email@example.com",
    description="A collection of tools for genomic data analysis and validation",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/genomic-tools",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    entry_points={
        "console_scripts": [
            "vcf-validator=genomic_tools.vcf_validator:main",
        ],
    },
)
