from setuptools import setup, find_packages

setup(
    name="crassify",
    version="0.1.0",
    description="Protein-based viral taxonomy tool",
    author="Your Name",
    url="https://github.com/linda5mith/crassify",
    packages=find_packages(),
    install_requires=[
        "snakemake",
        "pandas",
        "tqdm",
        "rich",
        "pyyaml",
    ],
    entry_points={
        "console_scripts": [
            "crassify=crassify.cli:main",
        ],
    },
    python_requires=">=3.8",
)
