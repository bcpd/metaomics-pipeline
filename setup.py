import io
from os.path import dirname, join
from setuptools import setup

__author__  = "Levy-Booth, Cardenas, Dimitriu"
__license__ = "BSD-3"

# read the contents of your README file
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()


setup(
    name="M5P",
    version=1.0,
    url="https://github.com/MicrobiomeInsights/metaomics-pipeline",
    license=__license__,
    author=__author__,
    zip_safe=False,
    description="Snakemake-based workflow for the taxonomic and functional annotation of metagenomics and metatranscriptomics datasets.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=["M5P"],
    package_data={
        "": [
            "M5P/*",
        ]
    },
    data_files=[(".", ["README.md"])],
    include_package_data=True,
    install_requires=[],
    # install via conda: atlas
    entry_points={"console_scripts": ["M5P = M5P.m5p:main"]},
    classifiers=["Topic :: Scientific/Engineering :: Bio-Informatics"],
)