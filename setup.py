from setuptools import setup, find_packages
import os

# Read the contents of your README file
with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

VERSION = '0.1.1'
DESCRIPTION = 'Preparing data for regression discontinuity design'

setup(
    name="geoRDDprep",
    version=VERSION,
    author="Shahir Shamim",
    author_email="shahir15314@gmail.com",
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/shahirshamim/geoRDDprep", # Placeholder, good practice to have
    packages=find_packages(),
    install_requires=[
        'geopandas', 
        'numpy', 
        'pandas', 
        'shapely',
        'scipy'
    ],
    keywords=['python', 'regression', 'discontinuity', 'geographic', 'spatial', 'gis'],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)