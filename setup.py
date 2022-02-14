from setuptools import setup, find_packages
import os

requires = [
    "pandas",
    "scikit-image",
    "numpy",
    "scipy",
    "scikit-learn",
    "tables",
    "umap-learn"
]

DESCRIPTION = "De-noising Spatial Transcriptomics Data Based on Position and Image Information"

setup(
    name="sprod",
    version=1.0,
    description=DESCRIPTION,
    url="https://github.com/yunguan-wang/QBRC-SPROD",
    author="Yunguan Wang",
    author_email="yunguan.wang@utsouthwestern.edu",
    license="MIT",
    packages=find_packages(),
    install_requires=requires,
    zip_safe=False,
)

