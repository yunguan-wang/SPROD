from setuptools import setup, find_packages
import os

requires = [
    "pandas",
    "scikit-image",
    "numpy",
    "scipy",
    "scikit-learn",
    "hdbscan",
    "tables"
]

DESCRIPTION = "De-noising of Spatial Expression Profiling Data Based on Latent Graph Learning of in situ Position and Image Data"

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

