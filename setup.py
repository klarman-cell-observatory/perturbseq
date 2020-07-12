import sys

from pathlib import Path

from setuptools import setup, find_packages


setup(
    name='perturbseq',
    description='Perturb-seq package',
    url='https://github.com/klarman-cell-observatory/perturbseq',
    author='Oana Ursu',
    author_email='oursu@broadinstitute.org',
    license='MIT-License',
    python_requires='>=3.6',
    packages=find_packages(),
    include_package_data=True,
)
