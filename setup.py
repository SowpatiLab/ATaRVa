import os
import sys
from setuptools import setup, find_packages
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/ATARVA')
from ATARVA import __version__
setup(
    name='atarva',
    version=__version__,
    description='ATaRVa',
    long_description='Analysis of Tandem Repeats Variations',
    url='https://github.com/sowpatilab/ATaRVa.git',
    author='Divya Tej Sowpati',
    author_email='tej@ccmb.res.in',
    license='',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Programming Language :: Python :: 3.9',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    packages=find_packages(),
    install_requires = [
        'python>=3.9.18'
        'pysam>=0.22.0',
        'numpy>=1.26.4',
        ],
   entry_points={
    'console_scripts': ['ATARVA=ATARVA.core:main']
   },
    )
