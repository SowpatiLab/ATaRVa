import os
import sys
from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext
import subprocess
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/ATARVA')

ssw_dir = "ATARVA/Complete-Striped-Smith-Waterman-Library"
git_repo = "https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library.git"

class CustomBuildExt(build_ext):
    def run(self):
        os.makedirs("ATARVA", exist_ok=True)

        if not os.path.exists(ssw_dir):
            subprocess.run(["git", "clone", git_repo, ssw_dir], check=True)

        subprocess.run(
            ["gcc", "-Wall", "-O3", "-pipe", "-fPIC", "-shared", "-rdynamic", "-o", "libssw.so", "ssw.c"],
            cwd=os.path.join(ssw_dir, "src"),
            check=True
        )


        os.rename(
            os.path.join(ssw_dir),
            os.path.join("ATARVA/Complete_Striped_Smith_Waterman_Library")
        )

        super().run()

from ATARVA import __version__
setup(
    name='ATaRVa',
    version=__version__,
    description='ATaRVa',
    long_description='Analysis of Tandem Repeats Variations',
    url='https://github.com/sowpatilab/ATaRVa.git',
    author='Divya Tej Sowpati',
    author_email='tej@ccmb.res.in',
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Programming Language :: Python :: 3.9',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    packages=find_packages(),
    install_requires = [
        'numpy>=1.26.4',
        'pyabpoa>=1.5.3',
        'pysam>=0.22.1',
        'scikit-learn>=1.4.1',
        'scipy>=1.12.0',
        'threadpoolctl>=3.3.0',
        ],
    python_requires='>=3.9.0',
    cmdclass={"build_ext":CustomBuildExt},
    entry_points={
    'console_scripts': ['ATARVA=ATARVA.core:main']
    },
    )