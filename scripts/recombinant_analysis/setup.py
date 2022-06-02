from setuptools import setup, find_packages
import sys

setup(
    name='ramifi',
    version= '0.1.0',
    author='Chienchi Lo',
    author_email='chienchi@lanl.gov',
    packages=find_packages(),
    scripts=['ramifi/ramifi'],
    url='https://github.com/chienchi/ramifi',
    license='LICENSE',
    package_data={'ramifi.data':['*']},
    include_package_data=True,
    description='Script to do recombinant read analysis',
    keywords="",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    install_requires=[
        "plotly >=5.5.0",
        "pandas >= 1.2.4",
        "pysam >= 0.16.0.1",
        "kaleido >= 0.2.1",
        'biopython == 1.78',
        'importlib-resources==5.7.1',
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bioinformatics",
        "Operating System :: OS Independent",
    ],
)
