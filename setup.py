import glob
from setuptools import setup, find_packages

setup(
    name='seroba',
    version='0.1.4',
    description='SEROBA: Serotyping for illumina reads',
    packages = find_packages(),
    author='Lennard Epping',
    author_email='path-help@sanger.ac.uk',
    url='',
    scripts=glob.glob('scripts/*'),
    test_suite='nose.collector',
    tests_require=['nose >= 1.3'],
    install_requires=[
        'ariba >= 2.9.1',
        'pymummer>=0.10.2',
        'PyYAML>=3.12',
        'biopython>=1.68',
        'pyfastaq>=3.15.0'
    ],
    license='GPLv3',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],
)
