# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

setup (
       name='funregulation',
       version='0.3',
       date='02/05/2020',
       packages=find_packages(),

       # Declare your packages' dependencies here, for eg:
       install_requires=['Biopython','Biopython','suds.jurko','natsort','pysqlite3','bcbio-gff','psutil'],

       # Fill in these to make your Egg ready for upload to
       # PyPI
       author='arlenz',
       author_email='arlenz@ucs.br',
       keywords = "bioinformatics",
       description='Gene regulatory networks (GRN) of Penicillium ucsensis 2HH and Penicillium oxalicum 114-2 inferred by a computational biology approach',

       #summary = 'Just another Python package for the cheese shop',
       url='https://github.com/alexandrelenz/funregulation.git',
        license='BSD 2-clause',
        long_description=read('README.md'),

        # could also include long_description, download_url, classifiers, etc.
        include_package_data=True,
        platforms='any',
        test_suite='tests',
        classifiers = [
        'Programming Language :: Python',
        'Natural Language :: Portuguese',
        'Intended Audience :: Bioinformatics Users',
        'License :: BSD 2-clause "Simplified" License',
        'Operating System :: OS Independent',
        "Topic :: Scientific/Engineering :: Bioinformatics",
        'Topic :: Software Development :: Libraries :: Python Modules',
        ],
  
       )
