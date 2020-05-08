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

       #summary = 'Just another Python package for the cheese shop',
       url='https://github.com/alexandrelenz/funregulation.git',
       license='BSD 2-Clause License',
       long_description='Long description of the package',

       # could also include long_description, download_url, classifiers, etc.

  
       )
