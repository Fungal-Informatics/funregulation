# -*- coding: utf-8 -*-
from __future__ import print_function
from setuptools.command.test import test as TestCommand
import io
import codecs
import os
import sys
from setuptools import setup, find_packages

here = os.path.abspath(os.path.dirname(__file__))

def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(filename, encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)

setup (
       name='funregulation',
       version='0.1',
       packages=find_packages(),

       # Declare your packages' dependencies here, for eg:
       install_requires=['foo>=3'],

        # Fill in these to make your Egg ready for upload to
        # PyPI
        author='Alexandre Lenz',
        author_email='arlenz@ucs.br',
        keywords = "bioinformatics",
        description='Pipeline to predict and annotate fungal promoter elements and transcription factor binding sites',

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