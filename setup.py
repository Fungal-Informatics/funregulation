from __future__ import print_function
from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand
import io
import codecs
import os
import sys

import funregulation

here = os.path.abspath(os.path.dirname(__file__))

def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(filename, encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)

long_description = read('README.md')

class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest
        errcode = pytest.main(self.test_args)
        sys.exit(errcode)

setup(
    name='funregulation',
    version=funregulation.__version__,
    url='https://github.com/alexandrelenz/funregulation.git',
    license='BSD 2-clause',
    author='Alexandre Lenz',
    tests_require=['pytest'],
    install_requires=['biopython>=1.65'],
    cmdclass={'test': PyTest},
    author_email='arlenz@ucs.br',
    keywords = "bioinformatics",
    description='Pipeline to predict and annotate fungal promoter elements and transcription factor binding sites',
    long_description=long_description,
    packages=['funregulation'],
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
    extras_require={
        'testing': ['pytest'],
    }
)
