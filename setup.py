import io
import os
import sys

from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand

import pyRLC

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
    name='pyRLC',
    version=pyRLC.__version__,
    url='http://github.com/axiom/pyRLC/',
    author='Jesse Lopez',
    tests_require=['pytest'],
    install_requires=[
        'numpy>=1.13.1',
        'scipy>=0.19.1',
        'matplotlib>=2.0.2',
        'Pillow>=4.2.1'
        ],
    cmdclass={'test': PyTest},
    author_email='jesse@axiomdatascience.com',
    description='Read and convert .rlc HF Radar files to images',
    long_description=long_description,
    packages=['pyRLC'],
    include_package_data=True,
    platforms='any',
    test_suite='sandman.test.test_sandman',
    classifiers=[
        'Programming Language :: Python',
        'Development Status :: 3 - Alpha',
        'Natural Language :: English',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Visualization',
        ],
    extras_require={
        'testing': ['pytest'],
    }
)
