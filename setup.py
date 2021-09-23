#!/usr/bin/env python
# -*- coding: utf-8 -*-

import io
import os
import re
from glob import glob
from setuptools import find_packages
from setuptools import setup


def read(*names, **kwargs):
    with io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get('encoding', 'utf8')
    ) as fh:
        return fh.read()


setup(
    name='md_davis',
    version='0.3.0',
    description='A tool for comparative analysis of molecular dynamics '
                'simulations of proteins.',
    long_description='%s' % (re.compile('^.. start-badges.*^.. end-badges', re.M | re.S).sub('', read('README.md'))),
    long_description_content_type='text/markdown',
    author='Dibyajyoti Maity',
    author_email='djdibs@gmail.com',
    url='https://github.com/djmaity/md-davis',
    packages=find_packages(),
    py_modules=[os.path.splitext(os.path.basename(path))[0] for path in glob('src/*.py')],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
    ],
    project_urls={
        'Documentation': 'https://md-davis.readthedocs.io/',
        'Changelog': 'https://md-davis.readthedocs.io/en/latest/changelog.html',
        'Issue Tracker': 'https://github.com/djmaity/md-davis/issues',
    },
    license='MIT license',
    keywords=[
        'analysis', 'data visualization', 'molecular dynamics', 'protein'
    ],
    python_requires='>=3.7',
    install_requires=['biopandas',
                      'biopython',  # TODO: Refactor dependency
                      'click',
                      'docopt',  # TODO: Refactor dependency
                      'h5py>=3.0',
                      'matplotlib',
                      'more_itertools',
                      'numpy',  # numpy should be before mdtraj
                      'pandas',
                      'plotly>=4.0',
                      'scikit-learn',
                      'scipy',
                      'toml',
                      'mdtraj',
                      # 'pymol', # Cannot Install automatically using pip
                      ],
    include_package_data=True,
    setup_requires=['flake8', 'pytest-runner'],
    test_suite='tests',
    tests_require=['pytest'],
    zip_safe=False,  # TODO: check if it can be made zip safe
    entry_points={
        'console_scripts': [
            'md_davis=md_davis.cli:main',
            'md-davis=md_davis.cli:main',
        ],
    },
)
