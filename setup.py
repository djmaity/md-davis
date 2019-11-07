#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['Click>=6.0',
                'biopython',
                'docopt',
                'matplotlib',
                'mdtraj',
                'plotly',
                'numpy',
                'pandas',
                'h5py',
]

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', ]

setup(
    author="Dibyajyoti Maity",
    author_email='djdibs@gmail.com',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        "Operating System :: OS Independent",
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="A package to analyze and visualize molecular dynamics simulations of proteins.",
    entry_points={
        'console_scripts': [
            'md_davis=md_davis.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    long_description_content_type="text/x-rst",
    include_package_data=True,
    keywords='analysis data visualization molecular dynamics protein',
    name='md_davis',
    packages=find_packages(include=['md_davis']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://djmaity.github.io/md_davis/',
    version='0.1.0',
    zip_safe=False,
    python_requires='>=3.6',
)
