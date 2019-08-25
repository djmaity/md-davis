import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="md_davis",
    version="0.0.2",
    author="Dibyajyoti Maity",
    author_email="djdibs@gmail.com",
    description="A package to analyze and visuazlize molecular dynamics simulations of proteins.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://djmaity.github.io/md_davis/",
    packages=setuptools.find_packages(),
    install_requires=['biopython','docopt','matplotlib','mdtraj','plotly','pymol'],
    entry_points={
        'console_scripts': [
            'md_davis=md_davis.cli:main'
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)