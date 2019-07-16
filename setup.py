import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="md_davis",
    version="0.0.1",
    author="Dibyajyoti Maity",
    author_email="djdibs@gmail.com",
    description="A package to analyze and visuazlize molecular dynamics simulations of proteins.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    # url="https://github.com/djmaity/md_davis",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)