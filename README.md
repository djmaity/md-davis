# MD Davis

## Installation
Clone this repository using:
```
git clone git@github.com:djmaity/md_davis.git
```
Optional: I recommend installing this package in a virtual environment. Create a virtual enviroment in Anaconda or Miniconda distribution using:
```
conda create --name md_davis python=3
```
Activate the environment before running the install command
```
conda activate md_davis
```
Install the dependencies using the following commands:
```
conda install scipy
conda install psutil
conda install -c schrodinger pymol
conda install -c conda-forge mdtraj
conda install -c plotly plotly-orca
conda update -c conda-forge h5py
```

Install this package using pip:
```
pip install path/to/md_davis
```
### Dependencies:
pip should automatically install the dependencies

In addtion the following external software are required as well:
[PyMOL] (https://pymol.org/2/)
[Delphi](http://compbio.clemson.edu/delphi)
[MSMS](http://mgltools.scripps.edu/downloads#msms)

This page is still under construction. Visit [https://djmaity.github.io/md_davis/](https://djmaity.github.io/md_davis/) to know more.


<!-- ## Plot .xvg file
To use this script type the following command in a terminal or command prompt and press 'Enter':
```
python plot_xvg.py <path/to/file.xvg>
```
Replace `<path/to/file.xvg>` with the location of your `.xvg` file.

## Plot DSSP file obtained from GROMACS
To obtain the input file for the script, run **do_dssp** command from GROMACS with the **-ssdump** option:
```
gmx do_dssp -f <trajectory> -s <structure> -o <ss.xpm> -ssdump <ssdump.dat>
```
Use the following script to count and plot the percentage secondary structure per residue throughout the trajectory:
```
python plot_do_dssp_per_residue.py <path/to/ssdump.dat>
```
Replace `<path/to/ssdump.dat>` depending on the location of your `.dat` file. This will show the plot on screen. To output the plot as an image instead of displaying it on screen use:
```
python plot_do_dssp_per_residue.py <path/to/ssdump.dat> -o image.png
```
### Getting Help
Providing **-h** option to each python script will print out its help message.
```
python <script_name.py> -h
```
NOTE: _Please replace the text wtih angular brackets < > by the respective filename or path._ 

# README #

These scripts extract dihedral angles from GROMACS (molecular dynamics) trajectory and analyse enregy landscapes.

### What is this repository for? ###

* Analysing GROMACS trajectories
* 0.1


* Database configuration
* How to run tests
* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Author: Dibyajyoti Maity (djmaity@ssl.serc.iisc.in)

-->