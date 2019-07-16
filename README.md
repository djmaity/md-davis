# gromacs_analysis_scripts

## Dependencies:
python 3 with numpy and matplotlib modules installed.

## Plot .xvg file
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