#! /usr/bin/bash

## To run this script:
## ./run_delphi.sh strcuture.pdb surface.pdb output_directory

# this_script_location="/home/djmaity/Documents/Repositories/protein_surface_electrostatics"
this_script_location="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

delphi_path="/home/djmaity/.opt/Delphi_Linux/executable/delphicpp_openmp"
# delphi_path="/home/djmaity/.opt/delphi8.4.1_sequential"

parameter_file() {
    ## Function input
    ## strcuture.pdb surface.pdb output_directory
    pdb_file=$1
    filename=${pdb_file##*/}
    output_filename=${filename%.*}
    surface=$2
    output_directory=$3

read -r -d '' parameter_file <<- EOM
in(pdb,file="${pdb_file}")
in(siz,file="${this_script_location}/charmm.siz")
in(crg,file="${this_script_location}/charmm.crg")
in(frc,file='${surface}')
acenter(0,0,0)
gsize=201
salt=0.10
exdi=80
linit=2000
maxc=0.0000000001
out(phi, file="${output_directory}/${output_filename}.cub", format="cube")
out(frc, file="${output_directory}/${output_filename}.pot")
site(Atom, Potential, Reaction, Coulomb, Field)

EOM
    echo -e "$parameter_file" > "tmp.prm"
}

mkdir $3

if [[ -f $1 ]]
then
    parameter_file $1 $2 $3
    "${delphi_path}" "tmp.prm"
fi

if [[ -d $1 ]]
then
    for pdb_file in $1/*.pdb
    do
        parameter_file $pdb_file $2 $3
        "${delphi_path}" "tmp.prm"
    done
fi

rm tmp.prm
