#! /usr/bin/bash

this_script_location="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
delphi_path="/home/djmaity/.opt/Delphi_Linux/executable/delphicpp_openmp"


parameter_file() {
    pdb_file=$1
    filename=${pdb_file##*/}
    output_filename=${filename%.*}
    output_directory=$2

read -r -d '' parameter_file <<- EOM
in(pdb,file="${pdb_file}")
in(siz,file="${this_script_location}/charmm.siz")
in(crg,file="${this_script_location}/charmm.crg")
salt=0.10
exdi=80
linit=2000
maxc=0.0000000001
out(phi, file="${output_directory}/${output_filename}.cub", format="cube")

EOM
    echo -e "$parameter_file" > "tmp.prm"
}

mkdir $2

if [[ -f $1 ]]
then
    parameter_file $1 $2
    "${delphi_path}" "tmp.prm"
fi

if [[ -d $1 ]]
then
    for pdb_file in $1/*.pdb
    do
        parameter_file $pdb_file $2
        "${delphi_path}" "tmp.prm"
    done
fi

rm tmp.prm
