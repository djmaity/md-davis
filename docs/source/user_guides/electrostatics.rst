How to calculate surface electrostatics using MD DaVis?
=======================================================

Delphi and MSMS must be installed for this to work. You may download these from: [Delphi](http://compbio.clemson.edu/delphi)  and [MSMS](http://mgl.scripps.edu/people/sanner/html/msms_home.html)

You will need the parameter files for delphi run, which can be downloaded from [paramter files](http://compbio.clemson.edu/downloadDir/delphi/parameters.tar.gz)

**Step 1:** Calculate the electrostatics using the following command

.. code-block:: bash
    
    md_davis electrostatics --surface_potential \
        -m ~/.opt/msms/ \
        -d ~/.opt/delphi8.4.2_sequential
        -c path/to/charge/file.crg \
        -r path/to/radius/file.siz \
        PDB_FILE.pdb \
        path/to/electrostatics/output/directory

The following charge and radius files are provided with the source code for MD&nbsp;DaVis.

.. code-block:: bash
    
    md_davis/md_davis/electrostatics/charmm.crg
    md_davis/md_davis/electrostatics/charmm.siz

If you receive a warning during Delphi run regarding missing charge or radius. Then the missing properties must be added to these files or whichever files you provide to `md_davis electrostatics`.

**Step 2:** Add the surface electrostatics data into the residue dataframe by modifying the command for creating residue dataframe to include the option `-d` and specifying the path to the directory containing the surface potential files. This will search for all  `.pot` files and include their mean and total in the residue dataframe.

.. code-block:: bash

    md_davis residue dataframe \
            --potential  path/to/electrostatics/output/directory \
            --annotation  annotations.json \
            --prefix  prefix \
            md_davis_collect_output_data.h5  \
            output_residue_wise_data.p

**Step 3:** Plot the residue dataframe with the usual command:

.. code-block:: bash
    
    md_davis plot residue output_residue_wise_data.p
