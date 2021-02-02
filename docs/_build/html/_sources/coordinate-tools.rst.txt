Coordinates and fragment tools
======================================


#############################
Modifying ASH fragment files
#############################

Convenient commands to modify the coordinates of the QM-region only can be found inside:
/path/to/ashdir/ash/scripts

Scripts are called: fragedit.py  and fragupdate.py

TO BE DOCUMENTED...


###########################
Working with PDB files
###########################

WARNING: PDB files are convenient for visualization purposes and for initial reading in of coordinates but are
generally not a file format to be used (one problem is the limited number of significant digits used
for coordinates).

**Reading in PDB file**

It is possible to read in coordinates from a PDB file to create an ASH fragment file.
This functionality is very basic, it will only read in the coordinates, not atom-types
or residue information. Atomtypes and residue information can be read-in via a PSF-file
by OpenMMTheory (see :doc:`MM-interfaces`).

This option should thus only be used to provide convenient starting coordinates.
.. code-block:: python

    pdbfrag = Fragment(pdbfile="mol.pdb")

**Writing out PDB file**

If you have an ASH fragment file created (loaded into memory), you can request to write out a PDB-file from it via the write_pdbfile function.

.. code-block:: python

    def write_pdbfile(fragment,outputname="ASHfragment", openmmobject=None, atomnames=None,
                    resnames=None,residlabels=None,segmentlabels=None):

An ASH fragment file needs to always be provided, and then optionally the outputname ("ASHfragment.pdb" will be created by default).


- Example 1 (dummy):

.. code-block:: python

    write_pdbfile(frag)

This will give you a PDB-file with the coordinates taken from inside the ASH fragment (here called frag) but without residue information (since none was provided).
All residues will be labelled 'DUM' and segments 'SEG', element information should be correct.

- Example 2 (manual correct specification):

.. code-block:: python

    openmmobject = OpenMMTheory(psffile=psffile, CHARMMfiles=True, charmmtopfile=topfile,charmmprmfile=parfile,
                    printlevel=1, platform='CPU' )
    write_pdbfile(frag, outputname="manual", atomnames=openmmobject.atomnames, resnames=openmmobject.resnames,
        residlabels=openmmobject.resids,segmentlabels=openmmobject.segmentnames)

Here the residue information is provided via keyword arguments and the information taken from an ASH OpenMMTheory object, previously created.
The residue information is present in openmmobject as it was read from the CHARMM PSF-file.
Could also be done completely manually if desired.

- Example 3 (simple and recommended way):

.. code-block:: python

    write_pdbfile(frag, outputname="simple",openmmobject=openmmobject)

Here an ASH openMMtheory object is provided to the function (defined like before) and the function will grab the information from it. It should then print a correct PDB-file with the residue, atom and segment information from the ASH OpenMM object. Note: all of this information is currently provided from the CHARMM PSF-file that is read into the ASH openMMtheory object
Note: the atomnames column differs from conventional CHARMM usage. Instead OpenMM atomnames are used. Should not matter too much.

Note: Only use PDB-files for basic visualization, when you want to be able to visualize the system and use the reside information etc in VMD to be able to select residues. PDB-file is not a good format for other things. We for example do not want to use it as a file format in general because the format only supports a limited number of decimal points for coordinates.