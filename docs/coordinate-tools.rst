Coordinates and fragment tools
======================================


############################################
Modifying coordinates in ASH calculations
############################################

One often needs to manually modify coordinates in QM or QM/MM calculations. While this is straightforward when working
with small molecules and for example XYZ-files (open the coordinates in a simple molecular builder and modify) it is more
of an issue when working with a large system (e.g. a protein or a molecular crystal cluster) and you want to modify only a few atoms buried in the center of a 100 000 atom file.

- If one prefers to work with XYZ-files then it might be possible to use a program like VMD to modify certain coordinates there.

- If one works with ASH fragment files (.ygg) or XYZ files then one can modify the coordinates of a group of atoms via the use of scripts. These scripts are located in: /path/to/ashdir/ash/scripts and are called: **fragedit.py**  and **fragupdate.py**

**Grab and visualize part of the fragfile (fragedit.py)**

If one wants to visualize or possibly modify the coordinates of a group of atoms (e.g. the QM active site of a protein) then one can use the fragedit.py script like this:

.. code-block:: shell

    fragedit.py systemfile.ygg atomlistfile
    fragedit.py systemfile.xyz atomlistfile


.. note:: The script needs to be in your PATH and might need to be made executable (chmod +x fragedit.py)

The script will then read the coordinate-file (.ygg or .xyz) and extract the coordinates corresponding to atom indices present
in the atomlistfile (e.g. named qmatoms or activeatoms). The atomlistfile should contain a list of atom indices in a single line : e.g. 1 2 3 4 5

This will create a file called fragment.xyz (coordinates in Å), containing only the part of the system (as defined by the atom indices).
This file can be visualized in a molecular builder (e.g. Chemcraft) and the coordinates can also be modified.

.. note:: If you are using 1-based atom indexing to manage your qmatoms and actatoms files, there is an option: index1, that will assume that the atomlistfile contains 1-based indexing instead of the default 0-based indexing.


**Update the part of the system (fragupdate.py)**

If the coordinates were modified in the molecular builder they could be copied back to the fragment.xyz file (careful not to modify the header) and use the same
unit (Å). The fragfile (containing coordinates of the full system) can then be updated using the modified coordinates in fragment.xyz.

.. code-block:: shell

    python3 fragupdate.py fragfile.ygg atomlistfile

This should update the coordinates of fragfile.ygg.


######################################################
**Define an active region**
######################################################

In QM/MM calculations in particular it is usually convenient or even necessary to divide a system into region that may e.g. be QM or MM, frozen or active etc.
In ASH this is done by defining a list of atomindices of the whole system (counting starts from zero), typically stored in a file 
which can be read into a Python list in a script like this:

.. code-block:: python

    #Creates Python list actatoms from file active_atoms
    #File active_atoms should contain a list of atom indices (counting from zero) in a single line
    actatoms = read_intlist_from_file("active_atoms")

Contents of active_atoms file:

.. code-block:: text

    716 717 718 719 720 721 722 723 724 725 726


**define_activeregion function:**

While defining a list of atoms can often be done manually, when selecting a large region (e.g. an active region of ~1000 atoms) it is usually more convenient
to automate this task by using the **actregiondefine** function which can select atoms based on distance and residue information of the MM system.

.. code-block:: python

    from ash import *

    #Defining fragment containing coordinates (can be read from XYZ-file, ASH fragment, PDB-file)
    lastpdbfile="final_MDfrag_laststep_imaged.pdb"
    fragment=Fragment(pdbfile=lastpdbfile)

    #Creating new OpenMM object from OpenMM XML files (built-in CHARMM36 and a user-defined one)
    omm = OpenMMTheory(xmlfiles=["charmm36.xml", "charmm36/water.xml", "./specialresidue.xml"], pdbfile=lastpdbfile, periodic=True,
                platform='CPU',  autoconstraints=None, rigidwater=False)


    #Defining active region as within X Å from originatom 755 (Fe)
    actregiondefine(mmtheory=omm, fragment=fragment, radius=12, originatom=755)

The script will create the following output:

.. code-block:: text

                      ###########################
                      #                         #
                    #     ActregionDefine     #
                      #                         #
                      ###########################


    Radius: 12
    Origin atom: 755 (Fe)
    Will find all atoms within 12 Å from atom: 755 (Fe)
    Will select all whole residues within region and export list
    Wrote list to file: active_atoms
    Active region size: 908
    Active-region indices written to file: active_atoms
    The active_atoms list  can be read-into Python script like this:	 actatoms = read_intlist_from_file("active_atoms")
    Wrote Active region XYZfile: ActiveRegion.xyz  (inspect with visualization program)


This active_atoms file just contains a list of atom indices indicating which atoms should be active (all others are frozen).
The file can be manually modified if required. The ActiveRegion.xyz file can be visualized to make sure that the active-region looks reasonable.


**VDM alternative**

An alternative to the actregiondefine function is to do the visualization in VMD which allows you to both 
visually create a suitable active-region and get a list of atom indices (VMD also counts from zero) that can be copy-pasted into ASH.

In the VMD-GUI you can creating a new representation in "Graphical representations" 
and test out different atom-selections using VMD-code such as:

.. code-block:: tcl
    
    same residue as within 11 of index 33138

Once you are happy with the selection you can get a list of atom indices by copy pasting a variant of the following code
into the VMD shell:

.. code-block:: tcl
    
    #VMD code to define active-region based on whole residues positioned X Å from a certain atom
    #Here all whole residues within 11 Å of atom 33138 are selected
    set mol [molinfo top] 
    set sel [atomselect $mol {same residue as within 11 of index 33138}]
    set num_sel [$sel num] 
    puts "Number of atoms in selection: $num_sel"
    puts $sel
    $sel list

The VMD shell will then output a list of atom indices that you can copy-paste into a file and read into ASH.



######################################################
**Adding/removing atoms of an MM system**
######################################################

If you need to add or remove atoms to your MM or QM/MM system this is a bit more involved than just modifying the coordinates. The reason is that both the coordinate and forcefield file needs to be updated and also: if you delete e.g. atom 4556 then all atom indices > 4556 change.
This requires updating of forcefield files, coordinate files as well as atom lists (qmatoms and active atoms) that reference atom indices of the system.

There are two options:

1. Go back to the original MM-system preparation and prepare a new MM model with the added/deleted atom(s). This is a safe option but inconvenient.

2. Modify the coordinate-file (XYZ-file, YGG-file, PDB-file), the forcefield file (e.g. PSF-file, topology file) and update atom-indices-files (e.g. active_atoms and qmatoms files).

    a. CHARMM files:
        The PSF-file has to be regenerated and the topology and parameter-files may also need modifications/additions.
        PSFgen is the best option for creating a new PSF-file.

        **Delete atoms (CHARMM)**

        Both the coordinate-deletion and PSF-file update can be performed with an ASH script like this:

        .. code-block:: python

            from ash import *

            #Path to dir containing PSFgen executable
            psfgendir="/home/bjornsson/QM-MM-Chemshell-scripts"

            #CHARMM Forcefield files
            topfile="top_all36_prot.rtf"
            psffile="newxplor.psf"

            #Reading coordinates into a fragment
            fragfile=Fragment(fragfile="Fragment-currentgeo.ygg")

            # Define qmatoms and actatoms lists
            qmatoms = read_intlist_from_file("qmatoms")
            actatoms = read_intlist_from_file("actatoms")

            #What atoms to delete
            deletionlist=[18840]

            #Delete atoms from system
            remove_atoms_from_system_CHARMM(atomindices=deletionlist, fragment=fragfile,psffile=psffile,topfile=topfile, 
                psfgendir=psfgendir, qmatoms=qmatoms, actatoms=actatoms)

        The script will delete the selected atoms (here 18840; note: ASH counts from zero) and create new fragmentfiles: 
        newfragment.xyz and newfragment.ygg
        and create the new PSF file named: newsystem_XPLOR.psf  . Also created is a PDB-file: new-system.pdb

        Remember that when you delete atoms from a system atom indices will have changed. 
        This means that you either have to update the qmatoms and actatoms list manually or do as in example above where the qmatoms and actatoms lists are provided to the remove_atoms_from_system_CHARMM function. These lists will then be updated.

    .. note:: If you are using 1-based atom indexing to manage your qmatoms and actatoms files, there is an option: offset_atom_indices=1, to remove_atoms_from_system_CHARMM  that will preserve the 1-based indexing.


    **Add atoms to system (CHARMM)**
            
    Both the coordinates and the PSF-file needs to be updated. 
    This can be performed with an ASH script like this:

    .. code-block:: python

        from ash import *

        #Path to dir containing PSFgen executable
        psfgendir="/home/bjornsson/QM-MM-Chemshell-scripts"

        #CHARMM Forcefield files
        topfile="top_all36_prot.rtf"
        psffile="newxplor.psf"

        #Reading coordinates into a fragment
        fragfile=Fragment(fragfile="Fragment-currentgeo.ygg")

        # Define qmatoms and actatoms lists
        qmatoms = read_intlist_from_file("qmatoms")
        actatoms = read_intlist_from_file("actatoms")

        #Defining the added coordinates as a string
        addition_string="""
        C        1.558526678      0.000000000     -0.800136464
        O        2.110366050     -0.126832008      0.222773815
        O        1.006687306      0.126832008     -1.823046743
        """
        #Name of resgroup to be added (this needs to be present in topfile!)
        resgroup='CO2'
        #Adding atoms
        add_atoms_to_system_CHARMM(fragment=fragfile, added_atoms_coordstring=addition_string, resgroup=resgroup, 
            psffile=psffile, topfile=topfile, psfgendir=psfgendir, qmatoms=qmatoms, actatoms=actatoms)

    The script will add the selected atom coordinates to the fragment (at the end) and create new fragmentfiles: 
    newfragment.xyz and newfragment.ygg
    and add the chosen resgroup to a PSF file named: newsystem_XPLOR.psf  . 
    Also created is a PDB-file: new-system.pdb

    Remember to add the new atom indices to QM-region and Active-Region definitions or provide the lists to the add_atoms_to_system_CHARMM function as above.

.. note:: If you are using 1-based atom indexing to manage your qmatoms and actatoms files, there is an option: offset_atom_indices=1, to add_atoms_to_system_CHARMM  that will preserve the 1-based indexing.


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