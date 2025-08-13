Coordinates and fragment tools
======================================

############################################
Read in XYZ-files from a directory
############################################

It can be convenient to read in multiple molecular structures at the same time to an ASH script, by pointing to a directory 
with multiple XYZ files present.
The *read_xyzfiles* function can be use for this purpose. One simply has to provide the path to the directory and ASH will then 
attempt to read all XYZ-files present in the directory (each file has to have a .xyz suffix), 
create an ASH fragment for each XYZ-file and then return a list of the ASH fragments

.. code-block:: python

    def read_xyzfiles(xyzdir,readchargemult=False, label_from_filename=True):

Additional options include the *readchargemult* option which if set to True then ASH will try to read charge/mult information (space-separated) from the 2nd-line of the XYZ-file (needs to be present).
Note that otherwise charge/mult information has to be provided in some other way when Fragment is used for a calculation (either by setting the charge/mult attribute to each Fragment or provide the info to the job-function).
The *label_from_filename* option is automatically set to True which ensures that each ASH fragment created will have a unique label (corresponding to original XYZ filename).
The label can be changed by changing the Fragment.label attribute if desired.

Example:

.. code-block:: python

    # Create ASH fragments from all XYZ-files in a directory.
    fragments = read_xyzfiles(xyzdir,readchargemult=False, label_from_filename=True)

    #Example: Loop over each fragment in the list of ASH fragments
    for frag in fragments:
        #Print info on each fragment
        print(frag)
        #Set charge/mult for each fragment 
        frag.charge=0; frag.mult=1
        #Do something with each fragment, e.g. run a single-point energy job
        #Note: theory object needs to first be defined.
        #charge and mult keywords can also be provided to Singlepoint if desired
        Singlepoint(fragment=frag, theory=theory)

.. note:: When reading the XYZ-files in the directory, ASH automatically sorts the filenames by `natural sort order <https://en.wikipedia.org/wiki/Natural_sort_order>`_. This ensures that XYZ-files labelled e.g. mol1.xyz, mol2.xyz, mol11.xyz are read in that order.


############################################
Read in a multi-XYZ or trajectory file
############################################

XMol-style XYZ-files can contain multiple geometries, which might come from a MD trajectory, optimization trajectory, NEB path or perhaps a collection
of lowest-energy conformers from e.g. a CREST calculation (see :doc:`crest-interface`). 
Such files might be created by ASH or by another program.

ASH can conveniently read-in such files using the *get_molecules_from_trajectory* function:

.. code-block:: python

    def get_molecules_from_trajectory(file, writexyz=False, skipindex=1, conncalc=False):

One simply has to provide the name of the multi-molecule XYZ-file but there are also options to write out each 
individual geometry as its own XYZ-file (i.e. split it) or skip geometries (e.g. every 2nd geometry by setting skipindex to 2).

Examples shown below:

.. code-block:: python

    # Read-in the MD trajectory from an ASH MD run (written in XYZ-format)
    Opt_traj_fragments = get_molecules_from_trajectory("OpenMMMD_traj.xyz")

    # Read-in the optimization trajectory from an ASH Optimization
    Opt_traj_fragments = get_molecules_from_trajectory("geometric_OPTtraj_optim.xyz")

    # Read-in the lowest-energy conformers predicted by CREST
    CREST_fragments = get_molecules_from_trajectory("crest_conformers.xyz")

    # Read-in the minimum energy path from an NEB calculation 
    NEB_fragments = get_molecules_from_trajectory("knarr_MEP.xyz")

    # Read in the Wigner ensemble 
    Wigner_fragments = get_molecules_from_trajectory("Wigner_traj.xyz")


While molecular dynamics trajectories are sometime written in XYZ format, it is inconvenient for large systems or long trajectories.
Instead it is more common to utilize a compressed format such as DCD for such trajectories.
To read in DCD files it is best to use the *mdtraj* library, either on its own or the ASH wrapper interface.
See :doc:`module_dynamics` for information on the *mdtraj* library info on how to read in DCD-files e.g. using the *MDtraj_slice* function.


############################################
Calculate RMSD
############################################

It can be useful to compare the similarity of 2 molecular geometries. Calculating the root-mean-square deviation of atomic positions is one suitable way of doing this.
This requires first superimposition or alignment of the 2 structures as they may be in different parts of Cartesian space. This is accomplished using the Kabsch algorithm.
The 2 structure generally must contain the same number of atoms and the atoms must be in the same order. However, a subset RMSD can be calculated by providing a subset of atom indices for both structures.
If the subsets for both fragments match then the RMSD will be calculated.

.. code-block:: python

    def calculate_RMSD(fragmentA, fragmentB, subset=None, heavyatomsonly=False, printlevel=2):


Examples on how to use the function:

.. code-block:: python

    #Calculate the RMSD (in Å) between 2 ASH fragments. All atoms used to calculate RMSD.
    rmsd_val = calculate_RMSD(reference_frag, frag)
    # Calculate the RMSD but using only heavy atoms (no hydrogens) included
    rmsd_val = calculate_RMSD(reference_frag, frag, heavyatomsonly=True)
    # Calculate the RMSD but using a subset of atom indices. Note: The 2 fragments must have exactly the same atom-order
    rmsd_val = calculate_RMSD(reference_frag, frag, subset=[5,6,7])
    # Calculate the RMSD but using a  list-of-lists definition of subset of atom indices (first list for first fragment etc.)
    rmsd_val = calculate_RMSD(reference_frag, frag, subset=[[5,6,7],[1,2,3]])

############################################
flexible_align
############################################

Sometimes it is useful to align a molecular geometry so that it is as similar as possible to another geometry. This is often performed for the purpose of calculating the RMSD (see above) but often
the purpose is the aligned geometry itself, the 2 structures might not fully match and one might even want the structure reoriented or even reorderered as much as possible for the purpose of maximum alignment.

The **flexible_align** function allows one to align a structure (fragmentA below) onto another fragment (fragmentB). One can choose to only allow rotation of the structure (rotate_only=True), 
only allow translation (translate_only=True) or allow both (default). One can also allow reordering which would use the Hungarian algorithm to reorder the atoms of fragmentA to match fragmentB as much as possible.
The subset option allows one to use in the comparison only a subset of atom indices (i.e. the atoms that fragmentA and fragmentB have in common). The resulting aligned fragmentA, however, will contain all atoms aligned.

.. code-block:: python

    def flexible_align(fragmentA, fragmentB, rotate_only=False, translate_only=False, reordering=False, reorder_method='brute', subset=None):

    # Versions that takes PDB-files or XYZ-files as input and outputs an aligned XYZ/PDB file
    def flexible_align_xyz(xyzfileA, xyzfileB, rotate_only=False, translate_only=False, reordering=False, reorder_method='brute', subset=None):

    def flexible_align_pdb(pdbfileA, pdbfileB, rotate_only=False, translate_only=False, reordering=False, reorder_method='brute', subset=None):


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


**actregiondefine function:**

While defining a list of atoms can often be done manually, when selecting a large region (e.g. an active region of ~1000 atoms) it is usually more convenient
to automate this task by using the **actregiondefine** function which can select atoms based on distance and residue information of the MM system. 
actregiondefine can either use residue information present in an OpenMMTheory object (created from CHARMM/Amber/XML forcefield-files)
or from a PDB-file.


*#Using the residue information from the PDB-file*

.. code-block:: python

    from ash import *

    #Defining fragment containing coordinates (can be read from XYZ-file, ASH fragment, PDB-file)
    pdbfile="final_MDfrag_laststep_imaged.pdb"
    fragment=Fragment(pdbfile=pdbfile)

    #Defining active region as within X Å from originatom 755 (Fe)
    actregiondefine(pdbfile=pdbfile, fragment=fragment, radius=12, originatom=755)


*#Using the residue information the OpenMMTheory object (there are cases where this fails)*

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
The file can be manually modified if required. The ActiveRegion.xyz file should be visualized to make sure that the active-region looks reasonable.

.. warning:: There are cases where an MM system might be set up in such a way that a residue definition can apply to multiple molecules/fragments in space.
    The actregiondefine function may not handle all such cases.

**VMD alternative**

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

If you need to add or remove atoms to your MM or QM/MM system this is a bit more involved than just modifying the coordinates. The reason is that both the coordinate and forcefield file needs to be updated and also: if you delete e.g. atom 4556 then all atom indices larger than 4556 change.
This requires updating of forcefield files, coordinate files as well as atom lists (qmatoms and active atoms) that reference atom indices of the system.

There are two options:

1. Go back to the original MM-system preparation and prepare a new MM model with the added/deleted atom(s). This is a safe option but inconvenient.

2. Modify the coordinate-file (XYZ-file, YGG-file, PDB-file), the forcefield file (e.g. PSF-file, topology file) and update atom-indices-files (e.g. active_atoms and qmatoms files).
The forceefield files are the tricky ones. For OpenMM-XML forcefield files it is relatively straightforward to modify the user-created XML-files. Amber prmtop files are much more difficult due to their format and it is instead recommended to go back to e.g. the Ambertools setup instead.
For CHARMM-files see below:

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

WARNING: PDB files are convenient for visualization purposes and for initial reading the initial set of coordinates but are
generally not a file format to be used (one problem is the limited number of significant digits used
for coordinates in the file).

----------------------
Reading in PDB file
----------------------

It is possible to read in coordinates from a PDB file to create an ASH fragment file.
This functionality is very basic, it will only read in the coordinates, not atom-types
or residue information. This option is thus only be used to provide convenient starting coordinates.

.. code-block:: python

    pdbfrag = Fragment(pdbfile="mol.pdb")

Note that OpenMMTheory objects (see :doc:`OpenMM-interface`) also have a pdbfile option, however, this 
option is primarily used for reading in topology information (residue information, atom types etc) and not for coordinates.

----------------------
Writing out PDB file
----------------------

ASH contains a few different options for writing out PDB-files which can be useful for visualization purposes etc.

**Fragment.write_pdbfile_openmm**: 

This writes out a PDB-file from an ASH fragment, using either topology and residue information that was read from original PDB-file.
If latter is not present (e.g. if an XYZ-file was read-in), a basic topology is automatically defined.
Routines from OpenMM library are used to read PDB-topology and write out the PDB-file.

.. code-block:: python
    
    #Initial fragment from a PDB-file
    frag = Fragment(pdbfile="initial.pdb")
    #Define theory
    theory = xTBTheory()
    #Geometry optimization, results in updated optimized coordinates in frag object
    Optimizer(theory=theory, fragment=frag)
    #Writing out PDB-file with optimized coordinates. Topology and residue information is reused (from initial.pdb)
    #Note: if a PDB-file was not used to create the fragment, basic topology and residue information will be guessed
    frag.write_pdbfile_openmm(filename="optimized.pdb")


**OpenMMTheory.write_pdbfile**: This is a method inside the OpenMMTheory object that writes out a PDB-file based on coordinates, residue and atom information present in the OpenMMTheory object.
Requires an OpenMMTheory object.

.. code-block:: python

    #omm is a predefined OpenMMTheory object
    omm.write_pdbfile(outputname="ASHfragment")

.. warning:: Make sure the OpenMMTheory object contains the desired coordinates.

**write_pdbfile_openMM**: 

Standalone function writing a PDB-file based on input OpenMM topology, positions and optionally connectivity information.
Uses OpenMM-library PDB-writing routines (usually pretty robust).

.. code-block:: python

    def write_pdbfile_openMM(topology, positions, filename, connectivity_dict=None):


**write_pdbfile**: 

This is a standalone flexible function that writes out a PDB-file based on an ASH fragment and other optional data.

.. code-block:: python

    def write_pdbfile(fragment,outputname="ASHfragment", openmmobject=None, atomnames=None,
                    resnames=None,residlabels=None,segmentlabels=None):


An ASH fragment file needs to always be provided.

.. code-block:: python

    #Example 1 (no residue information provided)
    #All residues will be labelled 'DUM' and segments 'SEG', element information should be correct.
    write_pdbfile(frag)
    #Example 2 (residue information provided manually, via information from OpenMMTheory object)
    openmmobject = OpenMMTheory(psffile=psffile, CHARMMfiles=True, charmmtopfile=topfile,charmmprmfile=parfile,
                    printlevel=1, platform='CPU' )
    write_pdbfile(frag, outputname="manual", atomnames=openmmobject.atomnames, resnames=openmmobject.resnames,
        residlabels=openmmobject.resids,segmentlabels=openmmobject.segmentnames)
    #Example 3: usually best way. Information taken from OpenMMTheoryobject
    #Note: the atomnames column differs from conventional CHARMM usage. Instead OpenMM atomnames are used. Should not matter too much.
    write_pdbfile(frag, outputname="simple",openmmobject=openmmobject)


.. warning:: While this function is flexible it does not always write out PDB-file that is compatible with all visualization programs. 




