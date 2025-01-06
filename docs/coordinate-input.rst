Coordinates and fragments
==========================

The ASH Fragment is the primary object-type that one works with in ASH.
Generally one creates a Fragment object, in one of several ways, and then passes it to a Job-function
along with a Theory object.

This page primarily describes how the ASH Fragment class works.
Also note that the page :doc:`coordinate-tools` contains descriptions of various helper-tools related to coordinates e.g.
on how to create a list of ASH fragments from an XYZ-directory (*read_xyzfiles*) or from a multi-structure XYZ-file or trajectory (*get_molecules_from_trajectory*).


######################################################
The Fragment class
######################################################

.. code-block:: python

    class Fragment:
        def __init__(self, coordsstring=None, fragfile=None, databasefile=None, xyzfile=None, pdbfile=None, grofile=None,
                    amber_inpcrdfile=None, amber_prmtopfile=None,
                    chemshellfile=None, coords=None, elems=None, connectivity=None, atom=None, diatomic=None, bondlength=None,
                    atomcharges=None, atomtypes=None, conncalc=False, scale=None, tol=None, printlevel=2, charge=None,
                    mult=None, label=None, readchargemult=False, use_atomnames_as_elements=False):

        def update_attributes(self):

        # Add coordinates from geometry string. Will replace.
        def add_coords_from_string(self, coordsstring, scale=None, tol=None, conncalc=False):

        # Replace coordinates by providing elems and coords lists. Optional: recalculate connectivity
        def replace_coords(self, elems, coords, conn=False, scale=None, tol=None):

        def delete_coords(self):

        def get_atomindices_except(self, excludelist):

        def get_nonH_atomindices(self):

        def get_atomindices_for_element(self, element):

        def get_atomindices_except_element(self, element):

        def get_XH_indices(self, conncode='julia'):

        def simple_get_water_constraints(self, starting_index=None):

        def delete_atom(self, atomindex):

        def add_coords(self, elems, coords, conn=True, scale=None, tol=None):

        def print_coords(self):

        def read_amberfile(self, inpcrdfile=None, prmtopfile=None, conncalc=False):

        def read_grofile(self, filename, conncalc=False, scale=None, tol=None):

        def read_chemshellfile(self, filename, conncalc=False, scale=None, tol=None):

        def read_pdbfile(self, filename, conncalc=True, scale=None, tol=None, use_atomnames_as_elements=False):

        def read_xyzfile(self, filename, scale=None, tol=None, readchargemult=False, conncalc=True):

        def set_energy(self, energy):

        def get_coordinate_center(self):

        # Get coordinates for specific atoms (from list of atom indices)
        def get_coords_for_atoms(self, atoms):

        # Calculate connectivity (list of lists) of coords
        def calc_connectivity(self, conndepth=99, scale=None, tol=None, codeversion=None):

        def update_atomcharges(self, charges):

        def update_atomtypes(self, types):

        # Adding fragment-type info (used by molcrys, identifies whether atom is mainfrag, counterfrag1 etc.)
        # This one is fast
        def add_fragment_type_info(self, fragmentobjects):

        def write_xyzfile(self, xyzfilename="Fragment-xyzfile.xyz", writemode='w', write_chargemult=True, write_energy=True):

        def write_XYZ_for_atoms(self,xyzfilename="Fragment-subset.xyz", atoms=None):

        # Print system-fragment information to file. Default name of file: "fragment.ygg
        def print_system(self, filename='fragment.ygg'):

        # Reading fragment from file. File created from Fragment.print_system
        def read_fragment_from_file(self, fragfile):


######################################################
Creating/modifying fragment objects
######################################################

Fragments in ASH are Python objects containing basic information about a molecule. You can create as many fragment objects
as you want. A fragment object will contain Cartesian coordinates about a molecule, elemental information and masses.
Sometimes additional information such as connectivity, constraints, charges and multiplicity information is present as well.
Fragments can be created in many different ways (from XYZ-file,PDB-file, coordinate string, lists or Numpy arrays etc.) but will behave the same after creation.

Fragments are created from the ASH *Fragment* object class above.

######################################################
Direct creation of an ASH fragment from coordinates
######################################################

*From string*

First define multi-line string (called fragcoords here) with element and coordinates (Å) separated by space:

.. code-block:: python

    fragcoords="""
    H 0.0 0.0 0.0
    F 0.0 0.0 1.0
    """

Then define object (here called **HF_frag**) of class *Fragment* by passing the coordinates to *coordsstring*, using coordinates from the string "fragcoords".
The *Fragment* class is an ASH class.

.. code-block:: python

    HF_frag=Fragment(coordsstring=fragcoords)



*From list*

Another way is if you have lists of coordinates and element information already available.

.. code-block:: python

    elems=['H', 'Cl']
    coords=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.9]]
    HCl_frag=Fragment(elems=elems,coords=coords)


*From external XYZ file*

Perhaps most convenient is to define the fragment directly from reading an XYZ-file (that exists in the same directory as the script):

.. code-block:: python

    HI_frag = Fragment(xyzfile="hi.xyz")

Note that the XYZ-file should be in XMol format meaning that a 2-line header must be present, containing the number of atoms (1st line) and a comment line (2nd line).

*From previous ASH fragment file*

ASH fragment files use the .ygg extension. They are typically not created manually but are often created automatically by ASH code and
can be created upon request. To read an old file from disk (here "previous.ygg") you do:

.. code-block:: python

    mol_frag = Fragment(fragfile="previous.ygg")


*From external PDB file*

It is also possible to read coordinates from a PDB file. This functionality is very rudimentary, only supporting read-in of
elements and coordinates, not atom-types or residue information.

.. code-block:: python

    pdbfrag = Fragment(pdbfile="mol.pdb")

*From XYZ-file in ASH database*

ASH contains an internal database of some small molecules.
These are XYZ-files that are present in the ASH repository ( see `ASH-code/databases/fragments directory <https://github.com/RagnarB83/ash/tree/master/databases/fragments>`_ )
Examples of available files: h2o.xyz, nh3.xyz, n2.xyz, butane.xyz, glycine.xyz  etc.

The coordinates in these files have been pre-optimized by some level of theory and are reasonable but should obviously be re-optimized for any serious calculations.

.. code-block:: python

    pdbfrag = Fragment(databasefile="h2o.xyz")


*From a SMILES-string*

ASH now also supports creating a Fragment using a `SMILES-string <https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system>`_ as input.
This feature requires OpenBabel to be installed in the same Python environment as ASH.
OpenBabel will parse the string, add H-atoms and guess the 3D-structure in Cartesian coordinates which is passed onto ASH.
If everything is successful the ASH Fragment can be used just like any other Fragment for further calculations.

.. code-block:: python

    #ASH fragment from a SMILES string for aspirine
    #From: https://en.wikipedia.org/wiki/Aspirin
    frag = Fragment(smiles="CC(=O)OC1=CC=CC=C1C(=O)O")
    #Write out XYZ file of fragment
    frag.write_xyzfile()

.. note:: OpenBabel can be installed in your ASH conda environment like this: conda install --yes -c conda-forge openbabel

######################################################
Adding coordinates to object
######################################################


*Add coordinates from string*



.. code-block:: python

    HCl_cluster = Fragment(xyzfile="hcl.xyz")

    fragcoords="""
    H 0.0 0.0 0.0
    Cl 0.0 0.0 1.0
    """
    HCl_frag.add_coords_from_string(fragcoords)


.. note:: This will append coordinates to fragment. If fragment already contains some coordinates the specified coordinates will be appended.

*Add coordinates from lists*

.. code-block:: python

    HCl_frag.add_coords(elems,coords)

where elems and coords are lists:

.. code-block:: python

    elems=['H', 'Cl']
    coords=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.9]]


.. note:: This will append coordinates to fragment. If fragment already contains some coordinates the added coordinates will follow.

*Add coordinates from XYZ file*

.. code-block:: python

    HF_frag.read_xyzfile("hcl.xyz")


.. note:: This will append coordinates to fragment. If fragment already contains some coordinates the added coordinates will follow.

######################################################
Replace coordinates of object
######################################################
If you want to replace coords and elems of a fragment object with new information this can be done conveniently through lists.

.. code-block:: python

    elems=['H', 'Cl']
    coords=[[0.0, 0.0, 0.0], [0.0, 0.0, 1.1]]
    HF_frag.replace_coords(elems,coords)

**TODO:** Add option here of replacing coords from XYZ file and string as well.


######################################################
Calculate connectivity of fragment object
######################################################

Connectivity between atoms can be an important attribute of a fragment object, especially for molecules.
The connectivity distinguishes atoms that are in close-contact (i.e. forming some kind of stable covalent bond) 
and atoms further apart and obviously not bonded. Correct connectivity is crucial for some ASH functionality (the Molcrys functionality in particular).
Connectivity is calculated based on a distance and covalent radii-based criterion.
Atoms A and B will be defined to be connected according to:

.. math::

    r(AtomA,AtomB) < scale*( covrad(AtomA) + covrad(AtomB) ) + tol

Thus, if the distance between atoms A and B is less than the sum of the elemental covalent radii
(which can be scaled by a parameter scale or shifted by a parameter tol) then the atoms are connected.
Using default parameters of the element radii (Alvarez 2008), the default scaling of 1.0 and a tolerance of 0.1
(global scale and tol parameters are defined in settings_ash file) works in many cases.

To calculate the connectivity table for a molecule:

.. code-block:: python

    mol_frag.calc_connectivity()

This creates a connectivity table which is a Python list of lists:
An example of a connectivity table would be: [[0,1,2],[3,4,5],[6,7,8,9,10]]
Atoms 0,1,2 are here bonded to each other as a sub-fragment (migh e.g. be an H2O molecule) and so are atoms 3,4,5 and also 6,7,8,9,10.
The connectivity table is available as:

.. code-block:: python

    mol_frag.connectivity


The connectivity table is calculated or recalculated automatically when needed. 
For large systems the connectivity is expensive to calculate and is thus not calculated by default (but only when needed).
For large systems, ASH will try to call a Julia routine for the calculation.

######################################################
Charge and Multiplicity
######################################################

Charge and spin multiplicity should usually be associated with the fragment.
One can also specify the charge and mult to the Job-function (e.g. **Singlepoint**).
When working with multiple fragment objects, however, it is highly convenient to associate a total charge and spin multiplicity with each fragment object.
Usually done when fragment is created like this:

.. code-block:: python

    NO_frag = Fragment(xyzfile="no.xyz", charge=0, mult=2)
    HF_frag=Fragment(coordsstring=fragcoords, charge=0, mult=1)

This can also be done afterwards:

.. code-block:: python

    NO_frag.charge = 0
    NO_frag.mult = 2

Yet another option is to read the charge and multiplicity information from the name/title line of the XYZ file.

.. code-block:: python

    NO_frag = Fragment(xyzfile="no.xyz", readchargemult=True)

This will only work if the 2nd-line of the XYZ file contains the charge and multiplicity, separated by a space as seen below:

.. code-block:: text

    2
    0 2
    N 0.0 0.0 0.0
    O 0.0 0.0 1.0

######################################################
Label
######################################################

If working with multiple fragment objects it can be useful to distinguish between them via a label-string.
The label can be added when fragment is first created:

.. code-block:: python

    benzene_frag = Fragment(xyzfile="c6h6.xyz", label='benzene')
    water_frag = Fragment(xyzfile="h2o.xyz", label='water')

or afterwards (by default, the label attribute is set to None).

.. code-block:: python

    benzene_frag.label='Benzene'


######################################################
Inspect defined fragment objects
######################################################

To inspect a defined fragment one can print out a Python dictionary of all defined attributes of the object.

.. code-block:: python

    print("HF_frag dict", HF_frag.info())

This will print out all defined attributes of the object including list of elements, coordinates, masses, connectivity etc.

One can also access individual attributes like accessing the pure coordinates only:

.. code-block:: python

    print("HF_frag.coords : ", HF_frag.coords)

For printing coordinates is may also be more convenient to use the print_coords function though (to print elems and coords):

.. code-block:: python

    HF_frag.print_coords()


Get coords and elems of specific atom indices:

.. code-block:: python

    specific_coords,specific_elems=HF_frag.get_coords_for_atoms([0,1,2])

Print connectivity:

.. code-block:: python

    conn = aspirine.connectivity
    print("conn:", conn)
    print("Number of subfragments in aspirine", len(conn))

Print number of atoms and number of connected atoms:

.. code-block:: python

    print("Number of atoms in aspirine", aspirine.numatoms)
    print("Number atoms in connectivity in aspirine", aspirine.connected_atoms_number)

Print various molecule attributes:

.. code-block:: python

    print("List of atom indices", frag.atomlist)
    print("Total mass of fragment", frag.mass)
    print("List of atom masses of fragment", frag.list_of_masses)
    print("Pretty elemental formula of fragment", frag.prettyformula)
    print("Elemental formula of fragment", frag.formula)
    print("Pretty elemental formula of fragment", frag.prettyformula)

The ASH fragment file can be printed conveniently to disk:

.. code-block:: python

    HF_frag.print_system(filename='fragment.ygg')

An XYZ file of coordinates can be printed out:

.. code-block:: python

    HF_frag.write_xyzfile(xyzfilename="Fragment-xyzfile.xyz")


Print charge and mult attributes (if not defined, then None will be outputted).

.. code-block:: python

    print(HF_frag.charge)
    print(HF_frag.mult)

##############################################################################
Calculate distances,angles and dihedral angles between atoms in a fragment
##############################################################################

Sometimes it is useful to get the distance, angle or dihedral angles defined by some atoms in a fragment.
This can be done using the following functions:

.. code-block:: python

    #Distance between atoms
    def distance_between_atoms(fragment=None, atoms=None):
    #Angle between atoms
    def angle_between_atoms(fragment=None, atoms=None):
    #Dihedral angle between atoms
    def dihedral_between_atoms(fragment=None, atoms=None):


Examples:

.. code-block:: python

    from ash import *

    #Defining an ethanol fragment
    ethanol = Fragment(databasefile="ethanol.xyz")
    ethanol.print_coords() #Printing coordinates

    #Distance between atoms 0 (C) and 7 (O)
    distance = distance_between_atoms(fragment=ethanol, atoms=[0,7])
    print(f"Distance between atoms 0 and 7 is {distance} Angstrom")
    
    #Angle between atoms 0, 7 and 8 (< C-O-H)
    angle = angle_between_atoms(fragment=ethanol, atoms=[0,7,8])
    print(f"Angle between atoms 0,7,8 is {angle} °")
    
    #Dihedral angle between atoms 3,0,7,8
    dihedral = dihedral_between_atoms(fragment=ethanol, atoms=[3,0,7,8])
    print(f"Dihedral between atoms 3,0,7,8 is {dihedral} °")

