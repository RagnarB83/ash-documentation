==========================
Coordinates and fragments
==========================


Creating/modifying fragment objects
***********************************

Fragments in Ash are Python objects containing basic information about a molecule. You can create as many fragment objects
as you want. A typical fragment will contain at least Cartesian coordinates about a molecule and the elemental information.
Fragments can be created in multiple ways but will behave the same after creation.

Fragments are Python objects created from the Ash *Fragment* object class.
See XXFragment-class-page-linkXX for an overview of all Fragment class attributes and functions.

Direct creation of fragment from coordinates
==============================================

*From string*

First define multi-line string (called fragcoords here) with element and coordinates (Å) separated by space:

.. code-block:: python

    fragcoords="""
    H 0.0 0.0 0.0
    F 0.0 0.0 1.0
    """

Then define object (here called **HF_frag**) of class *Fragment* by passing the coordinates to *coordsstring*, using coordinates from the string "fragcoords".
The *Fragment* class is an Ash class.

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

*From previous ASH fragment file*

ASH fragment files use the .ygg extension. They are typically not created manually but are often created automatically by ASH code and
can be created upon request. To read an old file from disk (here "previous.ygg") you do:

.. code-block:: python

    mol_frag = Fragment(fragfile="previous.ygg")


*From external PDB file*

Also possible to read coordinates from a PDB file. This functionality is very rudimentary, only supporting read-in of
elements and coordinates, not atom-types or residue information.

.. code-block:: python

    pdbfrag = Fragment(pdbfile="mol.pdb")


Adding coordinates to empty object
=====================================

An alternative to the direct way is to first create an empty fragment and then add the coordinate and element information later.
This can sometimes be useful and demonstrates here the built-in fragment object functions available (coords_from_string, add_coords, read_xyzfile)
First create empty fragment:

.. code-block:: python

    HCl_frag=Fragment()


*Add coordinates from string*


.. code-block:: python

    fragcoords="""
    H 0.0 0.0 0.0
    F 0.0 0.0 1.0
    """
    HCl_frag.add_coords_from_string(fragcoords)


**Note:** This will append coordinates to fragment. If fragment already contains some coordinates the specified coordinates
will be appended.

*Add coordinates from lists*

.. code-block:: python

    HCl_frag.add_coords(elems,coords)

where elems and coords are lists:

.. code-block:: python

    elems=['H', 'Cl']
    coords=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.9]]


**Note:** This will append coordinates to fragment. If fragment already contains some coordinates the added coordinates
will follow.

*Add coordinates from XYZ file*

.. code-block:: python

    HF_frag.read_xyzfile("hcl.xyz")


**Note:** This will append coordinates to fragment. If fragment already contains some coordinates the added coordinates
will follow.


Replace coordinates of object
==============================
If you want to replace coords and elems of a fragment object with new information this can be done conveniently through lists.

.. code-block:: python

    elems=['H', 'Cl']
    coords=[[0.0, 0.0, 0.0], [0.0, 0.0, 1.1]]
    HF_frag.replace_coords(elems,coords)

**TODO:** Add option here of replacing coords from XYZ file and string as well.


Delete coordinates of object
==============================
If you want to delete coordinates from object (both coords list and elems lists) then this is easily done.

.. code-block:: python

    HF_frag.delete_coords()


Calculate connectivity of fragment object
===========================================

Connectivity is an important aspect of the fragment as it distinguishes atoms that are in close-contact (i.e. forming some kind of stable covalent bond) and atoms further apart and obviously not bonded. Correct connectivity is crucial for some Ash functionality.
Connectivity is calculated based on a distance and covalen radii-based criterion.
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


The connectivity table is calculated or recalculated automatically when coordinates are added or when modified to the fragment.
It is typically unnecessary to request a calculation or recalculation.


Charge and Multiplicity
=================================

Charge and spin multiplicity can be associated with the fragment (either at creation or afterwards) but does not have to.
The QM theory level needs the charge and multiplicity information and it usually must be provided when the QMtheory object is created.
When working with multiple fragment objects, however, it is convenient to associate a total charge and spin multiplicity with each fragment object.
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

.. code-block:: shell

    2
    0 2
    N 0.0 0.0 0.0
    O 0.0 0.0 1.0

Label
=================================

If working with multiple fragment objects it can be useful to distinguish between them via a label-string.
The label can be added when fragment is first created:

.. code-block:: python

    benzene_frag = Fragment(xyzfile="c6h6.xyz", label='benzene')
    water_frag = Fragment(xyzfile="h2o.xyz", label='water')

or afterwards (by default, the label attribute is set to None).

.. code-block:: python

    benzene_frag.label='Benzene'



Inspect defined fragment objects
=================================

To inspect a defined fragment one can print out a Python dictionary of all defined attributes of the object.

.. code-block:: python

    print("HF_frag dict", HF_frag.__dict__)

One can also access individual attributes like accessing the pure coordinates only:

.. code-block:: python

    print("HF_frag.coords")

More conveniently would be to use the print_coords function though (to print elems and coords):

.. code-block:: python

    print("HF_frag.print_coords")


Get coords and elems of specific atom indices:

.. code-block:: python

    specific_coords,specific_elems=HF_frag.get_coords_for_atoms([0,1,2])

Print connectivity:

.. code-block:: python

    conn = FeFeH2ase.connectivity
    print("conn:", conn)
    print("Number of subfragments in FeFeH2ase", len(conn))

Print number of atoms and number of connected atoms:

.. code-block:: python

    print("Number of atoms in FeFeH2ase", FeFeH2ase.numatoms)
    print("Number atoms in connectivity in FeFeH2ase", FeFeH2ase.connected_atoms_number)

Print various molecule attributed:

.. code-block:: python

    print("Elemental formula of fragment", frag.formula)
    print("Pretty elemental formula of fragment", frag.prettyformula)
    print("Number atoms in connectivity in FeFeH2ase", FeFeH2ase.connected_atoms_number)


The Ash fragment file can be printed conveniently to disk:

.. code-block:: python

    HF_frag.print_system(filename='fragment.ygg')

An XYZ file of coordinates can be printed out:

.. code-block:: python

    HF_frag.write_xyzfile(xyzfilename="Fragment-xyzfile.xyz")


Print charge and mult attributes (if not defined, then None will be outputted).

.. code-block:: python

    print(HF_frag.charge)
    print(HF_frag.mult)