OpenBabel interface
======================================

`OpenBabel <https://openbabel.org>`_ is described as a chemical toolbox for speaking the many languages of chemical data.
It's an open-source project that has long been useful for converting between file formats and many other things.


ASH features a simple interface to using OpenBabel as a theory-level, currently only for the purpose of using the
built-in FF options in OpenBabel.
The interface hence allows one to utilize the simple MMFF94, UFF, GAFF and Ghemical forcefields.
These can be useful in ASH as e.g pre-optimizer or useful on their own if accurate enough.

Note that the interface can not be used as MM-theory in a hybrid QM/MM context.
However, it should work within a hybrid ONIOM scheme.

**OpenBabelTheory class:**

.. code-block:: python
    
    class OpenBabelTheory():
        def __init__(self, forcefield="UFF", chargemodel=None, label="OpenBabelTheory", 
                    printlevel=2, user_atomcharges=None):

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``printlevel``
     - integer
     - 2
     - Printlevel
   * - ``forcefield``
     - string
     - 'UFF'
     - Options: 'UFF', 'GAFF', 'MMFF94', 'Ghemical'

################################################################################
OpenBabel installation
################################################################################

Refer to installation instructions here:
https://openbabel.org/docs/Installation/install.html

A conda-forge installation (https://anaconda.org/conda-forge/openbabel) has worked for us:

.. code-block:: shell

    mamba install openbabel


################################################################################
Forcefield options
################################################################################

Please refer to the manual about the FF implementations
https://open-babel.readthedocs.io/en/latest/Forcefields/Overview.html


Note: 
While in principle different charge models can be used to modify the electrostatic interactions
that the FF may use, in practice this is not enabled yet.

################################################################################
Example usage. 
################################################################################



.. code-block:: python

    from ash import *

    # Define fragment
    frag = Fragment(xyzfile="cowley_full.xyz", charge=0, mult=1)
    # Define OpenBabelTheory object 
    # forcefield options : 'UFF', 'MMFF94', 'GAFF', 'Ghemical'
    theory = OpenBabelTheory(forcefield="UFF")

    # Optimize 
    Optimizer(theory=theory, fragment=frag)

################################################################################
Other OpenBabel features in ASH
################################################################################

OpenBabel is also occasionally used by ASH for other features.
The ability to define an ASH fragment via a SMILES string alone (see :doc:`coordinate-input`) 
is e.g. enabled by OpenBabel.

ASH has also wrapped some OpenBabel library functionality into special functions:

.. code-block:: python

    #Function to convert Mol file to PDB-file via OpenBabel
    def mol_to_pdb(file):

    #Function to convert SDF file to PDB-file via OpenBabel
    def sdf_to_pdb(file):

    #Function to read in PDB-file and write new one with CONECT lines (geometry needs to be sensible)
    def writepdb_with_connectivity(file):

    #Function to read in XYZ-file (small molecule) and create PDB-file with CONECT lines (geometry needs to be sensible)
    def xyz_to_pdb_with_connectivity(file, resname="UNL"):

    #Function to convert PDB-file to SMILES string
    def pdb_to_smiles(fname: str) -> str:

    #Function to convert SMILES string to coordinates
    def smiles_to_coords(smiles_string):