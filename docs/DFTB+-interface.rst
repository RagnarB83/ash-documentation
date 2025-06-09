DFTB+ interface
======================================

`DFTB+ <https://dftbplus.org>`_  is a density functional tight-binding (DFTB) program, intended for running semi-empirical DFT calculations.
Since DFTB can run at least 2-orders of magnitude faster than a DFT calculation it is an attractive option for 
DFTB+ can be used for both regular DFTB Hamiltonians as well as extended tightbinding Hamiltonians (xTB).
See also xTB interface for xTB calculations.

The ASH interface to DFTB+ is fairly new and not very extensive.
It supports energies and gradients and can thus be used within ASH for optimization, frequencies, dynamics etc.
It also supports pointcharge embedding and can thus be used for electrostatic embedding QM/MM.


**DFTBTheory class:**

.. code-block:: python
    
    class DFTBTheory():
        def __init__(self, dftbdir=None, hamiltonian="XTB", xtb_method="GFN2-xTB", printlevel=2, label="DFTB",
                    numcores=1, slaterkoster_dict=None, maxmom_dict=None, Gauss_blur_width=0.0,
                    SCC=True, ThirdOrderFull=False, ThirdOrder=False):

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``dftbdir``
     - string
     - None
     - Directory where DFTB+ binaries are.
   * - ``printlevel``
     - integer
     - 2
     - Printlevel
   * - ``numcores``
     - integer
     - 1
     - The number of CPU cores used.
   * - ``hamiltonian``
     - string
     - 'XTB'
     - Type of Hamiltonian to use. Options: 'DFTB', 'XTB'
   * - ``xtb_method``
     - string
     - 1
     - For XTB hamiltonian, what xTB method to use. Options: 'GFN1-xTB', 'GFN2-xTB', 'IPEA1-xTB'
   * - ``slaterkoster_dict``
     - dict
     - None
     - Dictionary pointing to Slaterkoster parameter files for each pair of elements.
   * - ``SCC``
     - Boolean
     - True
     - Whether self-consistence charge option should be enabled or not. 
   * - ``ThirdOrderFull``
     - Boolean
     - False
     - Whether full third-order term should be included or not.
   * - ``ThirdOrder``
     - Boolean
     - False
     - Whether on-site third-order term should be included or not. Typically not recommended.
   * - ``maxmom_dict``
     - dict
     - None
     - Dictionary pointing to max mom for each element.
   * - ``Gauss_blur_width``
     - float
     - 0.0
     - Gaussian blur width for pointcharge embedding 


################################################################################
DFTB+ installation
################################################################################

See DFTB+ download and installation page: https://dftbplus.org/download/index.html

In addition to the program,  DFTB parameters (Slater-Koster files) need to be downloaded for the specific parameterization.
Popular ones include e.g. '3ob', 'ob2', 'mio'
See Github repository
https://github.com/dftbparams


################################################################################
How to use 
################################################################################

The ASH interface to DFTB+ is basic and currently allows basic usage of the implemented DFTB and xTB methods.
ASH takes care of creating the DFTB+ inputfile.

DFTB and xTB methods inside DFTB+ are first differentiated by selecting *hamiltonian* and then by different options.

To use xTB methods within the DFTB program you select *hamiltonian* to be 'XTB' and then select the appropriate *xtb_method* string.

.. code-block:: python

    # GFN1-xTB
    DFTBTheory(hamiltonian="XTB", xtb_method="GFN1-xTB")

    # GFN2-xTB
    DFTBTheory(hamiltonian="XTB", xtb_method="GFN2-xTB")

To use DFTB methods you select hamiltonian to be 'DFTB' and you can then enable/disable SCC and ThirdOrderFull options to
get the desired DFTB Hamiltonian

.. code-block:: python

    # Original DFTB, a.k.a. DFTB1
    DFTBTheory(hamiltonian="DFTB", SCC=False)

    # DFTB2 a.k.a. SCC-DFTB.
    DFTBTheory(hamiltonian="DFTB", SCC=True)

    # DFTB3 a.k.a. SCC-DFTB with full third-order term
    DFTBTheory(hamiltonian="DFTB", SCC=True, ThirdOrderFull=True)

However, to use DFTB methods one must additionally specify the specific DFTB parameter set by defining a dictionary for all
element-pairs of the system and point to the individual Slater-Koster files.
The parameter sets must be downloaded separately.
Additionally third-order DFTB (DFTB3) requires the specification of atomic Hubbard derivatives,
and additionally, depending on the parameter set, a hydrogen bond correction may also have to be specified.


################################################################################
Examples
################################################################################

In the examples below, the DFTB parameter sets have already been downloaded from:
https://github.com/dftbparams

**DFTB2 with mio paremeter set**

.. code-block:: python

    from ash import *

    #H2O fragment
    frag = Fragment(databasefile="h2o.xyz")


    # DFTB2-mio
    skdir="/Users/rb269145/ash-tests/dftb_interface/dftb_for_website/DFTB2-mio/mio-main/skfiles"
    sldict_mio= {'O-O':f'{skdir}/O-O.skf', 'H-O':f'{skdir}/H-O.skf',
            'O-H':f'{skdir}/O-H.skf', 'H-H':f'{skdir}/H-H.skf'}

    # Defining DFTB2-mio Hamiltonian
    theory = DFTBTheory(hamiltonian="DFTB", SCC=True, slaterkoster_dict=sldict_mio)

    Singlepoint(theory=theory, fragment=frag)

**DFTB3 with 3ob paremeter set**

For DFTB3 calculations we have to enable the third-order term by *ThirdOrderFull* keyword.
We also have to use a compatible parameter set and here we use the `3ob set <https://github.com/dftbparams/3ob>`_ . 
Additionally, because of the third-order term we have to provide Hubbard derivatives for each element.
This information should be available with the parameter set (here in the README file of 3ob).
We pass the Hubbard derivatives as a dictionary with the *hubbard_derivs_dict* keyword.
Finally, we should also enable damping of the hydrogen interaction and here we 
set *hcorrection_zeta* to be 4.0 as recommended (see README of 3ob).

.. code-block:: python

    from ash import *

    #H2O fragment
    frag = Fragment(databasefile="h2o.xyz")


    # DFTB3-3ob
    skdir="/Users/rb269145/ash-tests/dftb_interface/dftb_for_website/DFTB3-3ob/3ob-main/skfiles"
    sldict_3ob = {'O-O':f'{skdir}/O-O.skf', 'H-O':f'{skdir}/H-O.skf',
            'O-H':f'{skdir}/O-H.skf', 'H-H':f'{skdir}/H-H.skf'}

    hubbard_derivs_dict={'O':-0.1575, 'H':-0.1857}
    zeta=4.0 #damping parameter for H

    # Defining DFTB3-3ob Hamiltonian
    theory = DFTBTheory(hamiltonian="DFTB", SCC=True, ThirdOrderFull=True,  slaterkoster_dict=sldict_3ob, 
       hubbard_derivs_dict=hubbard_derivs_dict, hcorrection_zeta=zeta)

    Singlepoint(theory=theory, fragment=frag)