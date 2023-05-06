NWChem interface
======================================

`NWChem <https://github.com/nwchemgit/nwchem>`_  is an open-source quantum chemistry code full of features.

ASH features a flexible interface to it that allows you to use DFT and WFT based methods within the program.
Energies and gradients are available in the interface so NWChemTheory in ASH can be used for single-point energies, geometry optimizations, 
numerical frequencies, surface scans, NEB and molecular dynamics. Pointcharge-support is also available so NWChemTheory can be used with QMMMTheory.
Currently, only the molecular code (not the plane-wave NWPW code) is supported.

**NWChemTheory class:**

.. code-block:: python
    
    class NWChemTheory:
        def __init__(self, nwchemdir=None, filename='nwchem', openshell=False, printlevel=2,
                    nwcheminput=None, method='scf', tce=False, numcores=1):

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``nwchemdir``
     - string
     - None
     - Directory where NWChem binaries are.
   * - ``filename``
     - string
     - 'pyscf'
     - Filename used for NWChem input/output files.
   * - ``nwcheminput``
     - Multiline-string
     - None
     - String containing NWChem block-input desired. Must contain basis set definition etc.
   * - ``openshell``
     - Boolean
     - False
     - Whether the calculation is open-shell or closed-shell: affects whether the restricted or unrestricted code is called.
   * - ``method``
     - string
     - 'scf'
     - The name of the NWChem method or task. Options: 'scf', 'dft', 'sodft', 'mp2, 'ccsd', 'ccsd(t)' and other NWChem tasks. See https://nwchemgit.github.io/TASK.html for options.
   * - ``tce``
     - Boolean
     - False
     - Whether the Tensor Contraction Engine (TCE) code is used for wavefunction methods or not. See: https://nwchemgit.github.io/TCE.html for options
   * - ``printlevel``
     - integer
     - 2
     - Printlevel
   * - ``numcores``
     - integer
     - 1
     - The number of CPU cores used for MPI parallelization.


################################################################################
NWChem installation
################################################################################

NWChem needs to be installed separately. 
See: `The NWChem Github repository <https://github.com/nwchemgit/nwchem>`_ and `Download options <https://nwchemgit.github.io/Download.html>`_ 
A convenient option is to install the program using conda in the same conda environment as perhaps was used to install OpenMM:

.. code-block:: text

    conda install -c conda-forge nwchem

Once NWChem has been installed you need to either :
i) use the nwchemdir keyword option to NWChemTheory to provide the location of the directory where the nwchem binary resides 
or 
ii) Make sure nwchem is in the PATH variable of your shell environment. ASH will try to find an executable 'nwchem' automatically when creating the NWChemTheory object.

################################################################################
Parallelization
################################################################################

NWChem can be parallelized using OpenMPI. Make sure OpenMPI is installed and that the mpirun command is in your environment.
Then provide the desired number of cores to NWChemTheory using the numcores keyword.
When ASH calls NWChem it will launch a command like this: mpirun -np nwchem file.nw


################################################################################
Examples
################################################################################

**DFT-SCF example:**

.. code-block:: python

    from ash import *

    #Variable defining number of cores
    numcores=4

    #Fragment
    frag = Fragment(databasefile="h2o.xyz", charge=0, mult=1)

    #Multi-line inputstring to define the main NWChem options
    #No geometry-input or TASK line should be provided here (handled by ASH)
    input_string=""""
    basis spherical 
      o library aug-cc-pvdz
      h library cc-pvdz
    end

    dft
      xc b3lyp
    end

    scf
      thresh 1e-10
    end
    """"

    nwchem = NWChemTheory(method='dft', nwcheminput=input_string, numcores=numcores)

    Singlepoint(theory=nwchem, fragment=frag)

**TCE CREOM-CCSD(T) example:**

.. code-block:: python

    from ash import *

    #Variable defining number of cores
    numcores=4

    #Fragment
    frag = Fragment(databasefile="h2o.xyz", charge=0, mult=1)

    #Multi-line inputstring to define the main NWChem options
    #No geometry-input or TASK line should be provided here (handled by ASH)
    input_string=""""
    basis spherical 
      o library aug-cc-pvdz
      h library cc-pvdz
    end

    tce
      freeze atomic
      creomccsd(t)
      tilesize 20
      2eorb
      2emet 13
      eomsol 2
    end

    """"

    nwchem = NWChemTheory(method='tce', tce=True, nwcheminput=input_string, numcores=numcores)

    Singlepoint(theory=nwchem, fragment=frag)
