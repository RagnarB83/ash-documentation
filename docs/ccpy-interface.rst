ccpy interface
======================================

`ccpy <https://github.com/piecuch-group/ccpy>`_ is a very recent Python-based coupled cluster code, developed by the research group of Piotr Piecuch and written by Karthik Gururangan.
The strength of the program lies in the availability of a large number of classic and modern coupled cluster methods, 
including the `Adaptive CC(P;Q) method <https://pubs.aip.org/aip/jcp/article/159/8/084108/2907794/Converging-high-level-coupled-cluster-energetics>`_ , CIPSI-driven CC(P;Q) or active-space CC(P;Q), not available anywhere else,
as well as being decoupled from the SCF and integral calculations and being Python-based (except for the core CC routines).
The program hence relies on SCF and integrals from another program or an FCIDUMP file.
It features a nice Python API that made it easy to integrate into ASH.

ASH features a basic interface to ccpy, primarily intended to make it easier to run some of the calculations in a more black-box way and allowing for integration
into other ASH workflows. ASH also features a way of running ccpy calculations using a reference calculated by ORCA.


**ccpyTheory class:**

.. code-block:: python


  class ccpyTheory:
      def __init__(self, pyscftheoryobject=None, orcatheoryobject=None, orca_gbwfile=None, orca_jsonformat="msgpack",
                  fcidumpfile=None, filename=None, printlevel=2, label="ccpy", delete_json=True,
                  frozencore=True, cc_tol=1e-8, numcores=1, dump_integrals=False,
                  cc_maxiter=300, cc_amp_convergence=1e-7, nact_occupied=None, nact_unoccupied=None, civecs_file=None, 
                  method=None, percentages=None, states=None, roots_per_irrep=None, EOM_guess_symmetry=False,
                  two_body_approx=False):

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``method``
     - string
     - None
     - Name of CC method to use. See below for all available methods.
   * - ``pyscftheoryobject``
     - PySCTheory
     - None
     - A PySCFTheory object that will be used to provide reference for ccpy
   * - ``orcatheoryobject``
     - ORCATheory
     - None
     - An ORCATheory object that will be used to provide reference for ccpy
   * - ``orca_gbwfile``
     - string
     - None
     - Name or path of a ORCA GBW file that will be used to provide reference for ccpy
   * - ``orca_jsonformat``
     - string
     - "msgpack"
     - Format of the intermediate ORCA JSON file. Can be "json", "bson" or "msgpack".
   * - ``delete_json``
     - Boolean
     - True
     - Whether to delete the ORCA JSON file after reading it.
   * - ``dump_integrals``
     - Boolean
     - False
     - Whether to dump processed integrals to file before ccpy starts.
   * - ``fcidumpfile``
     - string
     - None
     - Name or path of a FCIDUMP file that integrals will be read from.
   * - ``printlevel``
     - integer
     - 2
     - Printlevel
   * - ``numcores``
     - integer
     - 1
     - Number of CPU cores ccpy will use
   * - ``frozencore``
     - Boolean
     - True
     - Whether to use a frozen-core approximation or not. Will enforce ASH frozen core settings (based on ORCA settings).
   * - ``cc_tol``
     - float
     - 1e-8
     - Energy tolerance for convergence 
   * - ``cc_maxiter``
     - integer
     - 300
     - Max number of CC iterations
   * - ``cc_amp_convergence``
     - float
     - 1e-7
     - Convergence tolerance for CC amplitudes
   * - ``percentages``
     - list
     - None
     - Adaptive CC only: List of triples percentages to include during CC(P;Q). Example: [0.0, 0.5]

   * - ``nact_occupied``
     - integer
     - None
     - Active-space CC only: Number of occupied orbitals to use.
   * - ``nact_unoccupied``
     - integer
     - None
     - Active-space CC only: Number of unoccupied orbitals to use.
   * - ``civecs_file``
     - integer
     - None
     - CIPSI-driven CC only: File containing CI vectors to use.
   * - ``states``
     - list
     - None
     - EOM-CC only: List of states to calculate.
   * - ``roots_per_irrep``
     - dict
     - None
     - Number of EOM-CC roots per irrep to calculate. Example: {'A1':2, 'B1':1}

################################
Installation
################################

Follow the installation instructions at the `ccpy documentation page <https://piecuch-group.github.io/ccpy/installation.html>`_.

Briefly, ccpy requires a few libraries installable via conda and then compilation via *make* and installation within the Python environment (as a callable ccpy library).

As ccpy relies on Numpy with a MKL backend (instead of OpenBLAS) that may require some tweaking, it may be best to 
create a specific conda environment for ccpy and ASH. It is strongly recommended to install pyscf as well within the same environment, but after the installation of Numpy with MKL backend, to avoid pyscf installing another Numpy version.

####################################
Available methods in ccpy interface
####################################

The ASH interface to ccpy is designed to make it easier to use the methods of ccpy in a black-box way, i.e. call a method by its name and to use ccpy in ASH workflows.
This is a deliberate design-choice that has the downside that not every method of ccpy is yet available in the ASH interface.
It is best to use ccpy directly if more control over ccpy options are desired.

The following methods are available in the ASH interface (*method* name shown within quotes):

- **Regular CC methods:** "CCD", "CCSD", "CCSD(T)", "CCSDT", "CCSDTQ", "CC3" and "CC4".
- **Active-space CC methods:** "CCSDT1" (a.k.a. CCSDt), "CCT3" (a.k.a. CC(t;3)), "CCSDT_P" (CCSDt using CC(P) solver).
- **Completely renormalized CC methods:** "CRCC23", "CRCC24".
- **Adaptive CC(P;Q) method:** "adaptive-cc(p;q)""
- **CIPSI-driven methods:** "cipsi-cc(p;q)", "eccc23" , "eccc24" (a.k.a. CIPSI-CC(P;Q), ec-CC-II, ec-CC-II3)
- **EOM-CC methods:** "eomccsd", "eomcc3", "eomccsdt", "eomccsdt(a)_star", "creom23", "ipeom2", "ipeom3", "ipeomccsdta_star".


################################
Using the interface
################################

ccpy only solves the coupled-cluster problem and requires MOs and integrals from another program.
Currently, ccpy is limited to supporting RHF or ROHF reference wavefunctions (no UHF or broken-symmetry WFs).

ASH currently supports 3 ways of providing the necessary information to ccpy:

**via PySCFTheory**

By providing a :doc:`PySCF-interface` object when creating the ccpyTheory object, ccpyTheory will upon first run, request the run of the chosen SCF calculation (a HF/DFT RHF/ROHF SCF)
via pySCFTheory and then pass the resulting information to ccpy.
This is the easiest way as pyscf-installation requires only a quick: pip install pyscf in the same environment. The pyscf SCF-step will run quickly and then exchanges the MOs and integrals to ccpy in-memory.
Molecular symmetry information can be provided to the PySCFTheory object and this information will be picked up by ccpy.

.. code-block:: python

  from ash import *

  # RHF/cc-pVTZ PySCFTheory object
  pyscfobj = PySCFTheory(scf_type="RHF", basis="cc-pVTZ", symmetry="C1")
  # A CCSD ccpyTheory object using pyscfobj as input
  ccpy_theory = ccpyTheory(method="CCSD", pyscftheoryobject=pyscfobj, frozencore=True, 
              cc_tol=1e-10, numcores=1, cc_maxiter=300)

**via ORCATheory or ORCA GBW-file**

By providing a :doc:`ORCA-interface` object when creating the ccpyTheory object, an ORCA calculation will first be run (a HF/DFT RHF/ROHF SCF). 
Once the ORCA run is complete, the ORCA GBW file is automatically converted into a JSON-file, containing MOs and integrals. 
The JSON-file contents are processed and used to setup the ccpy calculation before running the ccpy coupled cluster job.
Because of the  GBW->JSON conversion and integral processing, using ORCATheory will take quite a bit longer than using PySCFTheory (above).
ASH supports reading ORCA-created JSON files in regular ASCII JSON, BSON and MSGPack formats. The MSGPack format is particularly fast and space-efficient.
As the ASCII JSON format is slow to process and can take up large space, a binary JSON-like format like MSGPack is strongly recommended.
The format is controlled by the *orca_jsonformat* keyword, the default value is : 'msgpack'. It does require installation of either the msgpack (pip install msgpack) 
or msgspec (pip install msgspec) library. The JSON or MSGPack file is by default deleted once the information has been read, this can be overridden by setting *delete_json*=False.


*Basic RHF example:*

.. code-block:: python

  from ash import *

  # RHF/def2-SVP ORCATheory object
  orcatheoryobj = ORCATheory(orcasimpleinput="! RHF def2-SVP tightscf")
  # A CCSD ccpyTheory object using orcatheoryobj as input
  ccpy_theory = ccpyTheory(method="CCSD", orcatheoryobject=orcatheoryobj, orca_jsonformat="msgpack",
              delete_json=True,
              frozencore=True, cc_tol=1e-10, numcores=1, cc_maxiter=300)

*An open-shell ROHF example:*

As ccpy only supports RHF and ROHF reference WFs, the ORCA reference calculation should be run as ROHF (rather than UHF) if the system is openshell. 
This also means that antiferromagnetic broken-symmetry UHF singlet reference WFs are not possible to calculate.

Due to a current bug in ORCA 6.0.0, the GBW->JSON conversion does not work for open-shell ROHF GBW files.
A current workaround is to non-iteratively convert the converged ROHF WF to a UHF WF. This can be accomplished automatically in ORCATheory by setting the *ROHF_UHF_swap* keyword to True.


.. code-block:: python

  from ash import *

  # ROHF/def2-SVP ORCATheory object using the ROHF_UHF_swap option
  orcatheoryobj = ORCATheory(orcasimpleinput="! ROHF def2-SVP tightscf", ROHF_UHF_swap=True)
  # A CCSD ccpyTheory object using orcatheoryobj as input
  ccpy_theory = ccpyTheory(method="CCSD", orcatheoryobject=orcatheoryobj, orca_jsonformat="msgpack",
              frozencore=True, cc_tol=1e-10, numcores=1, cc_maxiter=300)


*Using a natural orbital reference from a GBW-file*

An alternative to using the GBW-file created by running the ORCATheory object is to use a natural orbital reference WF created by ORCA.
The *ORCA_orbital_setup* function can be used to automatically create a natural orbital GBW-file from a correlated ORCA WF.
Such a GBW-file can be fed automatically to ccpyTheory.

.. code-block:: python

  from ash import *

  #Run ORCA calculation to get a natural orbital GBW file.
  newmofile, nat_occupations = ORCA_orbital_setup(orbitals_option="MP2", fragment=frag, basis="def2-SVP", MP2_density="unrelaxed",
                  charge=frag.charge, mult=frag.mult)

  # A CCSD ccpyTheory object using orcatheoryobj as input
  ccpy_theory = ccpyTheory(method="CCSD", orca_gbwfile=newmofile, orca_jsonformat="msgpack",
              frozencore=True, cc_tol=1e-10, numcores=1, cc_maxiter=300)


**via FCIDUMP-file**

An alternative way to use ccpyTheory involves providing a FCIDUMP file that contains the integrals.
An FCIDUMP file simply contains all the 1-electron and 2-electron integrals necessary for the CC calculation to be carried out.
Such a file can in principle be created by multiple QM programs.

Beware that FCIDUMP files can become very large and the FCIDUMP approach is only a viable option for small systems!
It is generally better to use the *pyscftheoryobject*, *orcatheoryobject* or *orca_gbwfile* options discussed above instead.

.. code-block:: python

  from ash import *

  theory = ccpyTheory(method="CCSD", fcidumpfile="FCIDUMP-file", frozencore=True, 
              cc_tol=1e-10, numcores=1, cc_maxiter=300)


As an example one, can create an FCIDUMP-file from PySCFTheory like this:

.. code-block:: python

  from ash import *

  pyscfobj = PySCFTheory(scf_type="RHF", basis="cc-pVTZ", symmetry="C1")
  Singlepoint(theory=pyscfobj, fragment=al2h2)
  pyscfobj.create_fcidump_file(filename="FCIDUMP.pyscf")

One can also create FCIDUMP-files from ORCA via **create_ORCA_FCIDUMP** using any GBW-file from ORCA.
Unfortunately the ORCA GBW->JSON conversion that **create_ORCA_FCIDUMP** uses, does not work for open-shell ROHF GBW-files (see above).

.. code-block:: python

  from ash import *

  orcatheory = ORCATheory(orcasimpleinput="! RHF def2-SVP")
  Singlepoint(theory=orcatheory, fragment=frag)
  create_ORCA_FCIDUMP(orcatheory.filename+'.gbw', header_format="FCIDUMP", filename="FCIDUMP_ORCA",
                          int_threshold=1e-16, mult=1)


It is also possible to use the ASH function **ORCA_orbital_setup** to conveniently create reference orbitals for ccpy from ORCA
that can then be used to create a FCIDUMP file.

.. code-block:: python

  from ash import *

  #Run ORCA calculation to get a natural orbital GBW file.
  newmofile, nat_occupations = ORCA_orbital_setup(orbitals_option="MP2", fragment=frag, basis="def2-SVP", MP2_density="unrelaxed",
                  charge=frag.charge, mult=frag.mult)

  #Create FCIDUMP file from MP2 natural orbital GBW-file
  create_ORCA_FCIDUMP(newmofile, header_format="FCIDUMP", filename="FCIDUMP_ORCA",
                          int_threshold=1e-16, mult=1)

################################
Calculate natural orbitals
################################

Natural orbitals can be requested for some of the methods in ccpy using the **get_natural_orbitals** method of a ccpyTheory object.
An input density is required which can be calculated by calling the *run_density* method of a ccpyTheory object.
A singlepoint energy calculation must have been carried out first and additionally the Hbar and LeftCC equations must have been solved before.

Example:

.. code-block:: python

  al2h2 = Fragment(databasefile="h2o.xyz", charge=0, mult=1)

  pyscfobj = PySCFTheory(scf_type="RHF", basis="cc-pVTZ", symmetry="C1")
  theory = ccpyTheory(pyscftheoryobject=pyscfobj, filename='input.dat', printlevel=2, label="ccpy",
                frozencore=True, cc_tol=1e-10, numcores=1,
                cc_maxiter=300, method="ccsd")
  #Run regular single-point energy
  Singlepoint(theory=theory, fragment=frag)
  # Run Hbar (similarity transformation) and LeftCC equations
  theory.driver.run_hbar(method="ccsd")
  theory.driver.run_leftcc(method="left_ccsd", state_index=[0])
  #Run density calculation
  rdm_matrix = theory.run_density()
  #Diagonalizes rdm to get natural orbitals (requires previous MO coefficients)
  natocc, natorb = theory.get_natural_orbitals(rdm_matrix, mo_coeffs=pyscfobj.mf.mf_coeff)
  #Write out Molden file of the natural orbitals
  theory.write_molden_file(natocc,natorb,mo_energies=None,label="CCSD")


################################
Examples
################################

**Regular CC calculations: CCSD(T)**

It is easy to use various regular coupled cluster methods such as: "CCD", "CCSD", "CCSD(T)", "CCSDT", "CCSDTQ", "CC3" and "CC4".
The method name simply needs to be provided to the ccpyTheory object as well as either a PySCFTheory object or a FCIDUMP file.

.. code-block:: python

  from ash import *

  frag = Fragment(databasefile="h2o.xyz"))

  #Reference WF via pySCF
  pyscfobj = PySCFTheory(scf_type="RHF", basis="def2-SVP", symmetry="C1")
  #CCpy object creation
  theory = ccpyTheory(method="CCSD(T)", pyscftheoryobject=pyscfobj, frozencore=True, 
              cc_tol=1e-10, numcores=1, cc_maxiter=300)

  result = Singlepoint(theory=theory, fragment=frag)

Example using ORCA as reference:

.. code-block:: python

  from ash import *

  frag = Fragment(databasefile="h2o.xyz"))

  #Reference WF via ORCA SCF
  orcatheory = ORCATheory(orcasimpleinput="! RHF def2-SVP")
  Singlepoint(theory=orcatheory, fragment=frag)

  #Create FCIDUMP file
  create_ORCA_FCIDUMP(orcatheory.filename+'.gbw', header_format="FCIDUMP", filename="FCIDUMP_ORCA",
                          int_threshold=1e-16, mult=1)
  #ccpy
  theory = ccpyTheory(method="CCSD(T)", fcidumpfile="FCIDUMP_ORCA", frozencore=True, 
              cc_tol=1e-10, numcores=1, cc_maxiter=300)

  result = Singlepoint(theory=theory, fragment=frag)

**Completely renormalized CC calculations: CR-CC(2,3)**

*Literature:*

`CR-CC(2,3) <https://pubs.aip.org/aip/jcp/article-abstract/123/22/224105/776282/Renormalized-coupled-cluster-methods-exploiting?redirectedFrom=fulltext>`_

Completely renormalized (CR) coupled cluster methods were designed to improve the performance of single-reference CC 
(in particular non-iterative corrections to CCSD like CCSD(T)) for electronically more complicated situations (bond stretching etc.).
CR-CC(2,3) in particular is a rigorously size-extensive formulation of CR-CCSD(T), 
based on the method of moments of CC equations, CC(P;Q), developed by P. Piechuch.
CR-CC(2,3) has often been found to be more reliable than CCSD(T) for electronically tricky situations.


.. code-block:: python

  from ash import *

  frag = Fragment(databasefile="h2o.xyz"))

  pyscfobj = PySCFTheory(scf_type="RHF", basis="def2-SVP", symmetry="C1")
  theory = ccpyTheory(method="CRCC23", pyscftheoryobject=pyscfobj, frozencore=True, 
              cc_tol=1e-10, numcores=1, cc_maxiter=300)

  result = Singlepoint(theory=theory, fragment=frag)

**Active-space coupled cluster methods**

*Literature:*

- `CCSDt <https://www.tandfonline.com/doi/abs/10.1080/00268976.2010.522608>`_
- `CC(t;3)-theory <https://www.sciencedirect.com/science/article/abs/pii/S0301010411005283?via%3Dihub>`_
- `CC(t;3)-benchmarks <https://pubs.aip.org/aip/jcp/article-abstract/136/14/144104/190989/Combining-active-space-coupled-cluster-methods>`_

Active-space coupled cluster (or active-orbital based) methods were designed to deal with electronically more complicated situations
where the perturbative correction (T) (such as in CCSD(T)) fails and explicit iterative triples may be required. Explicit quadruples may also be included in this way.
As full iterative CCSDT is too expensive, active-space CC allows selected explicit triples to be included via the definition of an active space. 

This has led to a hierarchy of methods in the literature and the following are available in ccpy:

- "CCSDT1" ; known as CCSDt in the literature. Dominant T3 components chosen by active space.
- "CCT3" ; the CC(t;3) method based on the CC(P;Q) moment expansion. Goes beyond CCSDt.
- "CCSDT_P" ; the CCSDt method but uses the general CC(P) solver.

These methods can be specified by the *method* keyword but also requires the definition of the active space
via the keywords *nact_occupied* and *nact_unoccupied*. This is similar to the definition of a CAS space in CASSCF.
If *nact_occupied* and *nact_unoccupied* are chosen to include the full orbital space, the methods should reduce to full CCSDT.

The drawback of these methods is that they are not as black-box as regular CC methods and require some insight into 
choosing the active space. The active-space aspect becomes even slightly problematic when calculating chemical reactions
where the active space will necessarily require redefinition for each species.

.. code-block:: python

  from ash import *

  frag = Fragment(databasefile="h2o.xyz"))

  pyscfobj = PySCFTheory(scf_type="RHF", basis="def2-SVP", symmetry="C1")
  theory = ccpyTheory(method="CCT3", pyscftheoryobject=pyscfobj,  frozencore=True, 
              cc_tol=1e-10, numcores=1, cc_maxiter=300)
  result = Singlepoint(theory=theory, fragment=frag)

**CIPSI-driven coupled cluster methods**

*Literature:*

- `CIPSI-driven CC(P;Q) <https://pubs.aip.org/aip/jcp/article/155/17/174114/565630>`_
- `EC-CC method analysis <https://pubs.acs.org/doi/full/10.1021/acs.jctc.1c00181>`_

ccpy also includes methods where the CC(P;Q) moment expansion equations are solved using information from a selected CI (CIPSI) calculation.
This allows the proper relaxation of T1 and T2 amplitudes in the presence of the most important T3 or T4 components,, according to information from the CI expansion.
This allows well-defined approximations to CCSDT or CCSDTQ to be calculated at a fraction of the cost of the full method and would have the advantage
of "only" requiring a CI-calculation. The precedingCI calculation could be performed with or without an active space.

The CIPSI-driven methods available in the ASH interface are:

- "cipsi-cc(p;q)", CIPSI-driven CC(P;Q).
-  "eccc23" , "eccc24" (a.k.a. ec-CC-II, ec-CC-II3)

These methods require the presence of a CI-vectors file, to be provided via the *civecs_file* keyword but do not require the definition of an active space 
(a possible active space definition would only apply to the preceding CI calculation).
An example of a civecs file can be found in the `ccpy repository <https://github.com/piecuch-group/ccpy/blob/main/tests/data/h2o/civecs-10k.dat>`_ 

Future versions of ASH will hopefully include an automatic way of generating the CI-vectors file using other programs.

.. code-block:: python

  from ash import *

  frag = Fragment(databasefile="h2o.xyz"))

  pyscfobj = PySCFTheory(scf_type="RHF", basis="def2-SVP", symmetry="C1")
  theory = ccpyTheory(method="cipsi-cc(p;q)", pyscftheoryobject=pyscfobj,  frozencore=True, 
              cc_tol=1e-10, numcores=1, cc_maxiter=300)
  result = Singlepoint(theory=theory, fragment=frag)

**Adaptive CC(P;Q) methods**

*Literature:*

- `Adaptive CC(P;Q) <https://pubs.aip.org/aip/jcp/article-abstract/159/8/084108/2907794/Converging-high-level-coupled-cluster-energetics>`_

The adaptive CC(P;Q) approach represents in some ways a considerable breakthrough for the approaches of methods previously discussed.
Similar to the active-space or CIPSI-based approaches above, explicit triples or quadruples amplitudes are included in a CC(P;Q) moment expansion
but instead of requiring an active-space (non-blackbox) or an external CI calculation for the selection,
adaptive CC(P;Q) automatically selects the necessary components based on the structure of the moment expansions.

It is easy to use the method with the ASH interface, as the method-name ("cipsi-cc(p;q)") simply needs to be chosen together with an estimate of the %-number of triples
to be included in the form of a list of percentages starting from 0 (e.g. [0.0, 0.5, 1.0, 2.0,3.0]). 
This results in an almost fully automated black-box method, requiring only a well-defined choice of the %-triples to include.

.. code-block:: python

  from ash import *

  frag = Fragment(databasefile="h2o.xyz"))

  pyscfobj = PySCFTheory(scf_type="RHF", basis="def2-SVP", symmetry="C1")
  theory = ccpyTheory(method="cipsi-cc(p;q)", percentages=[0.0, 0.5, 1.0], pyscftheoryobject=pyscfobj,  frozencore=True, 
              cc_tol=1e-10, numcores=1, cc_maxiter=300)
  result = Singlepoint(theory=theory, fragment=frag)

**EOM methods**

Equation of motion coupled cluster is a methodology to calculate excited states within the CC framework.
The ASH interface to ccpy allows for the calculation of various EOM-CC methods:

- "eomccsd", i.e. EOM-CCSD
- "eomcc3", i.e. IP-EOM-CC3
- "eomccsdt", i.e. EOM-CCSDT
- "eomccsdt(a)_star", i.e. EOM-CCSDT(a)*
- "creom23", i.e. CREOM(2,3)
- "ipeom2", i.e. IP-EOMCCSD(2h-1p)
- "ipeom3", i.e. IP-EOMCCSD(3h-2p)
- "ipeomccsdta_star", i.e. IP-EOM-CCSDT(a)*


To use, you specify the *method* keyword as one of the EOM-methods above.
Addtionally you need to specify the states to calculate via the *states* keyword and roots per irrep via the *roots_per_irrep* keyword.
The states are specified as a list of integers and the roots per irrep as a dictionary with the irrep as key and the number of roots as value.
For the case of no point-group symmetry (C1) the irrep is simply "A".

Finally, most of the EOM methods require the EOM-guess to be specified in more detail.
This requires definition of an active space using *nact_occupied* and *nact_unoccupied* keywords.
Using this active-space definition a specific diagonalization is carried out to obtain the EOM-guess.
The results may be sensitive to the quality of the guess.

.. code-block:: python

  from ash import *

  frag = Fragment(databasefile="h2o.xyz"))

  pyscfobj = PySCFTheory(scf_type="RHF", basis="cc-pVTZ", symmetry="C1")
  theory = ccpyTheory(pyscftheoryobject=pyscfobj, frozencore=True, cc_tol=1e-10, numcores=1,
                cc_maxiter=300, roots_per_irrep={'A':3},nact_occupied=3, nact_unoccupied=7,
                method="eomccsd", states=[0,1,2])
  result = Singlepoint(theory=theory, fragment=frag)