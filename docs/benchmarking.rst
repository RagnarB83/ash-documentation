Benchmarking in ASH
======================================

ASH contains convenient tools to do energy-benchmarking over test sets.
Available test sets are in ASH-dir/databases/Benchmarking-sets.

Currently available test sets are:

*Reaction energies*

-   MB08-165 (from `GMTKN30 <https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/GMTKN/gmtkn30>`_ , random reactions, reference-data: CCSD(T)/CBS)
-   `MOR41 <https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/mor41/metal-organic-reactions-mor>`_ (closed-shell organometallic reactions, reference-data: DLPNO-CCSD(T)/CBS)
-   S22 (from `GMTKN30 <https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/GMTKN/gmtkn30>`_, noncovalent interactions, reference-data: CCSD(T)/CBS)
-   SIE11 (from `GMTKN30 <https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/GMTKN/gmtkn30>`_, self-interaction error dominated reactions, reference-data: CCSD(T)/CBS)

*Electron affinities*

-   G21EA (from `GMTKN30 <https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/GMTKN/gmtkn30>`_, reference-data: back-corrected experiment)

*Ionization energies*

-   G21IP (from `GMTKN30 <https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/GMTKN/gmtkn30>`_, reference-data: back-corrected experiment)
-   IE-Pantazis (from `Isegawa et al. <https://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b00252>`_, reference-data: experiment)
-   IE-benzenes (reference-data: experiment)

#########################################
The module_benchmarking module
#########################################

The module_benchmarking module contains :

    - run_benchmark
    - run_geobenchmark (not ready)


#########################################
The run_benchmark function
#########################################

The run_benchmark function needs at minimum the set keyword argument and either a theory or workflow argument.

.. code-block:: python

    def run_benchmark(set=None, theory=None, workflow=None, orcadir=None, numcores=None, reuseorbs=False, corrections=None)

- set: Name of benchmark test set.
- theory: ASH Theory object
- workflow: name of valid function that can act as theory level
- orcadir : Only needed for workflow. Path to ORCA (string).
- numcores: Only needed for workflow. Number of cores to use.
- reuseorbs: Whether orbitals should be reused for each species in reaction. Only makes sense if geometries are similar (e.g. IE/EA reactions). Boolean True/False.
- corrections: Corrections to be applied to the calculated reaction energies, e.g. ZPE or thermal correction etc. List of floats. Can also be defined within testset.
- workflow_args: Optional workflow arguments if using workflow and the default settings need to be modified.

#########################################
Running a test set with a chosen QMtheory
#########################################

Running a test set is easy. First define a theory object (e.g. an ORCATheory object).
Then call the run_benchmark function choosing a specific testset as set keyword argument and the theory object as theory keyword argument.
The chosen set has to be a valid directory that is inside :
ASH-dir/databases/Benchmarking-sets
with a valid directory structure (see more info below).

.. code-block:: python

    from ash import *

    #Some variables for ORCA
    orcadir='/Applications/orca_4.2.1'
    numcores=4
    orcasimpleinput="! BP86 def2-SV(P) tightscf"
    orcablocks="%scf maxiter 200 end"

    #Define theory level for benchmarks
    ORCAcalc = ORCATheory(orcadir=orcadir, charge=0, mult=1, orcasimpleinput=orcasimpleinput, orcablocks=orcablocks, nprocs=numcores)

    #Running the benchmark
    run_benchmark(set="IE-benzenes", theory=ORCAcalc)


Running with added corrections and reusing orbitals within each reaction (convenient for IE reactions).
Note: Corrections can also be applied by providing a corrections.txt file to the test set (see below).

.. code-block:: python

    #ZPE corrections per reaction in eV. Theory level: XXX
    ZPE_corrections =[0.10, 0.01, 0.11, 0.12, 0.13]
    run_benchmark(set="IE-benzenes", theory=ORCAcalc, corrections=ZPE_corrections, reuseorbs=True)


Output:

.. code-block:: text

     ======================================================================
     FINAL RESULTS FOR TESTSET:   IE-benzenes
     ======================================================================
    Unit: eV
     Index   Reaction                                                 Ref.          Calc.         Calc.+corr.     Error
    ------------------------------------------------------------------------------------------------------------------------
     1       fluorobenzene-neut ⟶   fluorobenzene-ox                  9.2032        8.9576        9.0576         -0.1456
     2       benzene-neut ⟶   benzene-ox                              9.2438        9.1534        9.1634         -0.0803
     3       chlorobenzene-neut ⟶   chlorobenzene-ox                  9.0728        8.7748        8.8848         -0.1880
     4       bromobenzene-neut ⟶   bromobenzene-ox                    8.9975        8.6682        8.7882         -0.2093
     5       iodobenzene-neut ⟶   iodobenzene-ox                      8.7580        8.4440        8.5740         -0.1840
    ------------------------------------------------------------------------------------------------------------------------
     MAE               0.1614 eV
     ME               -0.1614 eV
     RMSE              0.1677 eV
     MaxError         -0.2093 eV


#########################################
Running a test set with a workflow
#########################################

The test set can also be run with a high-level workflow (multi-step theory).
The workflows are available inside the module: module_highlevel_workflows

- W1theory
- W1F12theory
- DLPNO_W1F12theory
- DLPNO_W1theory
- DLPNO_F12
- DLPNO_W2theory
- DLPNO_CC_CBS

When using an ORCA-based workflow the orcadir keyword argument and numcores argument needs to provided.

.. code-block:: python

    from ash import *

    orcadir='/Applications/orca_4.2.1'
    #Running the benchmark with a workflow
    run_benchmark(set="IE-benzenes", workflow=module_highlevel_workflows.DLPNO_W1theory, numcores=4, orcadir=orcadir)

If some of the default settings of each workflow needs to be modified this can be accomplished like this:

.. code-block:: python

    from ash import *

    orcadir='/Applications/orca_4.2.1'

    #Define a dictionary containing the arguments of the workflow to be modified
    DLPNO_CC_CBS_args = {'cardinals' : '2/3', "basisfamily" : "def2", 'stabilityanalysis' : True,
                        'memory' : 5112, 'extrablocks' : "%scf\ndirectresetfreq 1\nend\n", 'extrainputkeyword' : 'Slowconv'}

    #Running the benchmark with a workflow
    benchmarking.run_benchmark(set="IE-benzenes", workflow=module_highlevel_workflows.DLPNO_CC_CBS, workflow_args = DLPNO_CC_CBS_args, orcadir=orcadir, numcores=numcores)


################################
Creating or modifying a test set
################################

Each directory inside ASH-dir/databases/Benchmarking-sets is a separate benchmarking database for a group of molecular reactions.
Each testset-directory, e.g. "IE-benzenes" should contain a README file and a directory called data.
The README file should contain human-readable basic information about the dataset.
The data directory should contain XYZ-files for the dataset and a file: "Reference_data.txt" that contains definitions about the reactions.

*Example:*

::

    IE-benzenes
    ├── README
    └── data
       ├── benzene-neut.xyz
       ├── benzene-ox.xyz
       ├── etc.
       ├── Reference_data.txt
       └── corrections.txt (optional file)

**IMPORTANT**: Each XYZ-file should contain the charge and multiplicity in the title-line (2nd header-line of XYZ file format)


The **Reference_data.txt** contains information about the reactions in the following format:

    - The #TESTSET_INFO lines contain information on the number of reactions and the unit for the reference data.
      These special lines are read and parsed by ASH.
    - Other # lines are convenient comment-lines but are not read by ASH.
    - Each numbered line defines a reaction. The ASCII-string words (must contain a non-numeric character) in the line point to XYZ-files in the same dir
      while the integers indicate the stoichiometry of the reaction (negative number: reactant, positive number: product).
      The last floating point number is always the reference value (e.g. experimental value) in the unit indicated in the #TESTSET_INFO line.

If the **corrections.txt** file is present inside data dir (this is optional) then additive corrections per reaction will be read when
run_benchmark is run.
This correction can e.g. be ZPE, total enthalpy-correction, total free-energy correction etc.


**Reference_data.txt** example:

.. code-block:: text

    #TESTSET_INFO Numentries: 5
    #TESTSET_INFO Unit: eV
    #X-benzenes. Geometries: B3LYP-D3/def2-TZVP
    1 fluorobenzene-neut fluorobenzene-ox -1 1 9.2032
    2 benzene-neut benzene-ox -1 1 9.24378
    3 chlorobenzene-neut chlorobenzene-ox -1 1 9.0728
    4 bromobenzene-neut bromobenzene-ox -1 1 8.9975
    5 iodobenz


**corrections.txt** example:

.. code-block:: text


    ##################
    #TESTSET_INFO
    #TESTSET_INFO Numentries: 5
    #TESTSET_INFO Unit: eV
    #TESTSET_INFO Type: ZPE
    # ZPE corrections per reaction to be added to calculated reaction energies
    1 0.012
    2 0.013
    3 0.009
    4 0.010
    5 0.010

