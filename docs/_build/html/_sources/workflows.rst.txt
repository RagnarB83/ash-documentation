

Workflows in ASH
======================================

As an ASH-script is pure Python, this allows one to easily create advanced workflows in a single script.

For example, a geometry optimization of a structure in on QM-program can easily be combined with a subsequent frequency job and this
can be followed by a subsequent higher-level single-point energy job using another QM-program.

Simple for-loops can also be created to run multiple jobs with slightly different parameters.

##############################################################################
Example: Running multiple single-point energies with different functionals
##############################################################################


.. code-block:: python

    from ash import *
    import sys
    import PES
    settings_ash.init() #initialize

    h2string="""
    H 0 0 0
    H 0 0 0.7
    """

    h2=Fragment(coordsstring=h2string)

    #List of functional keywords (strings) to loop over
    functionals=['BP86', 'B3LYP', 'TPSS', 'TPSSh', 'PBE0', 'BHLYP', 'CAM-B3LYP']

    #Dictionary to keep track of results
    energies_dict={}

    for functional in functionals:
        print("FUNCTIONAL: ", functional)
        orcadir='/opt/orca_4.2.1'
        #Appending functional keyword to the string-variable that contains the ORCA inputline
        input="! def2-SVP Grid5 Finalgrid6 tightscf slowconv " + functional
        blocks="""
        %scf
        maxiter 200
        end
        """
        #Defining/redefining ORCA theory. Does not need charge/mult keywords.
        ORCAcalc = ORCATheory(orcadir=orcadir, orcasimpleinput=input, orcablocks=blocks, nprocs=8, charge=0, mult=1)

        # Run single-point job
        energy = Singlepoint(theory=ORCAcalc, fragment=h2)
        #Adding energy to dictionary
        energies_dict[functional] = energy

        #Cleaning up after each job (not always necessary)
        ORCAcalc.cleanup()
        print("=================================")

    print("Dictionary with results:", energies_dict)
    print("")
    #Pretty formatted printing:
    print("")
    print(" Functional   Energy (Eh)")
    print("----------------------------")
    for func, e in energies_dict.items():
        print("{:10} {:13.10f}".format(func,e))


Producing a nice table of results:

.. code-block:: shell

     Functional   Energy (Eh)
    ----------------------------
    BP86       -1.1689426849
    B3LYP      -1.1642632249
    TPSS       -1.1734355861
    TPSSh      -1.1729787552
    PBE0       -1.1610065506
    BHLYP      -1.1624650247
    CAM-B3LYP  -1.1625896338