Workflow examples in ASH
======================================

As an ASH-script is pure Python and the user has access to various ASH functionality for reading in coordinates, create fragments,
call QM and MM codes, calculate energy, minimize geometry, run dynamics etc. this allows one to easily create advanced workflows in a single script.
ASH already contains a lot of pre-made workflow to automize things but the user can easily write their own ASH-Python scripts to automize
more complex workflows, not yet available directly in ASH.

This page goes through several examples of workflows showing both the use of pre-coded workflows and how you can write your own.


##############################################################################
Example 1 : Optimization + Frequency + HighLevel-singlepoint
##############################################################################

Running a geometry optimization, frequency calculation (on optimized geometry) using e.g. a DFT protocol and then a single-point energy calculation
using a high-level theory (such as CCSD(T)) is a standard workflow in computational chemistry research, but is often performed manually as the QM-code is not always capable of performing all of the steps in a single job.
In ASH you can use the **thermochemprotocol_single** to automate such a workflow. 

See: :doc:`module_workflows` for details.

.. code-block:: python

    from ash import *

    numcores=2
    #Defining molecular fragment
    molstring="""
    N 0 0 0
    N 0 0 0.9
    """
    molecule=Fragment(coordsstring=molstring, charge=0, mult=1)
    
    #Defining theories object for optimization+frequency and High-level CCSD(T)/CBS
    ORCAcalc = ORCATheory(orcasimpleinput="! r2SCAN-3c tightscf", orcablocks="", numcores=numcores)
    HL = ORCA_CC_CBS_Theory(elements=["N"], cardinals = [2,3], basisfamily="cc", numcores=numcores)


    energy, components, thermochem = thermochemprotocol_single(fragment=molecule, numcores=numcores, Opt_theory=ORCAcalc, SP_theory=HL)
    print("Final HL energy:", energy, "Eh")
    print("ZPVE: ", thermochem['ZPVE'], "Eh")

But you could of course also write this kind of simple workflow manually, allowing you possibly more flexibility that better suits your needs:

.. code-block:: python

    from ash import *

    numcores=2
    #Defining molecular fragment
    molstring="""
    N 0 0 0
    N 0 0 0.9
    """
    molecule=Fragment(coordsstring=molstring, charge=0, mult=1)
    #Defining ORCA object for optimization
    ORCAcalc = ORCATheory(orcasimpleinput="! r2SCAN-3c tightscf", orcablocks="", numcores=numcores)

    #Geometry optimization of molecule and ORCAcalc theory object.
    geomeTRICOptimizer(theory=ORCAcalc,fragment=molecule)

    #Numfreq job of molecule (contains optimized coordinates). A 2-point Hessian is requested in runmode serial.
    thermochem = NumFreq(fragment=molecule, theory=ORCAcalc, npoint=2, runmode='serial')

    #Single-point HL job on optimized geometry
    cc = ORCA_CC_CBS_Theory(elements=["N"], cardinals = [2,3], basisfamily="cc", numcores=numcores)
    HLresult = Singlepoint(theory=HL, fragment=molecule)

    print("Final HL energy: ", HLresult.energy, "Eh")
    print("ZPVE: ", thermochem['ZPVE'], "Eh")


#######################################################################################################
Example 2a : Direct calculation of Reaction Energy:  N\ :sub:`2` \  + 3H\ :sub:`2`\  → 2NH\ :sub:`3`\
#######################################################################################################

Typically one is interested in more than a single energy of a single species and thus you might want to run calculations on multiple species
of a reaction in the same job and calculate the reaction energy directly.
In ASH you can easily define multiple fragments in a script, group the fragments in a list and then loop over the fragments with a simple Python for-loop
that calculates the single-point energy for each fragment.

.. code-block:: python

    from ash import *

    numcores=1
    #Defining all reaction species as ASH objects from XYZ-files
    N2=Fragment(xyzfile="n2.xyz", charge=0, mult=1)
    H2=Fragment(xyzfile="h2.xyz", charge=0, mult=1)
    NH3=Fragment(xyzfile="nh3.xyz", charge=0, mult=1)

    ##Defining reaction##
    # List of species from reactant to product
    specieslist=[N2, H2, NH3] #Use same order as stoichiometry

    #Defining ORCA theory object.
    ORCAcalc = ORCATheory(orcasimpleinput="! BP86 def2-SVP def2/J", orcablocks="", numcores=numcores)

    FinalEnergies=[] #Defining empty list to collect energies
    #Python for-loop that loops over each molecule in list specieslist
    for molecule in specieslist:
        result = Singlepoint(theory=ORCAcalc, fragment=molecule)
        #Adding energy to list. Note: Energy is also stored as part of fragment.
        FinalEnergies.append(result.energy)
        ORCAcalc.cleanup()

    print("Final list of energies:", FinalEnergies)
    reaction_energy = (2*FinalEnergies[2]-(1*FinalEnergies[0]+3*FinalEnergies[1]))*627.509
    print("Reaction energy:", reaction_energy, "kcal/mol")

The script above is verbose but the structure gives you a lot of flexibility that you can adapt to your needs.
Of course, ASH already contains functions to carry out such a job in a simpler way: **Singlepoint_fragments** and **ReactionEnergy**.
See :doc:`singlepoint` and :doc:`module_workflows` for more information.

.. code-block:: python

    from ash import *

    numcores=1
    #Haber-Bosch reaction: N2 + 3H2 => 2NH3
    N2=Fragment(diatomic="N2", diatomic_bondlength=1.0975, charge=0, mult=1) #Diatomic molecules can be defined like this also
    H2=Fragment(diatomic="H2", diatomic_bondlength=0.74, charge=0, mult=1) #Diatomic molecules can be defined like this also
    NH3=Fragment(xyzfile="nh3.xyz", charge=0, mult=1)
    specieslist=[N2, H2, NH3] #An ordered list of ASH fragments.
    stoichiometry=[-1, -3, 2] #Using same order as specieslist.
    ORCAcalc = ORCATheory(orcasimpleinput="! BP86 def2-SVP def2/J", orcablocks="", numcores=numcores)
    energies = Singlepoint_fragments(theory=ORCAcalc, fragments=specieslist) #Calculating list of energies

    #Calculating reaction-energy using list and stoichiometry
    reaction_energy, unused = ReactionEnergy(stoichiometry=stoichiometry, list_of_energies=energies, unit='kcal/mol', label='ΔE')


.. code-block:: text

      Reaction_energy: -37.157156917851935 kcal/permol

#######################################################################################################
Example 2b : Direct calculation of Reaction Energy with an Automatic Thermochemistry Protocol
#######################################################################################################

You can also combine the Opt+Freq+HL protocol from Example 1 with the multiple fragments-at-the-same-time approach from Example 2
and calculate the reaction energy directly at a high-level of theory together with thermochemical corrections from a frequency job.


.. code-block:: python

    from ash import *

    numcores=4

    N2=Fragment(xyzfile="n2.xyz", charge=0, mult=1)
    H2=Fragment(xyzfile="h2.xyz", charge=0, mult=1)
    NH3=Fragment(xyzfile="nh3.xyz", charge=0, mult=1)

    ##Defining reaction
    specieslist=[N2, H2, NH3] #Use same order as stoichiometry
    stoichiometry=[-1, -3, 2] #Use same order as specieslist

    #Define theories
    OptFreqtheory = ORCATheory(orcasimpleinput="! BP86 def2-SVP def2/J", orcablocks="", numcores=numcores)
    HL=ORCA_CC_CBS_Theory(elements=["N","H"], basisfamily="cc", cardinals=[3,4], numcores=numcores)
    #Thermochemistry protocol
    thermochemprotocol_reaction(fraglist=specieslist, stoichiometry=stoichiometry, Opt_theory=OptFreqtheory, SP_theory=HL, 
                numcores=numcores, memory=5000, analyticHessian=True, temp=298.15, pressure=1.0)


Final output:

.. code-block:: text

     Reaction_energy(ΔSCF):  -33.980155385058865
     Reaction_energy(ΔCCSD):  -6.937247193220541
     Reaction_energy(Δ(T)):  1.4333499904116154
     Reaction_energy(ΔCV+SR):  -0.07653672690344188
     Reaction_energy(ΔSO):  0.0
     Reaction_energy(ΔZPVE):  20.455727327700334
    ----------------------------------------------
     Reaction_energy(Total ΔE):  -19.104861987083417

The output shows the total reaction energy (0 K enthalpy) and the contribution from Hartree-Fock (SCF), singles-doubles excitations (ΔCCSD),
perturbative triples (Δ(T)), core-valence + scalar-relativistics (CV+SR), atomic spin-orbit coupling (ΔSO, here none), and zero-point
vibrational energy (ΔZPVE).
The agreement with experiment (-18.4 kcal/mol) is excellent.


############################################################################################
Example 3a : Running multiple single-point energies with different functionals (sequential)
############################################################################################

You might be interested in running multiple single-point energy calculations on a molecule with different functionals.
Such a job could be written directly like this:

.. code-block:: python

    from ash import *
    import os

    numcores=4
    h2string="""
    H 0 0 0
    H 0 0 0.7
    """

    h2=Fragment(coordsstring=h2string, charge=0, mult=1)

    #List of functional keywords (strings) to loop over. Need to be valid ORCA keywords.
    functionals=['BP86', 'B3LYP', 'TPSS', 'TPSSh', 'PBE0', 'BHLYP', 'CAM-B3LYP']

    #Dictionary to keep track of energies
    energies_dict={}
    for functional in functionals:
        print("FUNCTIONAL: ", functional)
        #Defining/redefining ORCA theory.
        #Appending functional keyword to the string-variable that contains the ORCA inputline
        input="! def2-SVP tightscf slowconv " + functional
        ORCAcalc = ORCATheory(orcadir=orcadir, orcasimpleinput=input, numcores=numcores)
        # Run single-point job
        result = Singlepoint(theory=ORCAcalc, fragment=h2)
        #Keep ORCA outputfile for each functional
        os.rename('orca-input.out', functional+'_orcajob.out')
        #Adding energy to dictionary
        energies_dict[functional] = result.energy
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

.. code-block:: text

     Functional   Energy (Eh)
    ----------------------------
    BP86       -1.1689426849
    B3LYP      -1.1642632249
    TPSS       -1.1734355861
    TPSSh      -1.1729787552
    PBE0       -1.1610065506
    BHLYP      -1.1624650247
    CAM-B3LYP  -1.1625896338


But could also be written a bit more succinctly using the **Singlepoint_theories** function, see :doc:`singlepoint` .

.. code-block:: python

    from ash import *

    numcores=4
    #Readomg h2.xyz from internal database
    H2=Fragment(databasefile="h2.xyz", charge=0, mult=1)

    #List of functional keywords (strings) to loop over. Need to be valid ORCA keywords.
    functionals=['BP86', 'B3LYP', 'TPSS', 'TPSSh', 'PBE0', 'BHLYP', 'CAM-B3LYP']
    theories=[]
    #Looping over strings to create a list of theories
    for functional in functionals:
        input="! def2-SVP tightscf slowconv " + functional
        ORCAcalc = ORCATheory(orcadir=orcadir, orcasimpleinput=input, numcores=numcores)
        theories.append(ORCAcalc)
    #Use Singlepoint_theories to run a SP calculation on fragment with each theory
    energies = Singlepoint_theories(theories=theories, fragment=H2)

with a final table being printed:

.. code-block:: text

    ======================================================================
    Singlepoint_theories: Table of energies of each theory:
    ======================================================================

    Theory class    Theory Label     Charge    Mult           Energy(Eh)
    ----------------------------------------------------------------------
    ORCATheory      BP86                  0       1        -1.1716903176
    ORCATheory      B3LYP                 0       1        -1.1668726382
    ORCATheory      TPSS                  0       1        -1.1756974430
    ORCATheory      TPSSh                 0       1        -1.1753091518
    ORCATheory      PBE0                  0       1        -1.1636152299
    ORCATheory      BHLYP                 0       1        -1.1646744530
    ORCATheory      CAM-B3LYP             0       1        -1.1653189937

############################################################################################
Example 3b : Running multiple single-point energies with different functionals (in parallel)
############################################################################################

The examples in 3a ran each job sequentially, one after the other, according to the list of functional strings defined.
While ORCA parallelization was utilized, it may be more economical to run the jobs simultaneously instead, especially if there are lot of jobs to go through.
This can be accomplished using the **Singlepoint_parallel** function inside ASH.
Here Python multiprocessing (pool.apply_async) is utilized. In this case ORCA parallelization is by default turned off, though it can be enabled if done carefully.
See :doc:`parallelization` for more information.

.. code-block:: python

    from ash import *
    numcores=4
    #Fragment
    H2=Fragment(xyzfile="h2.xyz", charge=0, mult=1)

    #Single-point job parallelization
    #Case: Multiple theories
    #Creating multiple ORCA objects and storing in list: orcaobjects
    #Important: use a label (here functional-name)for the created ORCA object to distinguish jobs
    orcaobjects=[]
    for functional in ['BP86', 'B3LYP', 'TPSS', 'TPSSh', 'PBE0', 'BHLYP', 'CAM-B3LYP']:
        ORCAcalc = ORCATheory(orcadir=orcadir, orcasimpleinput="! def2-SVP def2/J "+functional, orcablocks="", label=functional)
        orcaobjects.append(ORCAcalc)

    #Calling the Singlepoint_parallel function and providing list of fragments and list of theories:
    results = Singlepoint_parallel(fragments=[H2], theories=orcaobjects, numcores=numcores)

    #results is a dictionary of energies
    print("results :", results)

###########################################################################################
Example 4a : Running single-point energies on a collection of XYZ files (sequential)
###########################################################################################

At other times you are interested in using a single theory to run single-point energies on a collection of molecules.
This could again be accomplished using a straightforward for-loop where we first define the path to the XYZ-file directory, change directory to it (os.chdir), 
and then use glob to find all files with an ".xyz" file suffix.
Next, loop over those XYZ-files, define a fragment from each XYZ-file and then run a single-point calculation.

.. code-block:: python

    from ash import *
    import glob

    numcores=1
    #Directory of XYZ files. Can be full path or relative path.
    dir = '/path/to/xyz_files'
    #Changing to dir
    os.chdir(dir)

    energies=[]
    for file in glob.glob('*.xyz'):
        print("XYZ-file:", file)
        mol=Fragment(xyzfile=file, charge=0, mult=1) #Note: Here we have to assume that charge=0 and mult=1 for every molecule
        ORCAcalc = ORCATheory(orcasimpleinput="! BP86 def2-SVP def2/J", orcablocks="", numcores=numcores)
        result = Singlepoint(theory=ORCAcalc, fragment=mol)
        print("Energy of file {} : {} Eh".format(file, result.energy))
        ORCAcalc.cleanup()
        energies.append(energy)
        print("")
    #Pretty print
    print(" XYZ-file             Energy (Eh)")
    print("-"*50)
    for xyzfile, e in zip(glob.glob('*.xyz'),energies):
        print("{:20} {:>13.10f}".format(xyzfile,e))


Output:

.. code-block:: text

     XYZ-file             Energy (Eh)
    -------------------------------------------------
    h2.xyz               -1.1715257797
    h2o_MeOH.xyz         -192.0023991603
    O-O-dimer.xyz        -149.8555328055
    butane.xyz           -158.3248873844
    nh3.xyz              -56.5093301286
    n2.xyz               -109.4002969311
    hi.xyz               -298.3735362292
    h2o_strained.xyz     -76.2253312246

Again, we can simplify the script with the help of built-in ASH functionality: **read_xyzfiles** and **Singlepoint_fragments**
See: :doc:`coordinate-input` and :doc:`singlepoint` for more information.

.. code-block:: python

    from ash import *

    numcores=1
    #Directory of XYZ files. Can be full path or relative path.
    xyzdir = '/path/to/xyz_files'

    #This function reads in all XYZ-files from the chosen directory and returns a list of ASH fragments
    #Note: Each XYZ-file must have charge/mult defined in 2nd line of header for readchargemult=True to work
    fragments = read_xyzfiles(xyzdir, readchargemult=True, label_from_filename=True)
    #Define theory
    ORCAcalc = ORCATheory(orcasimpleinput="! BP86 def2-SVP def2/J", orcablocks="", numcores=numcores)

    #Call Singlepoint_fragments and get list of calculated energies at chosen theory
    energies = Singlepoint_fragments(theory=ORCAcalc, fragments=fragments)


Output:

.. code-block:: text

    ======================================================================
    Singlepoint_fragments: Table of energies of each fragment:
    ======================================================================
    Formula    Label                 Charge    Mult           Energy(Eh)
    ----------------------------------------------------------------------
    H1I1       hi.xyz                     0       1      -298.3737182333
    N1H3       nh3.xyz                    0       1       -56.5093324450
    H2         h2.xyz                     0       1        -1.1715262206
    H2O1       h2o_strained.xyz           0       1       -76.2253299452
    C4H10      butane.xyz                 0       1      -158.3249141864
    O2         O-O-dimer.xyz              0       1      -149.8555433766
    N2         n2.xyz                     0       1      -109.4002889693
    H6O2C1     h2o_MeOH.xyz               0       1      -192.0023967568

Such a protocol can be further simplified using the **calc_xyzfiles** function that even allows you to even run a thermochemistry workflow on each XYZ-file. 
See: :doc:`module_workflows`

.. code-block:: python

    from ash import *

    numcores=1
    #Directory of XYZ files. Can be full path or relative path.
    xyzdir = '/path/to/xyz_files'

    #Define theory
    ORCAcalc = ORCATheory(orcasimpleinput="! BP86 def2-SVP def2/J", orcablocks="", numcores=numcores)

    #Call calc_xyzfiles giving xyzdir and theory.
    #Geometry optimizations for each XYZ-file can be requested via Opt=True (default False, i.e. singlepoint)
    calc_xyzfiles(xyzdir=dir, theory=ORCAcalc)

############################################################################################
Example 4b : Running single-point energies on a collection of XYZ files (parallel)
############################################################################################
The examples in 4a had each job run sequentially, one job after the other, according to the list of XYZ-files available.
While ORCA parallelization was utilized, it may be more economical to run such embarrassingly parallel jobs simultaneously instead, especially if there are lot of XYZ-files.
This can be accomplished using the Singlepoint_parallel function inside ASH. This utilizes Python multiprocessing (pool.apply_async).
ORCA parallelization is here turned off.

See :doc:`parallelization` for more information.

.. code-block:: python

    from ash import *
    import glob
    
    numcores=4
    ORCAcalc = ORCATheory(orcadir=orcadir, orcasimpleinput="! BP86 def2-SVP def2/J", orcablocks="", numcores=1)
    #Directory of XYZ files. Can be full path or relative path.
    xyzdir = './xyz_files'

    molecules = read_xyzfiles(xyzdir,readchargemult=True, label_from_filename=True)

    #Calling the Singlepoint_parallel function and providing list of fragments and list of theories:
    results = Singlepoint_parallel(fragments=molecules, theories=[ORCAcalc], numcores=numcores)

    #results is a dictionary of energies
    print("results :", results)



###########################################################################################################
Example 5 : Calculate localized orbitals and create Cube files for multiple XYZ files or an XYZ-trajectory
###########################################################################################################

Analyzing electronic structure along a reaction path (e.g. a NEB or IRC path) or a trajectory (optimization or MD)
can be useful to understand the nature of the reaction. The code below shows how this can be accomplished in ASH
via a workflow involving single-point DFT, orbital localization and Cube-file creation (via orca_plot).

See :doc:`ORCA-interface` for information on **run_orca_plot**.

Using a collection of XYZ-files:

.. code-block:: python

    from ash import *
    import glob
    #
    numcores=1
    #Directory of XYZ files. Can be full path or relative path.
    dir = '/home/bjornsson/ASH-DEV_GIT/testsuite/localized-orbital-IRC-workflow/calcs/images'

    #Localization block in ORCA inputfile
    blockinput="""
    %loc
    LocMet IAOIBO
    end
    """

    #Looping over XYZ-files in directory, creating ASH fragments, running ORCA and calling orca_plot
    for file in sorted(glob.glob(dir+'/*.xyz')):
        basefile=os.path.basename(file)
        print("XYZ-file:", basefile)
        mol=Fragment(xyzfile=file)
        ORCAcalc = ORCATheory(orcasimpleinput="! BP86 def2-SVP def2/J", orcablocks=blockinput, numcores=numcores)
        result = Singlepoint(theory=ORCAcalc, fragment=mol, charge=-1, mult=1)
        print("Energy of file {} : {} Eh".format(basefile, result.energy))
        locfile=basefile.split('.')[0]+'_calc.loc'
        os.rename('orca-input.loc', locfile)
        #Call ORCA_plot and create Cube file for specific MO in locfile: here alpha-MOs 13 and 17
        run_orca_plot(orcadir, locfile, 'mo', gridvalue=30, mo_operator=0, mo_number=13)
        run_orca_plot(orcadir, locfile, 'mo', gridvalue=30, mo_operator=0, mo_number=17)

        ORCAcalc.cleanup()
        print("")


Using a multi-XYZ file containing multiple sets of geometries (could be a NEB path, MD/Opt trajectory, XYZ animation etc.)

.. code-block:: python

    from ash import *
    #
    numcores=1

    #Name of trajectory file containing multiple geometries (could be optimization traj, MD traj, NEB-path traj, Hessian XYZ animation etc.)
    #File should be in dir
    trajectoryfile="neb-ts_MEP_trj.xyz"

    blockinput="""
    %loc
    LocMet IAOIBO
    end
    """

    fraglist = get_molecules_from_trajectory(trajectoryfile)

    for index,frag in enumerate(fraglist):
        print("Frag :", index)
        ORCAcalc = ORCATheory(orcasimpleinput="! BP86 def2-SVP def2/J", orcablocks=blockinput, numcores=numcores)
        result = Singlepoint(theory=ORCAcalc, fragment=frag, charge=-1, mult=1)
        print("Energy of frag {} : {} Eh".format(index, result.energy))
        locfile='frag{}_calc.loc'.format(index)
        os.rename('orca-input.loc', locfile)
        #Call ORCA_plot and create Cube file for specific MO in locfile: here alpha-MOs 13 and 17
        run_orca_plot(orcadir, locfile, 'mo', gridvalue=30, mo_operator=0, mo_number=13)
        run_orca_plot(orcadir, locfile, 'mo', gridvalue=30, mo_operator=0, mo_number=17)

        ORCAcalc.cleanup()


###########################################################################################
Example 6 : Running conformer-sampling, geometry optimizations and High-level single-points
###########################################################################################

This example utilizes the interface to the powerful `crest <https://xtb-docs.readthedocs.io/en/latest/crest.html>`_ program to perform metadynamics-based conformational sampling from a starting geometry at a semi-empirical level of theory (GFN2-xTB).
From the conformational sampling we get a collection of low-energy conformers for that molecule (based on the GFN1-xTB or GFN2-xTB semi-empirical tightbinding Hamiltonian)
The conformational sampling is then followed by a DFT geometry optimization for each conformer.
Finally high-level coupled cluster single-point calculations (here DLPNO-CCSD(T)/CBS extrapolations) are performed for each conformer.

Such an example can be written in ASH like this in a rather verbose manner:

.. code-block:: python

    from ash import *

    numcores=4

    #0. Starting structure and charge and mult
    charge=0
    mult=1
    molecule = Fragment(xyzfile="ethanol.xyz", charge=charge, mult=mult)

    #1. Calling crest
    call_crest(fragment=molecule, xtbmethod='GFN2-xTB', numcores=numcores)

    #2. Grab low-lying conformers from crest_conformers.xyz as list of ASH fragments.
    list_conformer_frags, xtb_energies = get_crest_conformers()

    print("list_conformer_frags:", list_conformer_frags)
    print("")
    print("Crest Conformer Searches done. Found {} conformers".format(len(xtb_energies)))
    print("xTB energies: ", xtb_energies)

    #3. Run DFT geometry optimizations for each crest-conformer
    #ML Theory level. TODO: Run in ASH parallel instead of ORCA parallel?
    MLorcasimpleinput="! BP86 D3 def2-TZVP def2/J Grid5 Finalgrid6 tightscf"
    MLorcablocks="%scf maxiter 200 end"
    MLORCATheory = ORCATheory(orcasimpleinput=MLorcasimpleinput, orcablocks=MLorcablocks, numcores=numcores)

    DFT_energies=[]
    print("")
    for index,conformer in enumerate(list_conformer_frags):
        print("")
        print("Performing DFT Geometry Optimization for Conformer ", index)
        geomeTRICOptimizer(fragment=conformer, theory=MLORCATheory, coordsystem='tric', charge=charge, mult=mult)
        DFT_energies.append(conformer.energy)
        #Saving ASH fragment and XYZ file for each DFT-optimized conformer
        os.rename('Fragment-optimized.ygg', 'Conformer{}_DFT.ygg'.format(index))
        os.rename('Fragment-optimized.xyz', 'Conformer{}_DFT.xyz'.format(index))

    print("")
    print("DFT Geometry Optimization done")
    print("DFT_energies: ", DFT_energies)

    #4.Run high-level DLPNO-CCSD(T). Ash should now have optimized conformers
    HL_CC = ORCA_CC_CBS_Theory(elements=molecule.elems, cardinals = [2,3], basisfamily="cc", DLPNO=True, numcores=numcores)
    HL_energies=[]
    for index,conformer in enumerate(list_conformer_frags):
        print("")
        print("Performing High-level calculation for DFT-optimized Conformer ", index)
        HLresult = Singlepoint(theory=HL_CC, fragment=conformer, charge=charge, mult=mult)
        HL_energies.append(HLresult.energy)


    print("")
    print("=================")
    print("FINAL RESULTS")
    print("=================")

    #Printing total energies
    print("")
    print(" Conformer   xTB-energy    DFT-energy    HL-energy (Eh)")
    print("----------------------------------------------------------------")

    min_xtbenergy=min(xtb_energies)
    min_dftenergy=min(DFT_energies)
    min_HLenergy=min(HL_energies)

    for index,(xtb_en,dft_en,HL_en) in enumerate(zip(xtb_energies,DFT_energies, HL_energies)):
        print("{:10} {:13.10f} {:13.10f} {:13.10f}".format(index,xtb_en, dft_en, HL_en))

    print("")
    #Printing relative energies
    min_xtbenergy=min(xtb_energies)
    min_dftenergy=min(DFT_energies)
    min_HLenergy=min(HL_energies)
    harkcal = 627.50946900
    print(" Conformer   xTB-energy    DFT-energy    HL-energy (kcal/mol)")
    print("----------------------------------------------------------------")
    for index,(xtb_en,dft_en,HL_en) in enumerate(zip(xtb_energies,DFT_energies, HL_energies)):
        rel_xtb=(xtb_en-min_xtbenergy)*harkcal
        rel_dfT=(dft_en-min_dftenergy)*harkcal
        rel_HL=(HL_en-min_HLenergy)*harkcal
        print("{:10} {:13.10f} {:13.10f} {:13.10f}".format(index,rel_xtb, rel_dfT, rel_HL))

    print("")
    print("Workflow done!")


The manually defined workflow above is a bit verbose and can of course also be more conveniently run like below where we use the 
**confsampler_protocol** function (see :doc:`crest-interface`), that takes as input the ASH fragment and 2 levels of ASH theories to be used for geometry optimizations and high-level singlepoint energies.

.. code-block:: python

    from ash import *

    numcores=4
    #Fragment to define. Here taken from internal database
    molecule=Fragment(databasefile="ethanol.xyz")

    #Defining MLTheory: DFT optimization
    MLsimpleinput="! B3LYP D4 def2-TZVP TightSCF"
    MLblockinput="""
    %scf maxiter 200 end
    """
    ML_B3LYP = ORCATheory(orcasimpleinput=MLsimpleinput, orcablocks=MLblockinput, numcores=numcores)
    #Defining HLTheory: DLPNO-CCSD(T)/CBS
    HL_CC = ORCA_CC_CBS_Theory(elements=molecule.elems, cardinals = [2,3], basisfamily="cc", DLPNO=True, numcores=numcores)
    #Call confsampler_protocol
    confsampler_protocol(fragment=molecule, xtbmethod='GFN2-xTB', MLtheory=ML_B3LYP,
                            HLtheory=HL_CC, numcores=numcores)


Final result table of total and relative energies of calculated conformers at 3 different theory levels:

.. code-block:: text

    =================
    FINAL RESULTS
    =================

     Conformer   xTB-energy    DFT-energy    HL-energy (Eh)
    ----------------------------------------------------------------
             0 -25.8392205500 -346.2939482921 -345.2965932205
             1 -25.8377914500 -346.2884905132 -345.2911748671
             2 -25.8358803400 -346.2818766960 -345.2848279253
             3 -25.8313250600 -346.2788608396 -345.2815202116
             4 -25.8307377800 -346.2788662649 -345.2815419285
             5 -25.8303374700 -346.2775476223 -345.2792917601
             6 -25.8300128900 -346.2776089771 -345.2794648759

     Conformer   xTB-energy    DFT-energy    HL-energy (kcal/mol)
    ----------------------------------------------------------------
             0  0.0000000000  0.0000000000  0.0000000000
             1  0.8967737821  3.4248079602  3.4000680178
             2  2.0960134034  7.5750408530  7.3828340833
             3  4.9544947374  9.4675192805  9.4584557521
             4  5.3230184983  9.4641148891  9.4448282319
             5  5.5742168139 10.2915756050 10.8568301896
             6  5.7778938373 10.2530749008 10.7481984235
