Singlepoint
======================================

The Singlepoint jobtype is the most basic job to perform in ASH.
The main function is **Singlepoint** and special versions are **Singlepoint_fragments**, **Singlepoint_theories** and **Singlepoint_parallel**.


########################
Singlepoint function
########################

The base function is a simple Python function that only requires as input : an ASH fragment object and an ASH Theory object.
The ASHTHeory object can be a QMTheory from: :doc:`QM-interfaces`, or an
MMTheory (see :doc:`MM-interfaces`) or even a QM/MMTheory (see :doc:`module_QM-MM`).

.. code-block:: python

    def Singlepoint(fragment=None, theory=None, Grad=False, charge=None, mult=None):
        """Singlepoint function: runs a single-point energy calculation using ASH theory and ASH fragment.

        Args:
            fragment (ASH fragment): An ASH fragment. Defaults to None.
            theory (ASH theory): Any valid ASH theory. Defaults to None.
            Grad (Boolean): Do gradient or not. Defaults to False.
            charge (int, optional): Specify charge of system. Overrides fragment charge information.
            mult (int, optional): Specify mult of system. Overrides fragment charge information.            

        Returns:
            float: Energy
            or
            float,np.array : Energy and gradient array
        """

*Example*
In the example script below an ASH fragment is created from the XYZ-file "hf.xyz" that contains the Cartesian coordinates of hydrogen fluoride.
The charge and multiplicity is also defined as part of the Fragment object.
Next a theory-level object is defined, here an object is created from the ORCATheory class. 
For a single-point calculation one then simply passes the Theory object and the Fragment object to the **Singlepoint** function.

.. code-block:: python

    from ash import *

    #Defining fragment
    HF_frag=Fragment(xyzfile="hf.xyz", charge=0, mult=1)
    #ORCA
    orcasimpleinput="! BP86 def2-SVP def2/J tightscf"
    ORCAobject = ORCATheory(orcasimpleinput=orcasimpleinput, numcores=4)

    #Simple Energy SP calc. Energy will be printed to output and also returned as an energy variable
    Singlepoint(theory=ORCAobject, fragment=HF_frag)


Note that an alternative to defining the charge/mult as attributes of the Fragment is to instead pass the charge and multiplicity to the Singlepoint function.

.. code-block:: python

    #Simple Energy SP calc. Energy will be printed to output and also returned as an energy variable
    Singlepoint(theory=ORCAobject, fragment=HF_frag, charge=0, mult=1)

If this option is used, this will override any charge/mult information present in the fragment.

The energy calculated is printed to standard output by default.
The **Singlepoint** function returns an energy as a floating point number and usually you want to store the energy as a new variable so that you can do something more with the data.

.. code-block:: python

    #Simple Energy+Gradient SP calc
    # The function will return the energy that can be stored as a variable
    Energy = Singlepoint(theory=ORCAobject, fragment=HF_frag)
    print("Energy is", Energy)

It is also possible to request a gradient calculation in which case both the energy and gradient is returned:

.. code-block:: python

    #Simple Energy+Gradient SP calc
    Energy, Gradient = Singlepoint(theory=ORCAobject, fragment=HF_frag, Grad=True)
    print("Energy is", Energy)
    print("Gradient is:", Gradient)


By default, the files created by the Theory interface are not cleaned up. To have ORCA (in this example) clean up
temporary files (e.g. so they don't interfere with a future job), one can use the cleanup function:

.. code-block:: python

    #Clean up
    ORCAobject.cleanup()


The energy and gradient from the last Energy/Energy+Gradient run is also stored inside the Theory object and can be accessed:

.. code-block:: python

    print(ORCAobject.energy)
    print(ORCAobject.grad)

##################################
Singlepoint_fragments function
##################################

The **Singlepoint** function above is designed to be a simple function that does one job, returning 1 energy for the 1 theory level and the 1 fragment that was defined.
In a typical project, however, multiple calculations need to be performed. For example running the same single-point theory calculation on multiple fragments.

You could of course easily write a for-loop for this purpose in ASH, making sure to define first charge and multiplicity for each fragment first.

.. code-block:: python
    
    from ash import *
	
    #Species of the Haber-Bosch reaction: N2 + 3H2 => 2NH3
    N2=Fragment(diatomic="N2", diatomic_bondlength=1.0975, charge=0, mult=1)
    H2=Fragment(diatomic="H2", diatomic_bondlength=0.741, charge=0, mult=1)
    NH3=Fragment(xyzfile="nh3.xyz", charge=0, mult=1)

    fragment_list=[N2, H2, NH3] #An ordered list of ASH fragments.

    #Define theory
    xtbcalc=xTBTheory(xtbmethod='GFN1') # GFN1-xTB theory-level
    energies=[] #empty list to store energies

    #Iterating over fragments
    for fragment in fragment_list:
        energy = Singlepoint(theory=xtbcalc, fragment=fragment)
        energies.append(energy) #add energy to list

    print("List of energies:", energies)


More conveniently, however, you can instead use the **Singlepoint_fragments** function:

.. code-block:: python

    def Singlepoint_fragments(theory=None, fragments=None, stoichiometry=None):


that does the same thing:

.. code-block:: python

    from ash import *

    #Species of the Haber-Bosch reaction: N2 + 3H2 => 2NH3
    N2=Fragment(diatomic="N2", diatomic_bondlength=1.0975, charge=0, mult=1)
    H2=Fragment(diatomic="H2", diatomic_bondlength=0.741, charge=0, mult=1)
    NH3=Fragment(xyzfile="nh3.xyz", charge=0, mult=1)
    specieslist=[N2, H2, NH3] #An ordered list of ASH fragments.
    xtbcalc=xTBTheory(xtbmethod='GFN1') # GFN1-xTB theory-level

    #Call Singlepoint_fragments and get list of calculated energies
    energies = Singlepoint_fragments(theory=xtbcalc, fragments=specieslist)

In addition to returning a list of energies, a table is also printed in standard output:

.. code-block:: text

    ============================================================
    Singlepoint_fragments: Table of energies of each fragment:
    ============================================================
    Formula    Label       Charge    Mult           Energy(Eh)
    ------------------------------------------------------------
    N2         None             0       1        -6.3335016263
    H2         None             0       1        -1.0361629322
    N1H3       nh3              0       1        -4.8298958374

If you provide (optional) a stoichiometry (list order should match fragments list) to **Singlepoint_fragments** you will also get a print-out of the reaction energy.

.. code-block:: python

    energies = Singlepoint_fragments(theory=xtbcalc, fragments=specieslist, stoichiometry=[-1,-3,2])


.. code-block:: text

    Stoichiometry provided.
    Reaction_energy(ΔE):  -136.6723479900558 kcal/mol

##################################
Singlepoint_theories function
##################################

You might also have a single fragment that you want to run multiple single-point theory calculations on.
For this case you can use **Singlepoint_theories** instead.

.. code-block:: python

    #Define the fragment. Here providing charge/mult also.
    N2=Fragment(diatomic="N2", diatomic_bondlength=1.0975, charge=0, mult=1)

    #Defining theory levels. NOTE: For clearer printing it is recommended to add a label to each theory object.
    gfn1_xtbcalc=xTBTheory(xtbmethod='GFN1', label='GFN1-xTB') # GFN1-xTB theory-level
    gfn2_xtbcalc=xTBTheory(xtbmethod='GFN2', label='GFN2-xTB') # GFN2-xTB theory-level
    orca_r2scan=ORCATheory(orcasimpleinput='! r2SCAN-3c tightscf', label='ORCA-r2SCAN-3c') # ORCA r2SCAN-3c theory-level
    theories=[gfn1_xtbcalc,gfn2_xtbcalc,orca_r2scan] #Collecting all theories in a list

    energies = Singlepoint_theories(theories=theories, fragment=N2)

In addition to returning a list of energies, a table is also printed in standard output:

.. code-block:: text

    ======================================================================
    Singlepoint_theories: Table of energies of each theory:
    ======================================================================

    Theory class    Theory Label     Charge    Mult           Energy(Eh)
    ----------------------------------------------------------------------
    xTBTheory       GFN1-xTB              0       1        -6.3335016263
    xTBTheory       GFN2-xTB              0       1        -5.7639339581
    ORCATheory      ORCA-r2SCAN-3c        0       1      -109.5070425194


#############################################
Singlepoint_fragments_and_theories function
#############################################

You might even want to perform calculation on multiple fragments with multiple theories. For example calculating a reaction energy with multiple theory levels.
**Singlepoint_fragments_and_theories** makes this easy.

.. code-block:: python

    def Singlepoint_fragments_and_theories(theories=None, fragments=None, stoichiometry=None):


Example:

.. code-block:: python

    from ash import *

    #Haber-Bosch reaction: N2 + 3H2 => 2NH3
    N2=Fragment(diatomic="N2", diatomic_bondlength=1.0975, charge=0, mult=1)
    H2=Fragment(diatomic="H2", diatomic_bondlength=0.741, charge=0, mult=1)
    NH3=Fragment(xyzfile="nh3.xyz", charge=0, mult=1)
    specieslist=[N2, H2, NH3] #An ordered list of ASH fragments.
    stoichiometry=[-1, -3, 2] #Using same order as specieslist.
    xtbcalc=xTBTheory(xtbmethod='GFN1') # GFN1-xTB theory-level

    #Defining theories
    gfn1_xtbcalc=xTBTheory(xtbmethod='GFN1', label='GFN1-xTB') # GFN1-xTB theory-level
    gfn2_xtbcalc=xTBTheory(xtbmethod='GFN2', label='GFN2-xTB') # GFN2-xTB theory-level
    orca_r2scan=ORCATheory(orcasimpleinput='! r2SCAN-3c tightscf', label='ORCA-r2SCAN-3c') # ORCA r2SCAN-3c theory-level

    #All theories in a list
    theories=[gfn1_xtbcalc,gfn2_xtbcalc,orca_r2scan]

    # Running multiple fragments and theories
    results = Singlepoint_fragments_and_theories(theories=theories, fragments=specieslist, stoichiometry=stoichiometry)

This gives the output:

.. code-block:: text

    ============================================================
    Singlepoint_fragments_and_theories: FINAL RESULTS
    ============================================================

    Theory: xTBTheory
    Label: GFN1-xTB

    ============================================================
    Table of energies of each fragment:
    ============================================================
    Formula    Label       Charge    Mult           Energy(Eh)
    ------------------------------------------------------------
    N2         None             0       1      -109.5070425194
    H2         None             0       1        -1.1693814360
    H3N1       nh3              0       1       -56.5418434618

    Stoichiometry provided: [-1, -3, 2]
    Reaction_energy(GFN1-xTB):  -136.6723479900558 kcal/mol
    ____________________________________________________________

    Theory: xTBTheory
    Label: GFN2-xTB

    ============================================================
    Table of energies of each fragment:
    ============================================================
    Formula    Label       Charge    Mult           Energy(Eh)
    ------------------------------------------------------------
    N2         None             0       1      -109.5070425194
    H2         None             0       1        -1.1693814360
    H3N1       nh3              0       1       -56.5418434618

    Stoichiometry provided: [-1, -3, 2]
    Reaction_energy(GFN2-xTB):  -89.09909008527589 kcal/mol
    ____________________________________________________________

    Theory: ORCATheory
    Label: ORCA-r2SCAN-3c

    ============================================================
    Table of energies of each fragment:
    ============================================================
    Formula    Label       Charge    Mult           Energy(Eh)
    ------------------------------------------------------------
    N2         None             0       1      -109.5070425194
    H2         None             0       1        -1.1693814360
    H3N1       nh3              0       1       -56.5418434618

    Stoichiometry provided: [-1, -3, 2]
    Reaction_energy(ORCA-r2SCAN-3c):  -42.98445901864511 kcal/mol
    ____________________________________________________________

    Final list of list of total energies: [[-6.333501626274, -1.036162932168, -4.829895837389], 
        [-5.763933958102, -0.9820230341, -4.4259957498], [-109.507042519379, -1.16938143601, -56.541843461825]]
    Final reaction energies:
    Reaction_energy(GFN1-xTB):  -136.6723479900558 kcal/mol
    Reaction_energy(GFN2-xTB):  -89.09909008527589 kcal/mol
    Reaction_energy(ORCA-r2SCAN-3c):  -42.98445901864511 kcal/mol

A final list of lists of total energies is returned (each list containing the total energies of the fragment for each theory level )

##################################
Singlepoint_parallel function
##################################

The **Singlepoint_fragments** and **Singlepoint_theories** functions perform the calculations in a sequential fashion (via a for loop): i.e. one calculation after the other.
While convenient, the functions do not utilize the fact that each fragment-calculation (**Singlepoint_fragments**) or theory-calculation (**Singlepoint_theories**) is completely 
independent from each other and could thus run through the list of calculations (whether fragments or theories) in parallel on a multi-core CPU.
The **Singlepoint_parallel** function, however, allows you to do this.

See :doc:`parallelization` for information on using the **Singlepoint_parallel** function.