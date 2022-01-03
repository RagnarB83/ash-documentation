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

    def Singlepoint(fragment=None, theory=None, Grad=False):
        """Singlepoint function: runs a single-point energy calculation using ASH theory and ASH fragment.

        Args:
            fragment (ASH fragment): An ASH fragment. Defaults to None.
            theory (ASH theory): Any valid ASH theory. Defaults to None.
            Grad (Boolean): Do gradient or not. Defaults to False.

        Returns:
            float: Energy
            or
            float,np.array : Energy and gradient array
        """

*Example*
In the example script below an ASH fragment is created from the XYZ-file "hf.xyz" that contains the Cartesian coordinates of hydrogen fluoride.
Next a theory-level object is defined, here an object is created from the ORCATheory class. 

The charge and multiplicity is defined when creating the Theory object.
Note that an alternative is to define the charge/mult as attributes of the Fragment instead (see later).

For a single-point calculation one then simply passes the Theory object and the Fragment object to the **Singlepoint** function.

.. code-block:: python

    from ash import *

    #Defining fragment
    HF_frag=Fragment(xyzfile="hf.xyz")
    #ORCA
    orcasimpleinput="! BP86 def2-SVP def2/J tightscf"
    ORCAobject = ORCATheory(charge=0, mult=1, orcasimpleinput=orcasimpleinput, numcores=4)

    #Simple Energy SP calc. Energy will be printed to output and also returned as an energy variable
    Singlepoint(theory=ORCAobject, fragment=HF_frag)

The energy is printed to standard output by default.

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

The **Singlepoint** function is designed to be a simple function that does one job, returning 1 energy for the 1 theory level and the 1 fragment that was defined.
Usually, however, multiple calculations need to be performed. For example running the same single-point theory calculation on multiple fragments.

You could of course easily write a for-loop for this purpose, making sure to define first charge and multiplicity for each fragment first.
Here, defining the charge/mult for each fragment rather than the theory is desirable since the charge/mult is not always the same for all fragments.

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


More conveniently, however, you can instead use the Singlepoint_fragments function that does the same thing:

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

    Formula    Label       Charge    Mult           Energy(Eh)
    ----------------------------------------------------------------------
    N2         None             0       1        -6.3335016263
    H2         None             0       1        -1.0361629322
    N1H3       nh3              0       1        -4.8298958374

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

    Theory Label          Charge    Mult           Energy(Eh)
    ----------------------------------------------------------------------
    GFN1-xTB                   0       1        -6.3335016263
    GFN2-xTB                   0       1        -5.7639339581
    ORCA-r2SCAN-3c             0       1      -109.5070425194


##################################
Singlepoint_parallel function
##################################

The **Singlepoint_fragments** and **Singlepoint_theories** functions perform the calculations in a sequential fashion (via a for loop): i.e. one calculation after the other.
While convenient, the functions do not utilize the fact that each fragment-calculation (**Singlepoint_fragments**) or theory-calculation (**Singlepoint_theories**) is completely 
independent from each other and could thus run through the list of calculations (whether fragments or theories) in parallel on a multi-core CPU.
The **Singlepoint_parallel** function, however, allows you to do this.

See :doc:`parallelization` for information on using the **Singlepoint_parallel** function.