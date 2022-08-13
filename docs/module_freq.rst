Frequency calculation
======================================

A frequency calculation can be performed either numerically or analytically. The numerical approach is implemented directly in ASH and can always be applied to any QM, MM or QM/MM theory (requires the method to have an analytical gradient, however)
while the analytical approach is only available if the QM program contains the functionality and if the ASH interface to that theory supports it (see below)

For QM/MM Hamiltonians, the numerical approach is the only available option.

#########################################
Numerical frequencies
#########################################

.. code-block:: python

    def NumFreq(fragment=None, theory=None, npoint=2, displacement=0.005, hessatoms=None, numcores=1, runmode='serial', 
            temp=298.15, pressure=1.0, hessatoms_masses=None, printlevel=1, charge=None, mult=None):


**NumFreq** options:

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``fragment``
     - ASH Fragment
     - None
     - ASH Fragment object.
   * - ``theory``
     - ASHTheory object
     - None
     - ASH Theory object.
   * - ``runmode``
     - String
     - 'serial
     - | Whether to run the numerical frequency displacement calculations in serial (sequentially) 
       | or in parallel using Python multiprocessing.
   * - ``npoint``
     - Integer
     - 2
     - | Whether to do a 1-point or 2-point approximation for the Hessian. 1-point takes half as many
       | displacements as 2-point and is less accurate.
   * - ``displacement``
     - Float
     - 0.005
     - Size of displacement (in Å) to take.
   * - ``hessatoms``
     - list
     - None
     - | Optional list of atom indices that will be displaced and define a partial Hessian instead
       | of a full Hessian. If None, all atoms in Fragment will be displaced.
   * - ``temp``
     - Float
     - 298.15
     - Temperature (in K) to use in thermochemistry calculation. 
   * - ``pressure``
     - Float
     - 1.0
     - Pressure (in atm) to use in thermochemistry calculation.
   * - ``hessatoms_masses``
     - list
     - None
     - Optional list of masses for the Hessian atoms.
   * - ``charge``
     - integer
     - None
     - Optional charge. Will override charge attribute of ASH Fragment.
   * - ``mult``
     - integer
     - None
     - Optional spin multiplicity. Will override mult attribute of ASH Fragment.
   * - ``printlevel``
     - integer
     - 1
     - The printlevel.




Numerical frequencies can be performed with ASH using any QM, MM or QM/MM theory object.
Any method for which there is an analytical gradient (forces) available (in the external code) can be used (numerical 2nd derivative on top of numerical 1st derivative is not recommended).

Use the **NumFreq** function to request a numerical frequency job. The function requires a fragment object and a theory level at minimum.
The fragment object should typically contain a fragment with optimized coordinates at same level of theory (i.e. an already optimized minimum or saddlepoint).

*Type of Hessian*

Additionally you can select to do a 1-point Hessian or a 2-point Hessian by the *npoint* keyword (value of 1 or 2).
A 1-point Hessian makes a single displacement (+ direction) for each atom and each x,y and z-coordinate from the input geometry. This option is reasonably accurate and is the default.
A more accurate 2-point Hessian makes displacement in both + and - directions (for each x-, y- and z-coordinate of each atom), is twice as expensive (double the displacements)
but is more accurate.
The displacement step can be changed if required. The default recommended setting is: 0.0005 Å.

*Serial or parallel*

Two runmodes are available: 'serial' and 'parallel'. The 'serial' mode will run each displacement sequentially.
The Energy+Gradient step can still be run in parallel if e.g. the QM or QM/MM object has this information;
e.g. if an ORCA object has been defined with numcores=8 then ORCA will run each Energy+Gradient evaluation with 8 cores using the OpenMPI parallelization of ORCA.
For numerical frequencies, it is usually much more efficient, however, to run the displacement jobs simutaneously in parallel fashion.
This is accomplished using runmode='parallel' and the parallelization will be linear scaling (almost always recommended).
As there are almost always many more displacements available than CPUs, the parallelization of the QM or QM/MM object should usually be turned off and instead as many displacements
are run simultaneously as there are number of cores. For example, for a 30-atom system, there are 90 XYZ coordinates. For a 2-point Hessian, this means
that 180 displacements to be calculated. If 20 cores are available, then 20 displacements can be run simultaneously, fully utilizing all 20 cores.
This will require 9 runs in total (20*9=180).

.. warning:: Parallel runmode is currently not available for QM/MM calculations.

*Full or partial Hessian*

A partial Hessian can be easily performed instead of the full Hessian. This is an excellent approximation for vibrational modes with rather local character
and the quality of the approximation can be controlled. For a QM/MM model of a protein active site with an active region of a 1000 atoms, the full Hessian
of all 1000 atoms would typically not be a doable or recommended calculation; instead a partial Hessian job of the important atoms (e.g. the QM region) makes more sense.
A partial Hessian job is performed if a list of Hessian atoms (e.g. hessatoms=[0,1,2] ) is passed to the **NumFreq** function. In this case, the displacements
will only be calculated for the list of "hessatoms" and the result is a partial Hessian for the chosen atoms of the system

*Final output*

Once the displacements are complete, the gradients for all displacements are combined to give the full (or partial) Hessian.
The Hessian is then mass-weighted and diagonalized. (Limitation: translational and rotational modes are currently not projected out).
This gives the frequencies as eigenvalues and the normal mode eigenvectors.
A normal mode composition factor analysis is automatically performed as well as thermochemistry based on the rigid-rotor-harmonic-oscillator (RRHO) approximation.


**Examples:**

*Numerical frequencies in serial mode (QM-code parallelization instead used):*

.. code-block:: python

    from ash import *

    #the total number of CPU cores available to Ash (should match the job-script)
    numcores=8

    frag=Fragment(xyzfile="h2o.xyz", charge=0, mult=1)

    #ORCA theory object, ORCA parallelization turned off by not providing numcores keyword
    ORCAcalc = ORCATheory(orcasimpleinput="! r2SCAN-3c)

    #Serial Numfreq job (default):
    freqresult = NumFreq(fragment=Reactant, theory=ORCAcalc, npoint=2, runmode='serial')

    print("freqresult:", freqresult)


The resulting object from a NumFreq calculation is a dictionary (here called freqresult)
It contains the calculated frequencies and results from the Thermochemical analysis.
Individual items from the dictionary can be accessed by specifying the dictionary key:
Available keys: frequencies, ZPVE, vibenergy, transenergy, rotenergy, vibenergy, vibenergycorr

.. code-block:: python

  print("Frequencies : ", freqresult['frequencies'])
  print("ZPVE : ", freqresult['ZPVE'])



*Numerical frequencies in parallel mode (QM-code parallelization turned off):*

.. code-block:: python

    from ash import *

    #the total number of CPU cores available to Ash (should match the job-script)
    numcores=8

    frag=Fragment(xyzfile="h2o.xyz", charge=0, mult=1)

    #ORCA theory object, ORCA parallelization turned off by not providing numcores keyword
    ORCAcalc = ORCATheory(orcasimpleinput="! r2SCAN-3c)

    #Parallel mode: ASH will use the number of cores given to run same number of displacments simultaneously.
    freqresult = NumFreq(fragment=Reactant, theory=ORCAcalc, npoint=2, runmode='parallel', numcores=numcores)


#########################################
Analytical frequencies
#########################################

Some QM programs have analytical frequencies implemented and the ASH interface may support requesting the calculation of the analytical Hessian.
Currently analytical frequencies are supported in the QM codes: ORCATheory

An analytical Hessian calculation is requested by the AnFreq function that takes fragment and theory as necessary arguments:

.. code-block:: python

    def AnFreq(fragment=None, theory=None, numcores=1, temp=298.15, pressure=1.0)


Example:

.. code-block:: python

    HF_frag=Fragment(xyzfile="hf.xyz")
    ORCAcalc = ORCATheory(orcasimpleinput='BP86 def2-SVP def2/J', orcablocks="", numcores=1)
    thermochem_dict = AnFreq(theory=ORCAcalc, fragment=HF_frag)

    print("Thermochem properties dict:", thermochem_dict)
    print("Vibrational frequencies (cm**-1) : ", thermochem_dict['frequencies'])
    print("ZPVE (Eh) : ", thermochem_dict['ZPVE'])
    print("Gibbs energy corrections (Eh) : ", thermochem_dict['Gcorr'])


##############################################################################
thermochemistry corrections
##############################################################################

Thermochemistry corrections are automatically calculated when either a **Numfreq** or **Anfreq** job is requested.

.. code-block:: python

    thermochem_an = AnFreq(theory=ORCAcalc, fragment=HF_frag)
    thermochem_num = NumFreq(theory=ORCAcalc, fragment=HF_frag)

    print("Thermochem property dict:", thermochem_an)
    print("ZPVE (Eh) : ", thermochem_an['ZPVE'])

A dictionary containing various properties is returned (dictionary keys):
(frequencies, ZPVE, E_trans, E_rot, E_vib, E_tot, TS_trans, TS_rot, TS_vib, TS_el, vibenergycorr, Hcorr, Gcorr, TS_tot)

Alternatively, the thermochemcalc function can be called directly.

.. code-block:: python

    def thermochemcalc(vfreq,atoms,fragment, multiplicity, temp=298.18,pressure=1.0):

This function calculates the thermodynamic corrections from a list of available frequencies, number of atoms, ASH fragment object and spin multiplicity.
The temperature (default: 298.15 K) and pressure (default: 1.0 atm) can be specified.

.. code-block:: python

    h2o_frag = Fragment(xyzfile="h2o.xyz")
    #Manually defined frequencies for system
    frequencies=[1600.1, 2300.2, 2400.3]
    thermochemcalc(frequencies,3,h2o_frag, 1, temp=298.18, pressure=1.0)