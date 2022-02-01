Frequency calculation
======================================

A frequency calculation can be performed either numerically or analytically. The numerical approach is highly general and can always be applied to any QM, MM or QM/MM theory while the analytical approach
is only available if the QM program contains the functionality.
For QM/MM Hamiltonians, the numerical approach is currently the only available option.

#########################################
Numerical frequencies
#########################################

Numerical frequencies can be performed with Ash using any QM, MM or QM/MM theory object.
Any method for which there is an analytical gradient (forces) available can be used (numerical 2nd derivative on top of numerical 1st derivative is not recommended).

Use the **NumFreq** function to request a numerical frequency job. The function requires a fragment object and a theory level at minimum.
The fragment object should typically contain a fragment with optimized coordinates at same level of theory (i.e. an already optimized minimum or saddlepoint).

*Type of Hessian*
Additionally you can select to do a 1-point Hessian or a 2-point Hessian by the *npoint* keyword (value of 1 or 2).
A 1-point Hessian makes a single displacement (+ direction) for each atom and each x,y and z-coordinate from the input geometry. This option is reasonably accurate and is the default.
A more accurate 2-point Hessian makes displacement in both + and - directions (for each x-, y- and z-coordinate of each atom), is twice as expensive (double the displacements)
but is more accurate.
The displacement step can be chosen if wanted. The default setting is: 0.0005 Å.

*Serial or parallel*
Two runmodes are available: 'serial' and 'parallel'. The 'serial' mode will run each displacement sequentially.
The Energy+Gradient step can still be run in parallel if e.g. the QM or QM/MM object has this information;
e.g. if an ORCA object has been defined with numcores=8 then ORCA will run each Energy+Gradient evaluation with 8 cores using the OpenMPI parallelization of ORCA.
For numerical frequencies, it is usually much more efficient, however, to run the displacement jobs simutaneously in parallel fashion.
This is accomplished using runmode='parallel' and the parallelization will be linear scaling (almost always recommended).
As there are almost always many more displacements available than CPUs, the parallelization of the QM or QM/MM object is turned off and instead as many displacements
are run simultaneously as there are number of cores. For example, for a 30-atom system, there are 90 XYZ coordinates. For a 2-point Hessian, this means
that 180 displacements to be calculated. If 20 cores are available, then 20 displacements can be run simultaneously, fully utilizing all 20 cores.
This will require 9 runs in total (20*9=180).

*Full or partial Hessian*

A partial Hessian (NEEDS TO BE TESTED) can be easily performed instead of the full Hessian. This is an excellent approximation for vibrational modes with rather local character
and the quality of the approximation can be controlled. For a QM/MM model of a protein active site with an active region of a 1000 atoms, the full Hessian
of all 1000 atoms would typically not be doable; instead a partial Hessian job of the important atoms (e.g. the QM region) makes more sense.
A partial Hessian job is performed if a list of Hessian atoms (e.g. hessatoms=[0,1,2] ) is passed to the NumFreq function. In this case, the displacements
will only be calculated for the list of "hessatoms" and the result is a partial Hessian for the system.

*Final output*
Once the displacements are complete, the gradients for all displacements are combined to give the full (or partial) Hessian.
The Hessian is then mass-weighted and diagonalized. (Limitation: translational and rotational modes are currently not projected out).
This gives the frequencies as eigenvalues and the normal mode eigenvectors.
A normal mode composition factor analysis is automatically performed (NOT READY) as well as zero-point energy thermochemistry.


Example script below demonstrates a combined geometry optimization (using geomeTRIC).
The QM code used here is ORCA but any QM, MM or QM/MM object can be used.

.. code-block:: python

    from ash import *

    #the total number of CPU cores available to Ash (should match the job-script)
    numcores=8

    orcasimpleinput="! HF-3c "
    orcablocks="%scf maxiter 200 end"

    reactstring="""
       C  -2.66064921   -0.44148342    0.02830018
       H  -2.26377685   -1.23173358    0.68710920
       H  -2.29485851   -0.62084858   -0.99570465
       H  -2.27350346    0.53131334    0.37379014
       F  -4.03235214   -0.44462811    0.05296388
    """
    Reactant=Fragment(coordsstring=reactstring, charge=0, mult=1)

    #Calculator object without frag. numcores=8 is used here for parallelizing ORCA during optimization.
    ORCAcalc = ORCATheory(orcasimpleinput=orcasimpleinput, orcablocks=orcablocks, numcores=numcores)

    #Geometry optimization of Reactant object and ORCAcalc theory object.
    #Each Energy+Grad step is parallelized by ORCA.
    geomeTRICOptimizer(theory=ORCAcalc,fragment=Reactant)


    #Numfreq job. A 1-point or 2-point Hessian can be requested.
    # Either serial or parallell runmode can be used.
    # For parallel: Ash will use the number of cores given to run same number of displacments simultaneouslyu.
    #ORCA parallelization is turned off automatically.

    #Serial mode:
    freqresult = NumFreq(fragment=Reactant, theory=ORCAcalc, npoint=2, runmode='serial')
    #Parallel mode:
    freqresult = NumFreq(fragment=Reactant, theory=ORCAcalc, npoint=2, runmode='parallel', numcores=numcores)

    print("freqresult:", freqresult)
    #The resulting object from a NumFreq calculation is a dictionary (here called freqresult)
    # It contains the calculated frequencies and results from the Thermochemical analysis.
    #Individual items from the dictionary can be accessed by specifying the dictionary key:
    # Available keys: frequencies, ZPVE, vibenergy, transenergy, rotenergy, vibenergy, vibenergycorr
    # TO BE FINISHED...
    print("Frequencies : ", freqresult['frequencies])
    print("ZPVE : ", freqresult['ZPVE])



#########################################
Analytical frequencies
#########################################

Some QM programs have analytical frequencies implemented and the ASH interface may support requesting the calculation of the analytical Hessian.

Currently analytical frequencies are supported in: ORCATheory

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

A dictionary containing various properties is returned (dictionary keys) from an AnFreq job:
(frequencies, ZPVE, E_trans, E_rot, E_vib, E_tot, TS_trans, TS_rot, TS_vib, TS_el, vibenergycorr, Hcorr, Gcorr, TS_tot)