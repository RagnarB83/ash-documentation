Tutorial: Accelerating computational chemistry with ML potentials
=========================================================================

Machine-learning potentials are currently an area of much interest where a lot of progress is being made.
These (mostly) neural-network based models aim to essentially bypass the electronic structure problem of 
quantum chemistry, by predicting energy as a function of nuclear coordinates.
They can be at least 100-1000 times faster than DFT methods.
When trained directly by the user on a specific system they can essentially achieve arbitrary accuracy, depending strongly on the nature 
and amount of the input training data.
The aim of pre-trained foundational potentials, however, is to be used as black-box replacements for QM or MM methods.
The most recent ones have been trained on large data sets (100 million molecules or more) and are available for most elements of the periodic table.

While there is a lot buzz around these models for their cost-accuracy ratio, these models should not be trusted blindly. 
The models will not perform well on exotic systems with complicated electronic structure (charged and/or open-shell systems are still problematic).
and generally the models will fail the further the system is from the training set. 
Ultimately, being trained on DFT energy surfaces they can never surpass DFT either.

This tutorial will focus on how these pre-trained ML potentials can be used in a responsible way in ASH, serving the purpose
of accelerating computational chemistry jobs while not relying blindly on their accuracy.


For more general information on machine learning in ASH, please refer to :doc:`Machine_learning_in_ASH`
See also :doc:`Tutorial_delta_learning` .

#######################################################################################
1. Accelerating DFT-based geometry optimizations by ML pre-optimizations
#######################################################################################

Geometry optimizations using dispersion-corrected DFT methods in a triple-zeta basis set
can be considered close to the quantum chemistry state of the art. Such protocols 
come remarkably close to correlated wavefunction accuracy for geometries, though not for energies.
For large systems, and especially open-shell systems, DFT geometry optimizations can, however, be time-consuming, lasting for days or weeks even,
especially if started far away from the minimum. 
In a computational chemistry study where many geometry optimizations have to be performed in order to reliably 
explore the potential energy surface, this becomes time-consuming.

By pre-optimizing geometries with cheaper methods, one can reduce the number of expensive DFT optimization steps
required. This strategy is particularly useful for very large systems when far from equilibrium.
While semiempirical methods like GFN2-xTB often serve the purpose for such pre-optimizations,
new foundational ML potentials offer another alternative. 
These ML potentials have usually been trained to reproduce a certain dispersion-corrected DFT level of theory.
Models trained on the OMol25 training set contain e.g. 100 million molecular structures and associated energies and gradients at the wB97M-V/def2-TZVPD level of theory.
wB97M-V is a range-separate hybrid functional that together with the non-local VV10 dispersion functional and a close-to basis-set limit def2-TZVPD basis set,
is one of the most accurate DFT methods for maingroup and transition metal chemistry.

Below we show how an open-shell X-atom system can be pre-optimized with a MACE-OMOL-0 potential (a MACE neural-network potential trained on OMol25 data)
 before switching to the desired DFT protocol, in one ASH script.

.. code-block:: python

    from ash import *

    # Read in molecule XYZ into ASH Fragment
    frag = Fragment(xyzfile="largemolecule.xyz", charge=0, mult=1)

    # Define the ML-theory to be used for pre-optimization
    mltheory = MACETheory(model_name="mace_omol")
    # Define the DFT theory to be used for final optimization
    dfttheory = ORCATheory(orcasimpleinput="! PBE0 D4 def2-TZVP tightscf")

    # Pre-optimization with ML-theory
    Optimizer(theory=mltheory, fragment=frag)

    # Final optimization with the DFT theory 
    # note: frag object is automatically updated after each optimization
    Optimizer(theory=dfttheory, fragment=frag)

This is maybe the simplest way of reliably using ML potentials in computational chemistry without 
relying on the unknown accuracy of the ML potentials and potential for catastrophic failure.
The final geometry and energy of the system is still obtained at the DFT-level, the process is simply accelerated
by pre-optimizing using the ML potential.
Even if the ML potential happens to not be accurate enough for the system, it still may save some time.
Furthermore if the user finds by testing that the ML potential is generally accurate enough for predicting reliable structures,
another option is also to replace the final DFT optimization for a DFT single-point energy evaluation instead.

#######################################################################################
2. Accelerating saddlepoint searches and minimum energy paths
#######################################################################################

Saddlepoint (SP) searches can be one of the most time-consuming part of computational chemistry
An accurate saddle-point optimization needs to be performed by an eigenvector-following algorithm (usually a PRFO algorithm)
but requires as input: 

- i) a reasonably accurate guess for the SP geometry 
- ii) information about the curvature (Hessian information)

The process of obtaining the SP-geometry guess can often be time-consuming and laborious, requiring either
demanding minimum energy path calculations (such as NEB, growing/freezing string), time-consuming and/or user-demanding 
relaxed surface scans or considerable manual effort in setting up the SP-geometry guess in some other way.
Similarly, the computation of the Hessian can be highly demanding, especially if the system is large.

But such computations are actually primarily demanding because of the use of DFT methods in both the SP-geometry search
and in the computation of the Hessian. 
If one is able to replace the DFT-theory for a much cheaper (and sufficiently accurate) 
ML-potential theory, saddlepoint searches can become considerably cheaper.

Here we show several examples of how we can accelerate SP searches.

-------------------------------------------------------
Accelerating SP-search via relaxed surface scan
-------------------------------------------------------

While relaxed scans are not the most reliable way of finding saddlepoints,
they nonethess often work and additionally have the advantage of sometimes finding new minima.

In the script below we perform the 1D surface scan via a ML-potential by gradually changing a geometric variable,
serving the purpose of an approximate reaction coordinate. We then locate the highest energy point and use this 
geometry to start an eigenvector-following SP-search at a DFT-level of theory. 
Instead of computing a DFT Hessian we compute a cheaper Hessian at the ML-potential level of theory.

.. code-block:: python

    from ash import *

    # Read in molecule XYZ into ASH Fragment
    frag = Fragment(xyzfile="butane_CH10.xyz", charge=0, mult=1)

    # Define the ML-theory to be used for scan
    mltheory = MACETheory(model_name="mace_omol")
    # Define the DFT theory to be used for final optimization
    dfttheory = ORCATheory(orcasimpleinput="! PBE0 D4 def2-TZVP tightscf")

    # Run scan using ML
    result_scan = calc_surface(fragment=frag, theory=mltheory, scantype='relaxed',
        RC_list=[{'type': 'dihedral',  'indices': [[0,1,2,3]], 'range': [-180, 180, 10]}])
    surfacedictionary = result_scan.surfacepoints
    #Plot scan result
    reactionprofile_plot(surfacedictionary, finalunit="kcal/mol")

    # Find highest point of scan automatically by searching dictionary
    max_point = max(surfacedictionary, key=surfacedictionary.get)[0]
    print("max_point:", max_point)

    # Create SP geometry by reading the XYZ-file for the highest energy point
    SPguess = Fragment(xyzfile=f"surface_xyzfiles/RC1_{max_point}.xyz", charge=frag.charge, mult=frag.mult)
    print("SPguess:", SPguess)

    # Calculate numerical Hessian using MLtheory
    result_numfreq = NumFreq(theory=mltheory, fragment=SPguess)

    # Perform TS-Opt job using ML-calculated SP-guess and Hessian but at DFT-theory
    Optimizer(theory=dfttheory, fragment=SPguess, hessian=result_numfreq.hessian)

    # Perform final DFT NumFreq/Anfreq calculation (this will be most expensive step)
    NumFreq(theory=dfttheory, fragment=SPguess)

-------------------------------------------------
Accelerating SP-search via NEB
-------------------------------------------------

A Nudged elastic band (NEB) calculation minimizes a band of interconnected structures (images), initially obtained by interpolation between a reactant
and product structures.
An NEB calculation has the advantages of actually giving a minimum energy path between structures, of which the highest energy point
will correspond closely to the actual saddlepoint. The climbing-image variant of NEB (CI-NEB) is the most commonly employed NEB variant,
 which is often capable even of accurately locating the saddlepoint. The use of energy-weighted springs in the algorithm implemented in ASH (and in ORCA),
 is particularly useful for molecular reactions. By increasing the number of images, the path becomes more and more accurate
 and can describe even a complex reaction process with multiple intermediates and saddlepoints. Because the images are independent during each step,
 they can be run in parallel.
 
The disadvantage of NEB for the purpose of locating a single saddlepoint,
 is the need to compute many structures (images, 6-10 is a common choice) in each optimization step,
and overall slow convergence to the minimum energy path (so slow that one only rarely chooses to optimize a 
path to a tolerance where all images are converged).
Another disadvantage comes from the fact that because the algorithm only uses gradient-information, it does not always precisely locate the 1st-order saddlepoint 
 (characterized by 1 negative eigenvalue of the Hessian), even with CI-NEB.

Combining NEB with an eigenvector-following algorithm, as in the NEB-TS method,
provides a reliable way of both locating the area of the saddlepoint (and a partially converged minimum energy path), while also
precisely locating the saddlepoint.
This is accomplished by only partially converging a CI-NEB calculation, then using the highest energy climbing image as a SP-guess,
computing the Hessian at this geometry and then optimizing to the saddlepoint via an eigenvector-following algorithm.

It is possible to use ML potentials to accelerate these steps. If the ML potential is accurate enough it will be capable of predicting a minimum
energy path close to that predicted by DFT but at a much cheaper cost. We can modify the tolerances of the NEB job as we please:
the lower they are the better the SP-geometry is likely to be but they we do not have to worry about converging the images very tightly.
Because the ML potential is cheap we can also choose to increase the number of images which will both aid convergence and also provide more resolution 
in the saddlepoint region which will result in a more accurate guess.
Once the NEB job is over, we can take the highest energy structure and compute the Hessian at the ML-potential level of theory.

We now have all the ingredients for locating the saddlepoint precisely
and we can switch to the DFT-level of theory using the ML-located SP-guess and ML-computed Hessian.



NEB using a ML-potential:

.. code-block:: python

    from ash import *

    R = Fragment(xyzfile="reactant.xyz", charge=1, mult=1)
    P = Fragment(xyzfile="product.xyz", charge=1, mult=1)

    #Defining ML and DFT theory
    mltheory = MACETheory(model_name="mace_omol")
    dfttheory = ORCATheory(orcasimpleinput="! wB97X-V def2-TZVP tightscf")

    #Performing NEB using a decent number of images. Tolerances are increased by 10 to avoid convergence failure.
    NEB(theory=theory, reactant=R, product=P, images=12, tol_max_fci=0.26, tol_rms_fci=0.13)
    
    # Read in SP-geoemtry to new Fragment
    SPguess = Fragment(xyzfile="Saddlepoint-NEBCIapprox.xyz", charge=1, mult=1)

    # Compute Hessian on SP-geometry using ML-theory
    result_numfreq = NumFreq(theory=mltheory, fragment=SPguess)

    # OPT-TS using DFT-theory but using ML-located SP-guess and ML-computed Hessian
    Optimizer(theory=dfttheory, fragment=SPguess, TSOpt=True, hessian=result_numfreq.hessian)


#######################################################################################
3. Accelerating conformational sampling
#######################################################################################

Conformational sampling in computational chemistry is an important thing to consider for molecular systems
that feature a lot of flexibility. Due to the computational cost of DFT it is very difficult to perform conformational sampling at the DFT-level of theory.
Quantum-chemistry based conformational sampling has been popularized in recent years by the metadynamics-based CREST algorithm
by Grimme and coworkers that pairs well with the cheap semiempirical xTB methods from the same group.
CREST sampling jobs, however, are actually agnostic to the level of theory and a ML potential could be used instead of xTB for conformational sampling.
ASH offers a convenient way of performing CREST jobs at any level of theory (that is definable in ASH). See :doc:`crest-interface`.

The script below shows how one can combine CREST with an MACE-OMOL25 based ML potential to locate the lowest energy conformers,
prior to performing DFT-based geometry optimizations for the 

.. code-block:: python

    from ash import *

    molecule = Fragment(xyzfile="molecule.xyz", charge=0, mult=1)

    # Defining ML and DFT theory
    mltheory = MACETheory(model_name="mace_omol")
    dfttheory = ORCATheory(orcasimpleinput="! r2SCAN-3c tightscf")

    # Performing CREST imtd-gc conformational sampling job using the ML-theory
    list_conformers, list_energies = new_call_crest(fragment=molecule, theory=mltheory, runtype="imtd-gc")

    # Run single-point DFT calculations of each ML-located CREST conformer
    # Here we use Singlepoint_fragments for this purpose
    energies = Singlepoint_fragments(theory=xtbcalc, fragments=list_conformers)

    # Or run DFT geometry optimizations of each ML-located CREST conformer
    energies_dict={}
    for i,conf in enumerate(list_conformers):
        print(f"Running conformer {i} of {len(list_conformers})")
        optresult = Optimizer(theory=dfttheory, fragment=conf)
        os.rename("Fragment-optimized.xyz", f"conformer_{i}.xyz")
        # put in dict
        energies_dict[f"conformer_{i}"] = optresult.energy

    print(" Conformer   Energy (Eh)")
    print("----------------------------")
    for conf, e in energies_dict.items():
        print("{:10} {:13.10f}".format(conf,e))


#######################################################################################
4. Accelerating by finetuning or training bespoke potentials
#######################################################################################

TODO...