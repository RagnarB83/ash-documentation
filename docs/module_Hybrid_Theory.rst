Hybrid Theory
==========================

Hybrid Theories are multi-level ASH theories that combine multiple theory level objects to give some kind of combined theory description of the system.
The different theories might be used for different parts of the system (here called **Multilevel Theory methods** , e.g. QM/MM) or on the same part of the system (examples are **WrapTheory** and **DualTheory methods**)



######################################################
Multilevel Theory Methods
######################################################

The **QMMMTheory** class in ASH (:doc:`module_QM-MM`) is considered a special type of a multilevel method where a QMTheory is combined with an MMTheory to describe different parts of the system.
QM/MM is typically used when the system can naturally be divided up into a local important part (described by a QM method) and an environment part (described by a classical MM method).

ASH will at some point feature also ONIOM multilevel methods that are similar to QM/MM but have a different philosophy.


######################################################
WrapTheory methods
######################################################

**WrapTheory** is an ASH Theory class that wraps two different theory objects to get a combined theory description. It is assumed that the energy and gradient from each theory is completely additive.

.. code-block:: python

    class WrapTheory:
        """ASH WrapTheory theory.
        Combines 2 theories to give a modified energy and modified gradient
        """
        def __init__(self, theory1=None, theory2=None, printlevel=1, label=None):

        def run(self, current_coords=None, current_MM_coords=None, MMcharges=None, qm_elems=None, mm_elems=None,
        elems=None, Grad=False, PC=False, numcores=None, restart=False, label=None,
        charge=None, mult=None):


**WrapTheory** was created for the purpose of allowing one to combine a regular theory-level with a correction (both energy and gradient) from another source.
Specifically, it was created to allow one to easily add dispersion corrections using DFTD4Theory to a regular DFT calculation (without dispersion).

Example:

.. code-block:: python

    from ash import *

    #Glycine fragment from database
    frag = Fragment(databasefile="glycine.xyz")

    #PBE/def2-SVP via ORCA (no dispersion correction)
    orca = ORCATheory(orcasimpleinput="! PBE def2-SVP tightscf")
    #DFTD4 dispersion correction using DFTD4 library
    dftd4 = DFTD4Theory(functional="PBE")
    #Combining the two theories using WrapTheory
    dft_plus_dftd4_theory = WrapTheory(theory1=orca, theory2=dftd4)

    #Calling the Optimizer function using the WrapTheory object as theory 
    Optimizer(theory=dft_plus_dftd4_theory, fragment=frag)


WrapTheory could be used for many other purposes, one would simply have to make sure that the 2 theories are compatible and that the sum of the 2 theory-description does not result in double-counting of any similar physical energy terms.
A regular DFT calculation (barely describes dispersion) + an atom pairwise dispersion correction (DFT-D4) is a good example of this.

See DFTD4 section in :doc:`helper_programs` for more information on the DFTD4Theory object.


######################################################
DualTheory methods
######################################################

**DualTheory** is an experimental ASH Theory that combines two different theory objects, e.g. a low-level QM theory and a high-level QM theory in a specific way in order to speed up an otherwise expensive high-level calculation.
This only makes sense for an expensive multi-iteration job where the Theory object is called multiple times, e.g. a geometry optimization or NEB calculation (not a single-point calculation).

The idea is to approximate the accurate high-level potential energy surface description by a low-level potential energy surface desciption + a correction derived from the high-level theory.
If the correction is calculated in every step (of e.g. a geometry optimization) there is no advantage (in fact more expensive) to using a **DualTheory** description.
However, if the high-level correction is only occasionally calculated then it possible to cut down on the number of expensive high-level energy+gradient calculations required.

Both energy and the gradient (required for optimizations and NEB calculations) can be corrected.

Currently the only available correction option is: "Difference" which features a naive energy/gradient difference correction.
The update_freq keyword controls the interval between corrections.
To use a Dualtheory one needs to give valid ASH Theory objects to the theory1 and theory2 keywords where theory1 is assumed to be the low-level theory (called each time) while theory2 is the high-level theory (
called only when the high-level correction should be updated according to the value of *update_freq*).

.. code-block:: python

    class DualTheory:
        """ASH DualTheory theory.
        Combines two theory levels to give a modified energy and modified gradient
        """
        def __init__(self, theory1=None, theory2=None, printlevel=2, label=None, correctiontype="Difference", update_freq=5, numcores=1):




----------------------------------------------------------------------
Geometry optimization example using GFN1-xTB and DFT:
----------------------------------------------------------------------
.. code-block:: python

    from ash import *

    numcores=1
    frag=Fragment(xyzfile="react.xyz", charge=0, mult=1)

    #Defining theory levels
    xtb = xTBTheory(xtbmethod="GFN1", numcores=numcores)
    orca = ORCATheory(orcasimpleinput="!r2scan-3c tightscf CPCM", numcores=numcores)

    #Creating DualTheory object: 
    #theory1 is the cheaper low-level theory called in each step, theory2 is the less-called high-level theory
    dualcalc = DualTheory(theory1=xtb, theory2=orca, update_freq=15)

    #Calling the Optimizer function using the DualTheory object
    geomeTRICOptimizer(theory=dualcalc, fragment=frag, maxiter=250)


----------------------------------------------------------------------
A nudged elastic band job example using GFN1-xTB and DFT:
----------------------------------------------------------------------

.. code-block:: python

    from ash import *

    numcores=1

    #Fragment for an SN2 reaction
    Reactant=Fragment(xyzfile="react.xyz", charge=-1, mult=1)
    Product=Fragment(xyzfile="prod.xyz",charge=-1, mult=1)

    #Defining individual theory levels
    xtb = xTBTheory(numcores=numcores)
    orca = ORCATheory(orcasimpleinput="!r2scan-3c tightscf CPCM", numcores=numcores)

    #Creating DualTheory object: 
    #theory1 is the cheaper low-level theory called in each step, theory2 is the less-called high-level theory
    dualcalc = DualTheory(theory1=xtb, theory2=orca, update_freq=5)

    #Calling the NEB job function using the DualTheory object
    NEB(reactant=Reactant, product=Product, theory=dualcalc, images=12, printlevel=0, maxiter=200)
