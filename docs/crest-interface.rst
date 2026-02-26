CREST interface
======================================

ASH features an interface to the powerful conformational sampling program `crest <https://xtb-docs.readthedocs.io/en/latest/crest.html>`_ by the Grimme group.

The interface is evolving and currently contains 3 different ways of utilizing CREST together with ASH.

######################################################################################
new_call_crest  (general CREST interface, allowing any ASHTheory to be used)
######################################################################################

.. code-block:: python

    def new_call_crest(fragment=None, theory=None, crestdir=None, runtype="imtd-gc", 
                    energywindow=6.0, rthr=None, ethr=None, bthr=None,
                    shake=None, tstep=None, dump=None,length_ps=None,temp=None, hmass=None,
                    kpush=None, alpha=None, cvtype=None, dump_ps=None,
                    numcores=1, charge=None, mult=None, 
                    topocheck=True, constraints=None):

The **new_call_crest** function will call the CREST program to perform any CREST-runtype available, on a selected ASH fragment with an ASH theory.
The way the interface works is that *new_call_crest* will call CREST which will periodically run an ASH Python script to provide energies and gradients of the system.
In addition to ASH theories it is also possible to specify xtb-methods directly as valid string ('gfn1', 'gfn2', 'gfnff').

--------------------------
The CREST runtype options
--------------------------

- imtd-gc (the CREST MTD-based conformational sampler)
- nci-mtd (CREST sampling with a wall potential, NCI_MTD workflow)
- imtd-smtd (CREST sampling for calculating configurational entropy)
- mtd (Metadynamics as implemented in CREST)
- screen_ensemble
- ancopt_ensemble
- ancopt (the CREST ANC optimizer)
- numhess (numerical Hessian)
- imtd-gcimtd-gc

.. warning:: Not all of these runtypes have been tested 

See `crest documentation for details <https://crest-lab.github.io/crest-docs/page/documentation/inputfiles.html>`_

------------------------------------------------------------------------------
Example: CREST conformational sampling using the ORCATheory interface in ASH:
------------------------------------------------------------------------------

.. code-block:: python

    from ash import *
    frag = Fragment(xyzfile="molecule.xyz", charge=0, mult=1)
    theory  = ORCATheory(orcasimpleinput="! r2SCAN-3c tightscf")
    new_call_crest(fragment=frag, theory=theory, runtype="imtd-gc")


.. warning:: Unfortunately, while many other ASH theories will work for the example above, 
    the interface is currently limited to ASH Theory objects that can be serialized (pickled). Theory interfaces relying on Python libraries may not always work.

--------------------------
Modifying CREST parameters
--------------------------

See https://crest-lab.github.io/crest-docs/page/documentation/inputfiles.html#calculationconstraints-sub-blocks
for the CREST Inputfile documentation that ASH uses.

Currently some CREST options can be modified by the interface:

**CREGEN keywords**

- *ewin*
- *rthr*
- *ethr*
- *bthr

**Dynamics keywords**

- *shake*
- *tstep*
- *dump*
- *length_ps*
- *temp*
- *hmass*

**Metadynamics keywords**

- *kpush*
- *alpha*
- *cvtype*
- *dump_ps*


Topology check can be turned on/off  by the *topocheck* Boolean keyword.

--------------------------
Constraints
--------------------------
Finally constraints can also be defined by providing a dictionary of constraints via the *constraints* keyword.

A dictionary containing one or more of the following keys can be provided: 
'atoms', 'elements', 'distance', 'angle', 'dihedral', 'force', 'reference', 'metadyn_atoms'.
The corresponding values of the keys should be strings defining the constraints according to CREST's syntax.

*Example:*

.. code-block:: python

    from ash import *

    frag = Fragment(xyzfile="cowley_full.xyz", charge=0, mult=1)
    theory  = ORCATheory(orcasimpleinput="! r2SCAN-3c tightscf")
    # constraints
    constraints={'distance':'2, 66, 2.51630', 'bond':'1,2,auto'}
    new_call_crest(fragment=frag, theory=theory, runtype="imtd-gc", constraints=constraints)

.. warning:: The constraints defined here follow the CREST 1-based atom indexing in contrast to the general 0-based indexing that ASH uses. 

See CREST documentation for examples of constraints that can be applied:

- https://crest-lab.github.io/crest-docs/page/documentation/inputfiles_examples.html#constrained-geometry-optimization
- https://crest-lab.github.io/crest-docs/page/examples/qcg/example_3.html#constraining-the-solute
- https://crest-lab.github.io/crest-docs/page/examples/example_4.html


################################################################################
Using CREST to call ASH as a generic theory
################################################################################

It is also possible to use CREST as a driver and provide an ASH script as a theory-level to CREST. This is the generic option available in CREST.
This option has the advantage of being more flexible than the option above as it should in principle allow any ASH level of theory to be used within CREST.
Additionally all CREST options are more easily specified.

This option works like the following.  
CREST should be called by providing a CREST inputfile in TOML format (here called input.toml)

.. code-block:: shell

    crest --input input.toml

The input.toml file should look e.g. like this:

.. code-block:: text

    # CREST 3 input file
    input = "struc.xyz"
    runtype="ancopt"
    threads = 1

    [calculation]
    elog="energies.log"

    [[calculation.level]]
    method = "generic"
    binary = "python3 ../ash_input.py"
    gradfile = "genericinp.engrad"
    gradtype = "engrad"
    uhf = 0
    chrg = 0

where the system coordinates are provided in the form of an XYZ-file (struc.xyz), the runtype is chosen (here ancopt) as well as a few other options.
See `crest documentation <https://crest-lab.github.io/crest-docs/page/documentation/inputfiles.html>`_ for details.

The calculation level is next chosen to be "generic" and the syntax for running an ASH Python script is provided (don't modify).

The ASH script named *ash_input.py* should be created in the same directory as *input.toml* and should look like this:

.. code-block:: python

    from ash import *

    #Charge/mult settings
    charge=0
    mult=1
    #Definition of the ASH Theory that you want
    theory = ORCATheory(orcasimpleinput="! r2SCAN-3c tightscf")

    ###############################
    # No changing anything below !
    ###############################
    #ASH creation of fragment for CREST-generated XYZ-file (genericinp.xyz) in each CREST-step
    frag = Fragment(xyzfile="genericinp.xyz", charge=charge,mult=mult)
    #Singlepoint Energy+Gradient calculation
    result = Singlepoint(theory=theory, fragment=frag, Grad=True)
    #Print energy and gradient in the form of the ORCA-formatted engrad file (that CREST reads)
    print_gradient_in_ORCAformat(result.energy,result.gradient,"genericinp", extrabasename="")

This script is essentially just an ASH script for running a Singlepoint Energy+gradient calculation using CREST-created input coordinates (will be created/updated in each step in file genericinp.xyz)
and then writing the Energy and Gradient into a specifically formatted file (genericinp.engrad, note: in ORCA-format).

When CREST is running it will in each step create a new geometry, run the ASH script above and will then read the energy and gradient for each new geometry.

The advantage of this option is that the file *ash_input.py* can contain any valid ASH-level of theory (including hybrid theories)
and can in principle be customized if required.


################################################################################
call_crest  (simple xtb-based conformation sampling)
################################################################################

By providing an ASH fragment object to the **call_crest** function, the CREST xTB-metadynamics-based conformational sampling procedure is invoked.
This function can only perform xTB-based conformational sampling.

The output from crest is written to standard output. If successful, Crest will create a file crest_conformers.xyz
that can be directly read into ASH for further processing or further calculations.
This allows one to write a multi-step workflow of which the crest-procedure is one of many steps.

.. code-block:: python

    #Function to call crest
    def call_crest(fragment=None, xtbmethod=None, crestdir=None,charge=None, mult=None, solvent=None, energywindow=6, numcores=1,
                   constrained_atoms=None, forceconstant_constraint=0.5)
    #Function to grab conformers. Returns list of conformers and list of xtb energies
    def get_crest_conformers()


**call_crest** requires one to specify: an ASH fragment, xtbmethod (GFN1-xTB or GFN2-xTB ), location of crest directory, charge, multiplicity.

Optional keywords are: solvent, energywindow (default 6), numcores (default 1), constrained_atoms (list of integers) and the value of the force-consstant.

If you specify a list of constrained atoms then ASH will create an .xcontrol file that defines the constraints according to `crest Example Applications <https://xtb-docs.readthedocs.io/en/latest/crestxmpl.html>`_.




-----------------------------------------------------------------------------------
Example workflow 1. Call crest to get low-energy conformers as ASH fragments.
-----------------------------------------------------------------------------------
.. code-block:: python

    from ash import *

    crestdir='/opt/crest'
    numcores=24

    #0. Starting structure and charge and mult
    molecule = Fragment(xyzfile="ethanol.xyz", charge=0, mult=1)

    #1. Calling crest and getting list of conformer fragments and energies
    list_conformer_frags, xtb_energies = call_crest(fragment=molecule, xtbmethod='GFN2-xTB', crestdir=crestdir, 
        numcores=numcores)

    print("list_conformer_frags:", list_conformer_frags)
    print("")
    print("Crest Conformer Searches done. Found {} conformers".format(len(xtb_energies)))
    print("xTB energies: ", xtb_energies)


-----------------------------------------------------------------------------------
confsampler_protocol : Automatic Crest+DFTopt+DLPNO-CCSD(T) workflow
-----------------------------------------------------------------------------------

It is also possible to call the **confsampler_protocol** function that carries out an automatic multi-step workflow
at various levels of theory.

1. conformational sampling using crest and GFN-xTB (**low-level** theory).
2. Geometry optimizations for each low-energy conformer at a **medium-level** of theory (typically DFT using e.g. ORCATheory)
3. **High-level** single-point calculation (e.g. DLPNO-CCSD(T)/CBS using e.g. ORCA_CC_CBS_Theory)

.. code-block:: python

    def confsampler_protocol(fragment=None, crestdir=None, xtbmethod='GFN2-xTB', MLtheory=None, 
                            HLtheory=None, numcores=1, charge=None, mult=None, crestoptions=None,
                            optimizer_maxiter=200):


.. code-block:: python

    from ash import *

    #
    crestdir='/opt/crest'
    numcores=4
    #Fragment to define
    frag=Fragment(xyzfile="ethanol.xyz", charge=0, mult=1)

    #Defining MLTheory: DFT optimization
    MLsimpleinput="! B3LYP D3BJ def2-TZVP TightSCF "
    MLblockinput="""
    %scf maxiter 200 end
    """
    ML_B3LYP = ORCATheory(orcasimpleinput=MLsimpleinput, orcablocks=MLblockinput, numcores=numcores)
    #Defining HLTheory: DLPNO-CCSD(T)/CBS
    HL_CC = ORCA_CC_CBS_Theory(elements=frag.elems, cardinals = [2,3], basisfamily="def2", DLPNO=True, 
        pnosetting='extrapolation', pnoextrapolation=[1e-6,3.33e-7,2.38,'NormalPNO'], numcores=numcores)

    #Call confsampler_protocol
    confsampler_protocol(fragment=frag, crestdir=crestdir, xtbmethod='GFN2-xTB', MLtheory=ML_B3LYP,
                             HLtheory=HL_CC, orcadir=orcadir, numcores=numcores)

Final result table of calculated conformers at 3 different theory levels:

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