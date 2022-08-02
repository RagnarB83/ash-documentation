Highlevel workflows
======================================

These high-level workflows are multi-step singlepoint energy protocols can either be used on their own as a theory-level in Singlepoint calculations or used as a SP_theory in workflows such as **thermochemprotocol**, **calc_xyzfiles**, **confsampler_protocol** (see :doc:`module_workflows`) 
or as a theory in **run_benchmark** (see :doc:`module_benchmarking`) .
The ORCA_CC_CBS_Theory uses the ORCA quantum chemistry code for all steps of the workflows and gives a final 0 K electronic energy (no ZPVE). Gradients are not available and these can thus not be used in geometry optimizations or dynamics jobs.
The MRCC_CC_CBS_Theory uses the MRCC program (not yet available)


#########################################
ORCA_CC_CBS_Theory
#########################################

ORCA_CC_CBS_Theory is synonymous with CC_CBS_Theory.

This is an ASH Theory that carries out a multi-step single-point protocol to give a CCSD(T)/CBS estimated energy.
Multiple ORCA calculations for the given geometry are carried out and the SCF and correlation energies extrapolated to the CCSD(T)/CBS limit using either regular CCSD(T) theory or DLPNO-CCSD(T) theory.
This workflow is flexible and features multiple ways of approaching the complete basis set limit (CBS) or the complete PNO space limit (CPS).
Various options affecting the accuracy, efficiency and robustness of the protocol can be chosen.
Many basis set families can be chosen that are available for most of the periodic table.
Atomic spin-orbit coupling can be automatically included if system is an atom.

.. code-block:: python

  class ORCA_CC_CBS_Theory:
      def __init__(self, elements=None, scfsetting='TightSCF', extrainputkeyword='', extrablocks='', 
              guessmode='Cmatrix', memory=5000, numcores=1, 
              cardinals=None, basisfamily=None, SCFextrapolation=True, alpha=None, beta=None, 
              stabilityanalysis=False, CVSR=False, CVbasis="W1-mtsmall", F12=False, Openshellreference=None, 
              DFTreference=None, DFT_RI=False, auxbasis="autoaux-max",
              DLPNO=False, pnosetting='NormalPNO', pnoextrapolation=[1e-6,1e-7,1.5,'TightPNO'], FullLMP2Guess=False, 
              T1=False, T1correction=False, T1corrbasis_size='Large', T1corrpnosetting='NormalPNOreduced', 
              relativity=None, orcadir=None, FCI=False, atomicSOcorrection=False):

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``elements``
     - list of strings.
     - None
     - Required: List of all elements of the molecular system (reaction). Needed to set up basis set information. Duplicates are OK. fragment.elems is a valid list.
   * - ``cardinals``
     - list of integers
     - None
     - Required: List of cardinal numbers for basis-set extrapolation. Options: [2,3], [3,4], [4,5] or [5,6]. Single-item lists also valid: e.g. [4] (for a single QZ level calculation).
   * - ``basisfamily``
     - string
     - None
     - Required: Name of basis-set family to use. Various options. See table below. 
   * - ``relativity``
     - string
     - None
     - Scalar relativity treatment. Options: 'DKH', 'ZORA', 'NoRel', None. 
   * - ``orcadir``
     - string
     - None
     - Path to ORCA. Optional. 
   * - ``stabilityanalysis``
     - Boolean
     - False
     - Perform SCF stability analysis for each SCF calculation.
   * - ``numcores``
     - integer
     - 1
     - Number of cores to use in ORCA calculation.
   * - ``Openshellreference``
     - string
     - None
     - Use alternative reference WF in open-shell calculation. Options: 'UHF', 'QRO', None.
   * - ``DFTreference``
     - string
     - None
     - Use DFT reference WF (orbitals) in all CC calculations. Options: (any valid ORCA DFT keyword). Default: None
   * - ``DFT_RI``
     - Boolean
     - False
     - If using DFT-reference, if DFT_RI is True then RIJ/RIJCOSX with SARC/J and defgrid3 is used to calculate DFT orbitals.
   * - ``auxbasis``
     - string
     - "autoaux-max"
     - Auxiliary basis set for CC integrals (/C type). Options: 'autoaux, 'autoaux-max' Default:  "autoaux-max"
   * - ``memory``
     - integer
     - 5000
     - Memory in MBs to use by ORCA.
   * - ``scfsetting``
     - string
     - 'TightSCF'
     - SCF-convergence setting in ORCA. Options: 'NormalSCF', 'TightSCF', 'VeryTightSCF', 'ExtremeSCF'.
   * - ``DLPNO``
     - Boolean
     - False
     - Use of DLPNO approximation for coupled cluster calculations or not.
   * - ``T1``
     - Boolean
     - False
     - Option to use iterative triples, i.e. DLPNO-CCSD(T1) instead of the default DLPNO-CCSD(T0) in all steps.
   * - ``T1correction``
     - Boolean
     - False
     - Option to calculate T1 as a single-step correction instead.
   * - ``T1corrbasis_size``
     - string
     - 'Large'
     - Size of basis set in T1 correction. Options: 'Large' (larger cardinal basis), 'Small' (smaller cardinal basis).
   * - ``T1corrpnosetting``
     - string
     - 'NormalPNOreduced'
     - PNO setting for the T1  correction. Options: 'LoosePNO', 'NormalPNO', 'NormalPNOreduced' (TCutPNO=1e-6), 'TightPNO'.
   * - ``pnoextrapolation``
     - list of 2 integers and 1 string
     - [1e-6,1e-7,1.5,'TightPNO']
     - Parameters for PNO-extrapolation (X,Y,Z): X and Y being TCutPNO thresholds while Z signifies the setting for the other thresholds. 
   * - ``FullLMP2Guess``
     - Boolean
     - None
     - Whether to use Full-local MP2 guess in DLPNO calculations. Only use if all systems are closed-shell.
   * - ``alpha``
     - float
     - False
     - Manual alpha extrapolation parameter for SCF-energy extrapolation.
   * - ``beta``
     - float
     - None
     -  Manual beta extrapolation parameter for correlation-energy extrapolation.
   * - ``extrainputkeyword``
     - string
     - None
     - Optional extra simple-input-keyword to add in ORCA inputfile.
   * - ``extrablocks``
     - string
     - None
     - Optional extra ORCA block-input lines to add to ORCA inputfile.
   * - ``guessmode``
     - string
     - 'CMatrix'
     - What ORCA Guessmode to use when doing basis-set projections of orbitals. Options: 'CMatrix' (more robust), 'FMatrix' (cheaper).
   * - ``atomicSOcorrection``
     - Boolean
     - False
     - Whether to add the experimental atomic spin-orbit energy to system if the system is an atom.
   * - ``FCI``
     - Boolean
     - False
     - Whether to extrapolate the CCSD(T) calculation to the Full-CI limit by the Goodson formula.
   * - ``F12``
     - Boolean
     - False
     - Whether to do explicitly correlated CCSD(T)-F12 instead of CCSD(T)/CBS extrapolation. Use with basisfamily='cc-f12'.
   * - ``CVSR``
     - Boolean
     - False
     - Perform additional core-valence+scalar-relativistic correction.
   * - ``CVbasis``
     - string
     - "W1-mtsmall"
     - The core-valence basis set to use. The default "W1-mtsmall" is only available for elements H-Ar. Alternative: some other appropriate core-valence basis set.
   * - ``SCFextrapolation``
     - Boolean
     - True
     - Whether the SCF energies are extrapolated or not. If False then the largest SCF energy calculated will be used (e.g. the def2-QZVPP energy in a def2/[3,4] job).


**Basis-family options**

Appropriate all-electron or valence+ECP basis sets for each element with basis-families such as : cc, aug-cc, def2, ma-def2. 
If instead an all-electron relativistic approch is desired for all elements then basisfamily="cc-dk", "def2-zora", "def2-dkh" and relativity='DKH' or 'ZORA' can be chosen instead.


.. note:: - "def2" (Ahlrichs all-electron basis sets for H-Kr, valence basis+def2-ECP for K-Rn)
  - "ma-def2" (minimally augmented diffuse Ahlrichs basis sets)
  - "cc" (correlation consistent basis sets, cc-pVnZ for light elements and cc-pVnZ-PP (SK-MCDHF ECP) for heavy elements (Sr-Xe, Hf-Rn, Ba, Ru, U)). Note: not available for K.
  - "aug-cc" (augmented correlation consistent basis sets, cc-pVnZ for light elements and aug-cc-pVnZ-PP for heavy elements)
  - "cc-dk" (DKH-recontracted correlation consistent basis sets, cc-pVnZ-DK for light elements and cc-pVnZ-DK for heavy elements)
  - "aug-cc-dk" (DKH-recontracted aug correlation consistent basis sets, aug-cc-pVnZ-DK for light elements and aug-cc-pVnZ-DK for heavy elements)
  - "def2-zora" (ZORA-recontracted Ahlrichs basis sets or SARC-ZORA basis sets for heavy elements)
  - "ma-def2-zora" (minimally augmented ZORA-recontracted Ahlrichs basis sets or SARC-ZORA basis sets for heavy elements)
  - "def2-dkh" (DKH-recontracted Ahlrichs basis sets or SARC-DKH basis sets for heavy elements)
  - "def2-x2c" (All-electron X2C relativistic basis sets for H-Rn)
  - "ma-def2-dkh" (minimally augmented DKH-recontracted Ahlrichs basis sets or SARC-DKH basis sets for heavy elements)
  - "cc-CV" (Core-valence correlation consistent basis sets, cc-pwCVnZ)
  - "aug-cc-CV" (augmented core-valence correlation consistent basis sets, aug-cc-pwCVnZ)
  - "cc-CV-dk" (DKH-recontracted core-valence correlation consistent basis sets, cc-pwCVnZ-DK)
  - "aug-cc-CV-dk" (augmented DKH-recontracted core-valence correlation consistent basis sets, aug-cc-pwCVnZ-DK)
  - "cc-CV_3dTM-cc_L" (All-electron DKH protocol for 3d TM complexes. cc-pwCVnZ-DK on 3d transition metals, cc-pVNZ-DK on everything else.)
  - "aug-cc-CV_3dTM-cc_L" (Augmented all-electron DKH protocol for 3d TM complexes. cc-pwCVnZ-DK on 3d transition metals, aug-cc-pVNZ-DK on everything else.)
  - "cc-f12" (correlation consistent F12 basis sets for CCSD(T)-F12 theory.)


+---------------------+---------------------------------+------------------------------+----------------------------+
| Basis-family        | Basis-sets                      | Cardinals (n)                | ECP or relativity          |
+=====================+=================================+==============================+============================+
| def2                | Ahlrichs def2                   | - 2: def2-SVP                | def2-ECP                   |
|                     | on all atoms H-Rn               | - 3: def2-TZVPP              | on Rb-Rn                   |
|                     |                                 | - 4: def2-QZVPP              |                            |
+---------------------+---------------------------------+------------------------------+----------------------------+
| ma-def2             | Minimally augmented             | - 2: ma-def2-SVP             | def2-ECP                   |
|                     | diffuse def2                    | - 3: ma-def2-TZVPP           | on Rb-Rn                   |
|                     | on all atoms H-Rn               | - 4: ma-def2-QZVPP           |                            |
+---------------------+---------------------------------+------------------------------+----------------------------+
| def2-zora           | - H-Kr : ZORA-def2-TZVP         | - 2: (SARC-ZORA/def2)-SVP    | relativity='ZORA'          |
|                     | - Rb-Rn : SARC-ZORA-TZVP        | - 3: (SARC-ZORA/def2)-TZVPP  |                            |
|                     |                                 | - 4: (SARC-ZORA/def2)-QZVPP  |                            |
+---------------------+---------------------------------+------------------------------+----------------------------+
| ma-def2-zora        | - H-Kr : ma-ZORA-def2-TZVP      | - 2: (SARC-ZORA/def2)-SVP    | relativity='ZORA'          |
|                     | - Rb-Rn: (SARC-ZORA/def2)-TZVPP | - 3: (SARC-ZORA/def2)-TZVPP  |                            |
|                     |                                 | - 4: (SARC-ZORA/def2)-QZVPP  |                            |
+---------------------+---------------------------------+------------------------------+----------------------------+
| def2-dkh            | - H-Kr : DKH-def2-TZVP          | - 2: (SARC-DKH/def2)-SVP     | relativity='DKH'           |
|                     | - Rb-Rn : SARC-DKH-TZVP         | - 3: (SARC-DKH/def2)-TZVPP   |                            |
|                     |                                 | - 4: (SARC-DKH/def2)-QZVPP   |                            |
+---------------------+---------------------------------+------------------------------+----------------------------+
| ma-def2-dkh         | - H-Kr : ma-DKH-def2-TZVP       | - 2: (SARC-DKH/def2)-SVP     | relativity='DKH'           |
|                     | - Rb-Rn: (SARC-DKH/def2)-TZVPP  | - 3: (SARC-DKH/def2)-TZVPP   |                            |
|                     |                                 | - 4: (SARC-DKH/def2)-QZVPP   |                            |
+---------------------+---------------------------------+------------------------------+----------------------------+
| def2-x2c            | - H-Rn : x2c-nVP(P)all          | - 2: x2c-SVP-all             | relativity='DKH'           |
|                     |                                 | - 3: x2c-TZVPP-all           | ( later: relativity='X2C') |
|                     |                                 | - 4: x2c-QZVPP-all           |                            |
+---------------------+---------------------------------+------------------------------+----------------------------+
| cc                  | - H-Kr: cc-pVnZ                 | - 2: cc-pVDZ(-PP)            | SK-MCDHF-RSC               |
|                     | - Sr-Xe: cc-pVnZ-PP             | - 3: cc-pVTZ(-PP)            | on Sr-Xe, Hf-Rn,           |
|                     | - Hf-Rn: cc-pVnZ-PP             | - 4: cc-pVQZ(-PP)            | Ba,Ra,U                    |
|                     | - Ba,Ra,U: cc-pVnZ-PP           | - 5: cc-pV5Z(-PP)            |                            |
|                     |                                 | - 6: cc-pV6Z (H-Ar only)     |                            |
+---------------------+---------------------------------+------------------------------+----------------------------+
| cc-f12              | - H-Ar: cc-pVnZ-F12             | - 2: cc-pVDZ(-PP)-F12        | SK-MCDHF-RSC               |
|                     | - Ga-Kr: cc-pVnZ-PP-F12         | - 3: cc-pVTZ(-PP)-F12        | on Ga-Kr, In-Xe, Tl-Rn     |
|                     | - In-Xe: cc-pVnZ-PP-F12         | - 4: cc-pVQZ(-PP)-F12        |                            |
| (use with F12=True) | - Tl-Rn: cc-pVnZ-PP-F12         |                              |                            |
|                     |                                 |                              |                            |
+---------------------+---------------------------------+------------------------------+----------------------------+
| aug-cc              | - H-Kr: aug-cc-pVnZ,            | - 2: aug-cc-pVDZ(-PP)        | SK-MCDHF-RSC               |
|                     | - Sr-Xe: aug-cc-pVnZ-PP,        | - 3: aug-cc-pVTZ(-PP)        | on Sr-Xe, Hf-Rn,           |
|                     | - Hf-Rn: aug-cc-pVnZ-PP,        | - 4: aug-cc-pVQZ(-PP)        | Ba,Ra,U                    |
|                     | - Ba,Ra,U: aug-cc-pVnZ-PP       | - 5: aug-cc-pV5Z(-PP)        |                            |
|                     |                                 | - 6: aug-cc-pV6Z (H-Ar Only) |                            |
+---------------------+---------------------------------+------------------------------+----------------------------+
| cc-dk               | - cc-pVnZ-DK on H-Ar,           | - 2: cc-pVDZ-DK              |                            |
|                     | - Sc-Kr, Y-Xe, Hf-Rn,           | - 3: cc-pVTZ-DK              | relativity='DKH'           |
|                     | - 4: cc-pVQZ-DK                 | - 4: cc-pVQZ-DK              |                            |
|                     | - (missing QZ for Y-Cd)         | - 5: cc-pV5Z-DK              |                            |
+---------------------+---------------------------------+------------------------------+----------------------------+
| aug-cc-dk           | - cc-pVnZ-DK on H-Ar,           | - 2: aug-cc-pVDZ-DK          |                            |
|                     | - Sc-Kr, Y-Xe, Hf-Rn,           | - 3: aug-cc-pVTZ-DK          | relativity='DKH'           |
|                     | - 4: aug-cc-pVQZ-DK             | - 4: aug-cc-pVQZ-DK          |                            |
|                     | - (missing QZ for Y-Cd)         | - 5: aug-cc-pV5Z-DK          |                            |
+---------------------+---------------------------------+------------------------------+----------------------------+
| cc-CV               | - H-Kr: cc-pwCVnZ               | - 2: cc-pwCVDZ(-PP)          | SK-MCDHF-RSC               |
|                     | - Sr-Xe: cc-pwCVnZ-PP           | - 3: cc-pwCVTZ(-PP)          | on Sr-Xe, Hf-Rn,           |
|                     | - Hf-Rn: cc-pwCVnZ-PP           | - 4: cc-pwCVQZ(-PP)          | Ba,Ra,U                    |
|                     | - Ba,Ra,U: cc-pwCVnZ-PP         | - 5: cc-pWCV5Z(-PP)          |                            |
|                     |                                 |                              |                            |
+---------------------+---------------------------------+------------------------------+----------------------------+
| aug-cc-CV           | - H-Kr: aug-cc-pwCVnZ           | - 2: aug-cc-pwCVDZ(-PP)      | SK-MCDHF-RSC               |
|                     | - Sr-Xe: aug-cc-pwCVnZ-PP       | - 3: aug-cc-pwCVTZ(-PP)      | on Sr-Xe, Hf-Rn,           |
|                     | - Hf-Rn: aug-cc-pwCVnZ-PP       | - 4: aug-cc-pwCVQZ(-PP)      | Ba,Ra,U                    |
|                     | - Ba,Ra,U: aug-cc-pwCVnZ-PP     | - 5: aug-cc-pWCV5Z(-PP)      |                            |
|                     |                                 |                              |                            |
+---------------------+---------------------------------+------------------------------+----------------------------+
| cc-CV-dk            | - H-Be,Na-Mg: cc-pwCVnZ-DK      | - 2: cc-(pwC)VDZ-DK          |                            |
|                     | - B-Ne: cc-pVnZ-DK (!)          | - 3: cc-(pwC)VTZ-DK          | relativity='DKH'           |
|                     | - Al-Ar: cc-pVnZ-DK (!)         | - 4: cc-(pwC)VQZ-DK          |                            |
|                     | - Ca-Zn: cc-pVwCnZ-DK           | - 5: cc-(pwC)V5Z-DK          |                            |
|                     | - missing QZ for Y-Cd           |                              |                            |
+---------------------+---------------------------------+------------------------------+----------------------------+
| aug-cc-CV-dk        | - H-Be,Na-Mg: aug-cc-pwCVnZ-DK  | - 2: aug-cc-(pwC)VDZ-DK      |                            |
|                     | - B-Ne: aug-cc-pVnZ-DK (!)      | - 3: aug-cc-(pwC)VTZ-DK      | relativity='DKH'           |
|                     | - Al-Ar: aug-cc-pVnZ-DK (!)     | - 4: aug-cc-(pwC)VQZ-DK      |                            |
|                     | - Ca-Zn: aug-cc-pVwCnZ-DK       | - 5: aug-cc-(pwC)V5Z-DK      |                            |
|                     | - missing QZ for Y-Cd           |                              |                            |
+---------------------+---------------------------------+------------------------------+----------------------------+
| cc-CV_3dTM-cc_L     | - H-Kr: cc-pVnZ-DK              | - 2: cc-(pwC)VDZ-DK          |                            |
|                     | - Sc-Zn: cc-pwCVnZ-DK (!)       | - 3: cc-(pwC)VTZ-DK          | relativity='DKH'           |
|                     | - Ga-Rn: cc-pVnZ-DK             | - 4: cc-(pwC)VQZ-DK          |                            |
|                     |                                 | - 5: cc-(pwC)V5Z-DK          |                            |
|                     |                                 |                              |                            |
+---------------------+---------------------------------+------------------------------+----------------------------+
| aug-cc-CV_3dTM-cc_L | - H-Kr: aug-cc-pVnZ-DK          | - 2: (aug)-cc-(pwC)VDZ-DK    |                            |
|                     | - Sc-Zn: cc-pwCVnZ-DK (!)       | - 3: (aug)-cc-(pwC)VTZ-DK    | relativity='DKH'           |
|                     | - Ga-Rn: aug-cc-pVnZ-DK         | - 4: (aug)-cc-(pwC)VQZ-DK    |                            |
|                     |                                 | - 5: (aug)-cc-(pwC)V5Z-DK    |                            |
|                     |                                 |                              |                            |
+---------------------+---------------------------------+------------------------------+----------------------------+

.. note::  Note: often missing basis sets for K and Ca. Sometimes there are missing basis sets for specific elements and specific cardinals.


#########################################
ORCA_CC_CBS_Theory Examples
#########################################

**Basic examples**

.. code-block:: python
    
    N2=Fragment(xyzfile='n2.xyz')
    cc = ORCA_CC_CBS_Theory(elements=["N"], cardinals = [2,3], basisfamily="cc", numcores=1)
    Singlepoint(theory=cc, fragment=N2)


The example above defines an N2 fragment (from file n2.xyz) and runs a single-point calculation using the defined ORCA_CC_CBS_Theory object. 
Multiple CCSD(T) calculations are then carried out using the different basis sets specified by the basis-family and the cardinals.
Cardinals=[2,3] and basisfamily="cc" means that the cc-pVDZ and cc-pVTZ basis sets will be used.
Separate basis-set extrapolation of SCF and correlation energies is then performed. Appropriate extrapolation parameters for 2-point extrapolations with this basis set family are chosen.

.. code-block:: python

    ferrocene=Fragment(xyzfile='ferrocene.xyz')
    cc = ORCA_CC_CBS_Theory(elements=["Fe", "C", "H"], cardinals = [2,3], basisfamily="def2", numcores=1, 
        DLPNO=True, pnosetting="NormalPNO", T1=False)
    Singlepoint(theory=cc, fragment=ferrocene)

For a larger molecule like ferrocene, regular CCSD(T) is quite an expensive calculation and so here we invoke the DLPNO approximation via DLPNO=True.
We use the 'def2' basis family here with cardinals=[2,3] meaning that the def2-SVP and def2-TZVPP basis sets will be used.
The DLPNO approximation error can be controlled via threshold keywords ('LoosePNO', 'NormalPNO', 'TightPNO'), here we choose 'NormalPNO'. 
We also choose the regular triples approximation (DLPNO-CCSD(T0) by setting T1 to False.

.. code-block:: python

    ferrocene=Fragment(xyzfile='ferrocene.xyz')
    cc = ORCA_CC_CBS_Theory(elements=ferrocene.elems, cardinals = [3,4], basisfamily="cc-CV_3dTM-cc_L", relativity='DKH', numcores=1, 
        DLPNO=True, pnosetting="extrapolation", pnoextrapolation=[6,7] T1=True)
    Singlepoint(theory=cc, fragment=ferrocene)

Finally we crank up the accuracy even further by choosing cardinals=[3,4], switch to the basisfamily="cc-CV_3dTM-cc_L and activate the 'DKH' relativistic approximation.
This calculation will utilize a mixed metal-ligands basis set: cc-pwCVTZ-DK/cc-pwCVQZ-DK on Fe and cc-pVDZ-DK/cc-pVTZ-DK on C,H.
Instead of using a single DLPNO threshold we here calculate DLPNO-CCSD(T) energies using 2 PNO tresholds and extrapolate to the PNO-limit.
Finally we set T1 keyword to True which will tell ORCA to do a more accurate iterative triples DLPNO-CCSD(T1) approximation.

For additional examples on using ORCA_CC_CBS_Theory on real-world systems and showing real data see:  :doc:`Highlevel_CC_CBS_workflows`


##############################
Reaction_Highlevel_Analysis
##############################

In order to facilitate the analysis of basis-set and/or PNO convergence in CCSD(T) calculations, the **Reaction_Highlevel_Analysis** function can be used.
It will read in a list of ASH fragments and reaction stoichiometry and calculate the reaction energy with multiple levels of theory and plot the results using Matplotlib.
This allows one to easily see how well converged the results are.

CCSD(T) calculations are performed both with def2 (up to QZ level) and cc basis sets (up to 6Z level), explicitly correlated CCSD(T)-F12 calculations (up to QZ-F12) 
and complete basis set extrapolations are performed.
Note that the large-basis cc-pV5Z and cc-pV6Z calculations can not be carried out for all systems. Set highest_cardinal to a lower number if required.


.. warning:: The plots require the Matplotlib library to be installed. 

To be added: PNO-extrapolation options

.. code-block:: python

    def Reaction_Highlevel_Analysis(fraglist=None, stoichiometry=None, numcores=1, memory=7000, reactionlabel='Reactionlabel', energy_unit='kcal/mol',
                                    def2_family=True, cc_family=True, aug_cc_family=False, F12_family=True, DLPNO=False, extrapolation=True, highest_cardinal=6,
                                    plot=True ):
        """Function to perform high-level CCSD(T) calculations for a reaction with associated plots.
        Performs CCSD(T) with cc and def2 basis sets, CCSD(T)-F12 and CCSD(T)/CBS extrapolations

        Args:
            fragment ([type], optional): [description]. Defaults to None.
            fraglist ([type], optional): [description]. Defaults to None.
            stoichiometry ([type], optional): [description]. Defaults to None.
            numcores (int, optional): [description]. Defaults to 1.
            memory (int, optional): [description]. Defaults to 7000.
            reactionlabel (str, optional): [description]. Defaults to 'Reactionlabel'.
            energy_unit (str): Energy unit for ReactionEnergy. Options: 'kcal/mol', 'kJ/mol', 'eV', 'cm-1'. Default: 'kcal/mol'
            def2_family (bool, optional): [description]. Defaults to True.
            cc_family (bool, optional): [description]. Defaults to True.
            F12_family (bool, optional): [description]. Defaults to True.
            highest_cardinal (int, optional): [description]. Defaults to 5.
            plot (Boolean): whether to plot the results or not (requires Matplotlib). Defaults to True. 
        """

Example (Bond Dissociation Energy of N2): 

.. code-block:: python

    from ash import *

    #Define molecular fragments from XYZ-files or other
    N2=Fragment(xyzfile='n2.xyz', charge=0, mult=1, label='N2')
    N=Fragment(atom='N', charge=0, mult=4, label='N')

    #Create a list of fragments and define the stoichiometry
    specieslist=[N2, N]
    stoichiometry=[-1,2]
    reactionlabel='N2_BDE'

    # Call Reaction_Highlevel_Analysis
    Reaction_Highlevel_Analysis(fraglist=specieslist, stoichiometry=stoichiometry, numcores=1, memory=7000, reactionlabel=reactionlabel,
                                    def2_family=True, cc_family=True, F12_family=True, extrapolation=True, highest_cardinal=5 )

The outputfile will contain the CCSD(T) total energies and reaction energies for each species and basis set level.
Additionally energy vs. basis-cardinal plots are created for both the total energy for each species and the reaction energy.


.. image:: figures/N2_BDE.png
   :align: center
   :width: 700


.. image:: figures/N2_Energy.png
   :align: center
   :width: 700

.. image:: figures/N_Energy.png
   :align: center
   :width: 700


