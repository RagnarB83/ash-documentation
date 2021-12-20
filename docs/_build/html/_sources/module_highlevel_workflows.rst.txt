Highlevel workflows
======================================

These high-level workflows (singlepoint energy protocols) can either be used on their own as a theory-level or used as a SP_theory in **thermochemprotocol** (see :doc:`module_workflows`) or as a theory in **run_benchmark** (see :doc:`module_benchmarking`) .
All of these protocols use the ORCA quantum chemistry code and give the 0 K electronic energy (no ZPVE). Gradients are not available and these can thus not be used in geometry optimizations or dynamics jobs.


**CC_CBS_Theory**

This is an ASHTheory that carries out a multi-step single-point protocol to give a CCSD(T)/CBS estimated energy
Multiple ORCA calculations for the given geometry are carried out and the SCF and correlation energies extrapolated to the CCSD(T)/CBS limit using either regular CCSD(T) theory or DLPNO-CCSD(T) theory.
This workflow is flexible and features multiple ways of approaching the CBS limit and the PNO limit.
Various options affecting the accuracy, efficiency and robustness of the protocol can be chosen.
Many basis set families can be chosen that are available for most of the periodic table.
Atomic spin-orbit coupling can be automatically included if system is an atom.

.. code-block:: python

  class CC_CBS_Theory:
      def __init__(self, elements=None, cardinals = [2,3], basisfamily="def2", relativity=None, charge=None, mult=None, orcadir=None,
                stabilityanalysis=False, numcores=1, CVSR=False, CVbasis="W1-mtsmall", F12=False, Openshellreference=None, 
                DFTreference=None, DFT_RI=False, auxbasis="autoaux-max", memory=5000, scfsetting='TightSCF', DLPNO=False, T1=True, 
                pnosetting='extrapolation', pnoextrapolation=[6,7], FullLMP2Guess=False, alpha=None, beta=None, 
                extrainputkeyword='', extrablocks='', FCI=False, guessmode='Cmatrix', atomicSOcorrection=False):


.. code-block:: text

  Options:

  - elements (list of strings): A list of elements need to be provided when creating the object. This list must contain all the elements 
  found in the molecular system (duplicates are fine). This list of elements is used to select an appropriate basis-set member for the 
  chosen basis-set-family for each element. 
  - cardinals (list of ints): Cardinals can be set to [2,3], [3,4], [4,5] or [5,6] (though 5Z and 6Z basis-sets are not always available 
  for all elements or basis set families).
  - basisfamily (string): Can be set to various options. See table below. 
  - relativity (string): Options: 'DKH', 'ZORA', 'NoRel', None. Default: None
  - charge/mult (integer): Charge/multiplicity of system. Should be specified unless theory is used as part of run_benchmark or thermochemprotocol.
  - orcadir(string): Path-to-ORCA. Default: None
  - stabilityanalysis (Boolean): Perform SCF stability analysis. Default: families
  - numcores (integer): Number of cores to use in ORCA calculation.
  - Openshellreference (string): Use alternative reference WF in open-shell calculation. Options: 'UHF', 'QRO' Default: None
  - DFTreference (string): Use DFT reference WF (orbitals) in all CC calculations. Options: (any valid DFT keyword). Default: None
  - DFT_RI (Boolean): If using DFT-reference, if DFT_RI is True then RIJ/RIJCOSX with SARC/J and defgrid3 is used to 
  calculate DFT orbitals. Default: False
  - auxbasis (string): Auxiliary basis set for CC integrals (/C type). Options: 'autoaux, 'autoaux-max' Default:  "autoaux-max"
  - memory (integer): Memory in MBs to use by ORCA. Default: 5000
  - scfsetting (string): SCF-convergence setting. Options: 'NormalSCF', 'TightSCF', 'VeryTightSCF', 'ExtremeSCF'. Default: 'TightSCF'
  - DLPNO (Boolean): Use of DLPNO approximation or not. Default False.
  - T1 (Boolean): option to use iterative triples, i.e. DLPNO-CCSD(T1) instead of the default DLPNO-CCSD(T0) 
  (more accurate, more expensive). Default: True
  - pnosetting (string): Accuracy-threshold of the DLPNO approximation. Options: 'LoosePNO', 'NormalPNO', 'TightPNO', 'extrapolation'.
  Default: 'extrapolation'
  - pnoextrapolation (list of 2 ints, X and Y): If using PNO-extrapolation then 2 DLPNO-calculations will be performed with TCutPNO=1e-X and 
  TCutPNO=1e-Y and extrapolated to PNO limit. Default: [6,7] 
  - FullLMP2Guess (Boolean): Whether to use Full-local MP2 guess in DLPNO calculations. Only use if all systems are closed-shell. Default: False
  - alpha(float): Manual extrapolation parameter for SCF-extrapolation. Default: None
  - beta(float): Manual extrapolation parameter for correlation-energy extrapolation. Default: None
  - extrainputkeyword (string): Optional extra simple-input-keyword to add in ORCA inputfile. Default: None
  - extrablocks (string): Optional extra ORCA block-input lines to add to ORCA inputfile. Default: None
  - guessmode(string): What ORCA Guessmode to use when doing basis-set projections of orbitals. 
  Options: 'CMatrix' (more robust), 'FMatrix' (cheaper). Default: 'CMatrix'
  - atomicSOcorrection (Boolean). Whether to add the experimental atomic spin-orbit energy to system if the 
  system is an atom. Default: False
  - FCI (Boolean): Whether to extrapolate the CCSD(T) calculation to the Full-CI limit by the Goodson formula. (NOT ACTIVE)
  - CVSR (Boolean). Perform additional core-valence+scalar-relativistic correction. Default: False  (NOT ACTIVE)
  - CVbasis (string): The core-valence basis set to use. Default: "W1-mtsmall" (NOT ACTIVE)
  - F12 (Boolean): Use CCSD(T)-F12 approach. Default: False (NOT ACTIVE)


**Basisfamily options:**

Appropriate all-electron or valence+ECP basis sets for each element for basis-families such as : cc, aug-cc, def2, ma-def2. If instead an all-electron relativistic approch is desired for all elements then basisfamily="cc-dk" and relativity='DKH' can be chosen instead.


.. note:: - "def2" (Ahlrichs all-electron basis sets for H-Kr, valence basis+def2-ECP for K-Rn)
  - "ma-def2" (minimally augmented diffuse Ahlrichs basis sets)
  - "def2-zora" (ZORA-recontracted Ahlrichs basis sets or SARC-ZORA basis sets for heavy elements)
  - "ma-def2-zora" (minimally augmented ZORA-recontracted Ahlrichs basis sets or SARC-ZORA basis sets for heavy elements)
  - "def2-dkh" (DKH-recontracted Ahlrichs basis sets or SARC-DKH basis sets for heavy elements)
  - "ma-def2-dkh" (minimally augmented DKH-recontracted Ahlrichs basis sets or SARC-DKH basis sets for heavy elements)
  - "cc" (correlation consistent basis sets, cc-pVnZ for light elements and cc-pVnZ-PP (SK-MCDHF ECP) for heavy elements (Sr-Xe, Hf-Rn, Ba, Ru, U)). Note: not available for K.
  - "aug-cc" (augmented correlation consistent basis sets, cc-pVnZ for light elements and aug-cc-pVnZ-PP for heavy elements)
  - "cc-dk" (DKH-recontracted correlation consistent basis sets, cc-pVnZ-DK for light elements and cc-pVnZ-DK for heavy elements)
  - "aug-cc-dk" (DKH-recontracted aug correlation consistent basis sets, aug-cc-pVnZ-DK for light elements and aug-cc-pVnZ-DK for heavy elements)
  - "cc-CV" (Core-valence correlation consistent basis sets, cc-pwCVnZ)
  - "aug-cc-CV" (augmented core-valence correlation consistent basis sets, aug-cc-pwCVnZ)
  - "cc-CV-dk" (DKH-recontracted core-valence correlation consistent basis sets, cc-pwCVnZ-DK)
  - "aug-cc-CV-dk" (augmented DKH-recontracted core-valence correlation consistent basis sets, aug-cc-pwCVnZ-DK)


+--------------+-----------------------+---------------------+------------------+
| Basis-family | Basis-sets            | Cardinals (n)       | ECP on elems     |
+==============+=======================+=====================+==================+
| def2         | Ahlrichs def2         | - 2: def2-SVP       | def2-ECP         |
|              | on all atoms          | - 3: def2-TZVPP     | on Rb-Rn         |
|              |                       | - 4: def2-QZVPP     |                  |
+--------------+-----------------------+---------------------+------------------+
| cc           | - H-Kr: cc-pVnZ,      | - 2: cc-pVDZ        | SK-MCDHF-RSC     |
|              | - Sr-Xe: cc-pVnZ-PP,  | - 3: cc-pVTZ        | on Sr-Xe, Hf-Rn, |
|              | - Hf-Rn: cc-pVnZ-PP,  | - 4: cc-pVQZ        | Ba,Ra,U          |
|              | - Ba,Ra,U: cc-pVnZ-PP | - 5: cc-pV5Z        |                  |
|              |                       | - 6: cc-pV6Z (H-Ar) |                  |
+--------------+-----------------------+---------------------+------------------+
| cc           | - H-Kr: cc-pVnZ,      | - 2: cc-pVDZ        | SK-MCDHF-RSC     |
|              | - Sr-Xe: cc-pVnZ-PP,  | - 3: cc-pVTZ        | on Sr-Xe, Hf-Rn, |
|              | - Hf-Rn: cc-pVnZ-PP,  | - 4: cc-pVQZ        | Ba,Ra,U          |
|              | - Ba,Ra,U: cc-pVnZ-PP | - 5: cc-pV5Z        |                  |
|              |                       | - 6: cc-pV6Z (H-Ar) |                  |
+--------------+-----------------------+---------------------+------------------+

* Note: missing basis for K.
* 


**Basic example:**

.. code-block:: python

    N2=Fragment(xyzfile='n2.xyz')
    cc = CC_CBS_Theory(elements=["N"], cardinals = [2,3], basisfamily="cc", pnosetting='extrapolation', pnoextrapolation=[6,7], DLPNO=True, numcores=1)
    Singlepoint(theory=cc, fragment=N2)


The example above defines an N2 fragment (from n2.xyz) and runs multiple CCSD(T) calculations, utilizing basis-set extrapolation of SCF and correlation energies.
Cardinals=[2,3] and basisfamily="cc" means that the cc-pVDZ and cc-pVTZ basis sets will be used and extrapolated to the basis set limit. 
Appropriate extrapolation parameters for 2-point extrapolations with this basis set family are chosen.


**Example: DLPNO-CCSD(T1)/CBS with PNO extrapolation on a 4d-metal complex with the Ahlrichs def2-SVP/def2-TZVPP extrapolation:**

.. code-block:: python

    complex=Fragment(xyzfile='ru-phosphine-complex.xyz')
    cc = CC_CBS_Theory(elements=["Ru", "P", "H", "O", "N" ], cardinals = [2,3], basisfamily="def2", DLPNO=True, 
                  pnosetting='NormalPNO', numcores=1)
    Singlepoint(theory=cc, fragment=complex)

In this example of a large ruthenium metal complex we can not afford to do regular CCSD(T) calculations and utilize instead the powerful local-correlation DLPNO approximation.
Here we use the def2 basis family and a def2-ECP effective core-potential will be automatically selected for ruthenium. We choose cardinals=[2,3] here and this means that we do the relatively cheap def2-SVP/def2-TZVPP extrapolation.
The DLPNO approximation utilizes thresholds that determine the accuracy of the DLPNO approximation (compared to unapproximated CCSD(T)).
By setting pnosetting="NormalPNO" we get the default PNO settings that are reasonably accurate. Other options are: 'LoosePNO' (not recommended) and 'TightPNO' (more accurate, more expensive), and 'extrapolation' (see below).


**Example: DLPNO-CCSD(T)/CBS with PNO extrapolation on a 4d-metal complex with the cc-pVnZ and cc-pVnZ-PP (n=3,4) extrapolation:**

.. code-block:: python

    complex=Fragment(xyzfile='ru-phosphine-complex.xyz')
    #Note: here providing list of elements more conveniently from the defined fragment
    cc = CC_CBS_Theory(elements=complex.elems, cardinals = [3,4], basisfamily="cc", DLPNO=True, 
                  pnosetting='extrapolation', pnoextrapolation=[6,7], numcores=1)
    Singlepoint(theory=cc, fragment=complex)

For an even more accurate estimate of the coupled-cluster basis set limit the [3,4] extrapolation is much more reliable than [2,3] and here we also utilize the more accurate 
correlation-consistent basis set family ('cc'). For ruthenium, ASH tells ORCA to choose the cc-pVNZ-PP family for this heavy element and the 'SK-MCDHF' ECP.
To further reduce the error of the DLPNO approximation we use pnosetting="extrapolation" and pnoextrapolation=[6,7] which means that 2 DLPNO-CCSD(T) calculations will be performed
for each basis-set-cardinal calculation with different TCutPNO cutoffs (here TCutPNO=1e-6 and TCutPNO=1e-7). The results are then extrapolated to the PNO limit according to PNO extrapolation by Giovanni Bistoni and coworkers.
See these excellent papers: https://pubs.acs.org/doi/abs/10.1021/acs.jctc.0c00344 and https://pubs.acs.org/doi/abs/10.1021/acs.jpca.1c09106

**Example: DLPNO-CCSD(T)/CBS with PNO extrapolation on a 4d-metal complex with DKH relativistic approximation and cc-pwCVnZ-DK extrapolation:**

.. code-block:: python

    complex=Fragment(xyzfile='ru-phosphine-complex.xyz')
    #Note: here providing list of elements more conveniently from the defined fragment
    cc = CC_CBS_Theory(elements=complex.elems, cardinals = [3,4], basisfamily="cc-CV-dk", DLPNO=True, 
                  relativity='DKH', pnosetting='extrapolation', pnoextrapolation=[6,7], numcores=1)
    Singlepoint(theory=cc, fragment=complex)

While the cc-pVDZ-PP approach for ruthenium is affordable and accurate, for even greater accuracy one can opt for an all-electronic relativistic approach instead.
Here the Douglas-Kroll-Hess scalar relativistic Hamiltonian is used and this requires one to choose a basis-set family that has been recontracted for DKH Hamiltonians.
We could choose to use the 'cc-dk' but here we utilize the 'cc-CV-dk' family that in addition to being DKH-recontracted, features additional basis-functions typically used to describe core-valence 
correlation. The frozen-core approximation is still in use here, meaning that the extra basis functions instead serve to improve the valence-electron correlation problem instead.


https://pubs.acs.org/doi/abs/10.1021/acs.jctc.9b01109