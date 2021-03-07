Highlevel workflows module
======================================

These high-level workflows (singlepoint energy protocols) can either be called on their own (see below) or used as SP_theory keyword argument in thermochemprotocol (see above).
All of these protocols use the ORCA quantum chemistry code and give the 0 K electronic energy.


**DLPNO_CC_CBS_SP**

This workflow carries out multiple ORCA calculations for the given geometry and extrapolates to the DLPNO-CCSD(T)/CBS limit.
This workflow is flexible and features multiple ways of approaching the CBS limit and the PNO limit.
Various options affecting the accuracy, efficiency and robustness of the protocol can be chosen.
The basis set families: cc-pVnZ ('cc') and Ahlrichs def2 ('def2') can be chosen that are available for most of the periodic table (cc-pVnZ-PP for heavy elements).
A corevalence+scalar-relativistic correction (CVSR option) can be included upon request (important for high-accuracy atomization energies).
Atomic spin-orbit coupling is automatically included if system is an atom.

.. code-block:: python

    def DLPNO_CC_CBS(fragment=None, cardinals = [2,3], basisfamily="def2", charge=None, orcadir=None, mult=None, stabilityanalysis=False,
        numcores=1, CVSR=False, CVbasis="W1-mtsmall", memory=5000, pnosetting='NormalPNO', pnoextrapolation=[5,6], T1=False, scfsetting='TightSCF',
        extrainputkeyword='', extrablocks='', **kwargs):

Example:

.. code-block:: python

    N2=Fragment(xyzfile='n2.xyz')
    DLPNO_CC_CBS(fragment=N2, cardinals = [2,3], basisfamily="def2", charge=0, orcadir='/opt/orca_4.2.1', mult=1, stabilityanalysis=False,
    numcores=1, CVSR=False, memory=5000, pnosetting='extrapolation', pnoextrapolation=[5,6], T1=False, scfsetting='TightSCF')

The example above defines an N2 fragment (from n2.xyz) and runs multiple DLPNO-CCSD(T) calculations, utilizing basis-set and PNO extrapolation to give a final CCSD(T)/CBS estimate.
Cardinals=[2,3] and basisfamily="def2" means that the def2-SVP and def2-TZVPP basis sets will be used and extrapolated to the basis set limit.
pnosetting="extrapolation" and pnoextrapolation=[5,6] means that the DLPNO-calculations will be run using 2 different TCutPNO cutoffs and then extrapolated to the PNO limit.

- Cardinals can be set to [2,3] or [3,4].
- basisfamily can be set to "def2" (Ahlrichs basis sets) or "cc" (correlation consistent basis sets).
- If a heavy element is chosen (heavier than Kr) and basisfamily="cc", then the cc-pVnZ-PP and corresponding ECP will be used for the heavy element.
- T1 option utilizes iterative triples, i.e. DLPNO-CCSD(T1) instead (more accurate, more expensive).
- CVSR adds a Core-Valence-Scalar-Relativistic correction (more accurate, more expensive). The correction is performed at the DLPNO-CCSD(T) level (hardcoded to NormalPNO) using the W1-mtsmall basis set.
  Note: If "W1-mtsmall" is not available for the element involved you may have to provide an appropriate core-valence basis set keyword via CVbasis="basisname" .

TO BE DOCUMENTED:

- **W1theory**
- **W1F12theory**
- **DLPNO_W1theory**
- **DLPNO_W1F12theory**
- **DLPNO_F12**
- **DLPNO_W2theory**
