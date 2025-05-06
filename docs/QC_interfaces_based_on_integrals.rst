Quantum Chemistry Program Interfaces based on integrals
==========================================================

Some quantum chemistry programs allow calculation of correlated wavefunctions directly from integrals.
Typically this involves reading in the integrals from an FCIDUMP file.
This allows some freedom in performing quantum chemistry with integrals from different sources.

ASH has some basic support for facilitating this.
Currently supported programs are:

- MRCC : Arbitrary order CC, explicit correlation CC (F12), local natural orbital CC etc.
- pyscf : 
- ccpy : Specialized CC methods
- Block2: DMRG
- Dice: semi-stochastic heat-bath CI



######################################################
pySCF support
######################################################


For CC etc. Create script for both PySCFTheory and also simpler pyscf calls to MP2 and CC.