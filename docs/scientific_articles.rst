Scientific articles using ASH
================================

ASH is a relatively new code but has been used for a few different research projects that led to publications.

If you have used ASH in your research and you would like to mention it here, contact me!

###################################
MM and QM/MM protein studies
###################################

**QM/MM modelling of a CN-inhibited state of FeFe hydrogenase**

`2023 Chem. Sci. article <https://pubs.rsc.org/en/content/articlelanding/2023/sc/d2sc06098a>`_ 

This article used the QM/MM module of ASH together with the ORCA interface (:doc:`ORCA-interface`)
for the QM part and the OpenMM interface (:doc:`OpenMM-interface`) for the MM part. 
The OpenMMTheory interface used CHARMM-style forcefield files.



###################################
Highlevel WFT workflows
###################################

**Multistep DLPNO-CCSD(T)/CBS workflow for a transition metal complex**

`2023 PCCP article <https://pubs.rsc.org/en/content/articlelanding/2023/cp/d2cp04715b>`_ 

This article used ORCA_CC_CBS_Theory (:doc:`module_highlevel_workflows`) functionality in ASH.
Below is a script that describes a recommended DLPNO-CCSD(T)/CBS workflow that worked well for this class of metallocenes
and should be reasonably reliable in general (assuming coupled cluster is reliable).
It uses a CBS(3/4) basis set extrapolation using the cc-pVnZ-DK basis set family, together with BP86 reference orbitals, 
DKH scalar relativistic Hamiltonian, PNO extrapolation using the cheaper approach and the cheaper T1 correction described
in th article.

.. code-block:: python

    from ash import *
    numcores=24 #Number of cores reserved
    actualcores=16 #Number of cores used
    #Defining molecular fragments
    cpco0=Fragment(xyzfile="CpCo_0_gas.xyz", charge=0, mult=2)
    cpcoI=Fragment(xyzfile="CpCo_I_gas.xyz", charge=1, mult=1)
    # Defining species, stoichiometry and reaction specieslist=[cpco0,cpcoI]
    stoichiometry=[-1, 1]
    reaction = Reaction(fragments=specieslist, stoichiometry=stoichiometry)
    #Defining a ORCA_CC_CBS_Theory object
    cc = ORCA_CC_CBS_Theory(elements=cpco0.elems, cardinals=[3,4], basisfamily="cc-dk", DFTreference="BP86", 
        DLPNO=True, CVSR=False, T1correction=True, T1corrbasis_size='Small', T1corrpnosetting='NormalPNOreduced', 
        numcores=actualcores, pnosetting=" extrapolation", pnoextrapolation=[1e-6,3.33e-7,2.38,'NormalPNO'], 
        memory=20000, scfsetting="Verytightscf", relativity='DKH', SCFextrapolation=False)
    #Running reaction
    Singlepoint_reaction(theory=cc, reaction=reaction, unit='eV')



