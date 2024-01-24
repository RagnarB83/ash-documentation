Tutorial: Kohn-Sham potential inversion
====================================================================================================

Thanks to the `KS-pies library (Kohn-Sham inversion toolkit) <https://github.com/ssnam92/KSPies>`_ is is possible to easily perform 
Kohn-Sham potential inversion using the PySCF interface in ASH. 
This is the inverse Kohn Sham problem, where instead of solving the standard potential-> density problem, 
one solves the density->potential problem.
This in particular allows one to take a correlated density matrix and derive the Kohn-Sham potential (and Kohn-Sham MOs) 
corresponding to that density matrix. This in turn allows one to benchmark a density functional approximation (DFA) error by decomposing it into a density-driven error (DE) and a functional-driven error (FE).


`KS-pies <https://ssnam92.github.io/KSPies/>`_  is easily installed like this:

.. code-block:: shell
  
    pip install kspies
    pip install opt-einsum

If you use KS-pies make sure to cite the corresponding article: `S. Nam. J. Chem. Phys. 2021, 154, 124122 <https://pubs.aip.org/aip/jcp/article/154/12/124122/380751/KS-pies-Kohn-Sham-inversion-toolkit>`_

-------------------------------
density_potential_inversion
-------------------------------

Once kspies is installed on can in an ASH script use the ASH function **density_potential_inversion** to perform a Kohn-Sham inversion of a density matrix.

.. code-block:: python

  #Standalone density-potential inversion function: takes pyscfheoryobject and DM as input
  #Solves the inversion problem and returns MO coefficients, occupations,energies and new DM
  def density_potential_inversion(pyscftheoryobj, dm, method='WY', WY_method='trust-exact',
                                  ZMP_lambda=128, ZMP_levelshift=False, ZMP_cycles=400, DF=True):

**density_potential_inversion** requires a PySCFTheory object as input and a density matrix (dm) as input.
The output of the function are: new_MO_occupations, new_MO_energies, new_MO_coefficients, new_DM
which can be used for further calculations.

The input density matrix needs to have been calculated using pySCF but does not have been calculated by the same input pySCFTheory object (pyscftheoryobj).
A calculated dm object is often available as an attribute of a PySCFTheory object (e.g. pyscftheoryobj.dm) that has been run.

One can choose between the WY method or the ZMP method (*method* keyword). The accuracy of the ZMP method is controlled by 
*ZMP_lambda*. In case of problems with ZMP convergence it is best to enable level-shifting (*ZMP_levelshift* Boolean keyword) 
as well as increasing the number of ZMP cycles (*ZMP_cycles* keyword).
See also `KS-pies documentation <https://ssnam92.github.io/KSPies/userguide.html#userguide>`_ 


If problems with WY convergence (happens with large basis sets) then changing the optimization method (WY_method) from 'trust-exact' to either 'bfgs' or 'cg'
often works.

**Example: Inverting PBE density**

To make sure KS inversion works correctly it can be useful to first invert a PBE density matrix and then compare the resulting inverted KS potential to the PBE potential
by e.g. recalculating an energy using the new MOs.

Here we will use the H2 + F -> HHF reaction barrier example from `Kanungo et al. <https://pubs.acs.org/doi/full/10.1021/acs.jpclett.3c03088>`_
where KS-inversion was utilized to invert different densities.
Table S1 in the `SI <https://pubs.acs.org/doi/suppl/10.1021/acs.jpclett.3c03088/suppl_file/jz3c03088_si_001.pdf>`_ lists
barrier heights calculated at the PBE level vs. PBE\@PBE-inv level (inverted PBE density).
The PBE/aug-cc-pVQZ barrier height in the table (-12.68 kcal/mol) is well reproduced using PySCFTheory in ASH: -12.676 kcal/mol using the same geometries (see SI).

.. toggle::

  Regular PBE/aug-cc-pVQZ calculation

  .. code-block:: python

    from ash import *

    #Fragments
    hhf_cstring="""
    H 0.14656781 -0.24726460 0.00000000
    F 0.00000000 1.21154849 0.00000000
    H -0.14656781 -0.96428389 0.00000000
    """
    H2 = Fragment(diatomic="H2", bondlength=0.74187686, charge=0, mult=1)
    F = Fragment(atom="F", charge=0, mult=2)
    HHF = Fragment(coordsstring=hhf_cstring, charge=0, mult=2)
    fragments=[H2,F,HHF]
    reaction = Reaction(fragments=fragments, stoichiometry=[-1,-1,1], unit="kcal/mol")

    qm = PySCFTheory(scf_type="UKS", functional="PBE", basis="aug-cc-pVQZ", autostart=False)

    #Looping over fragment and running
    final_energies=[]
    for fragment in fragments:
        res = Singlepoint(fragment=fragment, theory=qm)
        final_energies.append(res.energy)

    print("final_energies:", final_energies)
    ReactionEnergy(list_of_energies=final_energies, stoichiometry=reaction.stoichiometry)


The KS-inverted PBE\@PBE-inv/aug-cc-pVQZ barrier height is also in the table (-12.48 kcal/mol) that reproduces the original the PBE value well. 
The inversion in the paper used a different algorithm than the WY or ZMP algorithms implemented in KS-Pies.
We can test how well we can reproduce the PBE-inversion using either the WY or ZMP algorithm.
One can run a WY-inversion like this:

.. toggle::

  WY-inverted PBE/aug-cc-pVQZ calculation

  .. code-block:: python

    from ash import *

    #Fragments
    hhf_cstring="""
    H 0.14656781 -0.24726460 0.00000000
    F 0.00000000 1.21154849 0.00000000
    H -0.14656781 -0.96428389 0.00000000
    """
    H2 = Fragment(diatomic="H2", bondlength=0.74187686, charge=0, mult=1)
    F = Fragment(atom="F", charge=0, mult=2)
    HHF = Fragment(coordsstring=hhf_cstring, charge=0, mult=2)
    fragments=[H2,F,HHF]
    reaction = Reaction(fragments=fragments, stoichiometry=[-1,-1,1], unit="kcal/mol")

    qm = PySCFTheory(scf_type="UKS", functional="PBE", basis="aug-cc-pVQZ", autostart=False)
    #Looping over fragment and running
    final_energies=[]

    for fragment in fragments:
        #Running DFA calc to get SCF-density
        dummy = Singlepoint(fragment=fragment, theory=qm)
        #Density->Potential inversion using above density
        new_mo_occ, new_mo_energy, new_mo_coeff, new_dm =density_potential_inversion(qm,qm.dm,method='WY')
        #Non-selfconsistent DFA calc using the new DM
        new = PySCFTheory(scf_type="UKS", functional="PBE", basis="aug-cc-pVQZ", autostart=False,scf_noiter=True,dm=new_dm)
        result = Singlepoint(fragment=fragment,theory=new)
        final_energies.append(result.energy)
    #Reaction energy
    ReactionEnergy(list_of_energies=final_energies, stoichiometry=reaction.stoichiometry,label="WY")

The WY algorithm in KS-Pies is found to here give a barrier height of -12.659 kcal/mol which compares really well to the original PBE-value (-12.676 kcal/mol).
And is even better than the inverted value of the paper.

A ZMP inversion is controlled by the *ZMP_lambda* parameter and is systematically improvable by increasing *ZMP_lambda*.
Here we use ZMP_lambda=512 with ZMP_levelshift enabled: 

.. toggle::

  ZMP-inverted PBE/aug-cc-pVQZ calculation

  .. code-block:: python

    from ash import *

    #Fragments
    hhf_cstring="""
    H 0.14656781 -0.24726460 0.00000000
    F 0.00000000 1.21154849 0.00000000
    H -0.14656781 -0.96428389 0.00000000
    """
    H2 = Fragment(diatomic="H2", bondlength=0.74187686, charge=0, mult=1)
    F = Fragment(atom="F", charge=0, mult=2)
    HHF = Fragment(coordsstring=hhf_cstring, charge=0, mult=2)
    fragments=[H2,F,HHF]
    reaction = Reaction(fragments=fragments, stoichiometry=[-1,-1,1], unit="kcal/mol")

    qm = PySCFTheory(scf_type="UKS", functional="PBE", basis="aug-cc-pVQZ", autostart=False)
    #Looping over fragment and running
    final_energies=[]

    for fragment in fragments:
        #Running DFA calc to get SCF-density
        dummy = Singlepoint(fragment=fragment, theory=qm)
        #Density->Potential inversion using above density
        new_mo_occ, new_mo_energy, new_mo_coeff, new_dm =density_potential_inversion(qm,qm.dm,method='ZMP', ZMP_lambda=512, ZMP_levelshift=True)
        #Non-selfconsistent DFA calc using the new DM
        new = PySCFTheory(scf_type="UKS", functional="PBE", basis="aug-cc-pVQZ", autostart=False,scf_noiter=True,dm=new_dm)
        result = Singlepoint(fragment=fragment,theory=new)
        final_energies.append(result.energy)
    #Reaction energy
    ReactionEnergy(list_of_energies=final_energies, stoichiometry=reaction.stoichiometry,label="WY")

Results using the ZMP algorithm depends on the ZMP_lambda parameter:

============  ====================
ZMP_lambda     Energy (kcal/mol)
============  ====================
32              -12.407
64              -12.615
128             -12.668
256             -12.677
512             -12.678
============  ====================

At ZMP_lambda=512 we have seemingly reached convergence at a value of -12.678 kcal/mol which is essentially perfect agreement with the exact result.
The ZMP inversion is generally more expensive than WY for similar accuracy but the advantage is that the result is systematically improvable by increasing ZMP_lambda.


**Example: Inverting a CCSD(T) density for use in a PBE[n_cc] calculation**

The main benefit of KS-inversion is that it allows one to calculate a KS-potential from a correlated density matrix.
If the correlated density matrix is a good approximation to the exact result, this effectively corresponds to solving 
the exact Kohn-Sham problem (in an expensive way).
This can give insight into errors of density functional approximations (DFA) by analyzing the corresponding near-exact MOs or KS-potential.
One can also use it to perform a non-selfconsistent DFA calculation on top of the near-exact KS density matrix (from the inversion). 

Below we show how we can use a similar script as above to do a PBE[n_cc] calculation using the inverted CCSD(T) (unrelaxed) density matrix.
Other WF or KS densities in pySCF could be used in the same way: e.g. from other-DFAs, MP2-unrelaxed, MP2-relaxed, CCSD-unrelaxed etc.
Any density matrix created by pySCF or from programs that use pySCF (Dice, Block2) could be used.

Example below uses a TZ basis instead of QZ to speed up the calculation.

.. toggle::

  WY-inverted CCSD(T) density as input to a PBE\@CCSD(T)-inv calculation

  .. code-block:: python

    from ash import *


    numcores=1
    #Fragments
    hhf_cstring="""
    H 0.14656781 -0.24726460 0.00000000
    F 0.00000000 1.21154849 0.00000000
    H -0.14656781 -0.96428389 0.00000000
    """
    H2 = Fragment(diatomic="H2", bondlength=0.74187686, charge=0, mult=1)
    F = Fragment(atom="F", charge=0, mult=2)
    HHF = Fragment(coordsstring=hhf_cstring, charge=0, mult=2)
    fragments=[H2,F,HHF]
    reaction = Reaction(fragments=fragments, stoichiometry=[-1,-1,1], unit="kcal/mol")

    #Looping over fragment and running
    final_energies=[]

    for fragment in fragments:
        #CCSD(T) density
        dens_theory = PySCFTheory(scf_type="UHF", basis="aug-cc-pVTZ", autostart=False, label="CCSD(T)", CC=True, CCmethod="CCSD(T)",
                CC_density=True,numcores=numcores)
        dummy = Singlepoint(fragment=fragment, theory=dens_theory)
        #Density->Potential inversion using above density
        new_mo_occ, new_mo_energy, new_mo_coeff, new_dm =density_potential_inversion(dens_theory,dens_theory.dm,method='WY')
        #Non-selfconsistent DFA calc using the new DM
        new = PySCFTheory(scf_type="UKS", functional="PBE", basis="aug-cc-pVTZ", autostart=False,
            scf_noiter=True,dm=new_dm,numcores=numcores)
        result = Singlepoint(fragment=fragment,theory=new)
        final_energies.append(result.energy)
    #Reaction energy
    ReactionEnergy(list_of_energies=final_energies, stoichiometry=reaction.stoichiometry,label="WY")

The resulting barrier height is -10.26 kcal/mol (aug-cc-pVQZ basis) which is comparable to the 
-11.78 kcal/mol error in the paper (aug-cc-pV5Z basis), see Table S2 (Forwards,PBE,CCSD(T)).

-----------------------------------------------------------------------
Decomposing DFA error into density-driven and functional-driven error
-----------------------------------------------------------------------

Kohn-Sham inversion of a reference density such as the CCSD(T) density above 
can be used to decompose the error of a DFA into a density-driven error (DE) and a functional-driven error (FE).

The **DFA_error_analysis** function in ASH can be used to easily do this analysis.
It takes as input an ASH fragment, a PySCFTheory object for the DFA calculation, a PySCFTheory object for the reference calculation (e.g. CCSD(T)),
the reference density matrix (e.g. CCSD(T) density matrix), the DFA density matrix (e.g. PBE density matrix), the reference energy (e.g. CCSD(T) energy) and
finally the DFA energy. It will then perform the inversion of the reference density matrix and calculate the corresponding KS-potential.

.. code-block:: python

  def DFA_error_analysis(fragment=None, DFA_obj=None, REF_obj=None, DFA_DM=None, REF_DM=None, REF_E=None, DFA_E=None,
                              inversion_method='WY', WY_method='trust-exact',
                              ZMP_lambda=128, ZMP_levelshift=False, ZMP_cycles=400, DF=True):


Example:

.. toggle::

  .. code-block:: python

    from ash import *

    numcores=1
    #Fragments
    hhf_cstring="""
    H 0.14656781 -0.24726460 0.00000000
    F 0.00000000 1.21154849 0.00000000
    H -0.14656781 -0.96428389 0.00000000
    """
    H2 = Fragment(diatomic="H2", bondlength=0.74187686, charge=0, mult=1)
    F = Fragment(atom="F", charge=0, mult=2)
    HHF = Fragment(coordsstring=hhf_cstring, charge=0, mult=2)
    fragments=[H2,F,HHF]
    reaction = Reaction(fragments=fragments, stoichiometry=[-1,-1,1], unit="kcal/mol")

    #Results lists
    final_energies=[]
    ref_energies=[]
    FEs=[];DEs=[];Errors_tot=[]

    #Looping over fragments
    for fragment in reaction.fragments:
        #WF calculation to get energy and density
        WF_theory = PySCFTheory(scf_type="UHF", basis="aug-cc-pVQZ", CC=True, CCmethod="CCSD(T)", CC_density=True,
                autostart=False, numcores=numcores)
        WF_ref_result = Singlepoint(fragment=fragment, theory=WF_theory)
        ref_energies.append(WF_ref_result.energy)
        #DFT-calc theory
        DFA_theory = PySCFTheory(scf_type="UKS", functional="PBE", basis="aug-cc-pVQZ", autostart=False,numcores=numcores)
        DFA_result = Singlepoint(theory=DFA_theory, fragment=fragment)
        final_energies.append(DFA_result.energy)
        #DFA_error_analysis
        FE,DE = DFA_error_analysis(fragment=fragment,DFA_obj=DFA_theory, REF_obj=WF_theory, DFA_DM=DFA_theory.dm,
                                REF_DM=WF_theory.dm, REF_E=WF_ref_result.energy, DFA_E=DFA_result.energy,
                                inversion_method='WY')
        FEs.append(FE);DEs.append(DE)
        Errors_tot.append(FE+DE)

    print("final_energies:", final_energies)
    print("ref_energies:", ref_energies)
    print("Errors_tot:", Errors_tot)
    print("FEs:", FEs)
    print("DEs:", DEs)

    #Printing reaction energy and respective errors
    ReactionEnergy(list_of_energies=ref_energies, stoichiometry=reaction.stoichiometry, label="deltaE-ref")
    ReactionEnergy(list_of_energies=final_energies, stoichiometry=reaction.stoichiometry, label="deltaE_DFA")
    ReactionEnergy(list_of_energies=Errors_tot, stoichiometry=reaction.stoichiometry, label="deltaError_DFA")
    ReactionEnergy(list_of_energies=FEs, stoichiometry=reaction.stoichiometry, label="deltaFE")
    ReactionEnergy(list_of_energies=DEs, stoichiometry=reaction.stoichiometry, label="deltaDE")


The script will print the CCSD(T) reaction energy (deltaE-ref), the PBE reaction energy (deltaE_DFA), 
the total error (deltaError_DFA), the functional-driven error (deltaFE) and the density-driven error (deltaDE) for the reaction energy.
The result compare well to the paper (Table S1).

=====================  ====================  ====================  ====================
Method                   BH (kcal/mol)          DE (kcal/mol)         FE (kcal/mol)
=====================  ====================  ====================  ====================
CCSD(T)/AQZ (paper)     1.29
PBE/AQZ (paper)         -12.68               -2.36                 -11.61
PBE/AQZ (ASH)           -12.X                 -2.X                  -11.X
=====================  ====================  ====================  ====================
