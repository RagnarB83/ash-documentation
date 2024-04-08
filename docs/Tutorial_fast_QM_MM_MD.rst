Tutorial: Running fast QM/MM MD simulations in ASH
=========================================================

In many research projects the question sometimes arises whether it might be feasible to
perform a long enough MD simulation (and perhaps enhanced sampling MD) of the system at some kind of quantum level of theory. 
Often such simulations are simply unfeasible at the direct QM-level, however, hybrid QM/MM methodology
can make such simulations possible, as one can focus the expensive QM calculation on a small important region of the system.

However, what is the best way to perform such a QM/MM MD simulation? 

- What program to use ? 
- What QM theory level? 
- Is fast communication between QM and MM program important ?
- Does the MM calculation need to be parallelized ?
- How to think about CPU parallelization and scaling?
- What about using GPUs instead of CPUs ?

This tutorial intends to give insight into how to choose a good QM/MM MD protocol
and to demonstrate how ASH is ideally suited for performing such simulations due to
the flexibility offered by the general QM/MM approch available and the many QM-code interfaces available.
The flexibility offered by ASH means that it is easy to switch from running classical MM (OpenMM), 
semi-empirical QM/MM (MNDO, xTB), Gaussian-basis DFT/MM (ORCA, pySCF, NWChem etc.), GPU-based codes (TeraChem, QUICK, pyscf), 
mixed Gaussian-planewave/MM (CP2K), and even post-HF/MM (ORCA,CFour, MRCC) etc.

This tutorial is not yet ready...