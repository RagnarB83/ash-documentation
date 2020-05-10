
==================================================
About Ash
==================================================

Ash is a Python-based computational chemistry and QM/MM environment, primarily for molecular calculations in the gas phase,
explicit solution, crystal or protein environment. Can do single-point calculations, geometry optimizations,
molecular dynamics (soon), numerical frequencies using a MM, QM or QM/MM Hamiltonian.
Interfaces to popular free-for-academic QM codes: ORCA, xTB, Psi4, PySCF.
Geometry optimizations and reaction

Requirements:

    - Python3 installation
    - Numpy
    - Julia (for faster versions of some routines)
    - Matplotlib (for some scientific plotting options).

Optional Python modules for specific functionality (can be installed via pip or conda):

    - PySCF (C++/Python quantum chemistry code)
    - Psi4 (C++/Python quantum chemistry code)
    - PyBerny (optimizer)
    - geomeTRIC (optimizer)
    - PyFrame (helper tool for Polarizable Embedding functionality)


We recommend Anaconda (https://www.anaconda.com/distribution/) for a good scientific Python distribution.
Contains Python3, Numpy, SciPy, Matplotlib.


#####################
Features
#####################

**Flexible coordinate input:**
    - coordinate string
    - XYZ file
    - CIF file
    - PDB file
    - Python lists


**Interfaces to various QM codes:**
    - ORCA (inputfile-based). Parallelization via OpenMPI. Flexible input, BS-DFT, pointcharge embedding.
    - xTB (inputfile-based). OpenMP parallelization
    - Psi4 (both as Python library and inputfile-based). OpenMP parallelization.
    - PySCF (as Python library). OpenMP parallelization.

**Parallelization :**
    - Parallelization via Python multiprocessing.
    - Support for simultaneous single-point jobs.
    - Support for simultaneous Numerical-Hessian displacement calculations.

**Single-point electrostic embedding QM/MM with ORCA, xTB and Psi4.**
    - **To do**: PySCF

**Polarizable embedding via Psi4 and CPPE library**
    - **Todo**: enable PySCF

**Molecular Mechanics (MM) via pointcharges and Lennard-Jones potentials**
    - Flexible definition of charges and Lennard-Jones potentials. Either via forcefield inputfile or in script.
    - Both energy and gradient available.
    - Slow Python version and fast Fortran version available (requires compilation of library).
    - Limitation: No bonded MM yet.

**Geometry optimization with multiple optimizers**
    - Knarr, Python LBFGS-optimizer in Cartesian coordinates (credit: Vilhjálmur Ásgeirsson). No internal coordinates but frozen atom support.
    - PyBerny optimizer interface with internal coordinates. Limitation: No frozen atoms or constraints. Todo: Manual frozen-atom feature to be done.
    - geomeTRIC interface: powerful optimizer supporting multiple internal coordinates (TRIC, HDLC, DLC etc.), frozen atoms, constraints.
    - **To do**: DL-FIND interface: powerful optimizer supporting DLC, HDLC internal coordinates, frozen atoms, constraints.

**(Nonbonded) QM/MM Geometry optimization:**
    - Possible with geomeTRIC optimizer currently, only.
    - **Todo**: Knarr-optimizer.

**Numerical frequencies: one-point (forward difference) and two-point (central difference)**
    - Partial Hessian possible
    - Full parallelization.
    - Support for any QM, MM or QM/MM Hamiltonian for which there is an Ash interface.
    - **Todo:** Request analytical Hessian from ORCA.

**Hessian analysis**
    - Diagonalization of Hessian (from Ash or ORCA). Print frequencies and normal modes.
    - **Todo:** projection of translation/rotational modes
    - Normal mode composition analysis in terms of individual atoms, elements and atom groups.
    - Print vibrational densities of states files (with linebroadening)
    - Mode mapping: compare normal modes of 2 Hessians (e.g. with isotope substitution) for similarity
    - Read/write ORCA-style Hessian files
    - Print XYZ-trajectory file for individual modes
    - Thermochemistry printing. **TODO:** finish
    - Write frequency output as pseudo ORCA-outputfile (enables visualization of modes in Chemcraft/Avogadro)

**Molecular dynamics**
    - **To be done**

**Submodules:**
    **molcrys: Automatic Molecular crystal QM/MM**

    - Read-in CIF-file, extract cell information and coordinates of asymmetric unit.
    - Fill-up coordinates of unitcell.
    - Expand unit cell.
    - Create spherical cluster from unitcell (with only whole molecules).
    - Near-automatic fragment indentification.
    - Intelligent reordering of fragments (supports inconsistently ordered CIF-files)
    - Automatic creation of nonbonded MM forcefield (charges and LJ potentials (**Todo**)).
    - Self-consistent QM/MM for charge definition of cluster.
    - QM/MM Geometry optimization of central fragment of cluster to capture solid-state geometrical effects.
    - **Todo:** QM/MM Numerical frequencies of central fragment of cluster.

**solvshell: Multi-shell solvation for redox potentials**

    - Reads snapshots from molecular dynamics trajectory and calculates VIE, VEA, redox pot. or other property.
    - Parallelization over snapshots. Averages over snapshots and finds representative snapshots of trajectory.
    - QM/MM single-points with/without increased QM-region.
    - Bulk correction for aqueous solutions.
    - Automatic procedure for accounting for short-range and long-range polarization effects.
    - Polarizable embedding via Psi4 or PySCF (soon available).


