
==================================================
About ASH
==================================================

ASH is a Python-based computational chemistry and multiscale modelling program. Can do single-point calculations, geometry optimizations, surface scans,
molecular dynamics and numerical frequencies using a MM, QM or QM/MM Hamiltonians as well as many complex workflows.
Interfaces available to popular free-for-academic QM codes: ORCA, xTB, Psi4, PySCF, Dalton, MRCC, CFour. Reaction profiles and saddlepoint optimizations
can be performed using the nudged elastic band method (NEB).

The code is open-source and available on `Github <https://github.com/RagnarB83/ash>`_

Strict requirements:

- Python3 installation (version 3.6 or higher)
- Numpy Python library

Strongly recommended (and mostly necessary):

- geomeTRIC (optimizer). Python library, easily installed via pip. Required for geometry optimizations.
- OpenMM (molecular mechanics library). Only needed for MM and QM/MM. Installed via conda.
- Julia installation. Only required for faster Julia versions of routines in molcrys-QM/MM code.
- Matplotlib. Only for plotting options.


Optional Python modules for specific functionality (can be installed via pip or conda):

- PySCF (C++/Python quantum chemistry code)
- Psi4 (C++/Python quantum chemistry code)
- PyBerny (optimizer)
- PyFrame (helper tool for Polarizable Embedding functionality)
- Scipy package
- mdtraj (MD trajectory analysis)

Optional external QM codes (separate installations):

- ORCA
- Dalton
- CFour
- MRCC
- xTB

#####################
Features
#####################

**Flexible coordinate input:**

- coordinate string
- XYZ file
- CIF file
- Fractional coordinate XTL file
- PDB file
- Amber CRD file
- Chemshell fragment file
- GROMACS .gro file
- Python lists
- ASH fragment file format (.ygg)


**Interfaces to various QM codes:**

- ORCA (inputfile-based). Parallelization via OpenMPI. Flexible input, BS-DFT, pointcharge embedding.
- xTB (both as Python library and inputfile-based). OpenMP parallelization
- Psi4 (both as Python library and inputfile-based). Threaded parallelization.
- PySCF (as Python library). OpenMP parallelization.
- CFour
- MRCC
- Dalton

**Parallelization :**

- Parallelization via Python multiprocessing: multiple jobs and numerical frequencies.
- QM code parallelization also available.
- Support for simultaneous single-point jobs.
- Support for simultaneous Numerical-Hessian displacement calculations.

**Single-point electrostic embedding QM/MM with ORCA, xTB and Psi4.**

- **To do**: PySCF

**Polarizable embedding via Psi4, PySCF and CPPE library**


**Nonbonded Molecular Mechanics (MM) via pointcharges and Lennard-Jones potentials**

- Flexible definition of charges and Lennard-Jones potentials. Either via forcefield inputfile or in script.
- Both energy and gradient available.
- Slow Python version and fast Julia version available.
- Limitation: No bonded MM yet.

**Full Molecular Mechanics (MM) via OpenMM interface**

- See :doc:`OpenMM-interface`

**Geometry optimization with internal coordinates**

- geomeTRIC interface: powerful optimizer supporting multiple internal coordinates (TRIC, HDLC, DLC etc.), frozen atoms, constraints.


**QM/MM Geometry optimization:**

- Possible with geomeTRIC optimizer currently, only. Support for HDLC coordinates.

**Numerical frequencies: one-point (forward difference) and two-point (central difference)**

- Partial Hessian possible
- Full parallelization.
- Support for any QM, MM or QM/MM Hamiltonian for which there is an ASH interface.
- Also possible to request analytical Hessian from ORCA.

**Hessian analysis**

- Diagonalization of Hessian (from ASH or ORCA). Print frequencies and normal modes.
- **Todo:** projection of translation/rotational modes
- Normal mode composition analysis in terms of individual atoms, elements and atom groups.
- Print vibrational densities of states files (with linebroadening)
- Mode mapping: compare normal modes of 2 Hessians (e.g. with isotope substitution) for similarity
- Read/write ORCA-style Hessian files
- Print XYZ-trajectory file for individual modes
- Thermochemistry according to RRHO approximation.
- Write frequency output as pseudo ORCA-outputfile (enables visualization of modes in Chemcraft/Avogadro)

**Molecular dynamics**

- via OpenMM library: MM-only, QM/MM and QM-only
- via ASE library

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

**solvshell: Multi-shell solvation protocol for redox potentials, spectroscopy etc.**

- Reads snapshots from molecular dynamics trajectory and calculates VIE, VEA, redox pot. or other property.
- Parallelization over snapshots. Averages over snapshots and finds representative snapshots of trajectory.
- QM/MM single-points with/without increased QM-region.
- Bulk correction for aqueous solutions.
- Automatic procedure for accounting for short-range and long-range polarization effects.
- Polarizable embedding via Psi4 or PySCF (soon available).


