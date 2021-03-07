xTB interface
======================================

The xTB interface comes in two forms, a shared-library interface and a file-based interface.
The shared-library interface is recommended as no disk I/O is required while running xTB. ASH and xTB communicate via a Python C-API.
As no files are written to disk, this makes the interface faster than the file-based interface, useful for e.g. fast MD.
The file-based interface writes an XYZ-file to disk, calls an xTB executable which reads the XYZ file, runs the job and writes the output to disk which is then read by ASH.
For regular jobs, e.g. geometry optimizations, the speed-difference between interfaces will probably not matter.

To use either interface is quite simple, when an xTB object is created, charge and multiplicity keywords should be provided
as well as the xtbmethod keyword argument that takes values: "GFN2", "GFN1" for the GFN2-xTB and GFN1-xTB Hamiltonians, respectively.
An optional fragment object can also be associated with the xTB-object (makes only sense for single-point jobs).
An optional runmode argument is also available: runmode='library' or runmode='inputfile'.

The runmode='library' option is used by default and requires the shell environment variable LD_LIBRARY_PATH to include the xtb library dir.
e.g. export LD_LIBRARY_PATH=/path/to/your/xtb_6_2_3/lib64:$LD_LIBRARY_PATH

The runmode='inputfile' option requires an additional xtbdir variable to be set that points to the dir containing the xtb executable, e.g. xtbdir=/path/to/xtb_6_2_3/bin .

.. code-block:: python

    #Create fragment object from XYZ-file
    HF_frag=Fragment(xyzfile='hf.xyz')
    xTBcalc = xTBTheory(charge=0, mult=1, xtbmethod='GFN2')
    #xTBcalc = xTBTheory(charge=0, mult=1, xtbmethod='GFN2', runmode='inputfile', xtbdir='/path/to/xtb_6_2_3/bin')

    #Run a single-point energy job on the fragment associated with the xtb-object
    Singlepoint(theory=xTBcalc, fragment=HF_frag)
    #An Energy+Gradient calculation running on 8 cores
    Singlepoint(theory=xTBcalc, fragment=HF_frag, Grad=True)


**Parallelization**

The xTB parallelization is OpenMP or MKL thread-based and can be controlled via the nprocs keyword.
Currently OMP threads are set equal to nprocs and MKL threads are set equal to 1.
Todo: confirm that this actually works

#################
Troubleshooting:
#################

- If the library-interface is not working, the reason is likely that something is missing from the LD_LIBRARY_PATH environment variable,  Make sure the lib64 dir of xtb is part of the LD_LIBRARY_PATH in the shell from which you are running (or in the jobscript you are submitting). e.g. export LD_LIBRARY_PATH=/path/to/your/xtb_6_2_3/lib64:$LD_LIBRARY_PATH
- Fortran libraries may also be missing for xTB. Make sure to load the necessary libraries (e.g. loading a module or  sourcing the Intel compilervars.sh script)
- If the problem is not resolved, try to load the Ash xtb-interface directly in a script:



.. code-block:: python

    import xtb_interface_library
    test = xtb_interface_library.XTBLibrary()

That should reveal what libraries are not found.