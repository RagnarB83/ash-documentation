Surface Scans
======================================

Potential Energy Surfaces can be conveniently scanned in ASH using the **calc_surface function** . The function uses the **geometric** optimization library.
Both unrelaxed and relaxed scans be calculated, using either 1 and 2 reaction coordinates.

.. code-block:: python

    def calc_surface(fragment=None, theory=None, scantype='Unrelaxed', resultfile='surface_results.txt', keepoutputfiles=True, keepmofiles=False,
                    runmode='serial', coordsystem='dlc', maxiter=50, extraconstraints=None, convergence_setting=None, 
                    ActiveRegion=False, actatoms=None, **kwargs):
        """Calculate 1D/2D surface

        Args:
            fragment (ASH fragment, optional): ASH fragment object. Defaults to None.
            theory (ASH theory, optional): ASH theory object. Defaults to None.
            scantype (str, optional): Type of scan: 'Unrelaxed' or 'Relaxed'. Defaults to 'Unrelaxed'.
            resultfile (str, optional): Name of resultfile. Defaults to 'surface_results.txt'.
            runmode (str, optional): Runmode: 'serial' or 'parallel. Defaults to 'serial'.
            coordsystem (str, optional): Coordinate system for geomeTRICOptimizer. Defaults to 'dlc'.
            maxiter (int, optional): Max number of Opt iterations. Defaults to 50.
            extraconstraints (dict, optional): Dictionary of additional constraints for geomeTRICOptimizer. Defaults to None.
            convergence_setting (str, optional): Convergence setting for geomeTRICOptimizer. Defaults to None.
            ActiveRegion (bool,optional): To use activeregion or not in optimization
            actatoms (list,optional): List of active atoms

        Returns:
            [type]: [description]
        """





The calc_surface function takes a fragment object and theory object as input. The type of scan is specified ('Unrelaxed' or 'Relaxed') and
then either 1 or 2 reaction coordinates are specified via keyword arguments: RC1_type, RC1_range and RC1_indices (and RC2 versions if using two reaction coordinates).

- The RC1_type/RC2_type keyword can be: 'bond', 'angle' or 'dihedral'.
- The RC1_indices/RC2_indices keyword defines the atom indices for the bond/angle/dihedral. Note: ASH counts from zero.
- The RC1_range/RC2_range keyword species the start coordinate, end coordinate and the stepsize (Å units for bonds, ° for angles/dihedrals).

The resultfile keyword should be used to specify the name of the file that contains the results of the scan ( format: coord1 coord2 energy).
This file can be used to restart an incomplete/failed scan. If ASH finds this file in the same dir as the script, it will read the data and skip unneeded calculations.
Default name : 'surface_results.txt'


**calc_surface** returns a dictionary of total energies for each surface point. The key is a tuple of coordinate value and the value is the energy, i.e.
(RC1value,RC2value) : energy

**1D scan:**

.. code-block:: python

    surfacedictionary = calc_surface(fragment=frag, theory=ORCAcalc, scantype='Unrelaxed', resultfile='surface_results.txt', 
    runmode='serial', RC1_range=[180,110,-10], RC1_type='angle', RC1_indices=[1,0,2], keepoutputfiles=True)

**2D scan:**

If both RC1 and RC2 keywords are provided then a 2D scan will be calculated.

.. code-block:: python

    surfacedictionary = calc_surface(fragment=frag, theory=ORCAcalc, scantype='Unrelaxed', resultfile='surface_results.txt', runmode='serial',
        RC1_type='bond', RC1_range=[2.0,2.2,0.01], RC1_indices=[[0,1],[0,2]], RC2_range=[180,110,-10], 
        RC2_type='angle', RC2_indices=[1,0,2], keepoutputfiles=True)

NOTE: It is possible to have each chosen reaction coordinate apply to multiple sets of atom indices by specifying a list of lists.
In the 2D scan example above, the RC1_indices keyword (a 'bond' reaction coordinate) will apply to both atoms [0,1] as well as [0,2].
This makes sense when preserving symmetry of a system e.g. the O-H bonds in H2O.

NOTE: Currently the runmode is serial which means that one surface point is run after the other and only the theory level can be parallelized.
A future parallel runmode will become available where X surfacepoints can be run simultaneously using X available cores.

Other options to calc_surface:

- coordsystem  (for geomeTRICOptimizer, default: 'dlc'. Other options: 'hdlc' and 'tric')
- maxiter (for geomeTRICOptimizer,default : 50)
- extraconstraints (for geomeTRICOptimizer, default : None. dictionary of additional constraints. Same syntax as constraints in **geomeTRICOptimizer**)
- convergence_setting (for geomeTRICOptimizer, same syntax as in **geomeTRICOptimizer**)
- keepoutputfiles  (Boolean, keep outputfiles for each point. Default is True. )
- keepmofiles (Boolean, keep MO files for each point in a directory. Default is False.)

Note: See :doc:`Geometry-optimization` for geomeTRICOptimizer-related features.

**Working with a previous scan from collection of XYZ files**

.. code-block:: python

    def calc_surface_fromXYZ(xyzdir=None, theory=None, dimension=None, resultfile=None, scantype='Unrelaxed',runmode='serial',
                            coordsystem='dlc', maxiter=50, extraconstraints=None, convergence_setting=None, numcores=None,
                            RC1_type=None, RC2_type=None, RC1_indices=None, RC2_indices=None, keepoutputfiles=True, keepmofiles=False,
                            read_mofiles=False, mofilesdir=None):
        """Calculate 1D/2D surface from XYZ files

        Args:
            xyzdir (str, optional): Path to directory with XYZ files. Defaults to None.
            theory (ASH theory, optional): ASH theory object. Defaults to None.
            dimension (int, optional): Dimension of surface. Defaults to None.
            resultfile (str, optional): Name of resultfile. Defaults to None.
            scantype (str, optional): Tyep of scan: 'Unrelaxed' or 'Relaxed' Defaults to 'Unrelaxed'.
            runmode (str, optional): Runmode: 'serial' or 'parallel'. Defaults to 'serial'.
            coordsystem (str, optional): Coordinate system for geomeTRICOptimizer. Defaults to 'dlc'.
            maxiter (int, optional): Max number of iterations for geomeTRICOptimizer. Defaults to 50.
            extraconstraints (dict, optional): Dictionary of constraints for geomeTRICOptimizer. Defaults to None.
            convergence_setting (str, optional): Convergence setting for geomeTRICOptimizer. Defaults to None.
            numcores (float, optional): Number of cores. Defaults to None.
            RC1_type (str, optional):  Reaction-coordinate type (bond,angle,dihedral). Defaults to None.
            RC2_type (str, optional): Reaction-coordinate type (bond,angle,dihedral). Defaults to None.
            RC1_indices (list, optional):  List of atom-indices involved for RC1. Defaults to None.
            RC2_indices (list, optional): List of atom-indices involved for RC2. Defaults to None.

        Returns:
            [type]: [description]
        """




If a surface scan has already been performed, it's possible to use the created XYZ-files and calculate single-point energies or optimizations for each surfacepoint with
any level of theory.

We can use the **calc_surface_fromXYZ** function to read in previous XYZ-files (named like this: RC1_2.0-RC2_180.0.xyz for a 2D scan and like this: RC1_2.0.xyz for a 1D scan).
These files should have been created from **calc_surface** already (present in surface_xyzfiles results directory).
By providing a theory level object we can then easily perform single-point calculations for each surface point or alternatively relax the structures employing constraints.
The results is a dictionary like before.

.. code-block:: python

    #Directory of XYZ files. Can be full path or relative path.
    surfacedir = '/users/home/ragnarbj/Fe2S2Cl4/PES/Relaxed-Scan-test1/SP-DLPNOCC/surface_xyzfiles'

    #Calculate surface from collection of XYZ files. Will read old surface-results.txt file if requested (resultfile="surface-results.txt")
    #Unrelaxed single-point job
    surfacedictionary = calc_surface_fromXYZ(xyzdir=surfacedir, scantype='Unrelaxed', theory=ORCAcalc, dimension=2, resultfile='surface_results.txt' )

    #Relaxed optimization job. A geometry optimization with constraints will be done for each point
    #The RC1_type and RC1_indices (and RC2_type and RC2_indices for a 2D scan) also need to be provided
    surfacedictionary = calc_surface_fromXYZ(xyzdir=surfacedir, scantype='Relaxed', theory=ORCAcalc, dimension=2, resultfile='surface_results.txt',
                        coordsystem='dlc', maxiter=50, extraconstraints=None, convergence_setting=None,
                        RC1_type='bond', RC1_indices=[[0,1],[0,2]], RC2_type='angle', RC2_indices=[1,0,2])


Other options:

- keepoutputfiles=True  (outputfile for each point is saved in a directory. Default True)
- keepmofiles=False (Boolean, MO-file for each point is saved in a directory. Default False)
- read_mofiles=False (Boolean: Read MO-files from directory if True. Default False.)
- mofilesdir=path   (Directory path containing MO-files (GBW files if ORCA) )
- ActiveRegion= True/False
- actatoms=list  (list of active atoms if doing relaxed scan)

**Plotting**

The final result of the scan is stored in a dictionary (named 'surfacedictionary' in the examples above) and can be easily
plotted by giving the dictionary as input to plotting functions (based on Matplotlib).
See :doc:`module_plotting`) page.

The dictionary has the format: (coord1,coord2) : energy  for a 2D scan  and (coord1) : energy for a 1D scan
where (coord1,coord2)/(coord1) is a tuple of floats and energy is the total energy as a float.

A dictionary using data from a previous job (stored e.g. in surface_results.txt) can be created via the **read_surfacedict_from_file** function:

.. code-block:: python

    surfacedictionary = read_surfacedict_from_file("surface_results.txt", dimension=1)