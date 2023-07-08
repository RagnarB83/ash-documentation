Electronic structure analysis
======================================

ASH contains some basic functionality for electronic structure analysis that can be useful.


NOTE: Some functions on this page requires you to import the functions_elstructure module:

.. code-block:: python

    from ash.functions.functions_elstructure import *

######################################################
Creating and modifying Gaussian Cube files
######################################################

Functions to read Gaussian Cube file, 

.. code-block:: python

    #Read Cube file
    def read_cube (cubefile):

    #Subtract one Cube-file from another
    def write_cube_diff(cubedict1,cubedict2, name="Default"):

    #Sum of 2 Cube-files
    def write_cube_sum(cubedict1,cubedict2, name="Default"):

    #Product of 2 Cube-files
    def write_cube_product(cubedict1,cubedict2, name="Default"):

    # CRead cubefile. Grabs coords. Calculates density if MO
    def create_density_from_orb (cubefile, denswrite=True, LargePrint=True):

Example of how to use the functions:

.. code-block:: python

    from ash import *

    #Read 2 Cube-files into dictionaries
    mo2=read_cube("hf.mo2a.cube")
    mo3=read_cube("hf.mo3a.cube")

    #Subtract, sum and multiply previously read-in Cube-file dictionaries
    write_cube_diff(mo2,mo3, name="MO_2-3-diff") #Creates file MO_2-3-diff.cube
    write_cube_sum(mo2,mo3, name="MO_2-3-sum") #Creates file MO_2-3-sum.cube
    write_cube_product(mo2,mo3, name="MO_2-3-prod") #Creates file MO_2-3-prod.cube

    # Create density from MO from a Cube-file.
    create_density_from_orb("hf.mo2a.cube") #Creates file hf.mo2a-dens.cube

######################################################
Various analysis tools
######################################################

CM5 charges can be calculated using the **calc_cm5** function. This function requires the atomic numbers (list), coordinates (numpy array) and Hirschfeld charges (list) of the system:

.. code-block:: python

    def calc_cm5(atomicNumbers, coords, hirschfeldcharges):

Functions to calculate J-couplings according to Yamaguchi, Bencini or Noodleman formulas.
All functions requires the energy of the high-spin and broken-symmetry energy.

.. code-block:: python

    #Yamaguchi equation also requires the <S^2> values of the high-spin and BS state.
    def Jcoupling_Yamaguchi(HSenergy,BSenergy,HS_S2,BS_S2):
    #The Bencini equation (strong-interaction limit, i.e. bond-formation) requires the maximum spin of the system.
    def Jcoupling_Bencini(HSenergy,BSenergy,smax):
    #The Noodleman equation (weak-interaction limit) also requires the maximum spin of the system.
    def Jcoupling_Noodleman(HSenergy,BSenergy,smax):



######################################################
NOCV analysis
######################################################

NOCV analysis can be performed in ASH in 2 different ways: **NOCV_density_ORCA** or **NOCV_Multiwfn**

**NOCV_density_ORCA** calls on ORCA to perform the NOCV and ETS-NOCV. 
It is unfortunately limited to closed-shell systems but the advantage is that the ETS-NOCV is performed exactly.

.. code-block:: python

    def NOCV_density_ORCA(fragment_AB=None, fragment_A=None, fragment_B=None, theory=None, griddensity=80,
                                NOCV=True, num_nocv_pairs=5, keep_all_orbital_cube_files=False,
                                make_cube_files=True):


The **NOCV_Multiwfn** function calls on Multiwfn to perform the NOCV and ETS-NOCV.
The advantage is that it can be used for open-shell systems but the disadvantage is that the energy decomposition analysis
is approximate as full ETS method is not performed.

.. code-block:: python

    def NOCV_Multiwfn(fragment_AB=None, fragment_A=None, fragment_B=None, theory=None, gridlevel=2, openshell=False,
                                num_nocv_pairs=5, make_cube_files=True, numcores=1, fockmatrix_approximation="ETS"):

######################################################
Various ORCA-specific analysis tools
######################################################

Read/write Fock matrix from/to ORCA outputfile.

.. code-block:: python

    # Convert Fock matrix into ORCA-format for printing. Returns string
    def get_Fock_matrix_ORCA_format(Fock):
    # Read Fock matrix from ORCA outputfile. Returns 2 numpy arrays (alpha and beta)
    def read_Fock_matrix_from_ORCA(file):
    # Write Fock matrix to disk as a dummy ORCA outputfile. Can be used by Multiwfn 
    def write_Fock_matrix_ORCA_format(outputfile, Fock_a=None,Fock_b=None, openshell=False):


Create difference density for 2 calculations differing in either fragment or theory-level.
Theory level has to be ORCATheory. Difference density is written to disk as a Cube-file.

.. code-block:: python

    #Create difference density for 2 calculations differing in either fragment or theory-level
    def difference_density_ORCA(fragment_A=None, fragment_B=None, theory_A=None, theory_B=None, 
        griddensity=80, cubefilename='difference_density'):