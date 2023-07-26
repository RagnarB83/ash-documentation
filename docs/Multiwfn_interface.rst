Multiwfn interface
======================================

Mulitiwfn is an open-source quantum chemistry code full of wavefunction and density analysis features.

ASH features a simple wrapper interface that allows you to quickly use some Multiwfn functionality in an
ASH script to e.g. create density or MO Cube-file.


.. code-block:: python

    def multiwfn_run(moldenfile, fchkfile=None, multiwfndir=None, option='density', 
    mrccoutputfile=None, mrccdensityfile=None,  
    grid=3, numcores=1, fragmentfiles=None, fockfile=None, openshell=False):

The **multiwfn_run**  function requires  a Molden file as input and can also read an FChk file.

The grid parameter controls the density grid used for the calculation. The default is 3, which is a fine grid.
The type of Multiwfn job to perform is controlled by the option parameter. 
The default option is 'density', for creating a Cubefile of the electron density.
Other options are: 
'nocv' for NOCV analysis (see NOCV_Multiwfn in :doc:`elstructure_analysis` for more information).
'hirshfeld' for Hirshfeld population analysis.




For reading a Molden-file produced by MRCC one should use the **mrccoutputfile** and **mrccdensityfile** options.