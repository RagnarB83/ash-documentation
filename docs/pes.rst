=================================================
PES: PhotoElectron/PhotoEmission Spectrum
=================================================

Workflow to calculate photoelectron/photoemission spectrum for a molecule using a TDDFT-Dysonnorm approach.
The workflow combines an SCF-calculation of the initial state, an SCF+TDDFT calculation of the ionized state to get
ionization energies. Dyson orbital norms are then calculated as approximate intensities.


- ORCA is the only supported QM-code for now.
- Requires `Wfoverlap <https://sharc-md.org/?page_id=309>`_ program
- Plotting option requires installation of Matplotlib.

TODO: Need to fix open-shell problem


######################################################
Example: H\ :sub:`2`\ O
######################################################

.. code-block:: python

    from ash import *
    import sys
    import PES
    settings_ash.init() #initialize

    h2ostring="""
    O        0.222646668      0.000000000     -0.752205128
    H        0.222646668      0.759337000     -0.156162128
    H        0.222646668     -0.759337000     -0.156162128
    """
    h2o=Fragment(coordsstring=h2ostring)

    orcadir='/opt/orca_4.2.1'
    input="! B3LYP def2-SVP Grid5 Finalgrid6 tightscf"
    blocks="""
    %scf
    maxiter 200
    end
    """

    #Define ORCA theory. Does not need charge/mult keywords.
    ORCAcalc = ORCATheory(orcadir=orcadir, orcasimpleinput=input, orcablocks=blocks, nprocs=1)

    #Calling PhotoElectronSpectrum to get IPs, dysonnorms and MO-spectrum
    IPs, dysonnorms, MOs_a, MOs_b = PES.PhotoElectronSpectrum(theory=ORCAcalc, fragment=h2o, InitialState_charge=0, Initialstate_mult=1,
                              Ionizedstate_charge=1, Ionizedstate_mult=2, numionstates=50,
                                path_wfoverlap="/home/bjornsson/sharc-master/bin/wfoverlap.x" )

    #Plotting TDDFT-IP spectrum with Dysonnorm-intensities as well as MO-spectrum.
    PES.plot_PES_Spectrum(IPs=IPs, dysonnorms=dysonnorms, mos_alpha=MOs_a, mos_beta=MOs_b,
                              start=0, finish=60, broadening=0.1, points=10000, plotname='H2O-B3LYP')