=================================================
PES: PhotoElectron/PhotoEmission Spectrum
=================================================

Workflow to calculate photoelectron/photoemission spectrum for a molecule using a TDDFT-Dysonnorm approach.
The workflow combines an SCF-calculation of the initial state, an SCF+TDDFT calculation of the ionized state to get
ionization energies. Dyson orbital norms are then calculated as approximate intensities.


- ORCA is the only supported QM-code for now.
- Requires `Wfoverlap <https://sharc-md.org/?page_id=309>`_ program
- Plotting option requires installation of Matplotlib.


Potential issues: Wfoverlap binary requires libblas.so.3 and liblapack.so.3 binaries. Make sure a directory containing these
libraries is in your LD_LIBRARY_PATH.

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




##########################################################################################################
Example: FeS\ :sub:`2` :sup:`-`\
##########################################################################################################
This example of the FeS\ :sub:`2` :sup:`-`\ - anion accounts for multiple Finalstate spin-multiplicities as we go from:

Initial state: FeS\ :sub:`2` :sup:`-`\ - S=5/2 to  Final state: FeS\ :sub:`2`\ S=2 and S=3

Also shown is how results with multiple functionals can be calculated at the same time.

.. code-block:: python

    from ash import *
    import sys
    import PES
    settings_ash.init() #initialize

    molecule=Fragment(xyzfile="FeS2-tpssh-opt.xyz")

    functionals=['BP86', 'BLYP', 'TPSS', 'TPSSh', 'B3LYP', 'PBE0', 'BHLYP', 'CAM-B3LYP', 'wB97M-D3BJ', 'HF']
    for functional in functionals:
        joblabel="FeS2min-"+functional
        orcadir='/opt/orca_4.2.1'
        input="! def2-TZVP RIJCOSX def2/J GridX5 Grid5 Finalgrid6 tightscf slowconv " + functional
        blocks="""
        %scf
        maxiter 1500
        directresetfreq 1
        diismaxeq 20
        end

        """

        #Define ORCA theory. Does not need charge/mult keywords.
        ORCAcalc = ORCATheory(orcadir=orcadir, orcasimpleinput=input, orcablocks=blocks, nprocs=4)

        #Calling PhotoElectronSpectrum to get IPs, dysonnorms and MO-spectrum
        IPs, dysonnorms, MOs_a, MOs_b = PES.PhotoElectronSpectrum(theory=ORCAcalc, fragment=molecule, InitialState_charge=-1, Initialstate_mult=6,
                              Ionizedstate_charge=0, Ionizedstate_mult=[5,7], numionstates=30,
                                path_wfoverlap="/home/bjornsson/sharc-master/bin/wfoverlap.x" )

        #Plotting TDDFT-IP spectrum with Dysonnorm-intensities as well as MO-spectrum.
        PES.plot_PES_Spectrum(IPs=IPs, dysonnorms=dysonnorms, mos_alpha=MOs_a, mos_beta=MOs_b,
                              start=-2, finish=10, broadening=0.1, points=10000, plotname=joblabel)

        os.rename("TDDFT-DOS.dat", joblabel+"_TDDFT-DOS.dat")
        os.rename("TDDFT-DOS.stk", joblabel+"_TDDFT-DOS.stk")
        PES.cleanup()
        print("=================================")


A nice table is printed out:

.. code-block:: shell

    -------------------------------------------------------------------------
    FINAL RESULTS
    -------------------------------------------------------------------------
    Initial state:
    State no.    Mult     TotalE (Eh)      State-type
        0       6    -2060.29687303000      SCF

    Final ionized states:
    State no.    Mult     TotalE (Eh)      IE (eV)  Dyson-norm State-type TDDFT Exc.E. (eV)
        0       5    -2060.17646751000      3.276    0.94885        SCF             0.000
        1       5    -2060.16669219030      3.542    0.93627        TDA             0.266
        2       5    -2060.15438116737      3.877    0.63286        TDA             0.601
        3       5    -2060.14129840868      4.233    0.00679        TDA             0.957
        4       5    -2060.14063692088      4.251    0.02222        TDA             0.975
        5       5    -2060.13957119054      4.280    0.61628        TDA             1.004
        6       5    -2060.13832171358      4.314    0.87886        TDA             1.038
        7       5    -2060.12435697115      4.694    0.00113        TDA             1.418
        8       5    -2060.12395272861      4.705    0.28032        TDA             1.429
        9       5    -2060.12185801725      4.762    0.01219        TDA             1.486
       10       5    -2060.11877107418      4.846    0.00003        TDA             1.570
       11       5    -2060.11634561892      4.912    0.01243        TDA             1.636
       12       5    -2060.11590462705      4.924    0.00225        TDA             1.648
       13       5    -2060.11583112841      4.926    0.05664        TDA             1.650
       14       5    -2060.11042897805      5.073    0.03065        TDA             1.797
       15       5    -2060.10917950110      5.107    0.00467        TDA             1.831
       16       5    -2060.10851801330      5.125    0.81624        TDA             1.849
       17       5    -2060.10238087649      5.292    0.05319        TDA             2.016
       18       5    -2060.10102115157      5.329    0.00405        TDA             2.053
       19       5    -2060.09738296868      5.428    0.00923        TDA             2.152
       20       5    -2060.09598649444      5.466    0.00326        TDA             2.190
       21       5    -2060.09367128714      5.529    0.00756        TDA             2.253
       22       5    -2060.09231156222      5.566    0.00653        TDA             2.290
       23       5    -2060.09080484001      5.607    0.00949        TDA             2.331
       24       5    -2060.09076809069      5.608    0.00402        TDA             2.332
       25       5    -2060.08507194575      5.763    0.01869        TDA             2.487
       26       5    -2060.08264649049      5.829    0.01427        TDA             2.553
       27       5    -2060.06949023315      6.187    0.01436        TDA             2.911
       28       5    -2060.06419833075      6.331    0.00118        TDA             3.055
       29       5    -2060.05736295683      6.517    0.07555        TDA             3.241
       30       7    -2060.17162372000      3.408    0.94915        SCF             0.000
       31       7    -2060.15927594775      3.744    0.93597        TDA             0.336
       32       7    -2060.14637693567      4.095    0.93261        TDA             0.687
       33       7    -2060.12476833423      4.683    0.26773        TDA             1.275
       34       7    -2060.12440084100      4.693    0.30968        TDA             1.285
       35       7    -2060.11852094946      4.853    0.61496        TDA             1.445
       36       7    -2060.11705097657      4.893    0.00015        TDA             1.485
       37       7    -2060.11525025978      4.942    0.30531        TDA             1.534
       38       7    -2060.11447852402      4.963    0.00146        TDA             1.555
       39       7    -2060.10429896177      5.240    0.00888        TDA             1.832
       40       7    -2060.10220425041      5.297    0.09174        TDA             1.889
       41       7    -2060.09805157700      5.410    0.00040        TDA             2.002
       42       7    -2060.09441339411      5.509    0.00172        TDA             2.101
       43       7    -2060.09224518410      5.568    0.02113        TDA             2.160
       44       7    -2060.08875399849      5.663    0.03280        TDA             2.255
       45       7    -2060.08787201476      5.687    0.49869        TDA             2.279
       46       7    -2060.08695328171      5.712    0.00422        TDA             2.304
       47       7    -2060.08151438203      5.860    0.02956        TDA             2.452
       48       7    -2060.07890518015      5.931    0.00197        TDA             2.523
       49       7    -2060.07677371946      5.989    0.03448        TDA             2.581
       50       7    -2060.07269454470      6.100    0.02572        TDA             2.692
       51       7    -2060.06953410300      6.186    0.37580        TDA             2.778
       52       7    -2060.06912986045      6.197    0.00396        TDA             2.789
       53       7    -2060.05487112345      6.585    0.03873        TDA             3.177
       54       7    -2060.05420963565      6.603    0.14670        TDA             3.195
       55       7    -2060.04469156121      6.862    0.00065        TDA             3.454
       56       7    -2060.03822368050      7.038    0.01066        TDA             3.630
       57       7    -2060.03579822524      7.104    0.00271        TDA             3.696
       58       7    -2060.01514510618      7.666    0.00638        TDA             4.258
       59       7    -2060.01429987177      7.689    0.00952        TDA             4.281