geomeTRIC interface
======================================

Options to **geometricOptimizer**:

- coordsystem (Default: 'hdlc', other options: 'tric', 'dlc', 'cart', 'prim')
- frozenatoms (default: None, provide list of atoms to be frozen in space)
- constraintsinputfile (default: None, provide name of constraints-inputfile according to geomeTRIC syntax.
- constraints (default: None, provide dictionary of constraint definitions, with or without the value of the constraint. Example: constraints = { 'bond' : [[0,1]]}
- constrainvalue (default: False, Boolean, whether constrain-value is provided in constraints dictionary or not)
- maxiter (default: 50, maximum number of iterations)
- ActiveRegion (default:False, whether to use an active region or not. Requires accompanying actatoms list.
- actatoms (default: None, list of atoms that are active during optimization, all others are frozen)
- convergence_setting (default: 'ORCA'. What type of convergence criteria to use. Valid options are: 'ORCA', 'Chemshell', 'ORCA_TIGHT', 'GAU', 'GAU_TIGHT', 'GAU_VERYTIGHT', 'SuperLoose'.

    - ORCA:    conv_criteria = {'convergence_energy' : 5e-6, 'convergence_grms' : 1e-4, 'convergence_gmax' : 3.0e-4, 'convergence_drms' : 2.0e-3, 'convergence_dmax' : 4.0e-3 }
    - Chemshell:    conv_criteria = {'convergence_energy' : 1e-6, 'convergence_grms' : 3e-4, 'convergence_gmax' : 4.5e-4, 'convergence_drms' : 1.2e-3, 'convergence_dmax' : 1.8e-3 }
    - ORCA_TIGHT:    conv_criteria = {'convergence_energy' : 1e-6, 'convergence_grms' : 3e-5, 'convergence_gmax' : 1.0e-4, 'convergence_drms' : 6.0e-4, 'convergence_dmax' : 1.0e-3 }
    - GAU:    conv_criteria = {'convergence_energy' : 1e-6, 'convergence_grms' : 3e-4, 'convergence_gmax' : 4.5e-4, 'convergence_drms' : 1.2e-3, 'convergence_dmax' : 1.8e-3 }
    - GAU_TIGHT:    conv_criteria = {'convergence_energy' : 1e-6, 'convergence_grms' : 1e-5, 'convergence_gmax' : 1.5e-5, 'convergence_drms' : 4.0e-5, 'convergence_dmax' : 6e-5 }
    - GAU_VERYTIGHT:    conv_criteria = {'convergence_energy' : 1e-6, 'convergence_grms' : 1e-6, 'convergence_gmax' : 2e-6, 'convergence_drms' : 4.0e-6, 'convergence_dmax' : 6e-6 }
    - SuperLoose:            conv_criteria = { 'convergence_energy' : 1e-1, 'convergence_grms' : 1e-1, 'convergence_gmax' : 1e-1, 'convergence_drms' : 1e-1, 'convergence_dmax' : 1e-1 }

- conv_criteria (Default: see ORCA setting above. Optionally provide your own dictionary of settings using syntax above.