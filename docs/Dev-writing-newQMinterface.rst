Writing a new Theory interface in ASH
======================================

ASH supports a number of interfaces to various quantum chemistry programs but obviously can not support
every single program out there.
If you are interested in writing an ASH-interface for a specific program then here is a guide on how to do this.

Thanks to the modular nature of ASH, it is very easy to add a new interface. 
If correctly done, the interface will work for all common job-functions in ASH as well as QM/MM coupling.

This page is a brief (not yet very detailed) documentation on how to create a new interface.
Contact the main author for details if you run into problems.

################################
Where should I add the code?
################################

Let's assume you want to write an interface for a QM-program named DummyQM.

- 1 new file needs to be created, it should be called: interface_DummyQM.py and should eventually be stored in ash/interfaces alongside all the other interfaces.
- 1 file should be updated: ash/__init__.py

In order for the new interface to be globally available you want to add a line:

.. code-block:: python

    from .interfaces.interface_DummyQM import DummyQMTheory

somewhere in the middle of the ash/__init__.py file.


#####################################
How should I write the interface?
#####################################

First create the file interface_DummyQM.py .
You then need to create a new class inside the file called DummyQMTheory.

The class should have this structure.

Here we create a subclass of the general Theory class in ASH (this makes sure that the new class will have all necessary attributes and methods to work).

.. code-block:: python

    class DummyQMTheory(Theory):
        def __init__(self, dummyqmdir=None, executable_name=None, filename='example', printlevel=2, label="DummyQM",numcores=1, 
                    input_option=None):

    # Run method. Takes coords, elems etc. arguments and computes E or E+G.
    def run(self, current_coords=None, current_MM_coords=None, MMcharges=None, qm_elems=None, mm_elems=None,
            elems=None, Grad=False, PC=False, numcores=None, restart=False, label=None, charge=None, mult=None):

            #Return energy, gradient and pcgradient

The input options to the **__init__** method of DummyQMTheory are flexible but they must contain keyword arguments *filename*, *printlevel*, *label*, *numcores*.
Other input options are up to you.

 **set_numcores** and **cleanup** methods are optional. **set_numcores** will be inherited from the general Theory class, while the **cleanup** method has the purpose of cleaning up temporary files that ASH may call in some workflows.

The **run** method is the most important method and should contain the options shown as ASH job-functions will try to call the runmethod, via **theory.run(...)**.

If you look inside ash/interfaces/interface_example.py of the ASH repository you can find an (inactive) example of how a new theory class should look like.

A complete ASH interface for a correlated wavefunction theory program (typically lacking gradients) needs only support the evaluation of the energy.
For a semiempirical or DFT program a complete interface should support energy and gradient so that it can be used for geometry optimizations, numerical frequencies, molecular dynamics etc. within ASH.

################################################
Properties beyond energy
################################################

If you would like the interface to also calculate and grab other properties (available after run) then you can do that.
These properties should typically be storied in the self.properties dictionary in this case.
See ORCATheory interface for an example (used by some workflows involving ORCATheory)


################################################
Support IR and Raman intensities
################################################

If the QM-code supports calculation of dipole moments and polarizabilities then it is possible for an ASH NumFreq job
to grab this information during the calculation of each displacement and use it to calculate IR intensities and Raman activities automatically.

In order to support this you need to add the following methods to your TheoryClass:
**get_dipole_moment** (for IR intensities)
**get_polarizability_tensor** (for Raman activities).

These 2 methods would be called by the NumFreq function after each theory.run call.
The **get_dipole_moment** and **get_polarizability_tensor** simply need to add code that grabs the information from the outputfile.
And also document the input options for making the dipole moment and polarizability available during a run.

See ORCATheory, CFourTheory and pySCFTheory for example interfaces that support IR intensites and Raman activities.


################################################
Supporting analytic Hessian
################################################

If the QM-program allows the calculation of an analytic Hessian then this can be added to the interface.
Analytical Hessian calculations use the **AnFreq** function.

If you want to add analytical-frequency option, the theory.run method should request the analytic Hessian to be calculated when the Hessian Boolean keyword (passed to theory.run method)
is True. 

See ORCATheory and CFourTheory as examples of interfaces.

################################################
How should I deal with QM/MM ?
################################################

It is best to first write the interface neglecting any QM/MM options and test it well without it first.
Next, you want to find out whether the QM program supports input of pointcharges and whether it is capable of calculating the pointcharge gradient.

Then write the necessary input-options for writing the pointcharge-coordinates and charges to disk (assuming a I/O based interface), e.g. directly to the same inputfile
or as a separate file read by the QM-program.
If the QM-program is not capable of calculating and printing the pointcharge gradient then your QM/MM interface can only support single-point QM/MM energies.

There are several interfaces in ASH that you can inspect to see how pointcharge handling was performed, e.g. ORCATheory, NWChemTheory, xTBTheory.

################################################
What if the QM-program has a Python API?
################################################

If the QM-program has a Python API then exchanging data via disk is probably not necessary and it is better to exchange data via the Python API instead.
Take a look at the ASH interfaces to pySCF, xtB (runmode='library') in this case as an example.

################################################
How do I make my new interface part of ASH?
################################################

Fork a version of ASH, add the file and necessary code changes (ideally nothing should change in any other files) and make a pull request.
Ask the main author to review.

