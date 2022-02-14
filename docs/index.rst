.. Ash documentation master file, created by
   sphinx-quickstart on Mon Jan 27 14:15:55 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Ash: Documentation
======================================

These are the documentation pages of ASH, a multiscale modelling program.

.. raw:: html

    <div align=center>
   <script id="asciicast-MUrhNGhDx9mAjdqomBppIGWsI" src="https://asciinema.org/a/MUrhNGhDx9mAjdqomBppIGWsI.js" async></script>
    </div>

.. toctree::
   :maxdepth: 2
   :caption: ASH

   About
   setup
   basics
   basic-examples
   ash_program_philosophy
   coordinate-input
   parallelization
   QM-interfaces
   MM-interfaces

.. toctree::
   :maxdepth: 2
   :caption: Jobtypes

   job-types
   singlepoint
   module_freq
   module_dynamics
   neb
   surfacescan

.. toctree::
   :maxdepth: 2
   :caption: Modules

   module_QM-MM
   module_benchmarking
   module_workflows
   module_highlevel_workflows
   module_molcrys
   module_PES
   module_plotting

.. toctree::
   :maxdepth: 2
   :caption: Interfaces

   ORCA-interface
   xTB-interface
   MRCC-interface
   CFour-interface
   Dalton-interface
   PySCF-interface
   Psi4-interface
   crest-interface
   geomeTRIC-interface
   knarr-interface
   OpenMM-interface

.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   Explicit-solvation
   Metalloprotein-I
   Metalloprotein-II
   QM-MM-protein
   workflows-examples
   Highlevel_CC_CBS_workflows
.. toctree::
   :maxdepth: 2
   :caption: Tools

   coordinate-tools

.. role:: underline
    :class: underline

.. warning:: This is Ash version 0.1. Use at your own risk!


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
