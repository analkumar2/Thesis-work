DDSP README file
================

David Sterratt, 4th July 2012
-----------------------------

How to run the simulations interactively
========================================

1. Compile all the ".mod" files in this directory using
   the NEURON utility "nrnivmodl" (on GNU/Linux) or mknrndll (on MS
   Windows).

2. Use NEURON to run "bpap-gui.hoc". A number of windows should appear
   on the screen. This is the interative version of the simulations.

3. Change parameters and options as desired, and then click on "BPAP
   Run". For example, try setting "AMPA synapses scaled", and changing
   the "Number of runs" to 5.

4. The simulation will run. Traces appear in the window at the right,
   and various analysis plots will appear in the centre window. The
   top row refers to voltage, and the bottom row to calcium. The left
   column is peak V or Ca, the middle column is integral of V or Ca
   and the right column is time to peak V or Ca.

5. To save the parameters and data, click on "Save Parameters and
   Data". Click on "Accept" in the window that pops up, and the data
   and parameters will be saved in a file in R format in the datastore
   directory.

6. To plot figures from this file, open up the R program, and type the
   following at the prompt:

   > source("figure-example.R")

   The first time this happens you may be prompted to download a
   package - this is intended.

7. To produce further simulations and plots, repeat steps 3 to 6. You
   will need to copy "figure-example.R" to a file with a new name and
   edit the "dataset" line to match the name of the dataset in the
   datastore directory.

How to reproduce the simulations with the standard morphology
=============================================================

1. Compile the mod files as in step (1) above.

2. Use NEURON to run "bpap-sims.hoc". This will take some time. Data
   files (R text files) will appear in the "datastore" directory.

3. To produced the figures, start R and type:

   > source("figures.R")

   Figures will appear in the "plots" directory.

How to reproduce the scaling simulations
========================================

1. Compile the mod files as above.

2. Use NEURON to run "bpap-sims-scaling.hoc". This will take some
   time. Data files (text files with a .dat suffix) will appear in the
   "datastore" directory.

3. Open MATLAB from the bpap directory and type "scaling_plots"

To run simulations with cells with alternative morphologies (Fig 3B,C)
======================================================================

1. Compile the mod files as above.

2. Use NEURON to run "bpap-sims-cell1.hoc" and
   "bpap-sims-cell2.hoc". This will take some time. Data files (R text
   files) will appear in the "datastore" directory.

3. These can be plotted using figure-example.R


