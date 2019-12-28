Parameter search with all Migliore kinetics
Error function is diff of 25pA injection protocols - sum of squares. Hence the errors cannot be compared to previous modesl. All models which have offset of more than 0 are stored.

The free parameters here:
Cm, Rm, Em, ENa, EK, Kinetics of Na, KDR, Gbars of Na, KDR
# The reason Gbar instead of gbars is use dis because then we have to calculate sm_len and sm_diam as CM=0.01

The fixed parameters which should be free and will be explored later:
1. Recheck all the kinetics if they perform correctly
3. Nerenst or GHK for each channel. It does not matter for Ca channel. So, only necessary for K and Na channels
4. Add R type Ca channel
5. There will be 4 channels which allows Ca ions. And, there are 5-6 channels which are modulated by Ca concentration. Thus, for each of these ca modulated channels, we make a 4 parameter function which will determine which channel affects which channel the most in terms of building the local Ca conc.


Make a machine learching algo since there would be ~60 variables here.


The current program will only store the parameters and the errors. At any point it will only store the top 1000 models. The pseudo code is:
1. Initialize the parameters in a dictionary
2. Check the output of the model against one of the experimental cell. In the current run - cell 4 of 61016. Calculate sm_len and sm_diam so that CM = 0.01.
3. Check if there are more than 1000 models already present in the top parameter list. If not then store the current one in the list, otherwise, check if the error is less than any of the previous models. If error is less, pop the last element and insert the current model (use bisect.insort function).
4. Every 1000runs, pickle the file
5. Currently it does not do scipy curve fit but that's on the horizon. Also, The cell size is fixed to 67e-6m len and diam. Which means that CM will not be 0.01. After the model is made, sm_len and sm_diam can be recalculated.

Scripts present:
1. The script plotexp.py is used to read experiment cell
2. RandomModel.py is used to generate a ranodm model.
