import numpy as np
import MOOSEModel_18 as mm
from glob import glob
from pprint import pprint
import matplotlib.pyplot as plt

filelist = glob('OutputModels/*.py')

best_error = 1000

# for i,file in enumerate(filelist):
#     print(i, end='\r')
#     exec(open(file).read())
#     for model in Models.keys():
#         if Models[model]['ErrorWT']<best_error:
#             best_error = Models[model]['ErrorWT']
#             best_no = model
#             best_file = file
#         mm.plotModel(Models[model])


for i,file in enumerate(filelist):
    print(i, end='\r')
    exec(open(file).read())
    for model in Models.keys():
        if abs(Models[model]['ScoresWT']['offset_1.5e-10'])<2 and abs(Models[model]['ScoresWT']['Absoffset_1.5e-10'])<2 and abs(Models[model]['ScoresWT']['freq_1.5e-10'])<2:
            # best_error = Models[model]['ErrorWT']
            # best_no = model
            # best_file = file
            print(file, model)
            pprint(Models[model])
            t,v,ca = mm.runModel(Models[model], 150e-12)
            plt.plot(t,v)
            t,v,ca = mm.runModel(Models[model], 300e-12)
            plt.plot(t,v)
            plt.xlim(0.8,1.8)
            plt.savefig(f'OutputModels_graphs/{file.split("/")[-1]}_{model}.png')
            plt.clf()


