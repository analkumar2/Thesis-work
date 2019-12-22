import numpy.random as nr
import time

############################
num_of_iterations = 200
errorthreshold = 6
seeed = nr.randint(2**32 - 1)
# seeed = 110043517
nr.seed(seeed)
print(seeed)
##############################
f = open("Modelparameters_temp/outputModels_dict_"+str(seeed)+".py","a+")
f.write("# exec(open('Modelparameters_temp/outputModels_dict_"+str(seeed)+".py').read())\n\n")
f.write("Models = {} \n\n")

ttt = time.time()
mnum = 1
for kk in range(500):
    f = open("Modelparameters_temp/outputModels_dict_"+str(seeed)+".py","a+")
    print(kk)
    exec(open('RandomModel_withsave.py').read())
    print(time.time() - ttt)
#
#
# # for model in Models.keys():
# #     MOOSEModel_2.plotModel(Models[model], 150e-12)


# ### So, parallelization does not work for MOOSE
# import concurrent.futures
# import time
# import numpy as np
# import numpy.random as nr
#
# def ransd(seeed):
#     num_of_iterations = 10
#     errorthreshold = 10
#     # seeed = nr.randint(2**32 - 1)
#     # seeed = 110043517
#     nr.seed(seeed)
#     print(seeed)
#     ##############################
#     f = open("Modelparameters_temp/outputModels_dict_"+str(seeed)+".py","a+")
#     f.write("# exec(open('Modelparameters_temp/outputModels_dict_"+str(seeed)+".py').read())\n\n")
#     f.write("Models = {} \n\n")
#
#     ttt = time.time()
#     mnum = 1
#     f = open("Modelparameters_temp/outputModels_dict_"+str(seeed)+".py","a+")
#     exec(open('RandomModel_withsave.py').read())
#     print(time.time() - ttt)
#
# ttt = time.time()
# print('Starting with ProcessPoolExecutor.....')
# with concurrent.futures.ProcessPoolExecutor() as executor:
#     print('In the order that it was started')
#     executor.map(ransd, nr.randint(2**32 - 1, size=5))
#
# print(time.time()-ttt)
