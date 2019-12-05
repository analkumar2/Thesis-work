import numpy.random as nr

############################
num_of_iterations = 100
errorthreshold = 6
seeed = nr.randint(2**32 - 1)
nr.seed(seeed)
print(seeed)
##############################
f = open("Modelparameters_temp/outputModels_dict_"+str(seeed)+".py","a+")
f.write("# exec(open('Modelparameters_temp/outputModels_dict_"+str(seeed)+".py').read())\n\n")
f.write("Models = {} \n\n")

mnum = 1
for kk in range(500):
    f = open("Modelparameters_temp/outputModels_dict_"+str(seeed)+".py","a+")
    print(kk)
    exec(open('RandomModel_withsave.py').read())
