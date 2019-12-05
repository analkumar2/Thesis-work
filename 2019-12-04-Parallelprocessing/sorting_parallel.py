# exec(open('sorting_parallel.py').read())
## Sort each row of a matrix. In this case too, prallelization is not recommended

import concurrent.futures
import time
import numpy as np

############## Prepare data ##############
np.random.RandomState(100)
arr = np.random.randint(0, 10, size=[4, 20000000])
data = arr.tolist()
#########################################

#########THe function#####################
def howmany_within_range(row, minimum, maximum):
    """Returns how many numbers lie within `maximum` and `minimum` in a given `row`"""
    count = 0
    for n in row:
        if minimum <= n <= maximum:
            count = count + 1
    return count

##########################################

#######without multiprocessing##########
ttt = time.time()
print('Starting without multiprocessing.....')
results = np.sort(data)
print(time.time()-ttt)
print(results[:5])
########################################

########## With ProcessPoolExecutor ############
def howmany_within_range_wrapper(asslist):
    return howmany_within_range(*asslist)

def ssort(row):
    print(time.time(), ' before')
    a= np.sort(row)
    print(time.time(), ' after')
    return a

ttt = time.time()
print('Starting with ProcessPoolExecutor.....')
with concurrent.futures.ProcessPoolExecutor() as executor:
    print('In the order that it was started')
    results = executor.map(ssort, data)

    # for result in results:
    #     print(result)
print(time.time()-ttt)
print(list(results)[:5])
##################################################
