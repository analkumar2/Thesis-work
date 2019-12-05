# exec(open('num_in_rows.py').read())
## To find out the number of elements between minimum and maximum in each row of a matrix
## Since, loading and uloading takes a lot of time, in this case parrallelization leads to a decline in performance. Lets try sorting in another program

import concurrent.futures
import time
import numpy as np

############## Prepare data ##############
np.random.RandomState(100)
arr = np.random.randint(0, 10, size=[20, 10000000])
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
results = []
for row in data:
    results.append(howmany_within_range(row, minimum=4, maximum=8))

print(time.time()-ttt)
print(results[:5])
########################################

########## With ProcessPoolExecutor ############
def howmany_within_range_wrapper(asslist):
    return howmany_within_range(*asslist)
ttt = time.time()
print('Starting with ProcessPoolExecutor.....')
with concurrent.futures.ProcessPoolExecutor() as executor:
    print('In the order that it was started')
    results = executor.map(howmany_within_range_wrapper, [[row, 4, 8] for row in data])

    # for result in results:
    #     print(result)
print(time.time()-ttt)
print(list(results)[:5])
##################################################
