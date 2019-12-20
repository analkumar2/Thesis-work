
import numpy as np
import time

# Prepare data
np.random.RandomState(100)
arr = np.random.randint(0, 10, size=[2000000, 5])
data = arr.tolist()
# print(data[:5])

###########################################################

# Solution Without Paralleization

def howmany_within_range(row, minimum, maximum):
    """Returns how many numbers lie within `maximum` and `minimum` in a given `row`"""
    count = 0
    for n in row:
        if minimum <= n <= maximum:
            count = count + 1
    return count

#####

ttt = time.time()
results = []
for row in data:
    results.append(howmany_within_range(row, minimum=4, maximum=8))

print(time.time()-ttt)
print(results[:10])
### > [3, 1, 4, 4, 4, 2, 1, 1, 3, 3]

############################################################

# # Parallelizing using Pool.apply()
#
import multiprocessing as mp

# Step 1: Init multiprocessing.Pool()
pool = mp.Pool(5)

# Step 2: `pool.apply` the `howmany_within_range()`
ttt = time.time()
results = [pool.apply(howmany_within_range, args=(row, 4, 8)) for row in data]
print(time.time()-ttt)

# Step 3: Don't forget to close
pool.close()

print(results[:10])
# #> [3, 1, 4, 4, 4, 2, 1, 1, 3, 3]

############################################################

# Parallel processing with Pool.apply_async()

import multiprocessing as mp
pool = mp.Pool(mp.cpu_count()-5)

results = []

# Step 1: Redefine, to accept `i`, the iteration number
def howmany_within_range2(i, row, minimum, maximum):
    """Returns how many numbers lie within `maximum` and `minimum` in a given `row`"""
    count = 0
    for n in row:
        if minimum <= n <= maximum:
            count = count + 1
    return (i, count)


# Step 2: Define callback function to collect the output in `results`
def collect_result(result):
    global results
    results.append(result)


# Step 3: Use loop to parallelize
ttt = time.time()
for i, row in enumerate(data):
    pool.apply_async(howmany_within_range2, args=(i, row, 4, 8), callback=collect_result)

# Step 4: Close Pool and let all the processes complete
# pool.close()
# pool.join()  # postpones the execution of next line of code until all processes in the queue are done.

# Step 5: Sort results [OPTIONAL]
# results.sort(key=lambda x: x[0])
# results_final = [r for i, r in results]
print(time.time()-ttt)

# print(results_final[:10])
#> [3, 1, 4, 4, 4, 2, 1, 1, 3, 3]
