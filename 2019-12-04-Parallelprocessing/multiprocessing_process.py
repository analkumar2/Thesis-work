# exec(open('multiprocessing_process.py').read())
## Returning is a headache here. Have to use Queue objects. Printing is fine.

import multiprocessing
import time

start = time.perf_counter()


def do_something(seconds):
    print(f'Sleeping {seconds} second(s)...')
    time.sleep(seconds)
    return f'Done Sleeping...{seconds}'

processes = []
for i in range(10):
    p = multiprocessing.Process(target=do_something, args=[i])
    p.start()
    processes.append(p)

for process in processes:
    process.join()

finish = time.perf_counter()

print(f'Finished in {round(finish-start, 2)} second(s)')
