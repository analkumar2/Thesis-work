# exec(open('concurrent_futures.py').read())
## executor map function returns in the order the functions were executed

import concurrent.futures
import time

start = time.perf_counter()


def do_something(seconds):
    print(f'Sleeping {seconds} second(s)...')
    time.sleep(seconds)
    return f'Done Sleeping...{seconds}'


with concurrent.futures.ProcessPoolExecutor() as executor:
    print('In the order that they were finished')
    secs = [5, 4, 3, 2, 1]
    results = [executor.submit(do_something,sec) for sec in secs]

    for f in concurrent.futures.as_completed(results):
        # to return as the jobs are finished
        print(f.results())


with concurrent.futures.ProcessPoolExecutor() as executor:
    print('In the order that it was started')
    secs = [5, 4, 3, 2, 1]
    results = executor.map(do_something, secs)

    for result in results:
        print(result)

finish = time.perf_counter()

print(f'Finished in {round(finish-start, 2)} second(s)')
