Random parameter search with All Gbars. Na KDR kinetics from Migliore2018. Except Ca_T is uncoupled from BK and SK
parallel -u 'sleep {1}; for i in $(seq 20); do python3 parametercrawl.py; echo $i; done > printoutputs/{1}.txt' ::: $(seq 20)
parallel -u 'sleep {1}; for i in $(seq 4); do python3 parametercrawl.py; echo $i; done' ::: $(seq 4)