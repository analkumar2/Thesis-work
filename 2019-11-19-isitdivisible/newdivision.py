n = 72763
n2 = 2*n
n3 = 3*n
k = 7276*3+1
N = 17066197835
for i in range(int(1e8)):
    while N != n and N != n2 and N != n3:
        N = int(str(N)[:-1]) + int(str(N)[-1])*k
