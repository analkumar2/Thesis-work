import sys

jarjar = ''
if len(sys.argv)==1:
	with open(sys.argv[0],'r') as file:
		jarjar = file.read()
elif len(sys.argv)==2:
	with open(sys.argv[1],'r') as file:
		jarjar = file.read()
print(jarjar)
print('larlar')
k=input("press close to exit") 
print(k)