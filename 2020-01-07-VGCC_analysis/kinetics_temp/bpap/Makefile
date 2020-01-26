all:
	nrnivmodl

clean: 
	rm -Rf i686 x86_64 *.eps *.o *.dll *.stackdump

check-units: 
	for f in *.mod ; do echo $$f ;  modlunit $$f ; done	
