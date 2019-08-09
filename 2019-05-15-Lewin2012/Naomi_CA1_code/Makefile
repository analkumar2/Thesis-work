all: mySim

mySim: main.o get_morphology.o vgated_channels.o synaptic_channels.o ion_transport.o update.o nr.o interneuron_update.o
	g++ -o mySim main.o get_morphology.o vgated_channels.o synaptic_channels.o ion_transport.o update.o nr.o interneuron_update.o 

main.o: main.cpp
	g++ -c main.cpp

get_morphology.o: get_morphology.cpp
	g++ -c get_morphology.cpp

vgated_channels.o: vgated_channels.cpp
	g++ -c vgated_channels.cpp

synaptic_channels.o: synaptic_channels.cpp
	g++ -c synaptic_channels.cpp

ion_transport.o: ion_transport.cpp
	g++ -c ion_transport.cpp

update.o: update.cpp
	g++ -c update.cpp

nr.o: nr.cpp
	g++ -c nr.cpp

interneuron_update.o: interneuron_update.cpp
	g++ -c interneuron_update.cpp

clean:
	rm -f *.o
