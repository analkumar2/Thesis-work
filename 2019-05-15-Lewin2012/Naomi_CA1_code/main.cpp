#include <iostream>
#include <math.h>
#include <stdio.h>
#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <unistd.h>
//#include <omp.h>
#include "nr.h"
#include "neuron_structures.h"
#include "get_morphology.h"
#include "vgated_channels.h"
#include "synaptic_channels.h"
#include "ion_transport.h"
#include "interneuron_update.h"
#include "update.h"

/*************************************************************************************************************************************************
*   Note: 
*   compile using Makefile; executable created with included Makefile is mySim 
*   nr.cpp contains matrix math functions from numerical recipes
*   neuron_structures.h contains the structure definition for the neuron
*   get_morphology.cpp contains the functions to read the morphology data for cell n123; requires n123_neighborslist.txt and n123_trim.txt
*   vgated_channels.cpp contains voltage gated channel functions
*   synaptic_channels.cpp contains the synaptic channel functions
*   ion_transport.cpp contains the membrane transporter functions
*   interneuron_update.cpp contains functions for a simple interneuron for small network simulations
*   update.cpp contains the function to call all the membrane channels, transporters, leaks update the current at each time point as well as the leak current functions;
*   update.cpp also contains the implicit euler function for updating voltage and the function to write data
*
*   this code can be run with parallel processing using intel omp if installed; makefile must be written accordingly
*****************************************************************************************************************************************************/



using namespace std;

int main(int argc, char * const argv[]){

    
    /**************************************************************************************
    * Data file structure and structure "Compartment" both defined in neuron_structures.h
    **************************************************************************************/
    Datafile filelist;
    Compartment old_neuron[183];   // this is the neuron that will be read from neuron n123 data; new neuron will be made depending on maximum compartment length specified
    Compartment neuron_copy[183];
    Compartment basket[1];
    
    /* if scaling of the gaba_a current is desired; these are passed to update_currents() and can be used there: */
    double A = 1.;  //soma gaba scale
    double B = 1.;  //prox apical gaba scale
    double C = 1.;  //distal apical gaba scale
    
    /* open data files for writing */
    filelist.vmptr = fopen("results/vmdata", "w");
    filelist.inaptr = fopen("results/inadata", "w");
    filelist.mnaptr = fopen("results/mnadata", "w");
    filelist.hnaptr = fopen("results/hnadata", "w");
    filelist.snaptr = fopen("results/snadata", "w");
    filelist.ikdrptr = fopen("results/ikdrdata", "w");
    filelist.mkdrptr = fopen("results/mkdrdata", "w");
    filelist.ikaptr = fopen("results/ikadata", "w");
    filelist.ihptr = fopen("results/ihdata", "w");
    filelist.mihptr = fopen("results/mihdata", "w");
    filelist.sahpptr = fopen("results/sahp", "w");
    filelist.mahpptr = fopen("results/mahp", "w");
    filelist.itotalptr = fopen("results/itotal", "w");
    filelist.ikmptr = fopen("results/ikmdata", "w");
    filelist.icatptr = fopen("results/icatdata", "w");
    filelist.ca_in_ptr = fopen("results/ca_in_data", "w");
    filelist.glutptr = fopen("results/glut_data", "w");
    filelist.nmdaptr = fopen("results/nmda_data", "w");
    filelist.ampaptr = fopen("results/ampa_data", "w");
    filelist.kiptr = fopen("results/ki", "w");
    filelist.koptr = fopen("results/ko", "w");
    filelist.naiptr = fopen("results/nai", "w");
    filelist.naoptr = fopen("results/nao", "w");
    filelist.cliptr= fopen("results/cli", "w");
    filelist.cloptr = fopen("results/clo", "w");
    filelist.kcc2ptr = fopen("results/kcc2", "w");
    filelist.nkcc1ptr = fopen("results/nkcc1", "w");
    filelist.clleak = fopen("results/clleak", "w");
    filelist.nakpump = fopen("results/nakpump", "w");
    filelist.capump = fopen("results/capump", "w");
    filelist.nacax = fopen("results/nacax", "w");
    filelist.gaba = fopen("results/gaba", "w");
    filelist.gaba_cl = fopen("results/gaba_cl", "w");
    filelist.gaba_open = fopen("results/gaba_open", "w");
    filelist.egaba = fopen("results/egaba", "w");
    filelist.gaba_conc = fopen("results/gaba_conc", "w");
    filelist.ecl = fopen("results/ecl", "w");
    filelist.iktotal = fopen("results/iktotal", "w");
	filelist.inatotal = fopen("results/inatotal", "w");
	filelist.icltotal = fopen("results/icltotal", "w");
	filelist.icatotal = fopen("results/icatotal", "w");
	filelist.icarptr = fopen("results/icar", "w");
	filelist.icalptr = fopen("results/ical", "w");
	filelist.bc_vm = fopen("results/bc_vm", "w");
	filelist.bc_ahp = fopen("results/bc_ahp", "w");
	filelist.bc_cai = fopen("results/bc_cai", "w");
	filelist.bc_Isyn = fopen("results/bc_Isyn", "w");
	filelist.bc_ina = fopen("results/bc_ina", "w");
	filelist.bc_ikdr = fopen("results/bc_ikdr", "w");
	filelist.inaleak = fopen("results/inaleak", "w");
	filelist.ikleak = fopen("results/ikleak", "w");
	filelist.na6_0 = fopen("results/na6_0", "w");
	filelist.na6_1 = fopen("results/na6_1", "w");
	filelist.na6_3 = fopen("results/na6_3", "w");
    filelist.na6_4 = fopen("results/na6_4", "w");
	filelist.na6_5 = fopen("results/na6_5", "w");
	filelist.dk_kcc2 = fopen("results/dk_kcc2", "w");
	
    FILE * hvst = fopen("results/hvst", "w");

    
  /* get morphological parameters for neuron n123 and calculate distance from soma for each compartment */
  read_dimensions(old_neuron);
  read_neighborlist(old_neuron);
  calc_distance(old_neuron);
	
	read_dimensions(neuron_copy);
	read_neighborlist(neuron_copy);
	calc_distance(neuron_copy);
	
	/* set max size of compartments: */
	double max_dx = 200.;  // max compartment length, in um
    int count_new = 0;
    int create_new[183];
    
    /* determine how many compartments necessary to represtent neuron, based on maximum size of compartments specified above */
    for(int seg = 0; seg<183; seg++){
        if(neuron_copy[seg].dx > max_dx){
            create_new[seg] = neuron_copy[seg].dx/max_dx;  // number of new compartments to create for each original compartment based on length of original
            //cout << "creating new compartments for seg " << seg << ": " << create_new << endl;
            count_new = count_new + create_new[seg];
        }
        else create_new[seg] = 0;
    }
    
    int new_total_segments = 183+count_new;
	count_new = new_total_segments;
	cout << "new total segments: " << new_total_segments << endl;
    Compartment neuron[new_total_segments];
	
    /* create the neuron, based on morphologic info from n123 with maximum compartment length specified above */
	reset_grid(neuron_copy, neuron, create_new);
	set_resistance(neuron, count_new);  // calculates the nonuniform membrane and axial resistances

	cout << "compartment 1 radius: " << neuron[1].radius << " and length: " << neuron[1].dx << endl;
    
    
  /* set simulation time */  
  double end_time =  2500.; 
  
  /* set time for GABAA stimulus to start */
  double stim_time = 10.;
    
  /* set time step size in ms */  
  double h = 2.e-2;   
    
    double duration = 100.; // this variable is passed to the update_currents function; not currently used but can be used to set different stimulus durations

  /* set interval for printing data */
  int printdata = 100000.;
  int printval = 10.;

  /* set initial conditions */
  for(int i = 0; i<count_new; i++){
	  neuron[i].na_i = 10.;
	  neuron[i].na_o = 140.;
	  neuron[i].k_i = 110.;
	  neuron[i].k_o = 3.;
	  neuron[i].cl_i = 4.5;
	  neuron[i].cl_o = 136.5; 
	  neuron[i].hco3_i = 12.; 
	  neuron[i].hco3_o = 20.; 
	  neuron[i].ca_i = 0.00006; 
      neuron[i].ca_o = 1.0;
	  neuron[i].Vm = -66.;
      neuron[i].cm = 1.5*1.0e-6;  // F/cm2
      neuron[i].lastrelease = -10000000.;
      neuron[i].glut_conc = 0.;
      neuron[i].Icl_gaba = 0.;
	  neuron[i].Icl_gaba_syn = 0.;
	  neuron[i].Ihco3_gaba_syn = 0.;
      neuron[i].Ihco3_gaba = 0.;
      neuron[i].na_i_diff= 0.;
      neuron[i].na_o_diff = 0.;
      neuron[i].k_i_diff = 0.;
      neuron[i].k_o_diff = 0.;
      neuron[i].cl_i_diff = 0.;
      neuron[i].cl_o_diff = 0.;
      neuron[i].ca_i_diff = 0.;
      neuron[i].ca_o_diff = 0.;
      /* set calcium buffer, using f = 100 mM-1 ms-1, b = .1 ms-1, total buffer = .03 mM */
      neuron[i].buffer = .1*(1.562)/(100.*neuron[i].ca_i + .1);  // changed for 30e-3
      /* set potassium buffer to steady state */
      double buffer_max = 50.0; // mM
	  double k1 =  50.*0.0002; // ms-1
	  double k2 =k1/(1+exp((neuron->k_o-15.0)/-1.15));
	  neuron[i].kbuffer = buffer_max*k1/(neuron[i].k_o*k2 + k1);    
  }
	
	/* basket cell initial conditions */
	basket[0].Vm=-67.;
	basket[0].cm = 1.e-6;
	basket[0].glut_conc=0.;
	basket[0].na_i = 10.;
	basket[0].na_o = 140.;
	basket[0].k_i = 110.;
	basket[0].k_o = 3.5;
	basket[0].ca_i = 5.e-5;
	basket[0].ca_o = 1.0;
	
	
	
	set_conductances(neuron, count_new);
	cout << "gmax_nmda, 0: " << neuron[0].gmax_nmda << endl;
	cout << "gmax_ampa, 0: " << neuron[0].gmax_ampa << endl;

  /* begin time loop */
  for(double time = 0.; time<end_time; time+=h){

    /* compartment loops - first calculates membrane currents, next calculates diffusion for each compartment */
    /* these loops can be parallelized */
      
    /* calculate membrane currents */
#pragma omp parallel for 
    for( int i = 0; i<count_new; i++){
      update_currents(&neuron[i], i, time, h, stim_time, duration, basket[0].Vm, A, B, C);  
	}  // end compartment loop 1


    /* if desire to calculate small network with a second cell named "basket" */
    //	basket_update(&basket[0], time, h, neuron[0].Vm);
	  
      
    /* calculate longitudinal diffusion of ions */  
#pragma omp parallel for 
	for(int i = 0; i<count_new; i++) {
		long_diffusion(neuron, i, h);
	}
	for(int i=0; i<count_new; i++) {
		neuron[i].na_i = neuron[i].na_i + neuron[i].na_i_diff;
		neuron[i].na_o = neuron[i].na_o + neuron[i].na_o_diff;
		neuron[i].k_i = neuron[i].k_i + neuron[i].k_i_diff;
		neuron[i].k_o = neuron[i].k_o + neuron[i].k_o_diff;
		neuron[i].cl_i = neuron[i].cl_i + neuron[i].cl_i_diff;
		neuron[i].cl_o = neuron[i].cl_o + neuron[i].cl_o_diff;
		neuron[i].ca_i = neuron[i].ca_i + neuron[i].ca_i_diff;
		neuron[i].ca_o = neuron[i].ca_o + neuron[i].ca_o_diff;

	}  // end compartment loop 2

      /*********************************************************
      * CALCULATE NEW VOLTAGE USING IMPLICIT EULER METHOD 
      * new voltage for each compartment is placed in Vm_new
      * *******************************************************/
	  cable_eqn_ie(neuron, h, count_new);

      
      /* UPDATE VOLTAGE */
	  for(int i = 0; i<count_new; i++) {
		  neuron[i].Vm = neuron[i].Vm_new;
	  }

	
      /* print to file at given intervals, defined by "printval" (printing at shorter intervals during action potentials) */
	  if(neuron[0].Vm > -57. || basket[0].Vm > -57) printval=10.;
	  else printval = 1000.;
      if(printdata >= printval){
          print_data(filelist, neuron, basket, time, count_new);
          printdata = 0;
          cout << "time is: " << time << endl;
      } printdata++;

    
  } // end time loop

  return 0;

}



