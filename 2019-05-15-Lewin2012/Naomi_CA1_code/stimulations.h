/********************************************************
 *
 * stimulations.h
 * containing various stimulation protocol functions
 *
 ***********************************************************/

void somatic_stim(Compartment * );

void somatic_stim(Compartment * neuron){

  double istim = -.38; // mA/cm2
  double duration = .1; // ms  (20 for ap train)

  double start_time = 0.0;
  double end_time = .1;
  double stim_time = 2.;
  double h = 0.00002;

  int printdata = 2000.;
  int printval = 100.;
  
  // time loop
  for( double time = start_time; time<end_time; time+=h){

    // compartment loop to update all ionic currents 
    #pragma omp parallel for
    for( int i = 0; i<183; i++ ){
      update_currents(neuron, i, time, h, stim_time, duration, istim);
      //cable_eqn_rk(neuron, i, h);  // update voltage
      neuron[i].Vm = -70.;
    } // end compartment loop

    /* put new voltages in Vm */
    for(int i = 0; i<183; i++) neuron[i].Vm= neuron[i].Vm_new;

    // print to file at intervals
    if(printdata >= printval){
      print_data(neuron, time);
      printdata = 0;
    } printdata++;

  } // end time loop
  
}
