/**********************************************************
 *
 * sodium.h
 * header file for sodium current functions
 * taken from Poirazi et al 2003 online supplement
 *
 * created by Naomi Lewin, Oct 3 2008
 * Clancy Lab, Weill Cornell Medical College
 * naomi.lewin@gmail.com
 *
 ***********************************************************/

void na(Compartment *neuron, double t, double dt); // updates gates and Ina
void  nap(double *m, double distance, double dt, double vm); // updates gates and Inap

void na(Compartment *neuron, double t, double dt){

  // Ina = gbar*m2*h*s*(V-Ena) or gbar*m2*h*s*GHK if gbar is converted to permeability
  /* this function returns gbar*m2*h*s; current must be calculated separately */

  double q = 96487./(8.315*(273.15+celsius));
  double tau_m = 0.05; // ms
  double tau_h = 0.5; // ms
  double tau_s = 0.00333*exp(0.0024*(neuron->Vm+60.)*q)/(1.0+exp(0.0012*(neuron->Vm+60.)*q));

  double natt;
  double maxdist = 2000.0; // um, initial estimate for sodium current attenuation
  if(neuron->distance > 0.0) natt = 1.0 - neuron->distance/maxdist;
  else natt = 1.0;

  double minf, hinf;
  if(neuron->distance < 0.1){   // soma
    minf = 1.0/(1.0+exp(-(neuron->Vm+40.)/3.));
    hinf = 1.0/(1.0+exp((neuron->Vm+45.)/3.));
  }
  else {
    minf = 1./(1.+exp(-(neuron->Vm + 44.)/3.));
    hinf = 1./(1.+exp((neuron->Vm + 49.)/3.5));
    tau_h = 1.0;
  }
  double sinf = (1.0+natt*exp((neuron->Vm+60.)/2.))/(1.0+exp((neuron->Vm+60.)/2.0));


  if( t < dt ){
    neuron->ina_m = minf;
    neuron->ina_h = hinf;
    neuron->ina_s = sinf;
  }
  else {
    neuron->ina_m = neuron->ina_m+(1.0-exp(-dt/tau_m))*(minf-neuron->ina_m);
    neuron->ina_h = neuron->ina_h+(1.0-exp(-dt/tau_h))*(hinf-neuron->ina_h);
    neuron->ina_s = neuron->ina_s+(1.0-exp(-dt/tau_s))*(sinf-neuron->ina_s);
  }

  /* sodium current densities:
     soma: 7 mS/cm2
     axon: 100 mS/cm2
     all dendrites: 7 mS/cm2
  */


  double gbar = 7.0*(1.e-3);  // S/cm2
  
  double gna = gbar*pow(neuron->ina_m, 2)*neuron->ina_h*neuron->ina_s; 
  // neuron->Ina = gna*ghk(neuron->na_i, neuron->na_o, neuron->Vm, 1.);
  neuron->Ina = gna*(neuron->Vm - 50.);

}





