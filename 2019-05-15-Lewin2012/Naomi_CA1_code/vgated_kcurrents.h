/***************************************************************************
 *
 * vgated_kcurrents.h
 * header file for voltage gated potassium currents
 * taken from Poirazi et al 2003 online supplement for CA1 neuron
 *
 * created by Naomi Lewin, Oct 29 2008
 * Clancy Lab, Weill Cornell Medical College
 * naomi.lewin@gmail.com
 *
 *****************************************************************************/

void kdr(Compartment *, double, double);
void ka(Compartment *, double, double);
void km(Compartment *, double, double);
void kdr_jun(Compartment *, double, double);

void km(Compartment * neuron, double t, double dt){
  
  double gbar = 60.*(1.e-3);  // should be double on oblique branches!
  double tadj = pow(2.3, (celsius-23.)/10.);
  double alpha = (1.e-3)*(neuron->Vm+30.)/(1.-exp(-(neuron->Vm+30.)/9.));
  double beta = -(1.e-3)*(neuron->Vm+30.)/(1.-exp((neuron->Vm+30.)/9.));
  double minf = alpha/(alpha+beta);
  double tau = 1./(alpha+beta);

  if(t < dt ) neuron->ikm_m = minf;
  else neuron->ikm_m = neuron->ikm_m + (1.-exp(-dt*tadj/tau))*(minf-neuron->ikm_m);

  neuron->Ikm = (1.e-4)*tadj*gbar*neuron->ikm_m*(neuron->Vm+80.);
  
}

void kdr(Compartment * neuron, double t, double dt){

  // Ikdr = gbar*m2*(Vm-Ek) or gbar*m2*GHK if gbar is converted to permeability
  
  double tau_m; // = 2.2; //ms
  double minf; // = 1./(1.+exp(-(neuron->Vm+42.)/2.));
  double gbar;

  /* somatic minf are modified: */
  if(neuron->distance < 0.1) { 
    minf = 1./(1.+exp(-(neuron->Vm+46.3)/3.));
    tau_m = 3.5;
    gbar = 1.4*(1.e-3); // S/cm2
  }
  else {
    gbar = 0.868*(1.e-3);
    minf = 1./(1.+exp(-(neuron->Vm + 42.)/2.));
    tau_m = 2.2;
  }

  if( t < dt ){
    neuron->ikdr_m = minf;
  }
  else {
    neuron->ikdr_m = neuron->ikdr_m + (1.-exp(-dt/tau_m))*(minf-neuron->ikdr_m);
  }

  /* delayed rectifier densities:
     soma: 1.4 mS/cm2
     axon: 20 mS/cm2
     other: 0.868 mS/cm2
  */
 

  double gkdr = gbar*pow(neuron->ikdr_m, 2);
  // neuron->Ikdr = gkdr*ghk(neuron->k_i, neuron->k_o, neuron->Vm, 1.);
  neuron->Ikdr = gkdr*(neuron->Vm - -80.0);
}

void ka(Compartment * neuron, double t, double dt){
  
  // fast inactivating A-type potassium channel
  // differing kinetics and densities from soma to apical trunk
  // Ia = ga * m * h * (V-Ek)

  // three regions:  proximal distance < 100; distal 100 < distance < 350; distal distance > 350
  double eps;
  double qt; 
  double alpha_m;
  double alpha_h;
  double minf;
  double beta_m;
  double beta_h;
  double hinf;
  double tau_m, tau_h;
  double q = F/(R*T);
  double gbar;

 
  qt = pow(5., (celsius-24.)/10.); 
  
  if(neuron->distance < 100.){  
    eps = -1.5 - 1./(1.+exp((neuron->Vm+40.)/5.));  // mV-1 
    alpha_m = exp((1.e-3)*eps*(neuron->Vm -11.)*q);
    beta_m = exp(0.00055*eps*(neuron->Vm - 11.)*q);
    alpha_h = exp(0.003*(neuron->Vm + 56.)*q);
    beta_h = alpha_h;
    tau_m = beta_m/(0.05*qt*(1.+alpha_m));
    if(tau_m<0.1) tau_m = 0.1;
    tau_h = 0.26*(neuron->Vm+50.);
    if(tau_h < 2.) tau_h = 2.;
    gbar = 7.5*(1.e-3); // S/cm2
  }
  else {
    eps = -1.8 - 1./(1.+exp((neuron->Vm+40.)/5.));
    alpha_m = exp((1.e-3)*eps*(neuron->Vm+1.)*q);
    beta_m = exp(0.00039*eps*(neuron->Vm+1.)*q);
    alpha_h = exp(0.003*(neuron->Vm+56.)*q);
    beta_h = alpha_h;
    tau_m = beta_m/(0.1*qt*(1.+alpha_m));
    if(tau_m < 0.1) tau_m = 0.1;
    tau_h = 0.26*(neuron->Vm+50.);
    if(tau_h < 2.) tau_h = 2.;
    if(neuron->distance <=350.){
      gbar = (1.e-3)*7.488*6.5*neuron->distance/350.;  // S/cm2
    }
    else gbar = (1.e-3)*7.488*6.5;
  } 

  minf = 1./(1.+alpha_m);
  hinf = 1./(1.+alpha_h);

  if( t < dt ) { //set initial steady-state conditions
    neuron->ika_m = minf;
    neuron->ika_h = hinf;
  }
  else {
    neuron->ika_m = neuron->ika_m + (1.-exp(-dt/tau_m))*(minf-neuron->ika_m);
    neuron->ika_h = neuron->ika_h + (1.-exp(-dt/tau_h))*(hinf-neuron->ika_h);
  }

  neuron->Ika = gbar*neuron->ika_m*neuron->ika_h*(neuron->Vm - (-80.));

}

void kdr_jun(Compartment * neuron, double t, double dt){
  // kdr function from jun's program ca3_hippocampus.cpp

  double an, bn, al, bl, tau_n, tau_l, n_inf, l_inf;

  if(t < dt){
    an = 0.03*exp(1.e-3*5.*0.4*(neuron->Vm+32.)*F/(R*T));
    bn = 0.03*exp(-1.e-3*5.*0.6*(neuron->Vm+32.)*F/(R*T));
    al = 0.001*exp(-1.e-3*2.*(neuron->Vm+61.)*F/(R*T));
    bl = 0.001;
    tau_n = 1./(an+bn);
    tau_l = 1./(al+bl);
    n_inf = an*tau_n;
    l_inf = al*tau_l;
    neuron->ikdr_n = n_inf;
    neuron->ikdr_l = l_inf;
  }
  else {
    an = 0.03*exp(1.e-3*5.*0.4*(neuron->Vm+32.)*F/(R*T));
    bn = 0.03*exp(-1.e-3*5.*0.6*(neuron->Vm+32.)*F/(R*T));
    al = 0.001*exp(-1.e-3*2.*(neuron->Vm+61.)*F/(R*T));
    bl = 0.001;
    tau_n = 1./(an+bn);
    tau_l = 1./(al+bl);
    n_inf = an*tau_n;
    l_inf = al*tau_l;

    neuron->ikdr_n = (neuron->ikdr_n + dt*n_inf/tau_n)/(1.+ dt/tau_n);
    neuron->ikdr_l = (neuron->ikdr_l + dt*l_inf/tau_l)/(1.+ dt/tau_l);
  }

  double gbar;
  if(neuron->distance<0.1) gbar = 1.4e-3;  // S/cm2
  else gbar = 0.868e-3;
  neuron->Ikdr = gbar*pow(neuron->ikdr_n, 3)*neuron->ikdr_l*(neuron->Vm + 80.);

}
