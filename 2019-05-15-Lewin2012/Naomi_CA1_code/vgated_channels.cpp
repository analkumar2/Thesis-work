/*
 *  vgated_channels.cpp
 *  
 *
 *  Created by Naomi Lewin on 11/14/08.
 *  Clancy lab, Weill Cornell Med
 *  naomi.lewin@gmail.com
 *
 */

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <cstdlib>
#include <cstdio>
#include <iomanip> 
#include "vgated_channels.h"
#include "synaptic_channels.h"
#include "neuron_structures.h"

using namespace std;

double ghk(double conc_in, double conc_out, double vm, double z){

  /* convert from S/cm2 to cm/s conductance to permeability
     assumes max conductance was calculated with no ion gradient Eion = 0 mV and 100 mM on each side
     based on typical experiments, so pmax = gmax*RT/(z2*F2*ion_conc)
  */
  vm = vm*1.e-3;  // convert mV to V
  double gtop = (1.0e3)*R*T/(pow(z,2)*pow(F,2)*0.1);
  double alpha = z*F/(R*T);
  double current = gtop*vm*z*alpha*F*(conc_in-conc_out*exp(-alpha*vm))/(1.-exp(-alpha*vm));  // nA/cm2
  /* !!! note this assumes gmax = 1 S/cm2; must be multiplied appropriately by true gmax in every function !!!! */

  return current;  // nA/cm2

}
 
void na(Compartment *neuron, double t, double dt, int i){

  // Ina = gbar*m2*h*s*(V-Ena) or gbar*m2*h*s*GHK if gbar is converted to permeability
  // THIS SODIUM CURRENT IS NOT USED -- 6 STATE FUNCTION IS USED INSTEAD


  double q = 96487./(8.315*(273.15+celsius));
  double tau_m = 0.05; // ms
	double tau_h =  0.5; // ms

  double natt;
	double maxdist = 350.0; // um, initial estimate for sodium current attenuation
   natt = 1.0 - neuron->distance/maxdist;
  if(natt < 0.) natt = 0.;
  //else natt = 1.0;

	/* dendritic slow, activity-dep inactivation, from Gasparini et al 2004 in Neuron ModelDB */
	double zs =12.;
	double vhalfs = -60.;
	double gms = 0.2;
	double a0s = 0.0003;
	double alphas = exp((1.e-3)*zs*(neuron->Vm-vhalfs)/rtf);
	double betas = exp((1.e-3)*zs*gms*(neuron->Vm-vhalfs)/rtf);
	double vvh = -62.; //-58.;
	double vvs = 2.;
	double c = 1./(1.+exp((neuron->Vm-vvh)/vvs));
	double tau_s = betas/(a0s*(1.+alphas));
	double tau_min = 10.;
	if(tau_s < tau_min) tau_s = tau_min;
	
  double minf, hinf;
	double am, bm, ah, bh;

	
	am = 0.32*(neuron->Vm+10.)/(1.-exp(-(neuron->Vm+10.)/30.)); //0.32*(neuron->Vm+47.)/(1.-exp(-0.1*(neuron->Vm+47.)));
	bm = 0.08*exp(-neuron->Vm/10.5);
	ah = 0.135*exp((65.+neuron->Vm)/-13.);   // vhalf = 82?
	bh = (3.56*exp(0.079*(neuron->Vm+0.))+exp(0.35*(neuron->Vm+0.))); // (3.56*exp(0.079*neuron->Vm)+exp(0.35*neuron->Vm));
	tau_m = 1./(am+bm);
	tau_h = 1./(ah+bh);
	minf = am*tau_m;
	hinf = ah*tau_h;
//    minf = 1.0/(1.0+exp(-(neuron->Vm+ 45.)/9.)); 
//	hinf = 1.0/(1.0+exp((neuron->Vm+62.)/6.9)); 

	double sinf = c + natt*(1.-c);


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


	double gbar; 
	if(neuron->distance < 1. ) gbar = .007; //7.0*(1.e-3); ///100.; //1.53; 
	else if( ( i<=24 && i>=23) ) gbar = .05; //.1; //.1;   high density of Na+ channels in AIS
	else if( i > 64 ) gbar = 0.007; //7; //natt*(8.*(1.e-3)); ///100.;  // S/cm2
	else  gbar = 0.; // 0.004; //.008;
	
			if(neuron->ina_h < 0.0) neuron->ina_h = 0.0;
  
	double gna = gbar*pow(neuron->ina_m, 3)*neuron->ina_h; //*neuron->ina_s; 
  // neuron->Ina = gna*ghk(neuron->na_i, neuron->na_o, neuron->Vm, 1.);
	neuron->Ina = gna*(neuron->Vm - neuron->erev_na);

}

void km(Compartment * neuron, double t, double dt, int i){
  
	
	// Storm et al estimate for activation time constant ~ 50 ms; faster in Chen and Johnston 2004
	double gbar = .003; //.0003; //60.*(1.e-3); ///100.;  // should be double on oblique branches!
  double tadj = pow(2.3, (celsius-23.)/10.);
	
	double alpha = exp(-(neuron->Vm+6.)/8.);
	if(neuron->distance > 150) alpha = exp(-(neuron->Vm+6.)/8.);
	double minf = 1./(1.+alpha);
	double tau = 30.; //50.; 

  if(t < dt ) neuron->ikm_m = minf;
  else neuron->ikm_m = neuron->ikm_m + (1.-exp(-dt*tadj/tau))*(minf-neuron->ikm_m);

	if(i >  64 ) gbar = gbar;
	if( (i <= 39 && i >= 35 ) || ( i<=27 && i>=23) ) gbar = gbar*3.;
  neuron->Ikm = tadj*gbar*neuron->ikm_m*(neuron->Vm - neuron->erev_k);
  
}

void kdr(Compartment * neuron, double t, double dt, int i){

  // Ikdr = gbar*m2*(Vm-Ek) or gbar*m2*GHK if gbar is converted to permeability
  
  double tau_m; // = 2.2; //ms
  double minf; // = 1./(1.+exp(-(neuron->Vm+42.)/2.));
  double gbar;


  /* somatic minf are modified: */
  if(neuron->distance < 10.) { 
	  minf = 1./(1.+exp(-(neuron->Vm+13)/8.)); //+46.3)/3.));
	  tau_m = 3.5;
	  gbar = .03; //1.6*.05; //.05; //50.*(1.e-3); ///100.; // S/cm2
  }
 else {
	 gbar =  .01; //3.*50.e-3; //20.e-3; //10.*(1.e-3); ///100.;
	  minf = 1./(1.+exp(-(neuron->Vm + 13)/8.)); //42.)/2.));
	 tau_m = 2.2;
  }
			
	/* axonal */
	if( (i <= 39 && i >= 35 ) || ( i<=27 && i>=23) ){
		gbar = .03; //1.6*.05; //.05; //.11;
		tau_m = 3.5;
		minf =  1./(1.+exp(-(neuron->Vm+13)/8.));
	}
  
	if( t < dt ){
    neuron->ikdr_m = minf;
  }
  else {
    neuron->ikdr_m = neuron->ikdr_m + (1.-exp(-dt/tau_m))*(minf-neuron->ikdr_m);
  }

	/* from Gasparini, Migliore, and Magee 2004: in NeuronModelDB */
//	double gbar = 0.008;
	double vhalf = 13.;
	double a0n = 0.02;
	double zn = -3.;
	double gmn = 0.7;
	double nmax = 2.;
	double q10 = 2.;
	
	double alphan = exp((1.e-3)*zn*(neuron->Vm-vhalf)/rtf);
	double betan = exp((1.e-3)*zn*gmn*(neuron->Vm-vhalf)/rtf);
	double qt = pow(q10, (celsius-24.)/10.);
	double ninf = 1./(1.+alphan);
	double taun = betan/(qt*a0n*(1.+alphan));
	if (taun < nmax) taun = nmax;
	
	if( t<dt) neuron->ikdr_m = ninf;
	else neuron->ikdr_m = neuron->ikdr_m + dt*(ninf - neuron->ikdr_m)/taun;
		
	
	
  /* delayed rectifier densities:
     soma: 1.4 mS/cm2
     axon: 20 mS/cm2
     other: 0.868 mS/cm2
  */
 

	double gkdr = gbar*neuron->ikdr_m*neuron->ikdr_m;
  // neuron->Ikdr = gkdr*ghk(neuron->k_i, neuron->k_o, neuron->Vm, 1.);
  neuron->Ikdr = gkdr*(neuron->Vm - neuron->erev_k);
}



void ka(Compartment * neuron, double t, double dt, int i){
  
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
  
  if(neuron->distance < 50.){  
    eps = -1.5 - 1./(1.+exp((neuron->Vm+40.)/5.));  // mV-1 
//    alpha_m = exp((1.e-3)*eps*(neuron->Vm -11.)*q);
	  alpha_m = exp((-11.-neuron->Vm)/8.);  // hoffman et al 1997
    beta_m = exp(0.00055*eps*(neuron->Vm - 11.)*q);
//    alpha_h = exp(0.003*(neuron->Vm + 56.)*q);
	  alpha_h = 1./exp((-64.-neuron->Vm)/9.);  // chen and johnston 2004
//	  alpha_h = 1./exp((-56.-neuron->Vm)/8.);
	beta_h = alpha_h;
    tau_m = beta_m/(0.05*qt*(1.+alpha_m));
    if(tau_m<0.01) tau_m = 0.01;
//    tau_h = 0.26*(neuron->Vm+50.);
//    if(tau_h < 2.) tau_h = 2.;
	  tau_h = 23.; // from chen and johnston 2004
	  gbar =.0001; //.5*.0008; //.0025; ///100.; // S/cm2  // gmax set to match current densities to somatic and dendritic patches from Hoffman et al. Nature 1997
  }
  else {
    eps = -1.8 - 1./(1.+exp((neuron->Vm+40.)/5.));
//    alpha_m = exp((1.e-3)*eps*(neuron->Vm+1.)*q);
	  alpha_m = exp((1.-neuron->Vm)/10.);  // hoffman et al 1997
	  beta_m = exp(0.00039*eps*(neuron->Vm+1.)*q);
//    alpha_h = exp(0.003*(neuron->Vm+56.)*q);
	  alpha_h = 1./exp((-64.-neuron->Vm)/9.);
//	  alpha_h = 1./exp((-56.-neuron->Vm)/8.);
    beta_h = alpha_h;
    tau_m = beta_m/(0.1*qt*(1.+alpha_m));
    if(tau_m < 0.01) tau_m = 0.01;
	  tau_m = 5.*tau_m; //3.5*tau_m;
//	  tau_h = 0.26*(neuron->Vm+50.);
//    if(tau_h < 2.) tau_h = 2.;
	  tau_h = 23.; //from chen and johnston 2004
	  gbar = 20.*.06*6.*neuron->distance/350.; //.004*6.*neuron->distance/350.;  // S/cm2

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

  if( (i <= 39 && i >= 35 ) || ( i<=27 && i>=23) ) gbar = 0.;
	neuron->Ika = gbar*neuron->ika_m*neuron->ika_m*neuron->ika_h*(neuron->Vm-neuron->erev_k);
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
  neuron->Ikdr = gbar*pow(neuron->ikdr_n, 3)*neuron->ikdr_l*(neuron->Vm - neuron->erev_k);

}


void sahp(Compartment * neuron, double t, double dt, int i){
  
  double gbar;
	double ginit = .0025; //2.*0.1*(1.e-3); ///100.; // S/cm2
	if(neuron->distance < 0.1) gbar =1.*ginit; // S/cm2
  else if(neuron->distance < 50.) gbar = 1.*ginit;
  else gbar = 0.; //ginit;
	   if( (i <= 39 && i >= 35 ) || ( i<=27 && i>=23) ) gbar = ginit;

	double cac = pow(neuron->ca_i/(.01), 2); //0.025), 2);  
	double tau = 250.; //1./(0.003*(1.+cac)*pow(3., (celsius-22.)/10.));
 // if (tau>0.5) tau = 0.5;

  double minf = cac/(1.+cac);
  if(t<dt) neuron->isahp_m = minf;
  else neuron->isahp_m = neuron->isahp_m + dt*(minf-neuron->isahp_m)/tau;

  neuron->Isahp = gbar*pow(neuron->isahp_m, 2)*(neuron->Vm - neuron->erev_k);
}

void sahp_markaki(Compartment * neuron, double t, double dt, int ){

	double gbar = 0.0025; // S/cm2
	//if(neuron->distance < 10.) gbar = 0.;
	double taumin = 150.; //100.; // ms
	double b = 0.005; // mM
	double a = neuron->ca_i/b;
	double minf = a/(a+1.);
	
	double tau = taumin + 1./(neuron->ca_i + b);
	if(t<dt) neuron->isahp_m = minf;
	else neuron->isahp_m = neuron->isahp_m + dt*(minf-neuron->isahp_m)/tau;
	
	neuron->Isahp = gbar*neuron->isahp_m*neuron->isahp_m*(neuron->Vm-neuron->erev_k);
}

void mahp(Compartment * neuron, double t, double dt, int i){

  double gbar;
	double ginit = .00165; // S/cm2
	if( neuron->distance < 0.1) gbar = 2.*ginit; //.9075*(1.e-3); ///100.;
  else if(neuron->distance < 50.) gbar = 2.*ginit;
  else gbar = 0.25*ginit;
		    if( (i <= 39 && i >= 35 ) || ( i<=27 && i>=23) ) gbar = 0.;

  double alpha = 0.48/(1.+(0.18/neuron->ca_i)*exp(-1.68*neuron->Vm*F/(R*T)));
  double beta = 0.28/(1.+(neuron->ca_i/(0.011*exp(-2.*neuron->Vm*F/(R*T)))));
  double tau = 1./(alpha+beta);
  double minf = alpha/tau;

  if(t<dt) neuron->imahp_m = minf;
  else neuron->imahp_m = neuron->imahp_m + (1.-exp(-dt/tau))*(minf-neuron->imahp_m);
  
  neuron->Imahp = gbar*neuron->imahp_m*(neuron->Vm - neuron->erev_k);
}

void update_ca_conc(Compartment * neuron, double t, double dt){

	double icalcium = neuron->Icat + neuron->Icar + neuron->Ical + neuron->icapump + neuron->Ica_nacax;
  
  /* poirazi calcium pumping method: */ 
  double f = 9.648e4; // Coulombs
  double dca;

  /* calculate how much is removed via the calcium pumping mechanism from poirazi et al */
  double drive;
  //double fe = (2.e4)*neuron->radius/((neuron->radius*neuron->radius)-((neuron->radius-.1)*(neuron->radius-.1)));  //10000./18.;
  double surface = 2.0*neuron->radius*PI*neuron->dx*1.e-8; // surface area in cm2
  double vol = PI*neuron->radius*neuron->radius*neuron->dx*1.e-15; // volume in liters
  


	dca = -dt*( (1.e-3)*icalcium*surface/(2.*f*vol) ); 
	double depth = 0.1;
	double taur = 200.; //200.; // 200.; //ms
	double ca_inf = 6.e-5; // mM
	if(icalcium < 0. ) drive = -(10000.*icalcium)/(2.*f*depth);
	else drive = 0.;
	
	neuron->ca_i = neuron->ca_i + dca/40.; //dca/40.; //dt*( (drive/1800.) );
	neuron->ca_i = neuron->ca_i + dt*(ca_inf-neuron->ca_i)/taur;
	//neuron->ca_i = neuron->ca_i + dca +dt*((6.e-5)-neuron->ca_i)/(5.);
 // neuron->ca_o = neuron->ca_o - dca/0.15;  // buffering doesn't affect external calcium directly
	if(t<dt) neuron->cai_Ltype = 0.;
	else{
		neuron->cai_Ltype = neuron->cai_Ltype +  -dt*( (1.e-3)*icalcium*surface/(2.*f*vol) ); 
		neuron->cai_Ltype = neuron->cai_Ltype + dt*(0.-neuron->cai_Ltype)/taur;
	}
}

void update_ca_kager(Compartment * neuron, double dt){

//  /* calculate current from ca-atpase pump */
//  double imax_pump = .5e-2; //estimate must be modified to determine time constant of calcium removal in soma
//  double kpump = 6.9e-3; // mM, from Kager 2007
//  neuron->icapump = imax_pump/(1.+(kpump/neuron->ca_i));
  
	double icalcium = neuron->Icat + neuron->Icar + neuron->Ical + neuron->Ica_nacax + neuron->icapump;
  double f = 9.648e4; // Coulomb/mol
  double dca;
  double surface = 2.0*neuron->radius*PI*neuron->dx*1.e-8; // surface area in cm2
  double vol = PI*neuron->radius*neuron->radius*neuron->dx*1.e-15; // volume in liters


  /* dca from ion movement across membrane: */
  dca = -dt*(1.e-3)*icalcium*surface/(2.*f*vol);  // this considers the entire volume of the segment!


	double ca_i_ss = 6.e-5; // mM
	double taur = 200.; // ms, rate of removal
  /* update calcium */
	neuron->ca_i = neuron->ca_i + dca; 
  neuron->ca_o = neuron->ca_o - dca/0.15;  // interstitial volume = 0.15*intracellular vol
 /* buffer calcium, from code by Kager et al 2008 in NEURON */
	double k1buff = 20.; // mM-1ms-1
	double k2buff = .5; // ms-1
	double kd = .008; // mM
	double totalbuffer = 1.562; // mM
	double BO = totalbuffer/(1.+kd*neuron->ca_i);
	double buffer = BO;
	double cabuffer = totalbuffer-BO;
	double catot = neuron->ca_i*(1.+(totalbuffer/(neuron->ca_i + kd)));
	double b = totalbuffer - catot + kd;
	double c = -kd*catot;
	double d = b*b - 4.*c;
//	neuron->ca_i = (-b+sqrt(d))/2.;
}

void car(Compartment * neuron, double t, double dt, int i){

  // HVAm (Rtype) calcium channel
  
  double gbar;
	
  // icar = gbar* m^3 * h * (v-eca);
  double alpha;
  double beta;
  double tau_m;
  double tau_h;
  if( neuron->distance < 0.1) {
	  
	  /* vhalf and k based on Metz et al 2005 */
    gbar = 2.*4.*(1.e-4); // S/cm2  40.*
	  alpha = 1./(1.+exp(-(neuron->Vm +22.)/7.4));
	  beta = 1./(1.+exp((neuron->Vm+59.)/13.));
	  tau_m = .7; //100.; // ms
	  tau_h = 700.; //7.; //5.; //200.; // ms
  }
  else{
	  gbar = 0.; //40.e-4;  // S/cm2
	  alpha = 1./(1.+exp(-(neuron->Vm +22.)/7.4));
	  beta = 1./(1.+exp((neuron->Vm + 59.)/13.));
    tau_m = .7; 
	  tau_h = 700.; //200.;
  }

  if(t<dt){
    neuron->icar_m = alpha;
    neuron->icar_h = beta;
  }
  else {
    neuron->icar_m = neuron->icar_m + (1.-exp(-dt/tau_m))*(alpha-neuron->icar_m);
    neuron->icar_h = neuron->icar_h + (1.-exp(-dt/tau_h))*(beta-neuron->icar_h);
  }

			     if( (i <= 39 && i >= 35 ) || ( i<=27 && i>=23)) gbar = 0.;
  double eca = (1.e3)*(R*T/(2.*F))*log(neuron->ca_o/neuron->ca_i);
  neuron->Icar = (gbar)*neuron->icar_m*neuron->icar_m*neuron->icar_m*neuron->icar_h*(neuron->Vm - eca);

}



void cal(Compartment * neuron, double t, double dt, int i){

	
	/* large conductance L-type dendritic current as described in Magee 1995 */
	/* no time-dep inactivation observed */
	double minf = 1./(1.+exp(-(neuron->Vm-9.)/6.));
	double tau_m = 1.5; // ms; slow time constant from Magee
	if(t<dt) neuron->ical_m = minf;
	else neuron->ical_m = neuron->ical_m + (1.-exp(-dt/tau_m))*(minf - neuron->ical_m);
	
	double gbar = .5*.001;
	 double eca = (1.e3)*(R*T/(2.*F))*log(neuron->ca_o/neuron->ca_i);
	neuron->Ical = gbar*neuron->ical_m*(neuron->Vm -eca);
	
 // // somatic HVA: similar to T-type
//  if(neuron->distance < 0.1){
//    double alpha_m = -0.055*(neuron->Vm + 27.01)/(exp(-(neuron->Vm + 27.01)/3.8)-1.);
//    double beta_m = 0.94*exp(-(neuron->Vm+63.01)/17.);
//    double tau_m = 1./(5.*(alpha_m+beta_m));
//    double minf = alpha_m/(alpha_m+beta_m);
//    
//    if(t<dt) {
//      neuron->ical_m = minf;
//    }
//    else {
//      neuron->ical_m = neuron->ical_m + (1.-exp(-dt/tau_m))*(minf-neuron->ical_m);
//    }
//    
//    double x = 0.0853*T/2.;
//    double z = neuron->Vm/x;
//    double f;
//    if(fabs(z)<0.0001) f = 1.-(z/2.);
//    else f = z/(exp(z)-1.);
//
//    double ghk = -x*(1.-(neuron->ca_i/neuron->ca_o)*exp(neuron->Vm/x))*f;
//	  double gbar = 1.e-3; //1.*1.e-4; // S/cm2, somatic conductance
//    neuron->Ical = 1.*(gbar)*neuron->ical_m*(0.001/(0.001+neuron->ca_i))*ghk;
//  }
//
//  // dendritic HVA similar to R-type
//  else {
//    double tau_m = 3.6; // ms
//    double tau_h = 29.; // ms
//    double alpha = 1./(1.+exp(-(neuron->Vm+37.)));
//    double beta = 1./(1.+exp((neuron->Vm+41.)/0.5));
//    
//    if(t<dt){
//      neuron->ical_m = alpha;
//      neuron->ical_h = beta;
//    }
//    else {
//      neuron->ical_m = neuron->ical_m + (1.-exp(-dt/tau_m))*(alpha - neuron->ical_m);
//      neuron->ical_h = neuron->ical_h + (1.-exp(-dt/tau_h))*(beta - neuron->ical_h);
//    }
//
//    double gbar;
//	  double ginit = 0.0116*1.e-2; //1.e-3;  // S/cm2
//    if(neuron->distance < 50.) gbar = 0.1*ginit;
//	  if(neuron->distance >= 50.) gbar = 0.; //4.6*ginit;
//	if( (i <= 39 && i >= 35 ) || ( i<=27 && i>=23) ) gbar = 0.;
//    double eca = (1.e3)*(R*T/(2.*F))*log(neuron->ca_o/neuron->ca_i);
//    neuron->Ical = (gbar)*neuron->ical_m*neuron->ical_m*neuron->ical_m*neuron->ical_h*(neuron->Vm - eca);
//  }

}

void cal_markaki(Compartment * neuron, double t, double dt, int i){


	double ki = 0.025; 
	double taumin = 100.;
	double vhalf = 9.; //-1.;
	double z = -4.6;
	double t0 = 1.5;
	double b = 0.01;
	
	double h2 = ki/(ki+ neuron->ca_i);
	double alpha = exp((1.e-3)*z*(neuron->Vm-vhalf)*F/(R*T));
	double minf = 1./(1.+alpha);
	double alpha2 = (neuron->ca_i/b)*(neuron->ca_i/b);
	double sinf = alpha2/(1.+alpha2);
	double tau_m = taumin + 1./(neuron->ca_i+b);
	
	if(t<dt){
		neuron->ical_m = minf;
		neuron->ical_h = sinf;
	}
	else {
		neuron->ical_m = neuron->ical_m + dt*(minf-neuron->ical_m)/t0;
		neuron->ical_h = neuron->ical_h + dt*(sinf-neuron->ical_h)/tau_m;
	}
	
	double gbar = .003;
	double eca = (1.e3)*(R*T/(2.*F))*log(neuron->ca_o/neuron->ca_i);
	neuron->Ical = gbar*(h2*neuron->ical_m*neuron->ical_m + 8.*neuron->ical_h*neuron->ical_h)*(neuron->Vm-eca);
}

void cat(Compartment * neuron, double t, double dt, int i){

  // LVA (T-type) calcium channel from Poirazi et al, adapted from Magee and Johnston 1995

  double x;
  double z;
  double gbar;
  
//  double alpha_m = -0.196*(neuron->Vm - 19.88)/(exp(-(neuron->Vm - 19.88)/10.) - 1.);
////	double alpha_m = -0.196*(neuron->Vm - 36)/(exp(-(neuron->Vm - 36)/10.) - 1.);
//  double beta_m = 0.046*exp(-(neuron->Vm/22.73));
////  double alpha_h = 0.00016*exp(-(neuron->Vm+57.)/19.);
//	double alpha_h = 0.00016*exp(-(neuron->Vm+57.)/19.);
//  double beta_h = 1./(1.+exp(-(neuron->Vm-15.)/10.));
//  double tau_m = 1./(alpha_m+beta_m);
//  double tau_h = 1./(0.68*(alpha_h+beta_h));
//
//  double minf = alpha_m/(alpha_m+beta_m);
//  double hinf = alpha_h/(alpha_h+beta_h);
	
	/* from seigelbaum -- tsay et al 2007  and magee and johnston 1995 */
	double alpha_m = exp(-(neuron->Vm + 30.)/5.);
	double minf = 1./(1.+alpha_m);
	double alpha_h = exp(-(neuron->Vm+67.)/-6.);
	double hinf = 1./(1+alpha_h);
	double tau_m = 1.1; //(2. + 32./(1.+exp(-(neuron->Vm+15.)/-9.)));   // time constants from text of magee and johnston 1995
	double tau_h = 50.; //60.+120./(1.+exp(-(neuron->Vm+40.)/-5.));

  if(t<dt) {
    neuron->icat_m = minf;
    neuron->icat_h = hinf;
  }
  else {
    neuron->icat_m = neuron->icat_m + (1.-exp(-dt/tau_m))*(minf-neuron->icat_m);
    neuron->icat_h = neuron->icat_h + (1.-exp(-dt/tau_h))*(hinf-neuron->icat_h);
  }

  /* conductance densities */
	double gfactor = 6.;
	double ginit = .5e-3; //.5E-3; // S/cm2  .01?
	double gsoma = 0.; //.0001e-3; //0.0001E-3; // S/cm2 (could be higher??)
  if( neuron->distance >= 150.0) gbar = ginit*gfactor*neuron->distance/350.;
  else if( neuron-> distance < 50. ) gbar = gsoma;
  else gbar = 0.0;

  if( (i <= 39 && i >= 35 ) || ( i<=27 && i>=23)) gbar = 0.;
  x = 0.0853*T/2.;
  z = neuron->Vm/x;
  double f;
  if(fabs(z)<0.0001) f = 1.-(z/2.);
  else f = z/(exp(z)-1.);

  double ghk = -x*(1.-(neuron->ca_i/neuron->ca_o)*exp(neuron->Vm/x))*f;
  neuron->Icat = 1.*(gbar)*pow(neuron->icat_m, 2)*neuron->icat_h*(0.001/(0.001+neuron->ca_i))*ghk;
	//neuron->Icat = gbar*neuron->icat_m*neuron->icat_h*(0.001/(0.001+neuron->ca_i))*ghk;
}


void ih(Compartment * neuron, double t, double dt, int i){

	double gbar_soma = 1.*.05*(1.e-3); //.05*(1.e-3); //.05*(1.e-3); //35.; // S/cm2; upper limit in Golding et al. is 1.7e-3 S/cm2
  double gend = 8.*gbar_soma;
  double dhalf = 250.; // um
  double steep = 50.; // um
  
  double gbar = gbar_soma + (gend-gbar_soma)/(1.+exp((dhalf-neuron->distance)/steep));

  /* Ih = gbar * m * (V-Eh) 
     in poirazi, assumed Eh = -10 mV */

  double tau;
  if (neuron->Vm > -30.) tau = 1.; // ms
  else tau = (2./(exp(-(neuron->Vm+145.)/17.5) + exp((neuron->Vm+16.8)/16.5)))+10.;

	double minf = 1.0-(1./(1.+exp(-(neuron->Vm+90.)/8.5)));
  if (t < dt ) neuron->ih_m = minf;
  else neuron->ih_m = neuron->ih_m + dt*(minf-neuron->ih_m)/tau;
  
  /* to test influence on input resistance: */
  //if(t>10.) gbar = gbar*.2;  // see figure s1 in poirazi et al
	 if( (i <= 39 && i >= 35 ) || ( i<=27 && i>=23) ) gbar = 0.;
  double erev_k = (1.e3)*(R*T/F)*log(neuron->k_o/neuron->k_i);
  double erev_na = (1.e3)*(R*T/F)*log(neuron->na_o/neuron->na_i);
  neuron->Ik_h = 0.5*(gbar)*neuron->ih_m*(neuron->Vm - neuron->erev_k);
  neuron->Ih_na = 0.5*(gbar)*neuron->ih_m*(neuron->Vm - neuron->erev_na);


}

void ih2(Compartment * neuron, double t, double dt){

  double gbar_soma = 18.72*(1.e-3); // S/cm2
  double gend = 9.*gbar_soma;
  double dhalf = 280.; // um
  double steep = 50.; // um
  double gbar = gbar_soma + (gend-gbar_soma)/(1.+exp((dhalf-neuron->distance)/steep));

  double Vhalfs = -82.; // from narayan and johnston, 2008 based on magee 1998
  double k = 8.0;
  double Vhalft = -75.0; 
  double sinf, taus;

  sinf = 1.0/(1.0+exp((neuron->Vm-Vhalfs)/k));
  taus = exp(0.033*(neuron->Vm-Vhalft))/(0.011*(1.0+exp(0.083*(neuron->Vm-Vhalft))));

  if(t < dt) neuron->ih_m = sinf;
  else neuron->ih_m = neuron->ih_m + dt*(sinf-neuron->ih_m)/taus;

  neuron->Ik_h = 0.6*gbar*neuron->ih_m*(neuron->Vm - neuron->erev_k);
  neuron->Ih_na = 0.4*gbar*neuron->ih_m*(neuron->Vm - neuron->erev_na);
  
//  if(t<dt) cout << "gh: " << gbar*neuron->ih_m << endl;

}

void na6state(Compartment * neuron, double t, double dt, int j){
	Vec_DP y(6);
	Vec_DP dydt(6);
	Vec_DP yout(6);
	Vec_DP dfdt(6);
	Mat_DP dfdy(6,6);
	
	Vec_DP y5(5);
	Vec_DP dydt5(5);
	Vec_DP yout5(5);
	Vec_DP dfdt5(5);
	Mat_DP dfdy5(5,5);
	
	for(int i=0; i<6; i++) yout[i] = neuron->na6_y[i];
	for(int i=0; i<5; i++) yout5[i] = neuron->na5_y[i];
	
	/* set initial conditions */
	if (t<dt) {
		
		if( neuron->Vm < -60. && neuron->Vm > -64.) {
		/* steady state values for vrest = -63 mV */
		yout[0] =  0.565492;
		yout[1] = 0.369513;
		yout[4] = 0.000534;
		yout[3] = 0.050157;
		yout[5] = 0.014035;
		yout[2] = 1.-yout[0]-yout[3]-yout[1]-yout[4]-yout[5];
		
		yout5[0] = 0.583754775036;
		yout5[1] = 0.415411661629;
		yout5[2] = 0.000279733193;
		yout5[4] = 0.000553830117;
		yout5[3] = 1.-yout5[0] - yout5[1] - yout5[2] - yout5[4];
		}
		
		else if( neuron->Vm <= -64. && neuron->Vm > -68.){
		/* steady state values for vrest = -66 mV */
		yout[0] =  .641218282554;
		yout[1] = 0.320073266685;
		yout[4] = 0.000201276546;
		yout[3] = 0.022151922814;
		yout[5] = 0.016224966276;
		yout[2] = 1.-yout[0]-yout[3]-yout[1]-yout[4]-yout[5];
		
		yout5[0] = .355386676164; //0.647544839544;
		yout5[1] = .644166906841; //0.352118452556;
		yout5[2] = .000072686938; //0.000132441802;
		yout5[4] = .000373730553; //0.000204290880;
		yout5[3] = 1.-yout5[0] - yout5[1] - yout5[2] - yout5[4];
		}
		
		else if( neuron->Vm <= -68.) {
		/* steady state values for vrest = -70 mV */
		yout[0] =  0.7222082720215;
		yout[1] = 0.251744902763;
		yout[4] = 0.000052288792;
		yout[3] = 0.007197388401;
		yout[2] = 0.000044918797;
		yout[5] = 1.-yout[0]-yout[3]-yout[1]-yout[4]-yout[2];
		
		yout5[0] = 0.724628078897;
		yout5[1] = 0.275271643557;
		yout5[2] = 0.000047627666;
		yout5[4] = 0.000052649852;
		yout5[3] = 1.-yout5[0] - yout5[1] - yout5[2] - yout5[4];
		}
		
		else {
		/* ss for -58; adjust 5-state */
		yout[0] =  .403367100012;
		yout[1] = 0.412886716500;
		yout[4] = 0.002386114588;
		yout[3] = 0.170871318832;
		yout[2] = 0.000794246208;
		yout[5] = 1.-yout[0]-yout[3]-yout[1]-yout[4]-yout[2];
		
		yout5[0] = .081947649813;
		yout5[1] = .91302286898; //0.352118452556;
		yout5[2] = .000162298285; //0.000132441802;
		yout5[4] = .00486718316; //0.000204290880;
		yout5[3] = 1.-yout5[0] - yout5[1] - yout5[2] - yout5[4];
		}
		
	}
	else {
		NR::jacobn_na6(t, yout, dfdt, dfdy, neuron->Vm);
		imp_euler(dfdy, yout, dt);
		NR::jacobn_na5(t, yout5, dfdt5, dfdy5, neuron->Vm);
		imp_euler(dfdy5, yout5, dt);
	}
	
	for(int i=0; i<6; i++) neuron->na6_y[i] = yout[i];
	for(int i=0; i<5; i++) neuron->na5_y[i]=yout5[i];
	
	double fraction = neuron->distance/250.; //1140.;
	if(fraction > 1.) fraction = 1.;
	
	double gbar = .007; // .08; //.012;
	if( ( j<27 && j>=23) ) {
		gbar = 0.1;
		fraction = 0.;
	}
	if( j <5 ) {
		gbar = 0.03; 
		fraction = 0.;
	}
	
	if (neuron->radius < .6) {
		gbar = 0.1*gbar;
		//if(t<dt) cout << "gbar = 0 for compartment " << j << endl;
	}
	
	gbar = gbar;
	
	neuron->Ina = (fraction*gbar*neuron->na6_y[2] + (1.-fraction)*gbar*neuron->na5_y[2])*(neuron->Vm - neuron->erev_na) ;
	
	/* persistent sodium current */
	double nap_m = 1./(1.+exp(-(neuron->Vm+50.4/4.5)));
	double gnap = 5.6*(1.e-7);
	double inap = gnap*nap_m*nap_m*nap_m*(neuron->Vm-neuron->erev_na);
	
	//if(neuron->distance > 350.) neuron->Ina = neuron->Ina + inap;
}

void NR::derivs_na6(const DP t,  Vec_I_DP &y, Vec_O_DP &dydt, double vm){

	double r02 = exp(5.218 + 0.1066*vm);
	double r12 = exp(2.187+0.04433*vm);
	double r14 = exp(6.863+0.22*vm);
	double r23 = exp(-11.53+0.03047*vm);
	double r25 = exp(0.5124+0.005264*vm);
	double r34 = exp(-2.802+0.053*vm);
	double r45 = exp(-3.671+0.04366*vm);
	double r20 = exp(-5.018+-0.1773*vm);
	double r21 = exp(-2.819+-0.1498*vm);
	double r41 = exp(-4.085+-0.05757*vm);
	double r32 = exp(-18.68+-0.0000025*vm);
	double r52 = exp(14.85+0.2956*vm);
	double r43 = exp(-1.599+0.*vm);
	double r54 = exp(16.61+0.4175*vm);
	
	dydt[0] = y[2]*r20 - y[0]*r02;
	dydt[1] = y[2]*r21 + y[4]*r41 - y[1]*(r12+r14);
	dydt[2] = y[3]*r32 + y[5]*r52 + y[0]*r02 + y[1]*r12 -y[2]*(r20+r21+r23+r25);
	dydt[3] = y[4]*r43 + y[2]*r23 - y[3]*(r32+r34);
	dydt[4] = y[5]*r54 + y[1]*r14 + y[3]*r34 - y[4]*(r41+r43+r45);
	dydt[5] = y[2]*r25 + y[4]*r45 - y[5]*(r52+r54);

}

void NR::jacobn_na6(const DP t, Vec_I_DP &y, Vec_O_DP &dfdt, Mat_O_DP &dfdy, double vm){

	double r02 = exp(5.218 + 0.1066*vm);
	double r12 = exp(2.187+0.04433*vm);
	double r14 = exp(6.863+0.22*vm);
	double r23 = exp(-11.53+0.03047*vm);
	double r25 = exp(0.5124+0.005264*vm);
	double r34 = exp(-2.802+0.053*vm);
	double r45 = exp(-3.671+0.04366*vm);
	double r20 = exp(-5.018+-0.1773*vm);
	double r21 = exp(-2.819+-0.1498*vm);
	double r41 = exp(-4.085+-0.05757*vm);
	double r32 = exp(-18.68+-0.0000025*vm);
	double r52 = exp(14.85+0.2956*vm);
	double r43 = exp(-1.599+0.*vm);
	double r54 = exp(16.61+0.4175*vm);
	
	int i; 
	int n=y.size();
	for (i=0; i<n; i++) dfdt[i]=0.0;
	dfdy[0][0] = -r02;
	dfdy[0][1] = 0.0;
	dfdy[0][2] = r20;
	dfdy[0][3] = 0.0;
	dfdy[0][4] = 0.0;
	dfdy[0][5] = 0.0;
	dfdy[1][0] = 0.0;
	dfdy[1][1] = -(r12+r14);
	dfdy[1][2] = r21;
	dfdy[1][3] = 0.0;
	dfdy[1][4] = r41;
	dfdy[1][5] = 0.0;
	dfdy[2][0] = r02;
	dfdy[2][1] = r12;
	dfdy[2][2] = -(r20+r21+r23+r25);
	dfdy[2][3] = r32;
	dfdy[2][4] = 0.0;
	dfdy[2][5] = r52;
	dfdy[3][0] = 0.0;
	dfdy[3][1] = 0.0;
	dfdy[3][2] = r23;
	dfdy[3][3] = -(r32+r34);
	dfdy[3][4] = r43;
	dfdy[3][5] = 0.0;
	dfdy[4][0] = 0.0;
	dfdy[4][1] = r14;
	dfdy[4][2] = 0.0;
	dfdy[4][3] = r34;
	dfdy[4][4] = -(r41+r43+r45);
	dfdy[4][5] = r54;
	dfdy[5][0] = 0.0;
	dfdy[5][1] = 0.0;
	dfdy[5][2] = r25;
	dfdy[5][3] = 0.0;
	dfdy[5][4] = r45;
	dfdy[5][5] = -(r52+r54);
	
}

void NR::jacobn_na5(const DP t, Vec_I_DP &y, Vec_O_DP &dfdt, Mat_O_DP &dfdy, double vm){
	
	double r02 = exp(5.218 + 0.1066*vm);
	double r12 = exp(2.147+0.04433*vm);
	double r14 = exp(5.498+0.22*vm);
	double r23 = exp(0.5124+0.004891*vm);
	double r34 = exp(16.61+0.4407*vm);
	double r20 = exp(-5.018+-0.1772*vm);
	double r21 = exp(-2.780+-0.1498*vm);
	double r41 = exp(-5.371+-0.05757*vm);
	double r32 = exp(14.85+-0.2960*vm);
	double r43 = exp(-3.671+0.06612*vm);
	
	int i; 
	int n=y.size();
	for (i=0; i<n; i++) dfdt[i]=0.0;
	dfdy[0][0] = -r02;
	dfdy[0][1] = 0.0;
	dfdy[0][2] = r20;
	dfdy[0][3] = 0.0;
	dfdy[0][4] = 0.0;
	dfdy[1][0] = 0.0;
	dfdy[1][1] = -(r12+r14);
	dfdy[1][2] = r21;
	dfdy[1][3] = 0.0;
	dfdy[1][4] = r41;
	dfdy[2][0] = r02;
	dfdy[2][1] = r12;
	dfdy[2][2] = -(r20+r21+r23);
	dfdy[2][3] = r32;
	dfdy[2][4] = 0.0;
	dfdy[3][0] = 0.0;
	dfdy[3][1] = 0.0;
	dfdy[3][2] = r23;
	dfdy[3][3] = -(r32+r34);
	dfdy[3][4] = r43;
	dfdy[4][0] = 0.0;
	dfdy[4][1] = r14;
	dfdy[4][2] = 0.0;
	dfdy[4][3] = r34;
	dfdy[4][4] = -(r41+r43);

	
}
