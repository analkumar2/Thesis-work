#include <iostream>
#include <math.h>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include "neuron_structures.h"
#include "ion_transport.h"

void nacax_old(Compartment * neuron, double t, double dt){
 
  double surface = 2.0*neuron->radius*PI*neuron->dx*1.e-8; // surface area in cm2
  double vol = PI*neuron->radius*neuron->radius*neuron->dx*1.e-15; // volume in liters 
  double ca_i_ss = 6.e-5; //mM
	/* calculate current from ca-atpase pump */
//  double vmax_pump = 2.*352.*5.e-7; //352.; // uM/s; from Somjen 2008
//  double kpump = 6.9e-3; // mM, from Kager 2007
//  if( (neuron->ca_i-ca_i_ss) > 0.) {
//	neuron->icapump =  vmax_pump*2.*F/(1.+(kpump/(neuron->ca_i-ca_i_ss)));
//  }
////  else if( (neuron->ca_i - ca_i_ss) < 0.) {
////	neuron->icapump = -vmax_pump*2.*F/(1.+(kpump/-(neuron->ca_i-ca_i_ss)));
////  }
//  else {
////	neuron->icapump = vmax_pump*2.*F*(neuron->ca_i-ca_i_ss)/(kpump);
//	neuron->icapump = 0.0;
//   }
//	neuron->icapump = 0.0;
 // pump equation from De Schutter and Smolen, in "Methods in Neuronal Modeling", p 219-220
 double vmax_pump = 2.88*9.e-11; // mol cm-2 ms-1                                      
 double kpump = 1.e-3; // mM                                                      
 neuron->icapump = (1.e3)*2.*F*vmax_pump*neuron->ca_i/(kpump + neuron->ca_i);
 
 /* calculating na/ca exchange */
  double kna = 87.5; // mM
  double kca = 1.38; // mM
  double gamma = 0.35; // voltage dependence parameter

  if(t < dt) { // need to set max current for nacax at steady state	
	double inacax_ss = -1.*(neuron->Icat + neuron->Icar + neuron->Ical) - neuron->icapump;  // mA/cm2
	neuron->imax_nacax =  inacax_ss*( (pow(kna, 3) + pow(neuron->na_o, 3))*(kca+neuron->ca_o)*(1.+0.1*exp((gamma-1.)*(neuron->Vm*1.e-3)*F/(R*T))) );
	neuron->imax_nacax = -.5*neuron->imax_nacax/( pow(neuron->na_i, 3)*neuron->ca_o*exp(gamma*(neuron->Vm*1.e-3)*F/(R*T)) - pow(neuron->na_o, 3)*neuron->ca_i*exp((gamma-1.)*(neuron->Vm*1.e-3)*F/(R*T)) );
  }

  double inacax = neuron->imax_nacax*( pow(neuron->na_i, 3)*neuron->ca_o*exp(gamma*(neuron->Vm*1.e-3)*F/(R*T)) - pow(neuron->na_o, 3)*neuron->ca_i*exp((gamma-1.)*(neuron->Vm*1.e-3)*F/(R*T)) );
  inacax = inacax/( (pow(kna, 3) + pow(neuron->na_o, 3))*(kca+neuron->ca_o)*(1.+0.1*exp((gamma-1.)*(neuron->Vm*1.e-3)*F/(R*T))) );
  
  neuron->Ina_nacax = 3.*inacax;
  neuron->Ica_nacax = -2.*inacax;
  
}

void ca_pump(Compartment * neuron, double t, double dt){
	
	// NOT USED ///////

 double kpump = 1.e-3; // mM                                                      
 double vpump = (1.e3)*2.*F*neuron->ca_i/(kpump + neuron->ca_i);

neuron->Ica_total = neuron->Icat + neuron->Icar + neuron->Ical + neuron->Ica_nacax;

	if(t<dt){
		neuron->imax_capump = -(neuron->Ica_total)/vpump;
		if(neuron->imax_capump < 0.0 ) cout << "capump wrong direction at rest!!" << endl;
	}
	
	neuron->icapump = neuron->imax_capump*vpump;
//	 cout << "neuron->Ica_total + neuron->icapump: " << neuron->Ica_total + neuron->icapump << endl;
}

void ca_pump_new(Compartment * neuron, double t, double dt){
	/* ca pump equation from Somjen et al, in Neuron modelDB */
	double surface = 2.0*neuron->radius*PI*neuron->dx*1.e-8;  // surface area in um2
	double vol = PI*pow(neuron->radius, 2.)*neuron->dx*1.e-15; // volume in um3
	
	double ratio = vol/surface;
	double km = 0.0069; // mM
	double hill = 1.;
	double scale = 1.e-4;

	double ca_i_ss = 6.e-5;
//	if(neuron->distance < 10.) ca_i_ss = 5.999e-5; // mM
//	else ca_i_ss = 6.e-5;
	double rate;
	
	double max_rate = .5*352.*1.e-2; // uM/s
	double imax = max_rate*2.*F*ratio*(1e-3); // mA/cm2

//	if( fabs(neuron->ca_i - ca_i_ss) < 1.e-7) neuron->icapump = imax*(neuron->ca_i-ca_i_ss)/km;
//	else 
	if (neuron->ca_i > ca_i_ss) {
		neuron->icapump = imax/(1.+km/(neuron->ca_i - ca_i_ss));
	}
	else neuron->icapump = 0.;
}

void nacax(Compartment * neuron, double t, double dt){
	
	double surface = 2.0*neuron->radius*PI*neuron->dx*1.e-8; // surface area in cm2
	double vol = PI*neuron->radius*neuron->radius*neuron->dx*1.e-15; // volume in liters 
	double ca_i_ss = 6.e-5; //mM
	
	 /* calculating na/ca exchange */
  double kna = 87.5; // mM
	double kca = 1.38; // mM
  double gamma = 0.35; // voltage dependence parameter

	/* from Zhu and Clancy 2007 */
	double vnacax = ( pow(neuron->na_i, 3)*neuron->ca_o*exp(gamma*(neuron->Vm*1.e-3)*F/(R*T)) - pow(neuron->na_o, 3)*neuron->ca_i*exp((gamma-1.)*(neuron->Vm*1.e-3)*F/(R*T)) );
	vnacax = vnacax/( (pow(kna, 3) + pow(neuron->na_o, 3))*(kca+neuron->ca_o)*(1.+0.1*exp((gamma-1.)*(neuron->Vm*1.e-3)*F/(R*T))) );
	neuron->Ina_total = neuron->Inaleak + neuron->na_atpase + neuron->na_nkcc1 + neuron->Ina + neuron->Ih_na + neuron->Ina_nmda + neuron->Ina_ampa;
	neuron->Ica_total = neuron->Icat + neuron->Icar + neuron->Ical + neuron->icapump;

	// nacax max transport rate
	if(t<dt){
//		neuron->imax_nacax = ( -neuron->Ina_total/(3.*vnacax) );
//		if(neuron->imax_nacax > 0.) cout << "nacaex wrong direction!" << endl;
		neuron->imax_nacax = (-neuron->Ica_total/(-2.*vnacax));
		if(neuron->imax_nacax < 0.) cout << "nacax transporting calcium in!!!" << endl;
	}

	neuron->Ina_nacax = neuron->imax_nacax*vnacax*3.;
	neuron->Ica_nacax = neuron->imax_nacax*vnacax*-2.;
	
}

void nacax_kager(Compartment * neuron, double t, double dt){
	
	double k = R*T/(F*1.e-3);
	double q10 = pow(3., (celsius-37.)/10.);
	double gamma = .35;  // voltage factor
	double kqa = exp(gamma*(neuron->Vm/k));
	double kb = exp((gamma-1.)*(neuron->Vm/k));
	
	double kna = 87.5; // mM
	double kca = 100.; //1.38;
	
	double rate = q10*(kqa*neuron->na_i*neuron->na_i*neuron->na_i*neuron->ca_o - kb*neuron->na_o*neuron->na_o*neuron->na_o*neuron->ca_i);
	rate = rate/((kna*kna*kna+neuron->na_o*neuron->na_o*neuron->na_o)*(kca+neuron->ca_o)*(1.+0.1*kb));
	
	neuron->Ica_total = neuron->Icat + neuron->Icar + neuron->Ical + neuron->icapump;
	if(t<dt){
		neuron->imax_nacax = (-neuron->Ica_total/(-2.*rate));
		if(neuron->imax_nacax < 0.) cout << "nacax transporting calcium in!!!" << endl;
	
	}
	
	neuron->Ina_nacax = neuron->imax_nacax*rate*3.;
	neuron->Ica_nacax = neuron->imax_nacax*rate*-2.;

}

void kcc2_chang(Compartment * neuron, double t, double dt, int i){

  double surface = 1.0*neuron->radius*PI*neuron->dx*1.e-8; // surface area in cm2
  double vol = PI*neuron->radius*neuron->radius*neuron->dx*1.e-15; // volume in liters
	double kk = 5.0; // mM   from staley/proctor
	double kcl = 9.0; //6.0;
	double vmax =1.*300.e-7; //300.e-7; // modified to match staley/proctor data
	if( i<5) vmax = 1.25*vmax;
  if( (i <= 39 && i >= 35 ) || ( i<=27 && i>=23) ) vmax = .06*1.25*vmax;
  	
  double numerator = (neuron->k_o*neuron->cl_o - neuron->k_i*neuron->cl_i)/(kk*kcl);
  double denom = 1.0+neuron->k_o*neuron->cl_o/(kk*kcl);
  denom = denom*(1.0+neuron->k_i/kk)*(1.+neuron->cl_i/kcl);
  denom = denom + (1.+neuron->k_i*neuron->cl_i/(kk*kcl))*(1.+neuron->k_o/kk)*(1.+neuron->cl_o/kcl);
  
//	fixed kcc2 k+:
/*
	double numerator = (3.*neuron->cl_o - 110.*neuron->cl_i)/(kk*kcl);
	double denom = 1.0+3.*neuron->cl_o/(kk*kcl);
	denom = denom*(1.0+110./kk)*(1.+neuron->cl_i/kcl);
	denom = denom + (1.+110.*neuron->cl_i/(kk*kcl))*(1.+3./kk)*(1.+neuron->cl_o/kcl);
*/
       double vkcc2 = vmax*numerator/denom;
//  if(fabs(vkcc2) < 1.e-15) vkcc2 = 0.0;
  neuron->k_kcc2 = -vkcc2*F;
  neuron->cl_kcc2 = vkcc2*F;
	if(i==0 && t<dt) cout << "Ikcc2,max: " << vmax*F << endl;

}

void atpase(Compartment * neuron, double t, double dt, int i){

	double km_k =  1.; //3; // mM, pump affinity for extracellular potassium
	double km_na = 10.; //8.; //10.0; // mM, pump affinity for intracellular sodium
  double Apump = pow((1. + (km_k/neuron->k_o)), -2.0)*pow((1.+(km_na/neuron->na_i)), -3.0);
//  double numerator = (pow(neuron->k_o,2.)*pow(neuron->na_i,3.) - pow(neuron->k_i,2.)*pow(neuron->na_o,3))/(pow(km_k,2)*pow(km_na,3));
//  double denom = 1.0+neuron->k_o*neuron->na_i/(km_k*km_na);
//  denom = denom*(1.0+neuron->k_i/km_k)*(1.+neuron->na_o/km_na);
//  denom = denom + (1.+neuron->k_i*neuron->na_o/(km_k*km_na))*(1.+neuron->k_o/km_k)*(1.+neuron->na_i/km_na);
//  double Apump = numerator/denom; 
  
	neuron->Ik_total = neuron->k_kcc2 + neuron->k_nkcc1 + neuron->Ikdr + neuron->Ika + neuron->Ik_h + neuron->Ikm + neuron->Isahp + neuron->Imahp + neuron->Ik_nmda + neuron->Ik_ampa;
	neuron->Ina_total =neuron->na_nkcc1 + neuron->Ina + neuron->Ih_na + neuron->Ina_nmda + neuron->Ina_ampa + neuron->Ina_nacax; 

  if(t<dt) {  // set max pump velocity at time 0 to achieve equilibrium
	  double a = -1.5*neuron->Ik_total/(-2.*Apump);  //-1.35*neuron->Ik_total/(-2.*Apump);
	  double b = -1.5*neuron->Ina_total/(3.*Apump); //-1.0001*neuron->Ina_total/(3.*Apump);
	  if ( a > b) neuron->imax_atpase = a;
	  else neuron->imax_atpase = b;
	  if(i==0) cout << "imax_atpase = " << -2.*neuron->imax_atpase*Apump << endl;
	}
  
  neuron->k_atpase = -2.*neuron->imax_atpase*Apump; // K+ pumped inward therefore negative current
  neuron->na_atpase = 3.*neuron->imax_atpase*Apump; // Na+ pumped outward therefore positive current

}

void LRatpase(Compartment * neuron, double t, double dt){

	double km_k = 1.5;
	double km_na = 10.;
	double sigma = (exp(neuron->na_o/67.3) - 1.0)/7.0;
	double fnak = 1./(1.+0.1245*exp(-0.1*(neuron->Vm*1.e-3)*F/(R*T))+0.0365*sigma*exp(-(neuron->Vm*1.e-3)*F/(R*T)));
	double Apump = fnak*(neuron->k_o/(neuron->k_o+km_k))/(1.+(km_na/neuron->na_i)*(km_na/neuron->na_i));
	
	neuron->Ik_total = neuron->Ikleak + neuron->k_kcc2 + neuron->k_nkcc1 + neuron->Ikdr + neuron->Ika + neuron->Ik_h + neuron->Ikm + neuron->Isahp + neuron->Imahp + neuron->Ik_nmda + neuron->Ik_ampa;
	
	if(t<dt) {  // set max pump velocity at time 0 to achieve equilibrium
		neuron->imax_atpase = -neuron->Ik_total/(-2.*Apump);
	}
	
	neuron->k_atpase = -2.*neuron->imax_atpase*Apump; // K+ pumped inward therefore negative current
	neuron->na_atpase = 3.*neuron->imax_atpase*Apump; // Na+ pumped outward therefore positive current
	
}

void nkcc1_2state(Compartment * neuron, double t, double dt, int i){

  
  double surface =2.0*neuron->radius*PI*neuron->dx*1.e-8;
  double vol = PI*neuron->radius*neuron->radius*neuron->dx*1.e-15;
  
  // rate constants, from table 1 of terashima et al 2006
  double kf_full = 3.065e3;
  double kb_full = 1.456e3;
  double kf_empty = 37.767e3; // s-1
  double kb_empty = kf_full*kf_empty/kb_full;
  double kna = 0.08445; // mM-1
  double kcl = 0.05735; // mM-1
  double kk = 1.16e-3; // mM-1
  double pE1, pE1full, pE2, pE2full, alpha_2state, beta_2state, factor;
  
  double na_o = neuron->na_o;
  double cl_o = neuron->cl_o;
  double k_o = neuron->k_o;
  double na_i = neuron->na_i;
  double k_i = neuron->k_i;
  double cl_i = neuron->cl_i;

  double ko1 = kna* na_o;
  double ko2 = kna*na_o*(kcl*cl_o);
  double ko3 = kna*na_o*(kcl*cl_o)*(kk*k_o);
  double ko4 = kna*na_o*(pow(kcl*cl_o,2))*(kk*k_o);


  double ki1 = kcl*cl_i;
  double ki2 = ki1*kk*k_i;
  double ki3 = ki2*kcl*cl_i;
	double ki4 = ki3*kna*na_i;


//   set initial conditions
  if(t == 0.0){
   

    // E1 state:
	  pE1 = 1.0/(1.0 + ko1 + ko2 + ko3 + ko4); 
	  pE1full = ko4*pE1; 
    //E2 state:
	  pE2 = 1.0/(1.0 + ki1 + ki2 +ki3 + ki4); 
	  pE2full = ki4*pE2; 

    alpha_2state = kf_full*pE1full + kb_empty*pE1;
    beta_2state = kb_full*pE2full + kf_empty*pE2;
    neuron->y_nkcc = beta_2state/(alpha_2state+beta_2state);
  }

  //after time = 0, use forward euler method to update y
  if (t>=dt){
    
    // E1 state:  external binding
	  pE1 = 1.0/( 1.0 + ko1 + ko2 + ko3 + ko4); 
	  pE1full = ko4*pE1; //ko10*pE1;
    // E2 state:  internal binding
	  pE2 = 1.0/( 1.0 + ki1 + ki2 +ki3 + ki4); 
	  pE2full = ki4*pE2; 

    alpha_2state = kf_full*pE1full + kb_empty*pE1;
    beta_2state = kb_full*pE2full + kf_empty*pE2;
//    factor = dt*( beta_2state*(1-neuron->y_nkcc) - alpha_2state*neuron->y_nkcc );
//    neuron->y_nkcc = neuron->y_nkcc + factor;
	//using a sort-of implicit euler scheme:
	neuron->y_nkcc = (neuron->y_nkcc + dt*beta_2state)/(1.+dt*beta_2state+dt*alpha_2state);
//  
  }

//	double numerator = (neuron->na_o*neuron->k_o*pow(neuron->cl_o,2) - neuron->na_i*neuron->k_i*neuron->cl_i*neuron->cl_i)/(kna*kk*kcl*kcl);
//	double denom = 1.0+neuron->k_o*neuron->cl_o*neuron->cl_o*neuron->na_o/(kna*kcl*kk*kcl);
//	denom = denom*(1.0+neuron->k_i/kk)*(1.+neuron->cl_i/kcl)*(1.+neuron->cl_i/kcl)*(1.+neuron->na_i/kna);
//	denom = denom + (1.+neuron->k_i*neuron->cl_i*neuron->na_i*neuron->cl_i/(kna*kcl*kk*kcl))*(1.+neuron->na_o/kna)*(1.+neuron->k_o/kk)*(1.+neuron->cl_o/kcl)*(1.+neuron->cl_o/kcl);

  double v_nkcc1;
  v_nkcc1 = (pE1full*neuron->y_nkcc*kf_full - pE2full*(1-neuron->y_nkcc)*kb_full);
//	v_nkcc1 = numerator/denom;

	if(t<dt) {
//		neuron->p_nkcc =20.e3; //.05*.3e-8; // .5e-4; //5.*(neuron->na_atpase + neuron->Ina + neuron->Ih_na + neuron->Ina_nmda + neuron->Ina_ampa + neuron->Ina_nacax +neuron->Inaleak)/(v_nkcc1*F);  // set max transport rate for equilibrium at initial conditions
//		if (neuron->distance > 1.) neuron->p_nkcc =.5e3; //0.1e-9; //.1*.1e-8;
//		if( neuron->distance > 300.) neuron->p_nkcc = .1e3; 

		neuron->p_nkcc = -neuron->cl_kcc2/(v_nkcc1*F);
	}
	v_nkcc1 = v_nkcc1*neuron->p_nkcc; 
	if( (i <= 39 && i >= 35 ) || ( i<=27 && i>=23) ) v_nkcc1 = v_nkcc1;
	//if(neuron->distance < 10.) v_nkcc1 = 0.;
  //for testing:
  //v_nkcc1[seg] = 0.0;
	
	if(t<dt && i==0) cout << "Inkcc1,max: " << neuron->p_nkcc*F << endl;

	neuron->na_nkcc1 = (-1./5.)*v_nkcc1*F;
	neuron->k_nkcc1 = (-4./5.)*v_nkcc1*F;
  neuron->cl_nkcc1 = v_nkcc1*F;


}
