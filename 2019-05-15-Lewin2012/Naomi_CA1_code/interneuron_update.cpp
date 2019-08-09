/*
 *  interneuron_update.cpp
 *  
 *
 *  Created by naomi on 12/9/09.
 *  Copyright 2009 Weill Medical College of Cornell University. All rights reserved.
 *
 */

#include "neuron_structures.h"
#include "synaptic_channels.h"
#include "nr.h"
#include "interneuron_update.h"

void bc_Ina(Compartment * neuron, double t, double dt){

	// Ina = gbar*m2*h*s*(V-Ena) or gbar*m2*h*s*GHK if gbar is converted to permeability
	
	
	double q = 96487./(8.315*(273.15+celsius));
	double tau_m = 0.05; // ms
	double tau_h = 0.5; // ms
	double tau_s = 0.00333*exp(0.0024*(neuron->Vm+60.)*q)/(1.0+exp(0.0012*(neuron->Vm+60.)*q));
	

	
	double minf, hinf;

	minf = 1./(1.+exp(-(neuron->Vm + 44.)/3.));
	hinf = 1./(1.+exp((neuron->Vm + 49.)/3.5));
	tau_h = 1.0;

	double sinf = (1.0+exp((neuron->Vm+60.)/2.))/(1.0+exp((neuron->Vm+60.)/2.0));
	
	
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
	

	
	double gbar; 
    gbar = .5*.004; //7.0*(1.e-3); ///100.; //1.53; 

	if(neuron->ina_h < 0.0) neuron->ina_h = 0.0;
	
	double gna = gbar*pow(neuron->ina_m, 2)*neuron->ina_h*neuron->ina_s; 
	neuron->Ina = gna*(neuron->Vm - neuron->erev_na);
		
}

void bc_Ikdr(Compartment * neuron, double t, double dt){
	// Ikdr = gbar*m2*(Vm-Ek) or gbar*m2*GHK if gbar is converted to permeability
	
	double tau_m; // = 2.2; //ms
	double minf; // = 1./(1.+exp(-(neuron->Vm+42.)/2.));
	double gbar;
	

	minf = 1./(1.+exp(-(neuron->Vm+46.3)/3.));
	tau_m = 3.5;
	gbar = .25*.002; //1.4*(1.e-3); ///100.; // S/cm2

	if( t < dt ){
		neuron->ikdr_m = minf;
	}
	else {
		neuron->ikdr_m = neuron->ikdr_m + (1.-exp(-dt/tau_m))*(minf-neuron->ikdr_m);
	}
	

	double gkdr = gbar*pow(neuron->ikdr_m, 2);
	neuron->Ikdr = gkdr*(neuron->Vm - neuron->erev_k);
}

void bc_Ikm(Compartment * neuron, double t, double dt){
	
	double gbar = .2*0.001;

	double tadj = pow(2.3, (celsius-23.)/10.);
	//  double alpha = (1.e-3)*(neuron->Vm+30.)/(1.-exp(-(neuron->Vm+30.)/9.));
	//  double beta = -(1.e-3)*(neuron->Vm+30.)/(1.-exp((neuron->Vm+30.)/9.));
	//  double minf = alpha/(alpha+beta);
	//  double tau = 1./(alpha+beta);
	
	double alpha = exp(-(neuron->Vm+6.)/8.);
	double minf = 1./(1.+alpha);
	double tau = 30.; //50.; 
	
	if(t < dt ) neuron->ikm_m = minf;
	else neuron->ikm_m = neuron->ikm_m + (1.-exp(-dt*tadj/tau))*(minf-neuron->ikm_m);
	
	neuron->Ikm = tadj*gbar*neuron->ikm_m*(neuron->Vm - neuron->erev_k);
}

void bc_Ica(Compartment * neuron, double t, double dt){
	
	/* interneuron v-gated current from Jeong and Gutkin 2007 */
	double gbar = 1.e-5;
	double minf = 1./( 1.+exp(-(neuron->Vm+25.)/2.5));
	double eca = (1.e3)*(R*T/(2.*F))*log(neuron->ca_o/neuron->ca_i);
	neuron->Icar = gbar*minf*(neuron->Vm-eca);
	
	/* update calcium concentration */
	double epsilon = 0.002; // uMcm2/msuA
	double tau = 80.; // ms
	double dca = dt*(-epsilon*neuron->Icar - (neuron->ca_i-5.e-5)/tau);
	neuron->ca_i = neuron->ca_i + dca;
	
}

void bc_Iahp(Compartment * neuron, double t, double dt){
	
	/* interneuron ahp current from Jeong and Gutkin 2007 */
	double gbar = 10.*.5e-5;
	neuron->Isahp = gbar*(neuron->ca_i/(1.+neuron->ca_i))*(neuron->Vm-neuron->erev_k);
}

void basket_update(Compartment * bc, double t, double dt, double pc_vm){

	pc_vm = -100.;
	
	bc->erev_k = (1.e3)*(R*T/F)*log(bc->k_o/bc->k_i);
	bc->erev_na = (1.e3)*(R*T/F)*log(bc->na_o/bc->na_i);	
	bc->gmax_nmda = .0005;
	bc->gmax_ampa = .0005;
	
	bc_Ina(bc, t, dt);
	bc_Ikdr(bc, t, dt);
	bc_Ikm(bc, t, dt);
	bc_Ica(bc, t, dt);
	bc_Iahp(bc, t, dt);
	nmda(bc, t, dt, pc_vm);
	
	/* generalized leak current */
	double gleak = .00005;
	double eleak = -67.;
	double Ileak = gleak*(bc->Vm - eleak);
	
	double Istim = 0.; 
	if(t>300. && t<700.) Istim=-0.007;
	//else if(t>10000. && t<10400.) Istim = -0.007;
	else Istim = 0.;
	bc->Iion_total = bc->Ina + bc->Ikdr + bc->Icar + bc->Isahp + bc->Ikm + Ileak + Istim + bc->Ik_nmda + bc->Ina_nmda + bc->Ik_ampa + bc->Ina_ampa;
	double vm_new = bc->Vm + (dt*1.e-3)*(-1.*bc->Iion_total/bc->cm);
	bc->Vm = vm_new;
}