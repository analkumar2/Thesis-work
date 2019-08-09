
#include "vgated_channels.h"
#include "synaptic_channels.h"
#include "ion_transport.h"
#include "neuron_structures.h"
#include "update.h"
#include "nr.h"

void calc_leak(Compartment * neuron, double time, double dt){  // this function calculates max na/katpase and leak currents for each ion

  double erevk = (1.e3)*(R*T/F)*log(neuron->k_o/neuron->k_i);
  double erevna = (1.e3)*(R*T/F)*log(neuron->na_o/neuron->na_i);
  double erevcl = (1.e3)*(R*T/F)*log(neuron->cl_i/neuron->cl_o);

	neuron->gkleak = 3.e-3; //3.e-1;  // S/cm2
  neuron->Ikleak = neuron->gkleak*(neuron->Vm - erevk);
	
  if(time < dt) {
    neuron->Ik_total = neuron->Ikleak + neuron->k_kcc2 + neuron->k_nkcc1 + neuron->Ikdr + neuron->Ika + neuron->Ik_h + neuron->Ikm + neuron->Isahp + neuron->Imahp + neuron->Ik_nmda + neuron->Ik_ampa;
    neuron->Ina_total = neuron->na_nkcc1 + neuron->Ina_nacax +  neuron->Ina + neuron->Ih_na + neuron->Ina_nmda + neuron->Ina_ampa;
	neuron->Icl_total = neuron->cl_kcc2 + neuron->cl_nkcc1;
		
    // calculate Imax for ATPase to balance potassium currents at steady state:
    double km_k = 3.5;  // mM, pump affinity for extracellular potassium
    double km_na = 10.0; // mM, pump affinity for intracellular sodium
    double Apump_ss = pow((1.+(km_k/neuron->k_o)), -2.)*pow((1.+(km_na/neuron->na_i)), -3.);
    neuron->imax_atpase = -neuron->Ik_total/(-2.*Apump_ss);
    neuron->k_atpase = -2.*neuron->imax_atpase*Apump_ss;
    neuron->na_atpase = 3.*neuron->imax_atpase*Apump_ss;

    // sodium leak conductance
    neuron->Ina_total = neuron->Ina_total + neuron->na_atpase;
    neuron->gnaleak = -neuron->Ina_total/(neuron->Vm - erevna);

    // calculate chloride leak conductance
	if(fabs(neuron->Icl_total) < 1.e-15) neuron->gclleak = 0.0;
	else neuron->gclleak = -neuron->Icl_total/(neuron->Vm - erevcl);
  }
  else {
   // atpase(neuron, time, dt, );
  }
  
  neuron->Inaleak = neuron->gnaleak*(neuron->Vm - erevna);
  neuron->Iclleak = neuron->gclleak*(neuron->Vm - erevcl);


}

void calc_leak2(Compartment * neuron, double time, double dt){

	neuron->gkleak = 1.e-4;
	neuron->Ikleak = neuron->gkleak*(neuron->Vm - neuron->erev_k);
	
	if(time < dt) {
    neuron->Ik_total = neuron->Ikleak + neuron->k_kcc2 + neuron->k_nkcc1 + neuron->Ikdr + neuron->Ika + neuron->Ik_h + neuron->Ikm + neuron->Isahp + neuron->Imahp + neuron->Ik_nmda + neuron->Ik_ampa;
    neuron->Ina_total = neuron->na_nkcc1 +  neuron->Ina + neuron->Ih_na + neuron->Ina_nmda + neuron->Ina_ampa;
    neuron->Icl_total = neuron->cl_kcc2 + neuron->cl_nkcc1;
		
    // calculate Imax for ATPase to balance potassium currents at steady state:
    double km_k = 3.5;  // mM, pump affinity for extracellular potassium
    double km_na = 10.0; // mM, pump affinity for intracellular sodium
    double Apump_ss = pow((1.+(km_k/neuron->k_o)), -2.)*pow((1.+(km_na/neuron->na_i)), -3.);
    neuron->imax_atpase = -neuron->Ik_total/(-2.*Apump_ss);
    neuron->k_atpase = -2.*neuron->imax_atpase*Apump_ss;
    neuron->na_atpase = 3.*neuron->imax_atpase*Apump_ss;

    // calculate sodium currents at steady state, assumes no sodium leak
    neuron->Ina_total = neuron->Ina_total + neuron->na_atpase;
	nacax(neuron, time, dt);  // sets Imax_nacax based on total sodium currents at steady state
	neuron->Ina_total = neuron->Ina_total + neuron->Ina_nacax;
	neuron->gnaleak = 0.0; //-neuron->Ina_total*(neuron->Vm - neuron->erev_na);

    // calculate chloride leak conductance
	neuron->gclleak = -neuron->Icl_total/(neuron->Vm - neuron->erev_cl);
  }
  else {
    //atpase(neuron, time, dt);
	nacax(neuron, time, dt);
  }

	neuron->Inaleak = neuron->gnaleak*(neuron->Vm - neuron->erev_na);
	neuron->Iclleak = neuron->gclleak*(neuron->Vm - neuron->erev_cl);

}

void calc_leak3(Compartment * neuron, double time, double dt){

	neuron->gkleak = .5/neuron->rm; //1.e-4;
	neuron->Ikleak = neuron->gkleak*(neuron->Vm - neuron->erev_k);
	
	if(time < dt) {
    neuron->Ik_total = neuron->Ikleak + neuron->k_kcc2 + neuron->k_nkcc1 + neuron->Ikdr + neuron->Ika + neuron->Ik_h + neuron->Ikm + neuron->Isahp + neuron->Imahp + neuron->Ik_nmda + neuron->Ik_ampa;
    neuron->Ina_total = neuron->na_nkcc1 +  neuron->Ina_nacax + neuron->Ina + neuron->Ih_na + neuron->Ina_nmda + neuron->Ina_ampa;
    neuron->Icl_total = neuron->cl_kcc2 + neuron->cl_nkcc1;
		
    // calculate Imax for ATPase to balance potassium currents at steady state:
    double km_k = 3.5;  // mM, pump affinity for extracellular potassium
    double km_na = 10.0; // mM, pump affinity for intracellular sodium
    //double Apump_ss = pow((1.+(km_k/neuron->k_o)), -2.)*pow((1.+(km_na/neuron->na_i)), -3.);
	double numerator = (neuron->k_o*neuron->na_i - neuron->k_i*neuron->na_o)/(km_k*km_na);
	double denom = 1.0+neuron->k_o*neuron->na_i/(km_k*km_na);
	denom = denom*(1.0+neuron->k_i/km_k)*(1.+neuron->na_o/km_na);
	denom = denom + (1.+neuron->k_i*neuron->na_o/(km_k*km_na))*(1.+neuron->k_o/km_k)*(1.+neuron->na_i/km_na);
	double Apump_ss = numerator/denom;
	neuron->imax_atpase = -neuron->Ik_total/(-2.*Apump_ss);
    neuron->k_atpase = -2.*neuron->imax_atpase*Apump_ss;
    neuron->na_atpase = 3.*neuron->imax_atpase*Apump_ss;

    // calculate sodium currents at steady state, assumes no sodium leak
    neuron->Ina_total = neuron->Ina_total + neuron->na_atpase;
	neuron->gnaleak = -neuron->Ina_total/(neuron->Vm - neuron->erev_na);
	
	neuron->gclleak = -neuron->Icl_total/(neuron->Vm - neuron->erev_cl);
	}
	else {
		//atpase(neuron, time, dt);
	}
	
	neuron->Inaleak = neuron->gnaleak*(neuron->Vm - neuron->erev_na);
	neuron->Iclleak = neuron->gclleak*(neuron->Vm - neuron->erev_cl);
}

void set_leaks(Compartment * neuron, double time, double h){

	if(time<h){
		// at beginning, calculate leak conductances
		double resist = neuron->rm; // - (1./neuron->gh);
		neuron->gkleak = .5/resist; //.5/neuron->rm;
		neuron->gclleak = (neuron->Vm - resist*neuron->gkleak*neuron->erev_k - (1.-resist*neuron->gkleak)*neuron->erev_na)/(resist*neuron->erev_cl - resist*neuron->erev_na);
		//neuron->gclleak = (neuron->Vm - neuron->rm*neuron->gkleak*neuron->erev_k - (1.-neuron->rm*neuron->gkleak)*neuron->erev_na)/(neuron->rm*neuron->erev_cl - neuron->rm*neuron->erev_na);
		neuron->gnaleak = 5.*( (1./resist) - neuron->gkleak - neuron->gclleak);
		//neuron->gnaleak = ( (1./neuron->rm) - neuron->gkleak - neuron->gclleak);
	}
	
	neuron->Ikleak = neuron->gkleak*(neuron->Vm - neuron->erev_k);
	neuron->Inaleak = neuron->gnaleak*(neuron->Vm - neuron->erev_na);
	neuron->Iclleak = neuron->gclleak*(neuron->Vm - neuron->erev_cl);
}

void set_kclleaks(Compartment * neuron, double time, double h){
	double surface = 2.0*neuron->radius*PI*neuron->dx*1.e-8;  // surface area in cm2

	if(time<h){
		neuron->Ik_total =  neuron->k_kcc2 + neuron->k_nkcc1 + neuron->Ikdr + neuron->Ika + neuron->Ik_h + neuron->Ikm + neuron->Isahp + neuron->Imahp + neuron->Ik_nmda + neuron->Ik_ampa + neuron->k_atpase;
		neuron->Ikleak = -neuron->Ik_total;
		neuron->gkleak = neuron->Ikleak/(neuron->Vm-neuron->erev_k);
		if(neuron->gkleak < 0.) cout << "error in k leaks! distance = " << neuron->distance << endl;
		if( neuron->distance <= 1.) cout << "gkleak: " << neuron->gkleak*surface << endl;
		}
	
	neuron->Ikleak = neuron->gkleak*(neuron->Vm - neuron->erev_k);

	
}

void set_naleaks(Compartment * neuron, double time, double h){

	if(time<h){
		neuron->Ina_total = neuron->na_atpase + neuron->na_nkcc1 + neuron->Ina + neuron->Ih_na + neuron->Ina_nmda + neuron->Ina_ampa + neuron->Ina_nacax; 
		neuron->Inaleak = -neuron->Ina_total;
		neuron->gnaleak = neuron->Inaleak/(neuron->Vm-neuron->erev_na);
		if(neuron->gnaleak < 0.) cout << "error in na leaks!  distance = " << neuron->distance << endl; //neuron->gnaleak = 0.; //.2/(neuron->rm); // neuron->Inaleak/(neuron->Vm-neuron->erev_na);

	}
	
	neuron->Inaleak = neuron->gnaleak*(neuron->Vm-neuron->erev_na);

}

void set_clleak(Compartment * neuron, double time, double h){
	if(time<h){
		neuron->Icl_total = neuron->cl_nkcc1 + neuron->cl_kcc2;
		neuron->Iclleak = -neuron->Icl_total;
		neuron->gclleak=0.;
	    neuron->gclleak = neuron->Iclleak/(neuron->Vm-neuron->erev_cl);
		if (neuron->gclleak < 0.) cout << "error in cl leaks! distance = " << neuron->distance << endl;
		//if(neuron->distance > 10.) neuron->gclleak = 1./(neuron->rm);
	}
	neuron->Iclleak=neuron->gclleak*(neuron->Vm-neuron->erev_cl);
}

void update_ions(Compartment * neuron, double dt, double time){  // only updates Na, K, and Cl

	double surface = 2.0*neuron->radius*PI*neuron->dx*1.e-8;  // surface area in cm2
	double vol = PI*pow(neuron->radius, 2.)*neuron->dx*1.e-15; // volume in Liters

	
	// d[ion]/dt = sum(ion currents)*surfacearea/(F*Vol) + iontransport
	//neuron->na_i = neuron->na_i + -1.0*(1.e-3)*neuron->Ina_total*surface*dt/(F*vol);  // 1e-3 for time unit conversion (ms to s)
	//neuron->k_i = neuron->k_i + -1.0*(1.e-3)*(neuron->Ik_total)*surface*dt/(F*vol);
	neuron->cl_i = neuron->cl_i + (1.e-3)*(neuron->Icl_total)*surface*dt/(F*vol);
	
	// calculating changes of EXTRACELLULAR concentrations:              
	// based on Kager et al, 2000
	// surface is same as for intracellular calcs
	// volume of intersititial compt is 15% of the neuronal compartment 

//	if(neuron-> distance < 10.){
	double tauko;
		if(neuron->distance <10.) {
		vol = 0.15*vol;  // L
			tauko = 40.; //5e3;
		}
		if(neuron->distance >= 10.) 
		{
			vol = vol*(.15);// + 200.*exp(-neuron->radius/.1));
			tauko = 40.;
		}
		//neuron->na_o = neuron->na_o + 1.0*(1.e-3)*neuron->Ina_total*surface*dt/(F*vol);
	
	/* potassium buffering */
	double buffer_max = 50.0; // mM
	double k1 = 50.*0.0002; // ms-1
	double k2 = k1/(1+exp((neuron->k_o-15.)/-1.15)); //-12.0)/-1.15));  // update forward rate constant
	double k_o_new; 
	//double tauko = 2.*5000.;
	k_o_new= neuron->k_o + 1.0*(1.e-3)*(neuron->Ik_total)*surface*dt/(F*vol); 
	if(neuron->distance < 70.) k_o_new = neuron->k_o + 5.0*(1.e-3)*(neuron->Ik_total)*surface*dt/(F*vol); 

	double dummy_Ik = neuron->Ik_total - neuron->Ikm - neuron->Ikdr - neuron->Ika - neuron->Ik_h - neuron->Isahp - neuron->Imahp - neuron->Ikleak;
	double konew_dummy;
	if(time<dt) neuron->k_o_dummy = neuron->k_o;
	konew_dummy = neuron->k_o_dummy + 1.0*(1.e-3)*(dummy_Ik)*surface*dt/(F*vol); 
	if(neuron->distance < 70.) konew_dummy = neuron->k_o_dummy + 5.0*(1.e-3)*(dummy_Ik)*surface*dt/(F*vol); 
	
	if(k_o_new > 3.){
		k_o_new = k_o_new + dt*((3.0)-neuron->k_o)/(tauko);  //200
		

		// buffer the new k concentration:
//	
//		double k_buff = dt*(k1*(buffer_max-neuron->kbuffer) - k_o_new*neuron->kbuffer*k2);
		//cout << "dk from buffering: " << k_buff << endl;
//		if ( ( neuron->kbuffer + k_buff) >= buffer_max){
//			k_buff = buffer_max-neuron->kbuffer;
//		}
//		if ( (neuron->kbuffer + k_buff) <= 0.0){
//			k_buff = 0.0-neuron->kbuffer;
//		}
//		neuron->kbuffer = neuron->kbuffer + k_buff;  // update buffer concentration based on new K+ concentration 
//		k_o_new = k_o_new + k_buff;
	
	neuron->k_o = k_o_new;
	}
	
	if(konew_dummy > 3.){
		konew_dummy = konew_dummy+dt*(3.-neuron->k_o_dummy)/tauko;
		neuron->k_o_dummy = konew_dummy;
	}							  
									
	neuron->cl_o = neuron->cl_o - (1.e-3)*neuron->Icl_total*surface*dt/(F*vol);
	//}
//	}
		
}

void long_diffusion(Compartment * neuron, int i, double dt){

	double dna = 1.33;  // um2/ms
	double dk = 1.96; 
	double dcl = 2.08; 
	double dca = 0.6;

	double dcdx1, dcdx2, dcdx3;  // concentration gradients
	double sa1, sa2, sa3; // surface areas
	int neighb1, neighb2, neighb3;
	
	double vol = PI*neuron[i].radius*neuron[i].radius*neuron[i].dx; // um3
	double vol_ratio = 0.15; // ratio of intracellular to extracellular volume
	if(neuron[i].distance > 100.) vol_ratio = 0.22; // + 200.*exp(-neuron[i].radius/.1);	
	if(neuron[i].num_neighbors == 3){
		neighb1 = neuron[i].neighbor_list[0];
		neighb2 = neuron[i].neighbor_list[1];
		neighb3 = neuron[i].neighbor_list[2];
		
		// na+:
		dcdx1 = (neuron[i].na_i - neuron[neighb1].na_i)/(.5*(neuron[i].dx+neuron[neighb1].dx));  // mM/um
		dcdx2 = (neuron[i].na_i - neuron[neighb2].na_i)/(.5*(neuron[i].dx+neuron[neighb2].dx));
		dcdx3 = (neuron[i].na_i - neuron[neighb3].na_i)/(.5*(neuron[i].dx+neuron[neighb3].dx));
		sa1 = PI*(neuron[i].radius*neuron[i].radius + neuron[neighb1].radius*neuron[neighb1].radius)/2.;  // um2
		sa2 = PI*(neuron[i].radius*neuron[i].radius + neuron[neighb2].radius*neuron[neighb2].radius)/2.;
		sa3 = PI*(neuron[i].radius*neuron[i].radius + neuron[neighb3].radius*neuron[neighb3].radius)/2.;
		neuron[i].na_i_diff = -dt*dna*(dcdx1*sa1 + dcdx2*sa2 + dcdx3*sa3)/vol;  // mM
		dcdx1 = (neuron[i].na_o - neuron[neighb1].na_o)/(.5*(neuron[i].dx+neuron[neighb1].dx));  // mM/um
		dcdx2 = (neuron[i].na_o - neuron[neighb2].na_o)/(.5*(neuron[i].dx+neuron[neighb2].dx));
		dcdx3 = (neuron[i].na_o - neuron[neighb3].na_o)/(.5*(neuron[i].dx+neuron[neighb3].dx));
		sa1 = vol_ratio*sa1;  // um2
		sa2 = vol_ratio*sa2;
		sa3 = vol_ratio*sa3;
		neuron[i].na_o_diff = -dt*dna*(dcdx1*sa1 + dcdx2*sa2 + dcdx3*sa3)/(.15*vol);
		
		// k+:
		dcdx1 = (neuron[i].k_i - neuron[neighb1].k_i)/(.5*(neuron[i].dx+neuron[neighb1].dx));  // mM/um
		dcdx2 = (neuron[i].k_i - neuron[neighb2].k_i)/(.5*(neuron[i].dx+neuron[neighb2].dx));
		dcdx3 = (neuron[i].k_i - neuron[neighb3].k_i)/(.5*(neuron[i].dx+neuron[neighb3].dx));
		sa1 = PI*(neuron[i].radius*neuron[i].radius + neuron[neighb1].radius*neuron[neighb1].radius)/2.;  // um2
		sa2 = PI*(neuron[i].radius*neuron[i].radius + neuron[neighb2].radius*neuron[neighb2].radius)/2.;
		sa3 = PI*(neuron[i].radius*neuron[i].radius + neuron[neighb3].radius*neuron[neighb3].radius)/2.;
		neuron[i].k_i_diff = -dt*dk*(dcdx1*sa1 + dcdx2*sa2 + dcdx3*sa3)/vol;  // mM
		dcdx1 = (neuron[i].k_o - neuron[neighb1].k_o)/(.5*(neuron[i].dx+neuron[neighb1].dx));  // mM/um
		dcdx2 = (neuron[i].k_o - neuron[neighb2].k_o)/(.5*(neuron[i].dx+neuron[neighb2].dx));
		dcdx3 = (neuron[i].k_o - neuron[neighb3].k_o)/(.5*(neuron[i].dx+neuron[neighb3].dx));
		sa1 = vol_ratio*sa1;  // um2
		sa2 = vol_ratio*sa2;
		sa3 = vol_ratio*sa3;
		neuron[i].k_o_diff = -dt*dk*(dcdx1*sa1 + dcdx2*sa2 + dcdx3*sa3)/(.15*vol);
		
		// cl+:
		dcdx1 = (neuron[i].cl_i - neuron[neighb1].cl_i)/(.5*(neuron[i].dx+neuron[neighb1].dx));  // mM/um
		dcdx2 = (neuron[i].cl_i - neuron[neighb2].cl_i)/(.5*(neuron[i].dx+neuron[neighb2].dx));
		dcdx3 = (neuron[i].cl_i - neuron[neighb3].cl_i)/(.5*(neuron[i].dx+neuron[neighb3].dx));
		sa1 = PI*(neuron[i].radius*neuron[i].radius + neuron[neighb1].radius*neuron[neighb1].radius)/2.;  // um2
		sa2 = PI*(neuron[i].radius*neuron[i].radius + neuron[neighb2].radius*neuron[neighb2].radius)/2.;
		sa3 = PI*(neuron[i].radius*neuron[i].radius + neuron[neighb3].radius*neuron[neighb3].radius)/2.;
		neuron[i].cl_i_diff = -dt*dcl*(dcdx1*sa1 + dcdx2*sa2 + dcdx3*sa3)/vol;  // mM
		dcdx1 = (neuron[i].cl_o - neuron[neighb1].cl_o)/(.5*(neuron[i].dx+neuron[neighb1].dx));  // mM/um
		dcdx2 = (neuron[i].cl_o - neuron[neighb2].cl_o)/(.5*(neuron[i].dx+neuron[neighb2].dx));
		dcdx3 = (neuron[i].cl_o - neuron[neighb3].cl_o)/(.5*(neuron[i].dx+neuron[neighb3].dx));
		sa1 = vol_ratio*sa1;  // um2
		sa2 = vol_ratio*sa2;
		sa3 = vol_ratio*sa3;
		neuron[i].cl_o_diff = -dt*dcl*(dcdx1*sa1 + dcdx2*sa2 + dcdx3*sa3)/(.15*vol);
		
		// ca+:
		dcdx1 = (neuron[i].ca_i - neuron[neighb1].ca_i)/(.5*(neuron[i].dx+neuron[neighb1].dx));  // mM/um
		dcdx2 = (neuron[i].ca_i - neuron[neighb2].ca_i)/(.5*(neuron[i].dx+neuron[neighb2].dx));
		dcdx3 = (neuron[i].ca_i - neuron[neighb3].ca_i)/(.5*(neuron[i].dx+neuron[neighb3].dx));
		sa1 = PI*(neuron[i].radius*neuron[i].radius + neuron[neighb1].radius*neuron[neighb1].radius)/2.;  // um2
		sa2 = PI*(neuron[i].radius*neuron[i].radius + neuron[neighb2].radius*neuron[neighb2].radius)/2.;
		sa3 = PI*(neuron[i].radius*neuron[i].radius + neuron[neighb3].radius*neuron[neighb3].radius)/2.;
		neuron[i].ca_i_diff = -dt*dca*(dcdx1*sa1 + dcdx2*sa2 + dcdx3*sa3)/vol;  // mM
		dcdx1 = (neuron[i].ca_o - neuron[neighb1].ca_o)/(.5*(neuron[i].dx+neuron[neighb1].dx));  // mM/um
		dcdx2 = (neuron[i].ca_o - neuron[neighb2].ca_o)/(.5*(neuron[i].dx+neuron[neighb2].dx));
		dcdx3 = (neuron[i].ca_o - neuron[neighb3].ca_o)/(.5*(neuron[i].dx+neuron[neighb3].dx));
		sa1 = vol_ratio*sa1;  // um2
		sa2 = vol_ratio*sa2;
		sa3 = vol_ratio*sa3;
		neuron[i].ca_o_diff = -dt*dca*(dcdx1*sa1 + dcdx2*sa2 + dcdx3*sa3)/(.15*vol);
	}
	else if(neuron[i].num_neighbors == 2){
		neighb1 = neuron[i].neighbor_list[0];
		neighb2 = neuron[i].neighbor_list[1];
		
		// na+:
		dcdx1 = (neuron[i].na_i - neuron[neighb1].na_i)/(.5*(neuron[i].dx+neuron[neighb1].dx));  // mM/um
		dcdx2 = (neuron[i].na_i - neuron[neighb2].na_i)/(.5*(neuron[i].dx+neuron[neighb2].dx));
		sa1 = PI*(neuron[i].radius*neuron[i].radius + neuron[neighb1].radius*neuron[neighb1].radius)/2.;  // um2
		sa2 = PI*(neuron[i].radius*neuron[i].radius + neuron[neighb2].radius*neuron[neighb2].radius)/2.;
		neuron[i].na_i_diff = -dt*dna*(dcdx1*sa1 + dcdx2*sa2)/vol;  // mM
		dcdx1 = (neuron[i].na_o - neuron[neighb1].na_o)/(.5*(neuron[i].dx+neuron[neighb1].dx));  // mM/um
		dcdx2 = (neuron[i].na_o - neuron[neighb2].na_o)/(.5*(neuron[i].dx+neuron[neighb2].dx));
		sa1 = vol_ratio*sa1;  // um2
		sa2 = vol_ratio*sa2;
		neuron[i].na_o_diff = -dt*dna*(dcdx1*sa1 + dcdx2*sa2)/(.15*vol);
		
		// k+:
		dcdx1 = (neuron[i].k_i - neuron[neighb1].k_i)/(.5*(neuron[i].dx+neuron[neighb1].dx));  // mM/um
		dcdx2 = (neuron[i].k_i - neuron[neighb2].k_i)/(.5*(neuron[i].dx+neuron[neighb2].dx));
		sa1 = PI*(neuron[i].radius*neuron[i].radius + neuron[neighb1].radius*neuron[neighb1].radius)/2.;  // um2
		sa2 = PI*(neuron[i].radius*neuron[i].radius + neuron[neighb2].radius*neuron[neighb2].radius)/2.;
		neuron[i].k_i_diff = -dt*dk*(dcdx1*sa1 + dcdx2*sa2)/vol;  // mM
		dcdx1 = (neuron[i].k_o - neuron[neighb1].k_o)/(.5*(neuron[i].dx+neuron[neighb1].dx));  // mM/um
		dcdx2 = (neuron[i].k_o - neuron[neighb2].k_o)/(.5*(neuron[i].dx+neuron[neighb2].dx));
		sa1 = vol_ratio*sa1;  // um2
		sa2 = vol_ratio*sa2;
		neuron[i].k_o_diff = -dt*dk*(dcdx1*sa1 + dcdx2*sa2)/(.15*vol);
		
		// cl+:
		dcdx1 = (neuron[i].cl_i - neuron[neighb1].cl_i)/(.5*(neuron[i].dx+neuron[neighb1].dx));  // mM/um
		dcdx2 = (neuron[i].cl_i - neuron[neighb2].cl_i)/(.5*(neuron[i].dx+neuron[neighb2].dx));
		sa1 = PI*(neuron[i].radius*neuron[i].radius + neuron[neighb1].radius*neuron[neighb1].radius)/2.;  // um2
		sa2 = PI*(neuron[i].radius*neuron[i].radius + neuron[neighb2].radius*neuron[neighb2].radius)/2.;
		neuron[i].cl_i_diff = -dt*dcl*(dcdx1*sa1 + dcdx2*sa2)/vol;  // mM
		dcdx1 = (neuron[i].cl_o - neuron[neighb1].cl_o)/(.5*(neuron[i].dx+neuron[neighb1].dx));  // mM/um
		dcdx2 = (neuron[i].cl_o - neuron[neighb2].cl_o)/(.5*(neuron[i].dx+neuron[neighb2].dx));
		sa1 = vol_ratio*sa1;  // um2
		sa2 = vol_ratio*sa2;
		neuron[i].cl_o_diff = -dt*dcl*(dcdx1*sa1 + dcdx2*sa2)/(.15*vol);
		
		// ca+:
		dcdx1 = (neuron[i].ca_i - neuron[neighb1].ca_i)/(.5*(neuron[i].dx+neuron[neighb1].dx));  // mM/um
		dcdx2 = (neuron[i].ca_i - neuron[neighb2].ca_i)/(.5*(neuron[i].dx+neuron[neighb2].dx));
		sa1 = PI*(neuron[i].radius*neuron[i].radius + neuron[neighb1].radius*neuron[neighb1].radius)/2.;  // um2
		sa2 = PI*(neuron[i].radius*neuron[i].radius + neuron[neighb2].radius*neuron[neighb2].radius)/2.;
		neuron[i].ca_i_diff = -dt*dca*(dcdx1*sa1 + dcdx2*sa2)/vol;  // mM
		dcdx1 = (neuron[i].ca_o - neuron[neighb1].ca_o)/(.5*(neuron[i].dx+neuron[neighb1].dx));  // mM/um
		dcdx2 = (neuron[i].ca_o - neuron[neighb2].ca_o)/(.5*(neuron[i].dx+neuron[neighb2].dx));
		sa1 = vol_ratio*sa1;  // um2
		sa2 = vol_ratio*sa2;
		neuron[i].ca_o_diff = -dt*dca*(dcdx1*sa1 + dcdx2*sa2)/(.15*vol);
	}
	else {
		neighb1 = neuron[i].neighbor_list[0];
		
		// na+:
		dcdx1 = (neuron[i].na_i - neuron[neighb1].na_i)/(.5*(neuron[i].dx+neuron[neighb1].dx));  // mM/um
		sa1 = PI*(neuron[i].radius*neuron[i].radius + neuron[neighb1].radius*neuron[neighb1].radius)/2.;  // um2
		neuron[i].na_i_diff = -dt*dna*(dcdx1*sa1)/vol;  // mM
		dcdx1 = (neuron[i].na_o - neuron[neighb1].na_o)/(.5*(neuron[i].dx+neuron[neighb1].dx));  // mM/um
		sa1 = vol_ratio*sa1;  // um2
		neuron[i].na_o_diff = -dt*dna*(dcdx1*sa1)/(.15*vol);
		
		// k+:
		dcdx1 = (neuron[i].k_i - neuron[neighb1].k_i)/(.5*(neuron[i].dx+neuron[neighb1].dx));  // mM/um
		sa1 = PI*(neuron[i].radius*neuron[i].radius + neuron[neighb1].radius*neuron[neighb1].radius)/2.;  // um2
		neuron[i].k_i_diff = -dt*dk*(dcdx1*sa1)/vol;  // mM
		dcdx1 = (neuron[i].k_o - neuron[neighb1].k_o)/(.5*(neuron[i].dx+neuron[neighb1].dx));  // mM/um
		sa1 = vol_ratio*sa1;  // um2
		neuron[i].k_o_diff = -dt*dk*(dcdx1*sa1)/(.15*vol);
		
		// cl+:
		dcdx1 = (neuron[i].cl_i - neuron[neighb1].cl_i)/(.5*(neuron[i].dx+neuron[neighb1].dx));  // mM/um
		sa1 = PI*(neuron[i].radius*neuron[i].radius + neuron[neighb1].radius*neuron[neighb1].radius)/2.;  // um2
		neuron[i].cl_i_diff = -dt*dcl*(dcdx1*sa1 )/vol;  // mM
		dcdx1 = (neuron[i].cl_o - neuron[neighb1].cl_o)/(.5*(neuron[i].dx+neuron[neighb1].dx));  // mM/um
		sa1 = vol_ratio*sa1;  // um2
		neuron[i].cl_o_diff = -dt*dcl*(dcdx1*sa1)/(.15*vol);
		
		// ca+:
		dcdx1 = (neuron[i].ca_i - neuron[neighb1].ca_i)/(.5*(neuron[i].dx+neuron[neighb1].dx));  // mM/um
		sa1 = PI*(neuron[i].radius*neuron[i].radius + neuron[neighb1].radius*neuron[neighb1].radius)/2.;  // um2
		neuron[i].ca_i_diff = -dt*dca*(dcdx1*sa1)/vol;  // mM
		dcdx1 = (neuron[i].ca_o - neuron[neighb1].ca_o)/(.5*(neuron[i].dx+neuron[neighb1].dx));  // mM/um
		sa1 = vol_ratio*sa1;  // um2
		neuron[i].ca_o_diff = -dt*dca*(dcdx1*sa1)/(.15*vol);
	}
	
	
	
}


void update_currents(Compartment * neuron, int i, double time, double h, double stim_time, double duration, double pre_Vm, double soma, double rad, double lm){


	neuron->erev_k = (1.e3)*(R*T/F)*log(neuron->k_o/neuron->k_i);
  neuron->erev_na = (1.e3)*(R*T/F)*log(neuron->na_o/neuron->na_i);
  neuron->erev_cl = (1.e3)*(R*T/F)*log(neuron->cl_i/neuron->cl_o);
  neuron->erev_ca = (1.e3)*(R*T/(2.*F))*log(neuron->ca_o/neuron->ca_i);
  neuron->erev_hco3 = (1.e3)*(R*T/F)*log(neuron->hco3_i/neuron->hco3_o);

	
	if(time<h && i==1) {
		cout << "Ek: " << neuron->erev_k << endl;
		cout << "Ena: " << neuron->erev_na << endl;
		cout << "Egaba: " << .8*neuron->erev_cl +.2*neuron->erev_hco3 << endl;
		cout << "Ehco3: " << neuron->erev_hco3 << endl;
		cout << "Eh: " << .6*neuron->erev_k + .4*neuron->erev_na << endl;
		cout << "Eca: " << neuron->erev_ca << endl;
	}
	
  double istim = 0.0;
	
  double gaba_time1 = 100.;
	double gaba_time2 = gaba_time1+8000.;
  neuron->istim = 0.;
  
    
    /* set stimulus */

  double pre = -100.;

	
	/* INJECTION CURRENTS */
//	if(time>5. && i == 5) neuron->istim = (50.e-9)/(2.*PI*neuron->radius*neuron->dx*1.e-8);
//	if(time > 5.) neuron->istim = -.0005; //(-2.e-9)/(2.*PI*neuron->radius*neuron->dx*1.e-8);	
//	if(time > 200) neuron->istim = (-3.e-9)/(2.*PI*neuron->radius*neuron->dx*1.e-8);
//	if(time > 20 && time < 1700. && i== 1) { //80) {
//		neuron->istim = (-160.e-9)/(2.*PI*neuron->radius*neuron->dx*1.e-8);
//		//neuron->k_o = 6.;
//	}

	
	/* ACTION POTENTIAL TRAIN */
/*
	int spike_num = 50.;
	int inter_spike = 100.; // ms
	int start = 50.;
	if(time<h) neuron->nspike = 0.;
	int n = neuron->nspike;
	if( i >=68 && i<= 78 ){
		if ( time > start && time < (start + spike_num*inter_spike) ){
			if ( time >= (start + n*inter_spike) && time < (start + n*inter_spike + 2.) ){
			  //neuron->istim = (-2.e-6)/(2.*PI*neuron->radius*neuron->dx*1.e-8);
			  pre = 10.;
			}
			else if ( time > (start + n*inter_spike) && time < (start + n*inter_spike + 2. + h) ) {
				neuron->nspike++;
				cout << "spike number: " << neuron->nspike << endl;
			}
		}
	}
*/

	
	/* HIGH FREQUENCY GABA STIMULUS at time gaba_time1 */
    gaba_6state(neuron, time, h, gaba_time1, 100., i);
  
	nmda(neuron, time, h, pre);
	
    /* call functions to calculate all membrane currents */
	ih(neuron, time, h, i); 
	ka(neuron, time, h, i);
	kdr(neuron, time, h, i);
	km(neuron, time, h, i);
	neuron->Ikm = .5*neuron->Ikm;
	neuron->Ikdr = neuron->Ikdr;
	sahp(neuron, time, h, i);

	cat(neuron, time, h, i);
	car(neuron, time, h, i);
	cal(neuron, time, h, i);
	na6state(neuron, time, h, i);

    neuron->Imahp = 0.;

	
  kcc2_chang(neuron, time, h, i); 
  nkcc1_2state(neuron, time, h, i);
  ca_pump_new(neuron, time, h); // this is set by max transport rate from kager et al 2008
  nacax(neuron, time, h); // this max transport rate determined by Ca balance
  atpase(neuron, time, h, i);  // this max transport rate determined by K leak 
	
	

	set_naleaks(neuron, time, h);
	set_clleak(neuron, time, h);
	set_kclleaks(neuron, time, h);
 

	update_ca_conc(neuron, time, h);


  /* define total membrane current */
	neuron->Ik_total = neuron->Ikleak + neuron->k_kcc2  + neuron->Ikdr + neuron->Ika + neuron->Ik_h + neuron->Ikm + neuron->Isahp + neuron->Imahp + neuron->Ik_nmda + neuron->Ik_ampa + neuron->k_atpase + neuron->k_nkcc1; 
	neuron->Ina_total = neuron->Inaleak + neuron->na_atpase + neuron->na_nkcc1 + neuron->Ina + neuron->Ih_na + neuron->Ina_nmda + neuron->Ina_ampa + neuron->Ina_nacax; 
	neuron->Icl_total = neuron->Iclleak + neuron->cl_kcc2 + neuron->Icl_gaba + neuron->Icl_gaba_syn + neuron->cl_nkcc1;
	neuron->Ica_total = neuron->Icat + neuron->Icar + neuron->Ical +  neuron->icapump + neuron->Ica_nacax;
  neuron->Iion_total = neuron->istim + neuron->Ik_total + neuron->Ina_total + neuron->Icl_total + neuron->Ica_total + neuron->Ihco3_gaba + neuron->Ihco3_gaba_syn;
  update_ions(neuron, h, time);

}

void update_noions(Compartment * neuron, int i, double time, double h, double stim_time, double duration, double stim){

  neuron->erev_k = (1.e3)*(R*T/F)*log(neuron->k_o/neuron->k_i);
  neuron->erev_na = (1.e3)*(R*T/F)*log(neuron->na_o/neuron->na_i);
  neuron->erev_cl = (1.e3)*(R*T/F)*log(neuron->cl_i/neuron->cl_o);
  neuron->erev_ca = (1.e3)*(R*T/(2.*F))*log(neuron->ca_o/neuron->ca_i);
  
  double istim = 0.0;

  /* set stimulus */
  if( i == 1 && (time > stim_time && time < stim_time+duration)) istim = (-.007e-3)/(PI*neuron->radius*neuron->dx*2.e-8);  // from jun
  else istim = 0.0;
  double pre = -1.;
//  if(i==0){
//    if( time>=stim_time && time < stim_time+h ){
//      cout << "setting pre to release glut!" << endl;
//      pre = 1.;
//    }
//  }
  nmda(neuron, time, h, pre);
  
  /* call functions to calculate all membrane currents */
  na(neuron, time, h, i);
  kdr(neuron, time, h, i);
  ka(neuron, time, h, i);
  ih2(neuron, time, h);  // based on narayanan and johnston 2008
  km(neuron, time, h, i);
  sahp(neuron, time, h, i);
  mahp(neuron, time, h, i);
  cat(neuron, time, h, i);
  car(neuron, time, h, i);
  cal(neuron, time, h, i);
  nacax(neuron, time, h);
  update_ca_kager(neuron, h);

  /* define total membrane current */
  neuron->Ik_total = neuron->Ikdr + neuron->Ika + neuron->Ik_h + neuron->Ikm + neuron->Isahp + neuron->Imahp + neuron->Ik_nmda + neuron->Ik_ampa;
  neuron->Ina_total = neuron->Ina + neuron->Ih_na + neuron->Ina_nmda + neuron->Ina_ampa + neuron->Ina_nacax;
  //neuron->Icl_total = neuron->Iclleak + neuron->cl_kcc2 + neuron->cl_nkcc1;
  neuron->Ica_total = neuron->Icat + neuron->Icar + neuron->Ical + neuron->Ica_nacax + neuron->icapump;

  neuron->Iion_total = neuron->Ik_total + neuron->Ina_total + neuron->Ica_total + istim;

  double gleak = 1./neuron->rm;
  // mA/cm2
  if( time < h) {
    neuron->eleak = (-neuron->Iion_total - gleak*neuron->Vm)/(-gleak);
  }
  //neuron->eleak = -80.0;

  neuron->Iion_total = neuron->Iion_total + gleak*(neuron->Vm - neuron->eleak);

}

void update_noca(Compartment * neuron, int i, double time, double h, double stim_time, double duration, double stim){

  neuron->erev_k = (1.e3)*(R*T/F)*log(neuron->k_o/neuron->k_i);
  neuron->erev_na = (1.e3)*(R*T/F)*log(neuron->na_o/neuron->na_i);
  neuron->erev_cl = (1.e3)*(R*T/F)*log(neuron->cl_i/neuron->cl_o);
  
  double istim = 0.0;

  /* set stimulus */
  if( i == 1 && (time > stim_time && time < stim_time+duration)) istim = (-.006e-3)/(PI*neuron->radius*neuron->dx*2.e-8);  // from jun
  else istim = 0.0;
  double pre = -1.;
//  if(i==0){
//    if( time>=stim_time && time < stim_time+h ){
//      cout << "setting pre to release glut!" << endl;
//      pre = 1.;
//    }
//  }
  nmda(neuron, time, h, pre);
  
  /* call functions to calculate all membrane currents */
  na(neuron, time, h, i);
  kdr(neuron, time, h, i);
  ka(neuron, time, h, i);
  ih2(neuron, time, h);  // based on narayanan and johnston 2008
  km(neuron, time, h, i);
  sahp(neuron, time, h, i);
  mahp(neuron, time, h,i);

  kcc2_chang(neuron, time, h, i);
  nkcc1_2state(neuron, time, h, i);
  calc_leak3(neuron, time, h);   // sets Imax for the Na/K pump and leak conductances for Na and Cl at time 0, otherwise calcs leak currents, and atpase


  /* define total membrane current */
  neuron->Ik_total = neuron->Ikleak + neuron->k_atpase + neuron->k_kcc2 + neuron->k_nkcc1 + neuron->Ikdr + neuron->Ika + neuron->Ik_h + neuron->Ikm + neuron->Isahp + neuron->Imahp + neuron->Ik_nmda + neuron->Ik_ampa;
  neuron->Ina_total = neuron->Inaleak + neuron->na_atpase + neuron->na_nkcc1 +  neuron->Ina + neuron->Ih_na + neuron->Ina_nmda + neuron->Ina_ampa;
  neuron->Icl_total = neuron->Iclleak + neuron->cl_kcc2 + neuron->cl_nkcc1;
  neuron->Iion_total = neuron->Ik_total + neuron->Ina_total + neuron->Icl_total + istim;
  update_ions(neuron, h, time);

  if(time<h){
	if(i==0) {
		cout << "Ina active + nkcc1 + atpase: " << neuron->Ina_total - neuron->Ina_nacax - neuron->Inaleak << endl;
	}

 }
}

void set_resistance(Compartment * neuron, int neuron_size){

  /* set membrane resistance; sigmoidally decreasing frm soma to apical trunk*/
	double rmsoma = 200.*(1.e3); // ohm*cm2  -- in golding et al, varies from 75-125 *10^3 ohm*cm2
	double rmend = 200.*(1.e3); //31.*(1.e3); //12.*(1.e3);   // in golding et al, varies from 10-31
  double rm_dhalf = 200.; // um
  double steep = 50.; // um

  /* set axial resistance; also sigmoidally decreasing from soma to apical trunk */
	double risoma = 200.; //2.*50.; // ohm*cm, in Kager this is 100 ohm*cm
	double riend = 200.; //2.*35.;
  double ri_dhalf = 210.; // um

  for(int i=0; i<neuron_size; i++){
	  neuron[i].rm = rmsoma + (rmend-rmsoma)/(1.+exp((rm_dhalf-neuron[i].distance)/steep));
	neuron[i].ra = (risoma + (riend-risoma)/(1.+exp((ri_dhalf-neuron[i].distance)/steep)))*(neuron[i].dx*(1.e-4)/(PI*neuron[i].radius*neuron[i].radius*(1.e-8))); // ohm
//    if( i == 64 || i == 68 || i == 70 || i ==  
//		neuron[i].ra = (1.e4)*risoma*neuron[i].dx/(PI*neuron[i].radius*neuron[i].radius) + (1.e4)*(riend-risoma)/(1.+exp((ri_dhalf-neuron[i].distance)/steep))*neuron[i].dx/(PI*neuron[i].radius*neuron[i].radius);  // ra is resistance, in ohm
    //cout << "seg: " << i << " rm: " << neuron[i].rm << " ra: " << neuron[i].ra << endl;
  }
}

void cable_eqn_ie(Compartment * neuron, double h, int neuron_size){
	
	h = h*1.e-3; // to convert from ms to s

	const int N = neuron_size; //183;
	Mat_DP a(N,N);
	Vec_DP b(N);
	Vec_INT indx(N);
	DP d;
	
	for(int i = 0; i<=N-1; i++)
		for(int j = 0; j<=N-1; j++) a[i][j] = 0.0;
	
	double cm, itotal, area;
	double vj_new;
	int neighb1, neighb2, neighb3;
	
	for(int i = 0; i<N; i++){ // set up matrix a (lhs) and b (rhs)
		
		area = PI*2.*neuron[i].radius*neuron[i].dx*(1.e-8); // surface area of compartment in cm2
		cm = neuron[i].cm*area;
		itotal = neuron[i].Iion_total*area;
		

		
		if(neuron[i].num_neighbors == 3){
			neighb1 = neuron[i].neighbor_list[0];
			neighb2 = neuron[i].neighbor_list[1];
			neighb3 = neuron[i].neighbor_list[2];
			a[i][i] = -cm-h*2./(neuron[i].ra+neuron[neighb1].ra)-h*2./(neuron[i].ra+neuron[neighb2].ra)-h*2./(neuron[i].ra+neuron[neighb3].ra);
			a[i][neighb1] = h*2./(neuron[i].ra+neuron[neighb1].ra);
			a[i][neighb2] = h*2./(neuron[i].ra+neuron[neighb2].ra);
			a[i][neighb3] = h*2./(neuron[i].ra+neuron[neighb3].ra);
			b[i] = -cm*neuron[i].Vm + h*itotal; 	
			neuron[i].Iaxial = (neuron[i].Vm - neuron[neighb1].Vm)*2./(neuron[i].ra+neuron[neighb1].ra) + (neuron[i].Vm-neuron[neighb2].Vm)*2./(neuron[i].ra+neuron[neighb2].ra) + (neuron[i].Vm-neuron[neighb3].Vm)*2./(neuron[i].ra+neuron[neighb3].ra);
		}
		if(neuron[i].num_neighbors == 2){
			neighb1 = neuron[i].neighbor_list[0];
			neighb2 = neuron[i].neighbor_list[1];
			a[i][i] = -cm-h*2./(neuron[i].ra+neuron[neighb1].ra)-h*2./(neuron[i].ra+neuron[neighb2].ra);
			a[i][neighb1] = h*2./(neuron[i].ra+neuron[neighb1].ra);
			a[i][neighb2] = h*2./(neuron[i].ra+neuron[neighb2].ra);
			b[i] = -cm*neuron[i].Vm + h*itotal;
			neuron[i].Iaxial = (neuron[i].Vm - neuron[neighb1].Vm)*2./(neuron[i].ra+neuron[neighb1].ra) + (neuron[i].Vm-neuron[neighb2].Vm)*2./(neuron[i].ra+neuron[neighb2].ra);
		}
		if(neuron[i].num_neighbors == 1){
			neighb1 = neuron[i].neighbor_list[0];
			a[i][i] = -cm-h*2./(neuron[i].ra+neuron[neighb1].ra);
			a[i][neighb1] = h*2./(neuron[i].ra+neuron[neighb1].ra);
			b[i] = -cm*neuron[i].Vm + h*itotal;
			neuron[i].Iaxial = (neuron[i].Vm - neuron[neighb1].Vm)*2./(neuron[i].ra+neuron[neighb1].ra);
		}
	}
	
	// sparse matrix solver, as in jun's ca3_hippocampus code
	const int ITOL=1,ITMAX=75;
	const DP TOL=1.0e-9;
    int iter;
    DP err;
    Vec_DP x(N);
    const int NMAX = 2*N*N+1;
    Vec_INT ija(NMAX);
	Vec_INT *ija_p;
    Vec_DP sa(NMAX);
	Vec_DP *sa_p;
//    NR::sprsin(a,0.0,sa,ija);  // indexes the sparse matrix, changed threshold from 1.e-9 to 0.
	const int nsize=sa.size();

    ija_p = &ija;
    sa_p = &sa;
	
	for(int i=0; i<N; i++) x[i] = neuron[i].Vm; // first guess is the vector of Vm's from previous time step
  //  NR::linbcg(b,x,ITOL,TOL,ITMAX,iter,err, ija_p, sa_p);

// regular matrix solver:
	NR::ludcmp(a,indx,d);
	NR::lubksb(a,indx,b);
	

	for(int i = 0; i<=N-1; i++)
	{   
		neuron[i].Vm_new = b[i];  // b when using ludcmp, x when using linbcg
	}   


}

void cable_eqn_fe(Compartment * neuron, int seg, double h){

  double area = PI*neuron[seg].radius*2.*neuron[seg].dx*(1.e-8);
  double itotal = neuron[seg].Iion_total;  // mA/cm2
  double cm = neuron[seg].cm;  // F/cm2
  double vj, vjp1, vjp2, vjm1, vj_new;
  double rj, rjp1, rjp2, rjm1;
  int neighb1, neighb2, neighb3;
  h = h*1.e-3; //converts ms to s

  /* update for segment with 2 neighbors */
  if( neuron[seg].num_neighbors == 2 ){
    neighb1 = neuron[seg].neighbor_list[0];
    neighb2 = neuron[seg].neighbor_list[1];
    vj = neuron[seg].Vm;
    vjp1 = neuron[neighb1].Vm;
    vjm1 = neuron[neighb2].Vm;
    rj = neuron[seg].ra;
    rjp1 = neuron[neighb1].ra; // ohm
    rjm1 = neuron[neighb2].ra;
    itotal = neuron[seg].Iion_total*area; // mA
    cm = neuron[seg].cm * area;  // F
    vj_new = vj + (h/cm)*(2./(rjp1+rj))*(vjp1-vj) + (h/cm)*(2./(rjm1+rj))*(vjm1-vj) - (h/cm)*itotal;
	neuron[seg].Vm_new = vj_new;
  }

  /* update for segment with 3 neighbors (branch point) */
  if( neuron[seg].num_neighbors == 3 ){
    neighb1 = neuron[seg].neighbor_list[0];
    neighb2 = neuron[seg].neighbor_list[1];
    neighb3 = neuron[seg].neighbor_list[2];
    vj = neuron[seg].Vm;
    vjp1 = neuron[neighb1].Vm;
    vjp2 = neuron[neighb2].Vm;
    vjm1 = neuron[neighb3].Vm;
    rj = neuron[seg].ra;
    rjp1 = neuron[neighb1].ra;
    rjp2 = neuron[neighb2].ra;
    rjm1 = neuron[neighb3].ra;
    itotal = neuron[seg].Iion_total * area;
    cm = neuron[seg].cm * area;
    vj_new = vj + (h/cm)*(2./(rjp1+rj))*(vjp1-vj) + (h/cm)*(2./(rjp2+rj))*(vjp2-vj) + (h/cm)*(2./(rjm1+rj))*(vjm1-vj) - (h/cm)*itotal;
	neuron[seg].Vm_new = vj_new;
  }

  /* update for segment with 1 neighbor (boundary) */
  if( neuron[seg].num_neighbors == 1 ){
    neighb1 = neuron[seg].neighbor_list[0];
    vj = neuron[seg].Vm;
    vjp1 = neuron[neighb1].Vm;
    rj = neuron[seg].ra;
    rjp1 = neuron[neighb1].ra;
    itotal = neuron[seg].Iion_total * area;
    cm = neuron[seg].cm * area;
    vj_new = vj + (h/cm)*(2./(rjp1+rj))*(vjp1-vj) - (h/cm)*itotal;
	neuron[seg].Vm_new = vj_new;
  }


}

void cable_eqn_rk(Compartment * neuron, int seg, double h, int neuron_size){

  double k1, k2, k3, k4;
  double vj_new, vj_orig, vj, vjp1, vjp2, vjm1;
  double rj, rjp1, rjp2, rjm1;
  int neighb1, neighb2, neighb3;
  double itotal, cm, area;
	int numCompts = neuron_size; //183;

  area = PI*neuron[seg].dx*neuron[seg].radius*2.e-8; // surface area in cm2

  /* note units:
     Vm in mV
     Iion_total in mA/cm2
     ra in ohm
     cm in F/cm2
     area in cm2
     h in ms
  */

  /* update for segment with 2 neighbors */
  if( neuron[seg].num_neighbors == 2 ){
    neighb1 = neuron[seg].neighbor_list[0];
    neighb2 = neuron[seg].neighbor_list[1];
    vj = neuron[seg].Vm;
    vjp1 = neuron[neighb1].Vm;
    vjm1 = neuron[neighb2].Vm;
    rj = neuron[seg].ra;
    rjp1 = neuron[neighb1].ra; // ohm
    rjm1 = neuron[neighb2].ra;
    itotal = neuron[seg].Iion_total*area; // mA
    cm = neuron[seg].cm * area;  // F
    k1 = h*( (vjm1-vj)*2./(rjm1+rj) - (vj-vjp1)*2./(rjp1+rj) - itotal )/cm;
    vj_orig = vj;
    vj = vj_orig + .5*k1;
    k2 = h*( (vjm1-vj)*2./(rjm1+rj) - (vj-vjp1)*2./(rjp1+rj) - itotal )/cm;
    vj = vj_orig + .5*k2;
    k3 = h*( (vjm1-vj)*2./(rjm1+rj) - (vj-vjp1)*2./(rjp1+rj) - itotal )/cm;
    vj = vj_orig + k3;
    k4 = h*( (vjm1-vj)*2./(rjm1+rj) - (vj-vjp1)*2./(rjp1+rj) - itotal )/cm;
    vj_new = vj_orig + (1.e-3)*(1./6.)*(k1 + 2.*k2 + 2.*k3 + k4);
    neuron[seg].Vm_new = vj_new;
  }

  /* update for segment with 3 neighbors (branch point) */
  if( neuron[seg].num_neighbors == 3 ){
    neighb1 = neuron[seg].neighbor_list[0];
    neighb2 = neuron[seg].neighbor_list[1];
    neighb3 = neuron[seg].neighbor_list[2];
    vj = neuron[seg].Vm;
    vjp1 = neuron[neighb1].Vm;
    vjp2 = neuron[neighb2].Vm;
    vjm1 = neuron[neighb3].Vm;
    rj = neuron[seg].ra;
    rjp1 = neuron[neighb1].ra;
    rjp2 = neuron[neighb2].ra;
    rjm1 = neuron[neighb3].ra;
    itotal = neuron[seg].Iion_total * area;
    cm = neuron[seg].cm * area;
    k1 = h*( (vjm1-vj)*2./(rjm1+rj) - (vj-vjp1)*2./(rjp1+rj) - (vj-vjp2)*2./(rjp2+rj) - itotal )/cm;
    vj_orig = vj;
    vj = vj_orig + .5*k1;
    k2 = h*( (vjm1-vj)*2./(rjm1+rj) - (vj-vjp1)*2./(rjp1+rj) - (vj-vjp2)*2./(rjp2+rj) - itotal )/cm;
    vj = vj_orig + .5*k2;
    k3 = h*( (vjm1-vj)*2./(rjm1+rj) - (vj-vjp1)*2./(rjp1+rj) - (vj-vjp2)*2./(rjp2+rj) - itotal )/cm;
    vj = vj_orig + k3;
    k4 = h*( (vjm1-vj)*2./(rjm1+rj) - (vj-vjp1)*2./(rjp1+rj) - (vj-vjp2)*2./(rjp2+rj) - itotal )/cm;
    vj_new = vj_orig + (1.e-3)*(1./6.)*(k1 + 2.*k2 + 2.*k3 + k4);
    neuron[seg].Vm_new = vj_new;
  }

  /* update for segment with 1 neighbor (boundary) */
  if( neuron[seg].num_neighbors == 1 ){
    neighb1 = neuron[seg].neighbor_list[0];
    vj = neuron[seg].Vm;
    vjp1 = neuron[neighb1].Vm;
    rj = neuron[seg].ra;
    rjp1 = neuron[neighb1].ra;
    itotal = neuron[seg].Iion_total * area;
    cm = neuron[seg].cm * area;
    k1 = h*( (vjp1-vj)*2./(rjp1+rj) - itotal )/cm;
    vj_orig = vj;
    vj = vj_orig + .5*k1;
    k2 = h*( (vjp1-vj)*2./(rjp1+rj) - itotal )/cm;
    vj = vj_orig + .5*k2;
    k3 = h*( (vjp1-vj)*2./(rjp1+rj) - itotal )/cm;
    vj = vj_orig + k3;
    k4 = h*( (vjp1-vj)*2./(rjp1+rj) - itotal )/cm;
    vj_new = vj_orig + (1.e-3)*(1./6.)*(k1 + 2.*k2 + 2.*k3 + k4);
    neuron[seg].Vm_new = vj_new;
  }
}

void print_data(Datafile filelist, Compartment * neuron, Compartment * basket, double t, int neuron_size){

  /* data will be written as follows:
     first column = time;
     then one column for each segment
  */
    

	int num_segments = 1; //neuron_size; //183;

  fprintf(filelist.vmptr, "%.6f\t", t);
  fprintf(filelist.inaptr, "%.6f\t", t);
  fprintf(filelist.mnaptr, "%.6f\t", t);
  fprintf(filelist.hnaptr, "%.6f\t", t);
  fprintf(filelist.snaptr, "%.6f\t", t);
  fprintf(filelist.ikdrptr, "%.6f\t", t);
  fprintf(filelist.mkdrptr, "%.6f\t", t);
  fprintf(filelist.ikaptr, "%.6f\t", t);
  fprintf(filelist.ihptr, "%.6f\t", t);
  fprintf(filelist.mihptr, "%.6f\t", t);
  fprintf(filelist.sahpptr, "%.6f\t", t);
  fprintf(filelist.mahpptr, "%.6f\t", t);
  fprintf(filelist.ikmptr, "%.6f\t", t);
  fprintf(filelist.itotalptr, "%.6f\t", t);
  fprintf(filelist.icatptr, "%.6f\t", t);
	fprintf(filelist.icarptr, "%.6f\t", t);
	fprintf(filelist.icalptr, "%.6f\t", t);
  fprintf(filelist.ca_in_ptr, "%.6f\t", t);
  fprintf(filelist.glutptr, "%.6f\t", t);
  fprintf(filelist.nmdaptr, "%.6f\t", t);
  fprintf(filelist.ampaptr, "%.6f\t", t);
  fprintf(filelist.kiptr, "%.6f\t", t);
  fprintf(filelist.koptr, "%.6f\t", t);
  fprintf(filelist.naiptr, "%.6f\t", t);
  fprintf(filelist.naoptr, "%.6f\t", t);
  fprintf(filelist.cliptr, "%.6f\t", t);
  fprintf(filelist.cloptr, "%.6f\t", t);
  fprintf(filelist.kcc2ptr, "%.6f\t", t);
  fprintf(filelist.nkcc1ptr, "%.6f\t", t);
  fprintf(filelist.clleak, "%.6f\t", t);
  fprintf(filelist.nakpump, "%.6f\t", t);
  fprintf(filelist.capump, "%.6f\t", t);
  fprintf(filelist.nacax, "%.6f\t", t);
  fprintf(filelist.gaba, "%.6f\t", t);
  fprintf(filelist.gaba_cl, "%.6f\t", t);
  fprintf(filelist.gaba_open, "%.6f\t", t);
  fprintf(filelist.egaba, "%.6f\t", t);
  fprintf(filelist.gaba_conc, "%.6f\t", t);
  fprintf(filelist.ecl, "%.6f\t", t);
  fprintf(filelist.iktotal, "%.6f\t", t);
  fprintf(filelist.inatotal, "%.6f\t", t);
	fprintf(filelist.icltotal, "%.6f\t", t);
	fprintf(filelist.icatotal, "%.6f\t", t);
	fprintf(filelist.inaleak, "%.6f\t", t);
	fprintf(filelist.ikleak, "%.6f\t", t);
	fprintf(filelist.bc_vm, "%.6f\t%.12f\n", t, basket[0].Vm);
	fprintf(filelist.bc_ina, "%.6f\t%.12f\n", t, basket[0].Ina);
	fprintf(filelist.bc_ikdr, "%.6f\t%.12f\n", t, basket[0].Ikdr);
	fprintf(filelist.bc_ahp, "%.6f\t%.12f\n", t, basket[0].Isahp);
	fprintf(filelist.bc_cai, "%.6f\t%.12f\n", t, basket[0].ca_i);
	fprintf(filelist.bc_Isyn, "%.6f\t%.12f\n", t, basket[0].Ik_nmda + basket[0].Ik_ampa + basket[0].Ina_nmda + basket[0].Ina_ampa);
	fprintf(filelist.na6_0, "%.6f\t", t);
	fprintf(filelist.na6_1, "%.6f\t", t);
	fprintf(filelist.na6_3, "%.6f\t", t);
	fprintf(filelist.na6_4, "%.6f\t", t);
	fprintf(filelist.na6_5, "%.6f\t", t);
	fprintf(filelist.dk_kcc2, "%.6f\t", t);
	
  for( int i=0; i<num_segments; i++){
	fprintf(filelist.ecl, "%.12f\t", neuron[i].erev_cl);
	  fprintf(filelist.iktotal, "%.12f\t", neuron[i].Ik_total-neuron[i].k_atpase-neuron[i].k_nkcc1);
	  fprintf(filelist.inatotal, "%.12f\t", neuron[i].Ina_total);
	  fprintf(filelist.icltotal, "%.12f\t", neuron[i].Icl_total);
	  fprintf(filelist.icatotal, "%.12f\t", neuron[i].cai_Ltype);
	  fprintf(filelist.inaleak, "%.12f\t", neuron[i].Inaleak);
	  fprintf(filelist.ikleak, "%.12f\t", neuron[i].Ikleak);
	  fprintf(filelist.gaba, "%.15f\t", (neuron[i].Icl_gaba+neuron[i].Ihco3_gaba)); //*2.*PI*neuron[i].radius*neuron[i].dx*(1.e-8)); //gaba6_y[3]);
	fprintf(filelist.gaba_cl, "%.12f\t", neuron[i].Icl_gaba);
	  fprintf(filelist.gaba_open, "%.12f\t", neuron[i].ggaba); //*2.*PI*neuron[i].radius*neuron[i].dx*(1.e-8));
	fprintf(filelist.gaba_conc, "%.4f\t", neuron[i].gaba_conc_syn+neuron[i].gaba_conc);
	fprintf(filelist.egaba, "%.12f\t", .8*neuron[i].erev_cl+.2*neuron[i].erev_hco3);
	fprintf(filelist.capump, "%.12f\t", neuron[i].icapump); //*2.*PI*neuron[i].radius*neuron[i].dx*1.e-8);
	fprintf(filelist.nacax, "%.12f\t", neuron[i].Ina_nacax); // + neuron[i].Ina_nacax); //*2.*PI*neuron[i].radius*neuron[i].dx*1.e-8);
	  fprintf(filelist.nakpump, "%.12f\t", neuron[i].na_atpase); // + neuron[i].na_atpase);
    fprintf(filelist.clleak, "%.12f\t", neuron[i].Iclleak);
    fprintf(filelist.nkcc1ptr, "%.12f\t", neuron[i].cl_nkcc1);
    fprintf(filelist.kcc2ptr, "%.12f\t", neuron[i].k_kcc2);
    fprintf(filelist.cloptr, "%.12f\t", neuron[i].cl_o);
    fprintf(filelist.cliptr, "%.12f\t", neuron[i].cl_i);
    fprintf(filelist.naoptr, "%.12f\t", neuron[i].na_o);
    fprintf(filelist.naiptr, "%.12f\t", neuron[i].na_i);
    fprintf(filelist.koptr, "%.12f\t", neuron[i].k_o);
    fprintf(filelist.kiptr, "%.12f\t", neuron[i].k_i);
    fprintf(filelist.ampaptr, "%.12f\t", neuron[i].ampa_r);
    fprintf(filelist.glutptr, "%.2f\t", neuron[i].glut_conc);
    fprintf(filelist.nmdaptr, "%.12f\t", neuron[i].nmda_r);
	fprintf(filelist.icatptr, "%.12f\t", neuron[i].Icat); // + neuron[i].Icar + neuron[i].Ical)); //*2.*PI*neuron[i].radius*neuron[i].dx*1.e-8);
	fprintf(filelist.icarptr, "%.12f\t", neuron[i].Icar);
	fprintf(filelist.icalptr, "%.12f\t", neuron[i].Ical);
	fprintf(filelist.ca_in_ptr, "%.12f\t", neuron[i].ca_i);
    fprintf(filelist.ikmptr, "%.12f\t", neuron[i].Ikm);
	  fprintf(filelist.itotalptr, "%.15f\t", (neuron[i].Iion_total - neuron[i].istim)*2.*PI*neuron[i].radius*neuron[i].dx*(1.e-8) + neuron[i].Iaxial);
    fprintf(filelist.vmptr, "%.18f\t", neuron[i].Vm);
	  fprintf(filelist.inaptr, "%.12f\t", neuron[i].Ina); //*2.*PI*neuron[i].radius*neuron[i].dx*1.e-8);
	  fprintf(filelist.mnaptr, "%.12f\t", neuron[i].na5_y[2]); //pow(neuron->ina_m, 3));
	  fprintf(filelist.hnaptr, "%.12f\t", neuron[i].na5_y[1]);
	  fprintf(filelist.snaptr, "%.12f\t", neuron[i].na5_y[0]);
	  fprintf(filelist.ikdrptr, "%.12f\t", neuron[i].Ikdr);
	  fprintf(filelist.mkdrptr, "%.12f\t", neuron[i].na5_y[4]);
	  fprintf(filelist.ikaptr, "%.12f\t", neuron[i].Ika);
	  fprintf(filelist.ihptr, "%.12f\t",neuron[i].Ik_h+neuron[i].Ih_na);
	  fprintf(filelist.mihptr, "%.12f\t", neuron[i].na5_y[3]); //ih_m);
    fprintf(filelist.sahpptr, "%.12f\t", neuron[i].Isahp);
    fprintf(filelist.mahpptr, "%.12f\t", neuron[i].Imahp);
	  fprintf(filelist.na6_0, "%.12f\t", neuron[i].na6_y[0]);
	  fprintf(filelist.na6_1, "%.12f\t", neuron[i].na6_y[1]);
	  fprintf(filelist.na6_3, "%.12f\t", neuron[i].na6_y[3]);
	  fprintf(filelist.na6_4, "%.12f\t", neuron[i].na6_y[4]);
	  fprintf(filelist.na6_5, "%.12f\t", neuron[i].na6_y[5]);
	  fprintf(filelist.dk_kcc2, "%.12f\t", neuron[i].k_o_dummy);
  }
  fprintf(filelist.ecl, "\n");
  fprintf(filelist.iktotal, "\n");
  fprintf(filelist.inatotal, "\n");
  fprintf(filelist.icltotal, "\n");
	fprintf(filelist.icatotal, "\n");
	fprintf(filelist.gaba, "\n");
  fprintf(filelist.gaba_cl, "\n");
  fprintf(filelist.gaba_conc, "\n");
  fprintf(filelist.gaba_open, "\n");
  fprintf(filelist.egaba, "\n");
  fprintf(filelist.capump, "\n");
  fprintf(filelist.nacax, "\n");
  fprintf(filelist.nakpump, "\n");
  fprintf(filelist.clleak, "\n");
  fprintf(filelist.nkcc1ptr,"\n");
  fprintf(filelist.kcc2ptr, "\n");
  fprintf(filelist.cloptr, "\n");
  fprintf(filelist.cliptr, "\n");
  fprintf(filelist.naoptr, "\n");
  fprintf(filelist.naiptr, "\n");
  fprintf(filelist.kiptr, "\n");
  fprintf(filelist.koptr, "\n");
  fprintf(filelist.ampaptr, "\n");
  fprintf(filelist.glutptr, "\n");
  fprintf(filelist.nmdaptr, "\n");
  fprintf(filelist.icatptr, "\n");
	fprintf(filelist.inaleak, "\n");
	fprintf(filelist.ikleak, "\n");
	fprintf(filelist.icarptr, "\n");
	fprintf(filelist.icalptr, "\n");
  fprintf(filelist.ca_in_ptr, "\n");
  fprintf(filelist.ikmptr, "\n");
  fprintf(filelist.itotalptr, "\n");
  fprintf(filelist.vmptr, "\n");
  fprintf(filelist.inaptr, "\n");
  fprintf(filelist.mnaptr, "\n");
  fprintf(filelist.hnaptr, "\n");
  fprintf(filelist.snaptr, "\n");
  fprintf(filelist.ikdrptr, "\n");
  fprintf(filelist.mkdrptr, "\n");
  fprintf(filelist.ikaptr, "\n");
  fprintf(filelist.ihptr, "\n");
  fprintf(filelist.mihptr, "\n");
  fprintf(filelist.sahpptr, "\n");
  fprintf(filelist.mahpptr, "\n");
	fprintf(filelist.na6_0, "\n");
	fprintf(filelist.na6_1, "\n");
	fprintf(filelist.na6_3, "\n");
	fprintf(filelist.na6_4, "\n");
	fprintf(filelist.na6_5, "\n");
	fprintf(filelist.dk_kcc2, "\n");
}

