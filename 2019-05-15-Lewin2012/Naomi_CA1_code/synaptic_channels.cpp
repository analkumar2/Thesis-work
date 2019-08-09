#include <iostream>
#include <math.h>
#include <stdio.h>
#include <cstdlib>
#include <cstdio>
#include "synaptic_channels.h"
#include "neuron_structures.h"

using namespace std;

void set_conductances(Compartment * neuron, int neuron_size){

	for(int i=0; i<neuron_size; i++){
		
		double gaba_scale;
		gaba_scale = 4.*neuron[i].distance/100.; 
		neuron[i].gmax_nmda = gaba_scale*3./2000.;
		neuron[i].gmax_ampa = 0.; // neuron[i].gmax_nmda/0.6; //this is true only for soma/apical trunk
		// scale gaba conductance densities based on pettit and augustine j neurophys 2000 with a distribution of synapses based on megias et al 2001
		double gsyn;
        gsyn = 1.0;
        if( i>=64 ) { // apical dendrites
			gaba_scale = 4.*neuron[i].distance/200.;  // G. Augustine et al. found 4*density at dendrites 200 um from soma
            if(gaba_scale > 4.) gaba_scale = 4.;
			if( neuron[i].distance <= 150. ) { // radiatum/thick/proximal; main branch
				neuron[i].gmax_gaba_a =gaba_scale* (1.7*gsyn)*(1.e-9)*neuron[i].dx/(2.*PI*neuron[i].radius*neuron[i].dx*(1.e-8)); // S/cm2; 
			}
			else if( neuron[i].distance > 150. && neuron[i].radius > 0.75 ) { //radiatum/thick/medial
				neuron[i].gmax_gaba_a = gaba_scale*(0.5*gsyn)*(1.e-9)*neuron[i].dx/(2.*PI*neuron[i].radius*neuron[i].dx*(1.e-8));
			}
			else if( neuron[i].distance < 350. ) { //radiatum/thick/distal and thin obliques (since all are distance > 150 and radius < .75)
				neuron[i].gmax_gaba_a = gaba_scale*(0.15*gsyn)*(1.e-9)*neuron[i].dx/(2.*PI*neuron[i].radius*neuron[i].dx*(1.e-8));
			}
			else { 
				neuron[i].gmax_gaba_a =gaba_scale* (0.12*gsyn)*(1.e-9)*neuron[i].dx/(2.*PI*neuron[i].radius*neuron[i].dx*(1.e-8)); 
			} // lacunosum-mol.
		}
		else if ( i>4 ) { // basal dendrites stratum oriens
			gaba_scale = 2.*neuron[i].distance/200.;
			if (neuron[i].distance > 350. ) { //oriens, distal
				neuron[i].gmax_gaba_a = (0.1*gsyn)*(1.e-9)*neuron[i].dx/(2.*PI*neuron[i].radius*neuron[i].dx*(1.e-8));
			}
			else 
				neuron[i].gmax_gaba_a = (0.61*gsyn)*(1.e-9)*neuron[i].dx/(2.*PI*neuron[i].radius*neuron[i].dx*(1.e-8));  // oriens, proximal
		}
		else neuron[i].gmax_gaba_a = (1.7*gsyn)*(1.e-9)*neuron[i].dx/(2.*PI*neuron[i].radius*neuron[i].dx*(1.e-8));  // S/cm2, somatic inhibitory input from megias et al 2001 

	}
}

void nmda(Compartment * neuron, double t, double dt, double pre){

  double cmax = 1.; // mM
  double cdur = 1.0; // ms
  double deadtime = 10.; // min time between release
  double thresh = -10.;
	double alpha = 10.; // ms-1 mM-1, forward binding 
	double alpha_ampa = 10.; //1.1;
	double beta = 5.*0.0125; // ms-1, backward unbinding
  double erev_k = (1.e3)*(R*T/F)*log(neuron->k_o/neuron->k_i); // epsp reversal potential, mV
  double erev_na = (1.e3)*(R*T/F)*log(neuron->na_o/neuron->na_i);
  double eta = 0.28; //0.33; // mM-1
  double mag = 1.; // mM
  double gamma = 0.062; // mV-1

  // to calculate ampa kinetics, slightly different binding/unbinding:
	double ampa_beta = 0.5; //0.19; //0.5;

  double g;
  double rinf = cmax*alpha/(cmax*alpha + beta);
  double rtau = 1./(alpha*cmax + beta);
  double ampa_rinf = cmax*alpha/(cmax*alpha + ampa_beta);
  double ampa_rtau = 1./(alpha*cmax + ampa_beta);

  if(t<dt){
    neuron->nmda_r = 0.;
    neuron->ampa_r = 0.;
    neuron->glut_conc = 0.;
  }  

  double q = ((t- neuron->lastrelease)-cdur); // time since last release ended
  if(q>deadtime){ // ready for another release
    if( pre > thresh ){
      neuron->glut_conc = cmax;
      neuron->nmda_r0 = neuron->nmda_r;
      neuron->ampa_r0 = neuron->ampa_r;
      neuron->lastrelease = t;
      cout << "last release set to: " << neuron->lastrelease << endl;
    }
    else if(q>0. && neuron->glut_conc==cmax){
      neuron->nmda_r1 = neuron->nmda_r;
      neuron->ampa_r1 = neuron->ampa_r;
      neuron->glut_conc = 0.;
      cout << "set glut conc back to zero, time: " << t << endl;
    }
    }
    else if(q>0) neuron->glut_conc = 0.;
  
  if( neuron->glut_conc>0.){ // transmitter is being released,receptor binding/activation 
    neuron->nmda_r = rinf + (neuron->nmda_r0-rinf)*exptable(-(t-neuron->lastrelease)/rtau);
    neuron->ampa_r = ampa_rinf + (neuron->ampa_r0-ampa_rinf)*exptable(-(t-neuron->lastrelease/ampa_rtau));
	  neuron->nmda_r1 = neuron->nmda_r;
      neuron->ampa_r1 = neuron->ampa_r; 
  }
  else { // no release, receptor unbinding
	  neuron->nmda_r = neuron->nmda_r1*exp(-beta*(t-(neuron->lastrelease+cdur)));
	  neuron->ampa_r = neuron->ampa_r1*exp(-ampa_beta*(t-(neuron->lastrelease+cdur)));
  }
 

  g = neuron->gmax_nmda*neuron->nmda_r/(1.+eta*mag*exp(-gamma*neuron->Vm));
  //if(g < 0.0000000000001 || g > 1.e10 ) g = 0.0;
  neuron->Ik_nmda = 0.43*g*(neuron->Vm - erev_k);
  neuron->Ina_nmda = 0.57*g*(neuron->Vm - erev_na);

  g = neuron->gmax_ampa*neuron->ampa_r;
  //if(g < 0.0000000000001 || g > 1.e10 ) g = 0.0;
  neuron->Ik_ampa = 0.43*g*(neuron->Vm - erev_k);
  neuron->Ina_ampa = 0.57*g*(neuron->Vm - erev_na);
  
}

double exptable(double x){

  double ans;
  if( x > -10. && x < 10.) ans = exp(x);
  else ans = 0.;

  return ans;
}

void gaba_6state(Compartment * neuron, double t, double dt, double t0, double duration, int j){

	Vec_DP y(6);
	Vec_DP dydt(6);
	Vec_DP yout(6);
	Vec_DP dfdt(6);
	Mat_DP dfdy(6,6);
	Vec_DP yE(6);
	Vec_DP dydtE(6);
	Vec_DP youtE(6);
	Vec_DP dfdtE(6);
	Mat_DP dfdyE(6,6);
	
	int spike_num;
	
	for(int i=0; i<6; i++) {
		yout[i]=neuron->gaba6_y[i];
		youtE[i]=neuron->gaba6_yE[i];
	}
	
	double cdur = 1.0;
	// define gaba concentration:
	// for a single pulse:
	neuron->gaba_conc=0.;
	neuron->gaba_concE=0.;
	double tau_syn = 0.1; // ms; tau for clearance of GABA from cleft; from Mozrzymas, Neuropharmacology 2004
	spike_num = (t-t0)/10; //10;
	if(spike_num > 40) spike_num = 40;
	if(t>t0 && t< (t0+400.)) {
		neuron->gaba_conc = 1.*exp(-(t-(t0 + 10.*spike_num))/tau_syn); //10.*spike_num))/tau_syn); // mM
		neuron->gaba_concE = 1.*exp(-(t-(t0 + 10.*spike_num))/tau_syn);
	}


	/* set initial conditions */
	if (t<dt) {
		yout[0] = 1.0;
		yout[1] = 0.0;
		yout[2] = 0.0;
		yout[3] = 0.0;
		yout[4] = 0.0;
		yout[5] = 0.0;
		youtE[0] = 1.0;
		youtE[1] = 0.0;
		youtE[2] = 0.0;
		youtE[3] = 0.0;
		youtE[4] = 0.0;
		youtE[5] = 0.0;
	}
	else  {
		
		NR::jacobn_gaba6(t, yout, dfdt, dfdy, neuron->gaba_conc, j);  // dfdy becomes the jacobian matrix
		imp_euler(dfdy, yout, dt);
	}
	for(int i=0; i<6; i++) {
		neuron->gaba6_y[i] = yout[i]; 
	}

	double open_ch;
	double open_ch_extra;
	
	if(j>=64){
		open_ch = neuron->gaba6_y[3]; 	}
	else open_ch = neuron->gaba6_y[3];

	
    
    neuron->ggaba = 0.5*neuron->gmax_gaba_a*(open_ch); // .5*      

	
    
    if( neuron->distance < 150.) neuron->ggaba = 0.5*neuron->gmax_gaba_a*open_ch;  // control: .2
	
    /* select for only apical and somatic gaba conductance; set basal and axon to zero */
	if( (j <= 39 && j >= 35 ) || ( j<=27 && j>=25) ) neuron->ggaba = 0.; //gmax_gaba_a = 0.;
	else if( j>=5 && j < 64) neuron->ggaba = 0.; //neuron->gmax_gaba_a = 0.; //neuron->gmax_gaba_a*open_ch;
	
    
    neuron->ggaba = 1.0*neuron->ggaba; 
	double pcl = 0.8*neuron->ggaba;
	double phco3 = 0.2*neuron->ggaba;
	neuron->Icl_gaba = pcl*(neuron->Vm - neuron->erev_cl);
	neuron->Ihco3_gaba = phco3*(neuron->Vm - neuron->erev_hco3);

	
}

void gaba_synapse(Compartment * neuron, double t, double dt, double pre, int j){
	
	Vec_DP y(6);
	Vec_DP dydt(6);
	Vec_DP yout(6);
	Vec_DP dfdt(6);
	Mat_DP dfdy(6,6);
	Vec_DP yE(6);
	Vec_DP dydtE(6);
	Vec_DP youtE(6);
	Vec_DP dfdtE(6);
	Mat_DP dfdyE(6,6);
	
	double tau_syn = 0.1;  // clearence of gaba from cleft
	double cdur = 1.;
	double deadtime = 5.;
	double cmax = 1.;
	double thresh = 0.;
	int spike_num;
	
	for(int i=0; i<6; i++) {
		yout[i]=neuron->gaba6_y[i];
		youtE[i]=neuron->gaba6_yE[i];
	}
	if(t<dt){
		neuron->gaba_conc_syn = 0.;
	}  
	
	double q = t- neuron->lastrelease; // time since last release ended
	if(q>deadtime){ // ready for another release
		if( pre > thresh ){
			neuron->gaba_conc_syn = cmax;
			neuron->lastrelease = t;
			cout << "last gaba release set to: " << neuron->lastrelease << endl;
		}
	}
	
	if( (t-neuron->lastrelease) <= cdur ) neuron->gaba_conc_syn = cmax;
	else neuron->gaba_conc_syn = 1.0*exp(-(t-neuron->lastrelease)/tau_syn);
	

	/* set initial conditions */
	if (t<dt) {
		yout[0] = 1.0;
		yout[1] = 0.0;
		yout[2] = 0.0;
		yout[3] = 0.0;
		yout[4] = 0.0;
		yout[5] = 0.0;
		youtE[0] = 1.0;
		youtE[1] = 0.0;
		youtE[2] = 0.0;
		youtE[3] = 0.0;
		youtE[4] = 0.0;
		youtE[5] = 0.0;
	}
	else {
		NR::jacobn_gaba6(t, yout, dfdt, dfdy, neuron->gaba_conc_syn, j);  // dfdy becomes the jacobian matrix
		imp_euler(dfdy, yout, dt);
		NR::jacobn_gaba6_extra(t, youtE, dfdtE, dfdyE, neuron->gaba_conc_syn, j);  // dfdy becomes the jacobian matrix
		imp_euler(dfdyE, youtE, dt);
	}
	
	for(int i=0; i<6; i++) {
		neuron->gaba6_y[i] = yout[i]; 
		neuron->gaba6_yE[i] = youtE[i];	
	}	
	double open_ch;

	
	open_ch = neuron->gaba6_y[3]+neuron->gaba6_yE[3];
	neuron->ggaba = .5*neuron->gmax_gaba_a*open_ch; 
	if( (j <= 39 && j >= 35 ) || ( j<=27 && j>=25) ) neuron->ggaba = 0.; //gmax_gaba_a = 0.;
	else if( j>=5 && j < 64) neuron->ggaba = 0.;
	else if(neuron->distance < 150. && neuron->distance > 10.) neuron->ggaba = .2*neuron->gmax_gaba_a*open_ch;
	else if(neuron->distance < 10.) neuron->ggaba = .1*neuron->gmax_gaba_a*neuron->gaba6_y[3]; 
	double pcl = 0.8*neuron->ggaba;
	double phco3 = 0.2*neuron->ggaba;
	neuron->Icl_gaba= pcl*(neuron->Vm - neuron->erev_cl);
	neuron->Ihco3_gaba = phco3*(neuron->Vm - neuron->erev_hco3);
	
	
}

void stylized_gaba(Compartment * neuron, double t, double stim_time ){
	

	double tau1 = 480.; // ms
	double tau2 = 600.; // ms
	double pcl = 0.;
	double phco3 = 0.;
	
	//neuron->ggaba = 30.* gmax*( -exp( -t/tau1 ) + exp( -t/tau2 ) );
	//if(neuron->radius > 0.35) neuron->gmax_gaba_a = neuron->gmax_gaba_a*0.9; //.5;
	if(t>stim_time){
		//if(neuron->distance < 100.) neuron->ggaba = 0.0001*neuron->gmax_gaba_a*( -exp( -(t-stim_time)/tau1 ) + exp( -(t-stim_time)/tau2 ) );
		//else 
		neuron->ggaba = neuron->gmax_gaba_a*( -exp( -(t-stim_time)/tau1 ) + exp( -(t-stim_time)/tau2 ) );
		double bicarb = 0.2;  // baseline = 0.2
		pcl = (1.-bicarb)*neuron->ggaba;
		phco3 = bicarb*neuron->ggaba;
	}
	neuron->Icl_gaba = pcl*(neuron->Vm - neuron->erev_cl);
	neuron->Ihco3_gaba = phco3*(neuron->Vm - neuron->erev_hco3);
}


// time derivsof 6 state GABAA receptor; one open state (y3) and 2 desensitized states (y4 and y5)
void NR::derivs_gaba6(const DP t, Vec_I_DP &y, Vec_O_DP &dydt, double gaba_conc){
	
	double kon = 1.e1; //0.01e6; // binding, mM-1 ms-1
	double koff = 0.103; //4.6e3; // unbinding, mM-1 ms-1
	double alph = 0.0613; //(1.0/53.)*1.0e3; // ms-1, 1/tau_decay
	double bet = 2.0833; //(1.0/4.0)*1.0e3; // ms-1, 1/tau_rise
	
	double rs = 0.004; 
	double rfast = 0.0952; 
	double ds = 0.0078;
	double df = 0.2024;
	
	dydt[0] = y[1]*koff - y[0]*2.*kon*gaba_conc;
	dydt[1] = y[0]*2.0*kon*gaba_conc + y[2]*2.*koff - y[1]*(koff+kon*gaba_conc);
	dydt[2] = y[1]*kon*gaba_conc + y[3]*alph + y[4]*rs + y[5]*rfast - y[2]*(2.*koff + ds + df + bet);
	dydt[3] = y[2]*bet  - y[3]*alph;
	dydt[4] = y[2]*ds - y[4]*rs;
	dydt[5] = y[2]*df - y[5]*rfast;
}

// jacobian matrix for 6-state gaba receptor
void NR::jacobn_gaba6(const DP t, Vec_I_DP &y, Vec_O_DP &dfdt, Mat_O_DP &dfdy, double gaba_conc, int j){

	double kon = 1.e1; //.2*1.e1; //0.01e6; // binding, mM-1 ms-1
	double koff = 0.103; //4.6e3; // unbinding, ms-1
	double alph = .0613; //.02; // ms-1, 1/tau_decay
	double bet = 2.0833; //(1.0/4.0)*1.0e3; // ms-1, 1/tau_rise
	
	double rs = 0.004; 
	double rfast = 0.0952; 
	double ds = 0.0078;
	double df = 0.2024;
	
	int i;
	int n=y.size();
	for (i=0; i<n; i++) dfdt[i]=0.0;
	dfdy[0][0] = -2.*kon*gaba_conc;
	dfdy[0][1] = koff;
	dfdy[0][2] = 0.0;
	dfdy[0][3] = 0.0;
	dfdy[0][4] = 0.0;
	dfdy[0][5] = 0.0;
	dfdy[1][0] = 2.*kon*gaba_conc;
	dfdy[1][1] = -(koff + kon*gaba_conc);
	dfdy[1][2] = 2.*koff;
	dfdy[1][3] = 0.0;
	dfdy[1][4] = 0.0;
	dfdy[1][5] = 0.0;
	dfdy[2][0] = 0.0;
	dfdy[2][1] = kon*gaba_conc;
	dfdy[2][2] = -(2.*koff + ds + df + bet);
	dfdy[2][3] = alph;
	dfdy[2][4] = rs;
	dfdy[2][5] = rfast;
	dfdy[3][0] = 0.0;
	dfdy[3][1] = 0.0;
	dfdy[3][2] = bet;
	dfdy[3][3] = -alph;
	dfdy[3][4] = 0.0;
	dfdy[3][5] = 0.0;
	dfdy[4][0] = 0.0;
	dfdy[4][1] = 0.0;
	dfdy[4][2] = ds;
	dfdy[4][3] = 0.0;
	dfdy[4][4] = -rs;
	dfdy[4][5] = 0.0;
	dfdy[5][0] = 0.0;
	dfdy[5][1] = 0.0;
	dfdy[5][2] = df;
	dfdy[5][3] = 0.0;
	dfdy[5][4] = 0.0;
	dfdy[5][5] = -rfast;
}

void NR::jacobn_gaba6_extra(const DP t, Vec_I_DP &y, Vec_O_DP &dfdt, Mat_O_DP &dfdy, double gaba_conc, int j){
	
	double kon = 1.e1; //.1*1.e1; //.2*1.e1; //0.01e6; // binding, mM-1 ms-1
	double koff = 0.103; //4.6e3; // unbinding, ms-1
	double alph = .01; // ms-1, 1/tau_decay
	double bet = .8; //(1.0/4.0)*1.0e3; // ms-1, 1/tau_rise
	
	double rs = 0.005; 
	double rfast = 0.; //0.0952; 
	double ds = 0.0029;
	double df = 0.; //0.2024;
	
	int i;
	int n=y.size();
	for (i=0; i<n; i++) dfdt[i]=0.0;
	dfdy[0][0] = -2.*kon*gaba_conc;
	dfdy[0][1] = koff;
	dfdy[0][2] = 0.0;
	dfdy[0][3] = 0.0;
	dfdy[0][4] = 0.0;
	dfdy[0][5] = 0.0;
	dfdy[1][0] = 2.*kon*gaba_conc;
	dfdy[1][1] = -(koff + kon*gaba_conc);
	dfdy[1][2] = 2.*koff;
	dfdy[1][3] = 0.0;
	dfdy[1][4] = 0.0;
	dfdy[1][5] = 0.0;
	dfdy[2][0] = 0.0;
	dfdy[2][1] = kon*gaba_conc;
	dfdy[2][2] = -(2.*koff + ds + df + bet);
	dfdy[2][3] = alph;
	dfdy[2][4] = rs;
	dfdy[2][5] = rfast;
	dfdy[3][0] = 0.0;
	dfdy[3][1] = 0.0;
	dfdy[3][2] = bet;
	dfdy[3][3] = -alph;
	dfdy[3][4] = 0.0;
	dfdy[3][5] = 0.0;
	dfdy[4][0] = 0.0;
	dfdy[4][1] = 0.0;
	dfdy[4][2] = ds;
	dfdy[4][3] = 0.0;
	dfdy[4][4] = -rs;
	dfdy[4][5] = 0.0;
	dfdy[5][0] = 0.0;
	dfdy[5][1] = 0.0;
	dfdy[5][2] = df;
	dfdy[5][3] = 0.0;
	dfdy[5][4] = 0.0;
	dfdy[5][5] = -rfast;
}

void imp_euler(Mat_O_DP &J, Vec_IO_DP &yt, const DP h){

	// given the Jacobian matrix J, use implicit euler method to solve for the vector yout
	// yt: rhs
	// yt is replaced by the vector at time t+h in the lubksb function
	
	int n = yt.size();
	Mat_DP matrix(n,n);
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			if(i==j) matrix[i][j] = (-J[i][j]*h)+1.; // diagonal elements
			else matrix[i][j] = -J[i][j]*h; // nondiagonal elements
		}
	}
	
	Vec_INT indx(n);
	DP d;
	NR::ludcmp(matrix, indx, d);
	NR::lubksb(matrix, indx, yt);

}

// LU decomposition matrix solver from
// Numerical Recipes in C++
/* "Given a matrix a[0...n-1][0...n-1], this routine replaces it by the LU 
decomposition of a rowwise permutation of itself.  a is input.  On output, 
it is arranged as in equation (2.3.14 [page 48]) . . ., indx[0...n-1] is an 
output vector that records the row permutation efected by the partial pivoting;
d is output as +/- 1 depending on whether the number of row interchanges was 
even or odd, respectively. This routine is used in combination with lubksb to solve linear
equations or invert a matrix." */
void NR::ludcmp(Mat_IO_DP &a, Vec_O_INT &indx, DP &d)
{
	const DP TINY=1.0e-20;  // a small number
	int i, imax, j, k;
	DP big, dum, sum, temp;
	
	int n=a.nrows();
	Vec_DP vv(n);			// vv stores the implicit scaling of each row
	d=1.0;					// no row interchanges yet
	for (i=0;i<n;i++) {		// loop over rows to get implicit scaling info
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}
	for (j=0;j<n;j++) {		// loop over columns of Crout's method
		for (i=0;i<j;i++) { // equation 2.3.12 except for i=j (p47)
			sum=a[i][j];
			for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big = 0.0;			// initialize for search for largest pivot element
		for (i=j; i<n; i++) { // this is i=j of 2.3.12 (p47) and i =j+1 . . . n-1 of 2.3.13
			sum=a[i][j];
			for (k=0;k<j;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ((dum=vv[i]*fabs(sum)) >= big) {  // is the figure of merit for the pivot better than the best so far ????
			big=dum;
			imax=i;
			}
		}
		if (j != imax) {	// do we need to interchange rows?
			for (k=0;k<n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			d = -d;			// change parity of d
			vv[imax]=vv[j]; // interchange scale factor
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;  // if pivot element is zero, matrix is singular; subs tiny for zero 
		if (j != n-1) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<n;i++) a[i][j] *= dum;
		}
	}
}

// backsubstitution matrix solver following LU decomposition
// from Numerical Recipes in C++
/* "Solves the set of n linear equations A.X = B.  Her a[0...n-1][0...n-1] is input, 
not as the matrix A but as its LU decomposition, determined by the routine ludcmp. 
indx[0...n-1]is input as the permutation vector returned by ludcmp.  b[0...n-1] is input
as the right-hand side vector B, and returns with the solution vector X.  a and indx 
are not modified by this routine and can be left in place for successive calls with 
different righthand sides b.  This routine takes into account the possibility that b
will begin with many zero elements, so it is efficient for use in matrix inversion." (p50) */
void NR::lubksb(Mat_I_DP &a, Vec_I_INT &indx, Vec_IO_DP &b)
{
	int i, ii=0, ip, j; 
	DP sum;
	
	int n=a.nrows();
	for (i=0; i<n; i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii != 0)
			for (j=ii-1; j<i; j++) sum -=a[i][j]*b[j];
		else if (sum != 0.0 )
			ii=i+1;
		b[i]=sum;
	}
	for (i=n-1; i>= 0; i--) {
		sum=b[i];
		for (j=i+1; j<n; j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}




void NR::rk4(Vec_I_DP &y, Vec_I_DP &dydx, const DP x, const DP h, Vec_O_DP &yout , void derivs_gaba6(const DP, Vec_I_DP &, Vec_O_DP &, double), double gaba_conc){

	int i;
	DP xh, hh, h6;
	int n = y.size();
	Vec_DP dym(n), dyt(n), yt(n);
	hh = h*0.5;
	h6 = h/6.0;
	xh=x+hh;
	for (i=0;i<n;i++) yt[i] = y[i]+hh*dydx[i];  // first step
	derivs_gaba6(xh,yt,dyt, gaba_conc); // second step
	for (i=0;i<n;i++) yt[i]=y[i]+hh*dyt[i]; 
	derivs_gaba6(xh,yt,dym, gaba_conc);  // third step
	for (i=0;i<n;i++) {
		yt[i]=y[i]+h*dym[i];
		dym[i] += dyt[i];
	}
	derivs_gaba6(x+h,yt,dyt, gaba_conc); // fourth step
	for (i=0;i<n;i++) yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);

}
