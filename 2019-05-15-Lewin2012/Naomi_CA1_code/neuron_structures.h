/*********************************************
 *
 * structures for CA1 model
 *
 *********************************************/

#ifndef STRUCTURES_H
#define STRUCTURES_H

#include "nr.h"

const double F = 96487.;
const double R = 8.314;
const double celsius = 34.;
const double T = celsius + 273.14;
const double PI = 3.14159;
const double rtf = R*T/F;


using namespace std;

struct list_of_files {
  /* sodium channel data files */
  FILE * vmptr;
  FILE * inaptr;
  FILE * mnaptr;
  FILE * hnaptr;
  FILE * snaptr;
	FILE * na6_0;
	FILE * na6_1;
	FILE * na6_3;
	FILE * na6_4;
	FILE * na6_5;

  FILE * itotalptr;

  /* kdr data files */
  FILE * ikdrptr;
  FILE *mkdrptr;

  /* ka data files */
  FILE * ikaptr;

  /* km data files */
  FILE * ikmptr;

  /* ih data files */
  FILE * ihptr;
  FILE * mihptr;

  /* Iahp data files */
  FILE * sahpptr;
  FILE * mahpptr;

  /* calcium data files */
  FILE * icatptr;
	FILE * icarptr;
	FILE * icalptr;
  FILE * ca_in_ptr;
  FILE * capump;
  FILE * nacax;
  
  /* synaptic data files */
  FILE * glutptr;
  FILE * nmdaptr;
  FILE * ampaptr;
  FILE * gaba;
  FILE * gaba_cl;
  FILE * gaba_open;
  FILE * egaba;
  FILE * gaba_conc;
  FILE * ecl;
  
  FILE * naiptr; 
  FILE * naoptr;
  FILE * kiptr;
  FILE * koptr;
  FILE * cliptr;
  FILE * cloptr;
	
	FILE * dk_kcc2;  /* change in ko if kcc2 is only source of outward potassium */
  
  
  /* transporter data files */
  FILE * kcc2ptr;
  FILE * nkcc1ptr;
  FILE * nakpump;
  
  /* leak currents */
  FILE * clleak;
	
	/* total currents */
	FILE * iktotal;
	FILE * inatotal;
	FILE * icltotal;
	FILE * icatotal;
	FILE * inaleak;
	FILE * ikleak;
	
	/* interneuron (basket cell) files */
	FILE * bc_vm;
	FILE * bc_Isyn;
	FILE * bc_ina;
	FILE * bc_ikdr;
	FILE * bc_ahp;
	FILE * bc_cai;


};

typedef struct list_of_files Datafile;

struct compartment {
  double t;
  double dt;
  double Vm;
  double Vm_new;
  double dV;
  double distance; // distance from soma
  double dx; // length
  double radius;
  double vol; // compartment volume
  int num_neighbors; // number of nearest neighbors, either 1 (boundary), 2 (unbranched), or 3 (branch point)
  int neighbor_list[3];  // lists the compartment number (identifier) of the nearest neighobors

  /* non uniform membrane resistance, based on poirazi ca1 model */
  double cm;
  double rm;
  double ra; // axial resistance for each compartment
  double eleak; // reversal potential for nonspecific leak current
  double Ikleak, Inaleak, Iclleak;
  double gkleak, gnaleak, gclleak;
  
  double erev_na, erev_k, erev_cl, erev_hco3, erev_ca;
  
  double istim;
	double nspike;

  /* active current gates */
  double Ina, ina_m, ina_h, ina_s;
  double Inap, inap_m; //persistent sodium, only in distal apical tuft
  double Ih, Ih_na, Ik_h, ih_m, gh;
  double Ikdr, ikdr_m, ikdr_n, ikdr_l;
  double Ika, ika_m, ika_h;
  double Ikm, ikm_m;
  double Isahp, isahp_m;
  double Imahp, imahp_m;
  double Icat, icat_m, icat_h;
  double Icar, icar_m, icar_h;
  double Ical, ical_m, ical_h;

	/* sodium markov model */
	double na6_y[6];
	double na5_y[5];
	
  /* synaptic currents */
  double gaba_conc;  // gaba concentration in cleft for HFS
	double gaba_conc_syn; //gaba conc for interneuron feedback
	double gaba_concE;  // for exrasyn receptors
  double Igabab, igabab_r, igabab_d, igabab_g;  // gaba(b) states
  double gaba6_y[6];
	double gaba6_yE[6];
	double gaba6_y_syn[6];
  double Icl_gaba, Ihco3_gaba;
  double Icl_gaba_syn, Ihco3_gaba_syn;
  double gmax_gaba_a;
  double ggaba;

  double glut_conc;
  double gmax_ampa, gmax_nmda;
  double ampa_r, ampa_r0, ampa_r1;
  double nmda_r, nmda_r0, nmda_r1;
  double Ik_nmda, Ik_ampa;
  double Ina_nmda, Ina_ampa, Ica_nmda;
  double lastrelease;

  double Iion_total; 
	double Iaxial;
  double Ik_total, Ina_total, Icl_total, Ica_total, Ihco3_total; 

  /* ion transport flux */
  double k_kcc2, cl_kcc2;
  double k_nkcc1, na_nkcc1, cl_nkcc1;
  double k_atpase, na_atpase, imax_atpase;
  double y_nkcc, p_nkcc;
  double imax_nacax, Ina_nacax, Ica_nacax;

  /* ion concentrations, intra and extra-cellular */
  double na_i, na_o;
  double k_i, k_o, k_o_dummy;
  double cl_i, cl_o;
	double ca_i, ca_o, cai_Ltype;
  double hco3_i, hco3_o;
  double kbuffer; // extracellular potassium buffer concentration
  
  /* longitudinal diffusion concentration changes */
  double na_i_diff, na_o_diff;
  double k_i_diff, k_o_diff;
  double cl_i_diff, cl_o_diff;
  double ca_i_diff, ca_o_diff;
  
  /* calcium related parameters */
  double buffer;
  double fe;
  double icapump, imax_capump;
  
};

typedef struct compartment Compartment;

#endif
