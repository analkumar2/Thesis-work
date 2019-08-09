#include "neuron_structures.h"


void set_conductances(Compartment *, int);
void nmda(Compartment *, double, double, double);
double exptable(double );
void gaba_6state(Compartment *, double, double, double, double, int);
void gaba_synapse(Compartment *, double, double, double, int);
void stylized_gaba(Compartment *, double, double );
void imp_euler(Mat_O_DP &, Vec_IO_DP &, const DP);