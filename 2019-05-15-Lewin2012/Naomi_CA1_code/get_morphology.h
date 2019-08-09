 /********************************************************************
 *
 * functions to read in length, radius, and nearest-neighbor list
 * for CA1 neuron n123
 *
 *******************************************************************/ 
#include "neuron_structures.h"

void read_dimensions(Compartment *);
void read_neighborlist(Compartment *);
void calc_distance(Compartment *);
void reset_grid(Compartment *, Compartment *, int *);

