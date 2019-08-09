/*
 *  poirazi_channels.h
 *  
 *
 *  Created by Naomi Lewin on 11/14/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */
#include "neuron_structures.h"

double ghk(double conc_in, double conc_out, double vm, double z);
void na(Compartment *neuron, double t, double dt, int i); 
void  nap(double *m, double distance, double dt,  int i); 
void kdr(Compartment *, double, double, int);
void ka(Compartment *, double, double, int);
void km(Compartment *, double, double, int );
void kdr_jun(Compartment *, double, double);
void sahp(Compartment *, double, double, int );
void sahp_markaki(Compartment *, double, double, int);
void mahp(Compartment *, double, double, int );
void cat(Compartment * , double , double, int );
void car(Compartment * , double, double, int );
void cal(Compartment * ,  double, double, int);
void cal_markaki(Compartment *, double, double, int);
void update_ca_conc(Compartment *, double, double);
void update_ca_kager(Compartment *, double);
void ih(Compartment * , double, double, int);
void ih2(Compartment *, double, double);
void na6state(Compartment * ,double, double, int);
