/*
 *  interneuron_update.h
 *  
 *
 *  Created by naomi on 12/9/09.
 *  Copyright 2009 Weill Medical College of Cornell University. All rights reserved.
 *
 */

void bc_Ina(Compartment *, double , double);
void bc_Ikdr(Compartment *, double, double);
void bc_Ikm(Compartment *, double, double);
void bc_Ica(Compartment *, double, double);
void bc_Iahp(Compartment *, double, double);
void bc_syn(Compartment *, double, double);
void bc_ahp(Compartment *, double, double);
void basket_update(Compartment * , double, double , double);