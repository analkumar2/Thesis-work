
void calc_leak(Compartment *, double, double);
void calc_leak2(Compartment *, double, double);
void calc_leak3(Compartment *, double, double); //no calcium
void set_leaks(Compartment *, double, double );
void set_kclleaks(Compartment *, double, double);
void set_naleaks(Compartment *, double, double);
void set_clleak(compartment *, double, double);
void update_ions(Compartment *, double, double );
void long_diffusion(Compartment *, int , double );
void update_currents(Compartment *, int, double, double, double, double, double, double, double, double);
void update_noions(Compartment *, int, double, double, double, double, double); // no ion flux
void update_noca(Compartment *, int, double, double, double, double, double); // no calcium currents
void set_resistance(Compartment * , int);
void cable_eqn_fe(Compartment *, int, double);
void cable_eqn_rk(Compartment * , int, double, int);
void cable_eqn_ie(Compartment *, double, int);
void print_data(Datafile, Compartment *, Compartment *,  double, int);

