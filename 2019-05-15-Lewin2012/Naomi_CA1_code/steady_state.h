/********************************************************
 * function to test steady state curves
 ******************************************************/

void steady(Compartment *);

void steady(Compartment * neuron){
  
  FILE * na_ss;
  FILE * na_m_ss;
  FILE * na_h_ss;
  FILE * na_s_ss;
  FILE * kdr_ss;
  FILE * kdr_m_ss;
  FILE * ka_ss;

  na_ss = fopen("na_ss", "w");
  kdr_ss = fopen("kdr_ss", "w");
  ka_ss = fopen("ka_ss", "w");

  double t = 0.0;
  double dt = 1.;

  for(double v = -80.; v < 40; v++){
    for(int i = 0; i<1; i++){
      neuron[i].Vm = v;
      na(neuron, t, dt);
      kdr(neuron, t, dt);
      ka(neuron, t, dt);
      if(i == 0){
	fprintf(na_ss, "%.1f\t%.12f\t%.12f\t%.12f\t%.12f\n", neuron[i].Vm, neuron[i].Ina, neuron[i].ina_m, neuron[i].ina_h, neuron[i].ina_s);
	fprintf(kdr_ss, "%.1f\t%.12f\t%.12f\n", neuron[i].Vm, neuron[i].Ikdr, neuron[i].ikdr_m);
	fprintf(ka_ss, "%.1f\t%.12f\t%.12f\t%.12f\n", neuron[i].Vm, neuron[i].Ika, neuron[i].ika_m, neuron[i].ika_h);
      } // end printing loop
    } // end compartment loop
  } // end voltage loop
  
  fclose(na_ss);
  fclose(kdr_ss);
  fclose(ka_ss);

}
