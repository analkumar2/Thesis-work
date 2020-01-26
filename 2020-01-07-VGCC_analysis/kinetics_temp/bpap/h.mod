TITLE  H-current that uses Na ions

NEURON {
	  SUFFIX h
    RANGE  gbar,vhalf, K, taun, ninf, g, gmax
	  USEION na WRITE ina      
    :	NONSPECIFIC_CURRENT i
}

UNITS {
	  (um) = (micrometer)
	  (mA) = (milliamp)
	  (uA) = (microamp)
	  (mV) = (millivolt)
	  (pmho) = (picomho)
	  (mmho) = (millimho)
}

PARAMETER { : parameters that can be entered when function is called in cell-setup
	  v              (mV)
    eh     = -10   (mV)
	  K      = 8.5   (mV)
	  gbar   = 0     (mho/cm2)          : initialize conductance to zero
	  vhalf  = -90   (mV)                 : half potential
}	


STATE {               : the unknown parameters to be solved in the DEs
	  n
}

ASSIGNED {                             : parameters needed to solve DE
	  ina  (mA/cm2)
	  ninf
	  taun (ms)
	  g    (mho/cm2)
    gmax (mho/cm2)
}

INITIAL {          : initialize the following parameter using states()
	  states(v)	
	  n = ninf
	  g = gbar*n
	  ina = g*(v-eh)
    gmax = g
}


BREAKPOINT {
	  SOLVE h METHOD cnexp
	  g = gbar*n
	  ina = g*(v-eh)  
    if (g > gmax) {
        gmax = g
    }
}

DERIVATIVE h {
	  states(v)
    n' = (ninf - n)/taun
}

PROCEDURE states(v(mV)) {  
 	  if (v > -30) {
	      taun = 1
	  } else {
        taun = 2(ms)*(1/(exp((v+145(mV))/-17.5(mV))+exp((v+16.8(mV))/16.5(mV))) + 5) :h activation tau
	  }  
    ninf = 1 - (1 / (1 + exp((vhalf - v)/K))) :steady state value
}



