TITLE Ca R-type channel with medium threshold for activation
: used in distal dendritic regions, together with calH.mod, to help
: the generation of Ca++ spikes in these regions
: uses channel conductance (not permeability)
: written by Yiota Poirazi on 11/13/00 poirazi@LNC.usc.edu

NEURON {
	  SUFFIX car
	  USEION ca READ eca WRITE ica
    RANGE gcabar, m, h, g, gmax
	  RANGE minf, mtau
}

UNITS {
	  (mA) = (milliamp)
	  (mV) = (millivolt)
}

PARAMETER {              : parameters that can be entered when function is called in cell-setup
    v             (mV)
    celsius = 34	(degC)
    gcabar = 0    (mho/cm2) : initialized conductance
	  eca = 140     (mV)      : Ca++ reversal potential
}  

STATE {	m h }            : unknown activation and inactivation parameters to be solved in the DEs  

ASSIGNED {               : parameters needed to solve DE
	  ica    (mA/cm2)
    minf[2]
	  mtau[2] (ms)
    g      (mho/cm2)
    gmax   (mho/cm2)
}

BREAKPOINT {
	  SOLVE states METHOD cnexp
    g = gcabar*m*m*m*h
	  ica = g*(v - eca)
    if (g > gmax) {
        gmax = g
    }
}

INITIAL {
    mhn(v)
    m = minf[0]
    h = minf[1]
    g = gcabar*m*m*m*h
    ica = g*(v - eca) : initial Ca++ current value
    gmax = g
}

DERIVATIVE states {
	  mhn(v)
	  m' =  (minf[0] - m)/mtau[0]
	  h' =  (minf[1] - h)/mtau[1]
}	

FUNCTION varss(v (mV), i) {
	  if (i==0) {
	      varss = 1 / (1 + exp((v+48.5(mV))/(-3(mV)))) : Ca activation
	  }
	  else if (i==1) {
        varss = 1/ (1 + exp((v+53(mV))/(1(mV))))    : Ca inactivation
	  }
}

FUNCTION varmtau(v (mV), i) (ms) {
	  if (i==0) {
        varmtau = 50  : activation variable time constant
    }
	  else if (i==1) {
        varmtau = 5   : inactivation variable time constant
    }
	  
}	

PROCEDURE mhn(v (mV)) {LOCAL a, b :rest = -70
    TABLE minf, mtau DEPEND celsius FROM -100 TO 100 WITH 200
  	FROM i=0 TO 1 {
	      mtau[i] = varmtau(v,i)
		    minf[i] = varss(v,i)
	  }
}















