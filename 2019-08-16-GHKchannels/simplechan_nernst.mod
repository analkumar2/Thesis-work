TITLE simp channel
: Simple calcium channel
: to test ghk

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)

	FARADAY = 96520 (coul)
	R = 8.3134 (joule/degC)
	KTOMV = .0853 (mV/degC)
}

PARAMETER {
	v (mV)
	celsius = 25	(degC)
	gcatbar=.003 (mho/cm2)
	cai = 50.e-6 (mM)
	cao = 2 (mM)
	q10 = 5
	mmin=0.2
	hmin=10
	a0h =0.015
	zetah = 3.5
	vhalfh = -75
	gmh=0.6
	a0m =0.04
	zetam = 2
	vhalfm = -28
	gmm=0.1
	eca =140 (mV)
}


NEURON {
	SUFFIX simpnernst
	USEION ca READ cai,cao WRITE ica
        RANGE gcatbar, ica, gcat
        GLOBAL minf,mtau
}

STATE {
	m
}

ASSIGNED {
	ica (mA/cm2)
        gcat (mho/cm2)
	minf
	mtau
}

INITIAL {
	rates(v)
	m = minf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gcat = gcatbar*m
	ica = gcat*(v-eca)

}

DERIVATIVE states {	: exact when v held constant
	rates(v)
	m' = (minf - m)/mtau
}

FUNCTION alpmt(v(mV)) {
  alpmt = exp(0.0378*zetam*(v-vhalfm))
}

FUNCTION betmt(v(mV)) {
  betmt = exp(0.0378*zetam*gmm*(v-vhalfm))
}

PROCEDURE rates(v (mV)) { :callable from hoc
	LOCAL a,b, qt
        qt=q10^((celsius-25)/10)

	a = 0.2*(-1.0*v+19.26)/(exp((-1.0*v+19.26)/10.0)-1.0)
	b = 0.009*exp(-v/22.03)
	minf = a/(a+b)
	mtau = betmt(v)/(qt*a0m*(1+alpmt(v)))
	if (mtau<mmin) {mtau=mmin}
}
