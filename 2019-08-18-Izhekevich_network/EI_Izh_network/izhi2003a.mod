: Izhikevich artificial neuron model from
: EM Izhikevich "Simple Model of Spiking Neurons"
: IEEE Transactions On Neural Networks, Vol. 14, No. 6, November 2003 pp 1569-1572
: V is the voltage analog, u controls
: see COMMENT below or izh.hoc for typical parameter values
: uncomment lines with dvv,du to graph derivatives

NEURON {
  POINT_PROCESS Izhi2003a
  RANGE a,b,c,d,f,g,thresh
  POINTER Ic, Iampar, Igabar
}

INITIAL {
  V=-65
  u=0.2*V
  net_send(0,1)
}

PARAMETER {
  a = 0.02
  b = 0.2
  c = -65
  d = 2
  f = 5
  g = 140
  Ic = 10
  Iampar = 10
  Igabar = 10
  taug = 1
  thresh=30
}

STATE { u V } : use V for voltage so don't interfere with built-in v of cell

ASSIGNED {
}

BREAKPOINT {
  SOLVE states METHOD derivimplicit
}

DERIVATIVE states {
  V' = 0.04*V*V + f*V + g - u + Ic + Iampar + Igabar
  u' = a*(b*V-u)
}

NET_RECEIVE (w) {
  if (flag == 1) {
    WATCH (V>thresh) 2
  } else if (flag == 2) {
    net_event(t)
    V = c
    u = u+d
  }
}
