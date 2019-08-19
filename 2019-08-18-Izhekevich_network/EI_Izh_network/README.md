Makes an Izhikevich model (2003a version), inserts GABAR and AMPAR (with modified Exp2Syn mod files),
one input which connects to AMPAR, and one which connects to GABAR.
As the in-built Exp2Syn and IClamp mechanisms use the cell's i and v, I have used POINTER so that they use the V and I used in Izhikevich mechanism.
This is a completed project. I just wanted to see how this can be done