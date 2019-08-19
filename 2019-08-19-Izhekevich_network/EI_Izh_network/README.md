Makes an Izhikevich model (2003a version), inserts GABAR and AMPAR (with modified Exp2Syn mod files),
one input which connects to AMPAR, and one which connects to GABAR.
As the in-built Exp2Syn and IClamp mechanisms use the cell's i and v, I have used USEION method to communicate
the izhikevich current to other mechanisms.
The difference between the POINTER method and USEION method is that in this method, amp, and weights have to be setup way too high. Don't know why.
This is a completed project. I just wanted to see how this can be done.  This may not be the most efficient method or the right way to do it but it works.