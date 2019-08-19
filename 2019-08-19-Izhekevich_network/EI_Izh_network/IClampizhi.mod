COMMENT
Since this is an electrode current, positive values of i depolarize the cell
and in the presence of the extracellular mechanism there will be a change
in vext since i is not a transmembrane current but a current injected
directly to the inside of the cell.
ENDCOMMENT

NEURON {
	POINT_PROCESS IClampizhi
	RANGE del, dur, amp
	USEION izh WRITE iizh
}
UNITS {
	(nA) = (nanoamp)
}

PARAMETER {
	del (ms)
	dur (ms)	<0,1e9>
	amp (nA)
}

STATE {
	iizh (nA)
}

ASSIGNED {
}

INITIAL {
 iizh = 0
}

BREAKPOINT {
	at_time(del)
	at_time(del+dur)

	if (t < del + dur && t >= del) {
		iizh = amp
	}else{
		iizh = 0
	}
}
