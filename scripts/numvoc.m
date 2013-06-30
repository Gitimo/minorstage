function y = numvoc (Jsc,J01,J02,x)
	k=1.3806488e-23;
	q=1.60217646e-19;
	T=293;	

	y=((Jsc)
	  -(J01)*(exp(q*x/(  k*T))-1)
	  -(J02)*(exp(q*x/(2*k*T))-1));

endfunction
