function y = numvoc (Jsc,J01,J02,Rsh,x)
	global k T q;
	if Rsh(1,1)== 0
		y=((Jsc)
		  -(J01)*(exp(q*x/(  k*T))-1)
		  -(J02)*(exp(q*x/(2*k*T))-1));
	else
		y=((Jsc)
		  -(J01)*(exp(q*x/(  k*T))-1)
		  -(J02)*(exp(q*x/(2*k*T))-1)
		  - x/Rsh(1,2));
	endif
endfunction
