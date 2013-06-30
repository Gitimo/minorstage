function [y]=single_diode(x,par)
	k=1.3806488e-23;
	q=1.60217646e-19;
	T=293;
	y=par(1)*(exp(q*x/(par(2)*k*T))-1);
endfunction
