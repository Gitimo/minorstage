function [y]=diode_single(V,par)
	global k T q;	
	y=par(1)*(exp(q*V/(par(2)*k*T))-1);
endfunction
