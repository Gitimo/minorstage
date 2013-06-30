function [y]=diode_double(x,par)
	global k T q;	
  	 y=par(1)*(exp(q*x/(k*T))-1) + par(2)*(exp(q*x/(2*k*T))-1);
end
