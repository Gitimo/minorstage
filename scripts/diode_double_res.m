function [y] = diode_double_res (Jm,V,par)
	global k T q;	
	y=((par(5)^2)
          -(par(1)^2)*(exp(q*(V-Jm*(par(3)^2))/(  k*T))-1) 
          -(par(2)^2)*(exp(q*(V-Jm*(par(3)^2))/(2*k*T))-1) 
          -(V+Jm*(par(3)^2))/(par(4)^2) 
          -Jm);
endfunction
