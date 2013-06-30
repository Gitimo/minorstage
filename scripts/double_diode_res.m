function y = double_diode_res (V,Jm,x)
	k=1.3806488e-23;
	q=1.60217646e-19;
	T=293;
	y=((x(5)^2)
         -(x(1)^2)*(exp(q*(V-Jm*(x(3)^2))/(  k*T))-1) 
         -(x(2)^2)*(exp(q*(V-Jm*(x(3)^2))/(2*k*T))-1) 
         -(V+Jm*(x(3)^2))/(x(4)^2) 
         -Jm);
endfunction
