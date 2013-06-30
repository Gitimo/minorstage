function y = ddresfit (Jm,V,x)
# ddresfit is a function used to convert a measured I-V curve into a set off
# nonlinear equations according to the double-diode model. 
# It will be used in a call-by-handle to fsolve to minimise this set.
# Parameters are physically constrained to positive values; however,
# fsolve does not support constraints. Solution: square parameters!
# Usage example: 
# J=measuredJ; V=measuredV; xguess=sqrt([J01,J02,Rs,Rsh,Jsc]);
# [par,X,info]=fsolve(@(x) ddresfit(J,V,x),x,optimset(...))

	k=1.3806488e-23;
	q=1.60217646e-19;
	T=293;
	y=zeros(size(V));
	for m=1:max(size(V)) 
		y(m)=((x(5)^2)
                     -(x(1)^2)*(exp(q*(V(m)-Jm(m)*(x(3)^2))/(  k*T))-1) 
                     -(x(2)^2)*(exp(q*(V(m)-Jm(m)*(x(3)^2))/(2*k*T))-1) 
                     -(V(m)+Jm(m)*(x(3)^2))/(x(4)^2) 
                     -Jm(m) );
	endfor
endfunction
