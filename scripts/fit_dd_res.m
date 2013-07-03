function [y] = fit_dd_res (par,Jm,V)
# ddresfit is a function used to convert a measured I-V curve into a set off
# nonlinear equations according to the double-diode model. 
# It will be used in a call-by-handle to fsolve to minimise this set.
# Parameters are physically constrained to positive values; however,
# fsolve does not support constraints. Solution: square parameters!
# Usage example: 
# J=measuredJ; V=measuredV; parguess=sqrt([J01,J02,Rs,Rsh,Jsc]);
# [par,X,info]=fsolve(@(x) ddresfit(J,V,x),x,optimset(...))
	y=zeros(size(V));
	y=diode_double_res(Jm,V,par);
endfunction
