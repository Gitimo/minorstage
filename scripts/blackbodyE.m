function [y] = blackbodyE(dE,T)
	# Black-body photon density emission spectrum as function of energy in J or power in W ...
	global h; global c; global k;
	y = 2*pi*dE.^2/(h^3*c^2).*1./(exp(dE/(k *T))-1) ; 
	# y = 2*pi*dE.^3/ (h^3*c^2).*1./(exp(dE/(k *T))-1) ;  from own algebra
end
