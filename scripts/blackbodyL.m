function [y] = blackbodyL(dL,T)
	#spectral photon density as function of wavelength ... checked, gives same results as blackbodyE
	global h; global c; global k;
	y = 2*pi ./(h*dL.^2) .*1./ (exp(h*c./(dL*k*T))-1);
end
