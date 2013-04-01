function [y] = blackbody(lam,Fs,T)
	# Black-body energy density emission spectrum as function of wavelength in m
	k=1.3806488e-23;
	h=6.62606957e-34;
	c=299792458;
	y = 8*Fs*h*c./(lam.^5) ./(exp(h*c./(lam*k *T)-1) ); # not sure about the 8
end