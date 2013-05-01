function [y] = blackbody(lam,Fs,T)
	# Black-body energy density emission spectrum as function of wavelength in m
	global h; global c; global k;
	y = 2*Fs*(h*c^2)./(lam.^5).*1./(exp(h*c./(lam*k *T))-1 );
end
