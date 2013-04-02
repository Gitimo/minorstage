function [E]=wavetoen(lam)
	k=1.3806488e-23;
	h=6.62606957e-34;
	c=299792458;
	E=h*c./lam;
end
