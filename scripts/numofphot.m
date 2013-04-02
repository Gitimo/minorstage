function [out]=numofphot(wavel,energyd)
	h=6.62606957e-34;
	c=299792458;
	out=energyd.*wavel/(h*c);
end
