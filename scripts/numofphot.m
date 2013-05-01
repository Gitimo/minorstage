function [out]=numofphot(wavel,energyd)
	global h; global c;
	out=energyd.*wavel/(h*c);
end
