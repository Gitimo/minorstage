function y = ddresfitlsq (Jm,V,pin)

	Y=zeros(size(Jm));
	=leasqr(V,Y,pin,@(x) double_diode_res(V,Jm,x) )	

endfunction
