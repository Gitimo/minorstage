function [J01,J02,Jret,r2]=fit_dd(Jm,V,J0fromsd)
	[Jret,par,cvg,iter,corp,covp,covr,stdresid,Z,r2]=(
	leasqr(V,Jm,[J0fromsd*0.001,J0fromsd*0.999],"diode_double"));
	J01=par(1);
	J02=par(2);
endfunction
