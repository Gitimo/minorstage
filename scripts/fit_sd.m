function [nonlinN,nonlinJ0,Jret,r2]=fit_sd(Jm,V,J0guess,nguess)
	[Jret,par,cvg,iter,corp,covp,covr,stdresid,Z,r2]=(
	leasqr(V,Jm,[J0guess,nguess],'diode_single'));
	nonlinN=par(2);
	nonlinJ0=par(1);
end
