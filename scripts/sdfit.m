function [nonlinN,nonlinJ0]=sdfit(J,V,J0guess,nguess)
	[fA,par,cvg,iter,corp,covp,covr,stdresid,Z,r2A]=leasqr(V,J,[J0guess,nguess],'single_diode');
	disp(r2A)
	nonlinN=par(2);
	nonlinJ0=par(1);
end
