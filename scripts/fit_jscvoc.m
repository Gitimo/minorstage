function [n,J0,r]=fit_jscvoc(Jm,V)
	global k T q;
	[a,r]=fit_lin(log(Jm),V);
	n=a(1)*q/(k*T);
	J0=exp(-a(2)*q/(n*k*T));
end
