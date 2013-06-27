function [n,J0,r]=jscvoc(J,V)
	T=298;
	k=1.3806488e-23;
	q=1.60217646e-19;

	[a,r]=linfit(log(J),V);
	n=a(1)*q/(k*T);
	J0=exp(-a(2)*q/(n*k*T));
end
