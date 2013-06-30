function [Voc]=analvoc(Jsc,J01,J02)
	k=1.3806488e-23;
	q=1.60217646e-19;
	T=293;

	a=(J01/J02)^2;
	b=-2*Jsc*J01/J02^2;
	c=(Jsc/J02)^2-b;
	x=(-b+sqrt(b^2-4*a*c))/(2*a);
	if x<=0
		x=(-b-sqrt(b^2-4*a*c))/(2*a);
	end

	Voc=k*T/q*log(x);
end
