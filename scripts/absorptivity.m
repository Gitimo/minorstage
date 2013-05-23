function [y]=absorptivity(E,Eg,L)
	format long;
	q=1.60217657e-19;	#elementary charge
	a_0=800000;		#absorption-coefficient for GaAs 8000/cm -> 8000 00/m
	E_0=q*6.7e-3;		#Urbach energy: 6.7 meV
	E_d=q*140e-3;		#
	
	g(1:max(size(E)))=0;

	for m=1:max(size(E))
		if E(m)<=Eg
			g(m)=1-exp(-2*a_0*exp((E(m)-Eg)/E_0)*L);
		else
			g(m)=1-exp(-2*a_0*(1+(E(m)-Eg)/E_d)*L);
		endif
	endfor
	y=g;
end
