%Refractive indices
na=1;	 		%air
nag=0.1630+i*5.95; 	%silver at 886nm
nau=0.21015+i*5.8835;	%gold at 886 nm

n1=3.57189;		%AlGaAs (9.9% AL) at 886nm
				%n1=3.7203+i*2.6089e-1;	%InGaP at 886 nm
n2=3.61354+i*1.7188e-3; %GaAs at 886nm 
la0 = 886;		%target wavelength (nm)
la=linspace(350,1000,301);

%material thickness in nm
L1=[1000];
L2=[10];
%optical thickness
L=[L1/n1,L2/n2];
L=L/la0;
clear L;
L=[];
Gag=abs(multidiel([na,nag],L,la/la0)).^2;
Gau=abs(multidiel([na,nau],L,la/la0)).^2;
plot(la,Gag,la,Gau,'rx')
xlabel "Wavelength (nm)"
ylabel "Reflectance"
legend("Silver","Gold","location","sout")
