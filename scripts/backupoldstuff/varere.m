clear all; clf; format long;

# Define some physical constants
global k=1.3806488e-23;
global q=1.60217646e-19;
global h=6.62606957e-34;
global c=299792458;
global Troom=293;
global Tsun=6000;

# Black-body energy density emission spectrum as function of wavelength in m
function [y] = bb(lam,Fs,T)
global k q h c;
	y = 8*Fs*h*c./(lam.^5) ./(exp(h*c./(lam*k *T)-1) ); # not sure about the 8
end
# Fuction to compute the number of photons
function [y] = numofphot(lam,Fs,T)
global k q h c;
	y = bb(lam,Fs,T).*lam/(h*c);
end
#Utility-functions
function [y] = wavetoen(lambda)
global k q h c;
	y = h*c ./ lambda ;
end

function [y] =entowave(E);
global k q h c;
	y = h*c ./ E ;
end


nofcells=input("How many cells?");
Voc(1:nofcells)=0;
names(1,:)=input("Path to first data-file?","s");
Voc(1)=input("First cells Voc?");
fcell=dlmread(names(1,:));
npoints=size(fcell)(1);
wave=fcell(:,1);
eqearray=zeros(npoints,nofcells);
eqearray(:,1)=fcell(:,2);
for m=2:nofcells
	names(m,:)=input("Path to next file?","s");
	data=dlmread(names(m,:));
	eqearray(:,m)=data(:,2);
	Voc(m)=input("This cells Voc?");
end
plot(wave,eqearray)
xlabel "Wavelength (nm)"
ylabel "External Quantum Efficiency"
legend(names,"location","outside")
print -dpng eqes.png
print -deps eqes.eps


longwave=(10:10:920000)*10^(-9); #wavelength-vector in m	
longeqearray=zeros(size(longwave)(2),nofcells);
for n=1:nofcells
	for m=1:npoints
		longeqearray(wave(m)/10,n)=eqearray(m,n);	#first row of EQE: lambda =350 assigned to 35th
	end
end
