clear all; clf;

# Define some physical constants
k=1.3806488e-23;
q=1.60217646e-19;
T=293;

# fitting to find idealities and currents

funcion [y]=n1fit(x,par)
	k=1.3806488e-23;
	q=1.60217646e-19;
	T=293;
	y=q*x/(k*T) +par;
end
funcion [y]=n2fit(x,par)
	k=1.3806488e-23;
	q=1.60217646e-19;
	T=293;
	y=q*x/(2*k*T) +par;
end

function [y]=nonlincurrentfitfix(x,par)
	k=1.3806488e-23;
	q=1.60217646e-19;
	T=293;
	y=par(1)*(exp(q*x/(k*T))-1)+par(2)*(exp(q*x/(2*k*T))-1);
end

function [y]=nonlincurrentfit(x,par)
	k=1.3806488e-23;
	q=1.60217646e-19;
	T=293;
	y=par(1)*(exp(q*x/(par(2)*k*T))-1)+par(3)*(exp(q*x/(par(4)*k*T))-1);
end

# Function that will be fit initially to follow flash-curve
function [y]=double_exp(x,par)
  y=par(1)*exp(-par(2)*x)-par(1)*exp(-par(3)*x);
end

# Function to save
function savefig()
	s=input("Save the plot?[y/n]","s");
	if (s=="y")
		clear plotname;
		plotname=input("Filename for plot?","s");
		print(gcf,plotname,'-deps');
	endif
	clf;
end
# Sum of squares-function
function y=sos(a,b)
	y=sum(abs((a-b)));
end

#Load data and assign vectors
fhalf="../data/" ; shalf=input("Filename?","s") ; fname=strcat(fhalf,shalf);
m=load(fname);
t=m(:,1) ; ir1=m(:,2) ; ir2=m(:,3) ; it1=m(:,4) ; vt=m(:,5) ;

# Calculate illumination in "suns" and convert voltages to currents-densities (ref cell area:1cm^2)
tsurf=input("Test-cel surface area in cm^2?");
Rref=0.04248 ; equivref=13.7e-3 ;
Rtest=0.0631 ; 
ir1=ir1/Rref; ir2=ir2/Rref; 
it1=it1/(Rtest*tsurf);

# Generic weights matrix
for m=1:8000
	if (ir1(m)!=0)
	wir1(m)=0.001/ir1(m);
	else
	wir1(m)=0.001/(0.5*(ir1(m-1)+ir1(m+1)));
	end
	if (ir2(m)!=0)
	wir2(m)=0.001/ir2(m);
	else
	wir2(m)=0.001/(0.5*(ir2(m-1)+ir2(m+1)));
	end
	if (it1(m)!=0)
	wit1(m)=0.001/it1(m);
	else
	wit1(m)=0.001/(0.5*(it1(m-1)+it1(m+1)));
	end
end

 
# Perform the fits
j=5e-02;
t1= 2.6032e+02;
t2= 4.2686e+02;
pin=[j,t1,t2];
wght=ones(size(it1));
[fir1,pir1,cvg,iter,corp,covp]=leasqr(t,ir1,pin,"double_exp",.0001,20,wght);
[fir2,pir2,cvg,iter,corp,covp]=leasqr(t,ir2,pin,"double_exp",.0001,20,wght);
[fit1,pit1,cvg,iter,corp,covp]=leasqr(t,it1,pin,"double_exp",.0001,20,wght);
illu=fir2/(equivref);

# Calculate reference cell's I_sc - is that correct?
it2=fit1.*fir2./fir1;

# Some more fitting
[nlinJf,pnlinf,cvg,iter,corp,covp]=leasqr(vt,it2,[1e-21,1e-11],"nonlincurrentfitfix",.0001,20,ones(size(it2))*0.001./it2);
pnlinf
[nlinJ,pnlin,cvg,iter,corp,covp]=leasqr(vt,it2,[1e-11,1,1e-21,2],"nonlincurrentfit",.0001,20,ones(size(it2))*0.001./it2);
pnlin

sos_nlin_fix=sos(nlinJf,it2)
sos_nlin=sos(nlinJ,it2)
