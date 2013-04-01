clear all; clf;

# Define some physical constants
global k=1.3806488e-23;
global q=1.60217646e-19;
global T=293;

# fitting to find idealities and currents

function [y]=n1fit(x,par)
global k;
	q=1.60217646e-19;
	T=293;
	y=q*x/(k*T) +par;
end
function [y]=n2fit(x,par)
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

# Perform the fits
j=5e-02;
t1= 2.6032e+02;
t2= 4.2686e+02;
pin=[j,t1,t2];
wght=ones(size(it1(1:2000)));
[fir1,pir1,cvg,iter,corp,covp]=leasqr(t(1:2000),ir1(1:2000),pin,"double_exp",.0001,20,wght);
[fir2,pir2,cvg,iter,corp,covp]=leasqr(t(1:2000),ir2(1:2000),pin,"double_exp",.0001,20,wght);
[fit1,pit1,cvg,iter,corp,covp]=leasqr(t(1:2000),it1(1:2000),pin,"double_exp",.0001,20,wght);
fir1=double_exp(t,pir1);
fir2=double_exp(t,pir2);
fit1=double_exp(t,pit1);
illu=ir2/(equivref);

# Calculate reference cell's I_sc - is that correct?
it2=it1.*ir2./ir1;

# Some more fitting
[nlinJf,pnlinf,cvg,iter,corp,covp]=leasqr(vt(100:8000),it2(100:8000),[1e-21,1e-11],"nonlincurrentfitfix",.0001,20,ones(size(it2(100:8000))));#*0.001./it2);
J01=pnlinf(1);
J02=pnlinf(2);
V=linspace(0.6,1.4,2000);
J=J01*(exp(q*V/(k*T))-1)+J02*(exp(q*V/(2*k*T))-1);
wght=ones(size(V(1:500)));
#[fn1,pn1,cvg,iter,corp,covp]=leasqr(V(1501:2000),log(J(1501:2000)),1,"n1fit",.0001,20,wght);
#[fn2,pn2,cvg,iter,corp,covp]=leasqr(V(1:500),log(J(1:500)),1,"n2fit",.0001,20,wght);

plot(V,log(J))
axis([min(vt(30:8000)),max(vt),min(log(it2)),max(log(it2))])
title "base e"
hold on
plot(V,q*V/(k*T)+pn1)
plot(V,q*V/(2*k*T)+pn2)
#plot(vt(100:8000),log(it2(100:8000)))
input("wait");
clf

# Sum of squares-function
function y=sos(a,b)
	y=sum(abs((a-b)));
end

# Try simpler model
function [y]=singlefree(x,par)
	k=1.3806488e-23;
	q=1.60217646e-19;
	T=293;
	y=par(1)*(exp(q*x/(par(2)*k*T))-1);
end

[sf,psf,cvg,iter,corp,covp]=leasqr(vt(1:4000),it2(1:4000),[1e-18,1],"singlefree",.0001,20,ones(size(it2)));#*0.001./it2);
psf
Jsf=psf(1)*(exp(q*V/(psf(2)*k*T))-1);


#simpler still
p1=polyfit(vt(100:7000),log(it2(100:7000)),1);
n=p1(1)*k*T/q
J0=exp(p1(2))
#plot(V,log(Jsf))
plot(vt(100:8000),log(it2(100:8000)))
axis([min(vt(30:8000)),max(vt),min(log(it2)),max(log(it2))])
title "base e"
hold on

plot(V,V*p1(1)+p1(2))




