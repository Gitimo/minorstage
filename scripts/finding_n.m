clear all; clf;

# Define some physical constants
global k=1.3806488e-23;
global q=1.60217646e-19;
global T=293;

# Function that will be fit to initial data
function [y]=double_exp(x,par)
  y=par(1)*exp(-par(2)*x)-par(1)*exp(-par(3)*x);
end

function [y]=n1fit(x,par)
global q k T;
	y=q*x/(k*T) +par;
end
function [y]=n2fit(x,par)
global q k T;
	y=q*x/(2*k*T) +par;
end

# Function to save
function savefig()
s=input("Save the plot?[y/n]","s");
if (s=="y")
	plotname=input("Filename for plot?","s");
	print(gcf,plotname,'-deps');
endif
clf;
end

# Generate a double exponent
j=5e-02;
t1= 2.6032e+02;
t2= 4.2686e+02;

#Load data and assign vectors
fhalf="../data/" ; shalf=input("Filename?","s") ; fname=strcat(fhalf,shalf);
m=load(fname);
t=m(:,1) ; ir1=m(:,2) ; ir2=m(:,3) ; it1=m(:,4) ; vt=m(:,5) ;
tsurf=input("Test-cel surface area in cm^2?");
Rref=0.04248 ; equivref=13.7e-3 ;
Rtest=0.0631 ; equivref_t=20.7e-3;
ir1=ir1/Rref; ir2=ir2/Rref; 
it1=it1/(Rtest*tsurf);
nfit2=it1.*ir2./ir1;

# Generic weights matrix
weights=ones(size(ir1));
 
# Perform the fits
pin=[j,t1,t2];
[fir1,pir1,cvg,iter,corp,covp]=leasqr(t,ir1,pin,"double_exp",.0001,20,weights);
[fir2,pir2,cvg,iter,corp,covp]=leasqr(t,ir2,pin,"double_exp",.0001,20,weights);
[fit1,pit1,cvg,iter,corp,covp]=leasqr(t,it1,pin,"double_exp",.0001,20,weights);
illu=fir2/(equivref);
fit2=fit1.*fir2./fir1;
logJ=log(fit2);
# Function to fit J-V

function [J]=currentvoltage(V,par)
global q k T;
	J=par(1)*(exp(q*V/(k*T))-1)+par(2)*(exp(q*V/(2*k*T))-1);
end


[J,parJ,cvg,iter,corp,covp]=leasqr(vt(100:5000),nfit2(100:5000),[1e-21,1e-11],"currentvoltage",.0001,20,ones(size(vt(100:5000))));

V=linspace(1,1.2,8000);
J=currentvoltage(V,parJ);
numofpoints=10;
for m=1:numofpoints
	st=(m-1)*(size(vt)(1)/numofpoints)+1;
	sp=m*(size(vt)(1)/numofpoints);
	ldata(m,:)=polyfit(vt(st:sp),log(real(fit2(st:sp))),1);
	lfit(m,:)=polyfit(V(st:sp),log(real(J(st:sp))),1);
#	plot(vt(st:sp),log(real(fit2(st:sp))));
#	input("wait...");
	avgV(m)=sum(V(st:sp))/(size(V)(2)/numofpoints);
	I(m)=sum(illu(st:sp))/(size(illu)(1)/numofpoints);
	Avgerror(m)=(sum(abs(nfit2(st:sp)-fit2(st:sp))./fit2(st:sp)))/(size(illu)(1)/numofpoints);
end

n=1./lfit(:,1)*q/(k*T);
hg = errorbar(avgV,n,Avgerror);
set(hg,"marker","+");
set(hg,"linestyle","none");

#plot(V,n,"+");
axis([0.99*min(V),1.01*max(V),0,1.01*max(n)])
title "n vs. V"
input("wait...");
#plot(I,n,"+");
#axis([0.99*min(I),1.01*max(I),0,1.01*max(n)])
#input("wait...");
#plot(Avgerror,n,"+");
#title "n vs. avg. rel-deviation"
