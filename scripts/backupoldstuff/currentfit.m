clear all; clf;

# Define some physical constants
global k=1.3806488e-23;
global q=1.60217646e-19;
global T=293;
global dVdt;

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
Rtest=0.0631 ; 
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

# Function to fit J-V

function [J]=currentvoltage(V,par)
global q k T dVdt;
	J=par(1)*(exp(q*V/(k*T))-1)+par(2)*(exp(q*V/(2*k*T))-1) + par(3)*dVdt;
end

[J,param_it2,cvg,iter,corp,covp]=leasqr(vt(100:4000),nfit2(100:4000),[1e-21,1e-11],"currentvoltage",.0001,20,ones(size(vt(100:4000))));
[J,param_fit2,cvg,iter,corp,covp]=leasqr(vt(100:4000),fit2(100:4000),[1e-21,1e-11],"currentvoltage",.0001,20,ones(size(vt(100:4000))));
V=linspace(1,1.2,8000);
Jit2=currentvoltage(V,param_it2);
Jfit2=currentvoltage(V,param_fit2);
#Plotting
tcut(1:floor(8000/100)+1)=0;
Jcut(1:floor(8000/100)+1)=0;
[n1,pn1,cvg,iter,corp,covp]=leasqr(V(7000:8000),log(Jit2(7000:8000)),-46,"n1fit",.0001,20,ones(size(log(Jit2(7000:8000)))));
[n2,pn2,cvg,iter,corp,covp]=leasqr(V(1:1000),log(Jit2(1:1000)),-23,"n2fit",.0001,20,ones(size(log(Jit2(1:1000)))));


plot(vt,log(real(nfit2)),".g");
hold on
plot(V,log(Jit2),"-b");
axis([1,1.2,min(log(Jit2)),max(log(Jit2))]);
plot(V,V*q/(k*T)+pn1,".r",V,V*q/(2*k*T)+pn2,".r");
xlabel "Voc in V";
ylabel "J sc";
secondj=num2str(param_fit2(2));
firstj=num2str(param_fit2(1));
fhalf=strcat(firstj,"*(exp(qV/kT)-1) + ");
shalf=strcat(fhalf,secondj); 
tot=strcat(shalf,"*(exp(qV/2kT)-1)");
tit=strcat("Fit: ",tot);
title(tit);
savefig();
hold off;



