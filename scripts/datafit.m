clear all; clf;

# Define some physical constants
k=1.3806488e-23;
q=1.60217646e-19;
T=293;

# Function that will be fit
function [y]=double_exp(x,par)
  y=par(1)*exp(-par(2)*x)-par(1)*exp(-par(3)*x);
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

# Generic weights matrix
weights=ones(size(ir1));
 
# Perform the fits
pin=[j,t1,t2];
[fir1,pir1,cvg,iter,corp,covp]=leasqr(t,ir1,pin,"double_exp",.0001,20,weights);
[fir2,pir2,cvg,iter,corp,covp]=leasqr(t,ir2,pin,"double_exp",.0001,20,weights);
[fit1,pit1,cvg,iter,corp,covp]=leasqr(t,it1,pin,"double_exp",.0001,20,weights);

# Calculate reference cell's I_sc - is that correct?
it2=fit1.*fir2./fir1;

# Calculate illumination in "suns"
r=0.04248 ; equiv=13.7e-3 ; illu=fir2/(r*equiv);

# Plotting ... saving does not work properly
tcut(1:floor(8000/100)+1)=0;
ir1(1:floor(8000/100)+1)=0;
ir2cut(1:floor(8000/100)+1)=0;

for m=1:8000
	if (m==30)
		tcut(m)=t(m); ir1cut(m)=ir1(m); ir2cut(m)=ir2(m);
	elseif (mod(m,100)==0)
		tcut(m)=t(m); ir1cut(m)=ir1(m); ir2cut(m)=ir2(m);
	end
end

hold on;
hg = errorbar(tcut,ir1cut,0.0005);
set(hg,"marker","+");
set(hg,"linestyle","none");
axis([0,max(t),min(fir1),max(fir1)]);
plot(t,fir1);
hold off;
savefig();

hold on;
hg = errorbar(tcut,ir2cut,0.0005);
set(hg,"marker","+");
set(hg,"linestyle","none");
axis([0,max(t),min(fir2),max(fir2)]);
plot(t,fir2);
hold off;
savefig();

plot(t(100:8000),illu(100:8000));
axis([0,max(t),1.01*min(illu(100:8000)),1.01*max(illu(100:8000))]);
savefig();

plot(vt(100:8000),log(illu(100:8000)));
axis([min(vt(100:8000)),max(vt(100:8000)),1.01*min(log(illu(100:8000))),1.01*max(log(illu(100:8000)))]);
savefig();

plot(vt(100:8000),log(it2(100:8000)));
axis([min(vt(100:8000)),max(vt(100:8000)),1.01*min(log(it2(100:8000))),1.01*max(log(it2(100:8000)))]);
savefig();38
