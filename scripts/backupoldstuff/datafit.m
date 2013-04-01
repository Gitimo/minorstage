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
	clear plotname;
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

# Calculate illumination in "suns" and convert voltages to currents
r=0.04248 ; equiv=13.7e-3 ; illu=fir2/(r*equiv);
ir1=ir1/r; ir2=ir2/r;
it1=it1/0.05; it2=it2/0.05; #Have to find out Resistance precisely!

# fitting to find idealities and currents
function [y]=currentfit_n1(x,par)
	k=1.3806488e-23;
	q=1.60217646e-19;
	T=293;
	y=q*x/(k*T)+par;
end
function [y]=currentfit_n2(x,par)
	k=1.3806488e-23;
	q=1.60217646e-19;
	T=293;
	y=q*x/(2*k*T)+par;
end
weights=ones(size(vt(300:8000)));
[nis1,pnis1,cvg,iter,corp,covp]=leasqr(vt(300:8000),log(it2(300:8000)),-60,"currentfit_n1",.0001,20,weights);
weights=ones(size(vt(100:300)));
[nis2,pnis2,cvg,iter,corp,covp]=leasqr(vt(100:300),log(it2(100:300)),-30,"currentfit_n2",.0001,20,weights);

# Plotting ... saving does not work properly
#tcut(1:floor(8000/100)+1)=0;
#ir1(1:floor(8000/100)+1)=0;
#ir2cut(1:floor(8000/100)+1)=0;

#for m=1:8000
#	if (m==30)
#		tcut(m)=t(m); ir1cut(m)=ir1(m); ir2cut(m)=ir2(m);
#	elseif (mod(m,100)==0)
#		tcut(m)=t(m); ir1cut(m)=ir1(m); ir2cut(m)=ir2(m);
#	end
#end

#hold on;
#hg = errorbar(tcut,ir1cut,0.0005);
#set(hg,"marker","+");
#set(hg,"linestyle","none");
#axis([0,max(t),min(fir1),max(fir1)]);
#plot(t,fir1);
#hold off;
#savefig();

#hold on;
#hg = errorbar(tcut,ir2cut,0.0005);
#set(hg,"marker","+");
#set(hg,"linestyle","none");
#axis([0,max(t),min(fir2),max(fir2)]);
#plot(t,fir2);
#hold off;
#savefig();

plot(t(100:8000),illu(100:8000));
axis([0,max(t),1.01*min(illu(100:8000)),1.01*max(illu(100:8000))]);
savefig();

plot(vt(100:8000),log(illu(100:8000)));
axis([min(vt(100:8000)),max(vt(100:8000)),1.01*min(log(illu(100:8000))),1.01*max(log(illu(100:8000)))]);
savefig();

plot(vt(100:8000),log(it2(100:8000)));
axis([min(vt(100:8000)),max(vt(100:8000)),1.01*min(log(it2(100:8000))),1.01*max(log(it2(100:8000)))]);
hold on;
plot(vt(100:8000),vt(100:8000)*q/(k*T)+pnis1)
plot(vt(100:8000),vt(100:8000)*q/(2*k*T)+pnis2)
xlabel "Voc in V";
ylabel "ln(I sc)";
title "Semilog plot, quality factor given, J01 and J02 fitted";
fh="J01=";sh=num2str(exp(pnis1));tot=strcat(fh,sh);
text(1.12,-4,tot);
fh="J02=";sh=num2str(exp(pnis2));tot=strcat(fh,sh);
text(1.12,-4.5,tot);
savefig();


