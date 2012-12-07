clear all; clf;

# Define some physical constants
global k=1.3806488e-23;
global q=1.60217646e-19;
global T=293;
global dVdt;

# Function that will be fit to initial data
function [y]=diodcap(x,par)
global k q T dVdt;
	y=par(1).*(exp(q.*x./(k.*T))-1)+par(2).*(exp(q.*x/(2.*k.*T))-1)+par(3)./(sqrt(par(4)-x)).*(1./(2.*(par(3)-x)).*x.*dVdt+dVdt)+par(5).*exp(par(6).*q.*x./(k.*T)).*(dVdt.*(par(6).*q./(k.*T).*x+1));

end

# Function to approximate series by lines

function [out]=divideup(x,y,pieces)
	steps=size(x)(1)/pieces;
	if (abs(mod(steps,1))==0)
	clf; hold on;
	#plot(x,y,".g")
	for m=1:pieces
		st=(m-1)*steps+1;
		sp=m*steps;
		out(m,:)=polyfit(x(st:sp),y(st:sp),1);
		plot(x(st:sp),(x(st:sp)*out(m,1)+out(m,2)-y(st:sp))./y(st:sp));
	end

	else
	disp("Number of data points in vector cannot evenly be divided")
	end

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


#Load data and assign vectors
fhalf="../data/" ; shalf=input("Filename?","s") ; fname=strcat(fhalf,shalf);
m=load(fname);
t=m(:,1) ; ir1=m(:,2) ; ir2=m(:,3) ; it1=m(:,4) ; vt=m(:,5) ;
tsurf=input("Test-cel surface area in cm^2?");
Rref=0.04248 ; equivref=13.7e-3 ;
Rtest=0.0631 ; equivref_t=20.7e-3;
ir1=ir1/Rref; ir2=ir2/Rref; 
it1=it1/(Rtest*tsurf);
it2=it1.*ir2./ir1;

pieces=input("Number of lines to be used?");
dVdt=divideup(t,vt,pieces);
tpieces=min(t):max(t)/pieces:max(t);
dVdt=interp1(tpieces,dVdt,t);
dVdt=dVdt(:,1);


# Generic weights matrix
wt=ones(size(it2));

pin=[1e-20,1e-10,0.08,5,1e-12,1];
# Perform the fits
[model,parameters,cvg,iter,corp,covp]=leasqr(vt,it2,pin,"diodcap",.0001,20,wt);

