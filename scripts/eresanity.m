clear all; clf; format long;
# Define some physical constants
k=1.3806488e-23;
q=1.60217646e-19;
h=6.62606957e-34;
c=299792458;
Troom=293;
Tsun=6000;
am15=load("/home/timo/minorstage/data/am15short");

names=cellstr('');
nofcells=input("How many cells?");
Voc(1:nofcells)=0;
scov(1:nofcells)=0;
names{1,1}=input("Path to first data-file?","s");
Voc(1)=input("First cells Voc?");
scov(1)=input("First cells Surface-Coverage?");
fcell=dlmread(names{1,1});
npoints=size(fcell)(1);
wave=fcell(:,1);
eqearray=zeros(npoints,nofcells);
eqearray(:,1)=fcell(:,2)*scov(1);
for m=2:nofcells
	names{1,m}=input("Path to next file?","s");
	scov(m)=input("This cells Surface-Coverage?");
	data=dlmread(names{1,m});
	eqearray(:,m)=data(:,2)*scov(m);
	Voc(m)=input("This cells Voc?");
end
pden=numofphot(wave,am15(:,2));
clear data fcell am15 scov;

plot(wave,eqearray)
xlabel "Wavelength (nm)"
ylabel "External Quantum Efficiency"
legend(names,"location","south")

print(gen_pname(names,'_eqe','.png'),'-dpng');
print(gen_pname(names,'_eqe','.eps'),'-deps');
clf;

#longwave=(10:10:920000)*10^(-9); #wavelength-vector in m	
#longeqearray=zeros(size(longwave)(2),nofcells);
#for n=1:nofcells
#	for m=1:npoints
#		longeqearray(wave(m)/10,n)=eqearray(m,n);	#first row of EQE: lambda =350 assigned to 35th
#	end
#end


F_sun=2.16e-5*pi/2;			#geometric factors, angular dependence of bb-law
F_bb=pi;

J_sc=zeros(1,nofcells);
J_rad=zeros(1,nofcells);
Q_led=zeros(1,nofcells);
Voc_lim=zeros(1,nofcells);
lum=zeros(nofcells,size(wave)(1));
#eqearray=eqearray;
[fid,msg]=fopen(gen_pname(names,'','.txt'),'w');
fprintf(fid,'cell\tVoc(V)\tVlim(V)\tJ_sc(A/m^2)\tJ_rad(A/m^2)\tQ_led\n')
for m=1:nofcells
	J_sc(1,m)=q*trapz(flipud(wavetoen(wave*10^-9)),flipud(pden).*flipud(eqearray(:,m)));
	J_rad(1,m)=q*trapz(flipud(wavetoen(wave*10^-9)),blackbodyE(flipud(wavetoen(wave*10^-9)),F_bb,Troom).*flipud(eqearray(:,m)));
	Q_led(1,m)=J_rad(m)*exp(q*Voc(m)/(k*Troom))/J_sc(m);
	Voc_lim(1,m)=k*Troom/q * log(J_sc(1,m) / J_rad(1,m) +1);
	lum(m,:)=exp(q*Voc(m)/(k*Troom))*flipud(pden).*flipud(eqearray(:,m));
	fprintf(fid,'%s\t%.4g\t%g\t%g\t%g\t%g\n',names{1,m},Voc(m),Voc_lim(1,m),J_sc(1,m),J_rad(1,m),Q_led(1,m));
end
fclose(fid);


xmin=min(wave);
xmax=max(wave);
ymin=0;
ymax=1.1*max(max(lum));
plot(flipud(wave),lum);
axis([xmin,xmax,ymin,ymax]);
xlabel "Wavelength (nm)"
ylabel "Calculated spectral luminescence (A/m^2)"
legend(names,"location","northwest")

print(gen_pname(names,'_lum','.png'),'-dpng');
print(gen_pname(names,'_lum','.eps'),'-deps');
clf;
