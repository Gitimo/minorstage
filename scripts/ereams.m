clear all; clf; format long;
############## Define some physical constants ##################
global k=1.3806488e-23;
global q=1.60217646e-19;
global h=6.62606957e-34;
global c=299792458;
global Troom=300;
global Tsun=5777;
############## Load spectra ####################################

load("am15full");
load("am153nm");
load("am1510nm")
# Renormalise shortened spec to correct energy-density: 1kW/m^2
# Do I have to normalise? W/o normalising: J_sc=28mA/cm^2; with: J_sc=42mA/cm^2
## As is: normalised to 683.56W/m^2, i.e. int 350 to 920 in am15full
wave_3=am153nm(:,1); #wavelength in m (useful for plotting, finding, integrating)
wave_10=am1510nm(:,1); 
wavem_3=wave_3*10^-9;
wavem_10=wave_10*10^-9;
energy_3=flipud(wavetoen(wavem_3));
energy_10=flipud(wavetoen(wavem_10));
start=find(am15full==(min(wave_3)));
stop=find(am15full==(max(wave_3)));
shouldpower(1)=trapz(am15full(start:stop,1),am15full(start:stop,3));
start=find(am15full==(min(wave_10)));
stop=find(am15full==(max(wave_10)));
shouldpower(2)=trapz(am15full(start:stop,1),am15full(start:stop,3));
specpower_3=am153nm(:,3)*shouldpower(1)/trapz(am153nm(:,1),am153nm(:,3));
specpower_10=am1510nm(:,3)*shouldpower(2)/trapz(am1510nm(:,1),am1510nm(:,3));
clear am15full am153nm am1510nm
specphot_3= specpower_3.*wavem_3/(h*c);
specphot_10= specpower_10.*wavem_10/(h*c);
#specphot=numofphot(wave,specpower)*10^9; # to convert to SI, i.e. photons per m^2 *per m* instead of ...per nm
############## Get cell data ###################################


names=cellstr('');
nofcells=input("How many cells?");
nofpoints_3=size(wave_3)(1);
eqearray_3=zeros(nofpoints_3,nofcells);
nofpoints_10=size(wave_10)(1);
eqearray_10=zeros(nofpoints_10,nofcells);

Voc=zeros(1,nofcells);
res=zeros(1,nofcells);
Sco=0;
for m=1:nofcells
	names{1,m}=input("Path to cell-eqe Data? ","s");
	Voc(1,m)=input("Voc? (V)");
	%Sco=input("Surface-Coverage? ");
	data=dlmread(names{1,m});
	res(1,m)=data(2,1)-data(1,1);
	if res(1,m)==3
		eqearray_3(:,m)=data(:,2);
	elseif res==10
		eqearray_10(:,m)=data(:,2);
	endif
endfor

clear data
############# Plot Eqe-data ####################################

ans=input("Plot eqe?y/n","s")
if ans=="y"
	plot(wave_3,eqearray)
	xlabel "Wavelength (nm)"
	ylabel "External Quantum Efficiency, surface coverage cleared"
	legend(names,"location","south")

	print(gen_pname(names,'_eqe','.png'),'-dpng');
	print(gen_pname(names,'_eqe','.eps'),'-deps');
	clf;
endif
########################## Now do the real calculations #################


# Preallocate vectors and matrices for calculated parameters
J_sc=zeros(1,nofcells);
J_0_r=zeros(1,nofcells);
J_0_nr=zeros(1,nofcells);
Q_led=zeros(1,nofcells);
Voc_lim=zeros(1,nofcells);
ERE=zeros(1,nofcells);
lum_3=zeros(nofpoints_3,nofcells);
lum_10=zeros(nofpoints_10,nofcells);
# Open txt file for saving parameters
[fid,msg]=fopen(gen_pname(names,'','.txt'),'w');
fprintf(fid,'cell\tVoc(V)\tVlim(V)\tJ_sc(mA/cm^2)\tJ_0_r\tJ_0_nr\tERE(%%)\n');

for m=1:nofcells
	if res(1,m)==3
		J_sc(1,m)=q*trapz(wave_3,specphot_3.*eqearray_3(:,m))/10;
		J_0_r(1,m)=q*trapz(energy_3,blackbodyE(energy_3,Troom).*flipud(eqearray_3(:,m)))/10;
		lum_3(:,m)=q*exp(q*Voc(m)/(k*Troom))*blackbodyL(wavem_3,Troom).*eqearray_3(:,m); #check units!!! WHAT THE Fuck?!?
	elseif res(1,m)==10
		J_sc(1,m)=q*trapz(wave_10,specphot_10.*eqearray_10(:,m))/10;
		J_0_r(1,m)=q*trapz(energy_10,blackbodyE(energy_10,Troom).*flipud(eqearray_10(:,m)))/10;
		lum_10(:,m)=q*exp(q*Voc(m)/(k*Troom))*blackbodyL(wavem_10,Troom).*eqearray_10(:,m); #check units!!! WHAT THE Fuck?!?
	else
		disp("Error with resolution");
		break;
	endif
	
	Voc_lim(1,m)=k*Troom/q * log(J_sc(1,m) / J_0_r(1,m) +1);
	ERE(1,m)=J_0_r(1,m)*exp(q*Voc(m)/(k*Troom))/J_sc(1,m)*100;
	J_0_nr(1,m)=J_sc(1,m)/(exp(Voc(1,m)*q/(k*Troom))-1)-J_0_r(1,m);
	#dVoc=-k*Troom/q*log(ERE(1,m)/100); can be used as sanity check: Voc+dVoc=Voc_lim -> works!
	fprintf(fid,'%s\t%.4g\t%g\t%g\t%g\t%g\t%g\n',names{1,m},Voc(m),Voc_lim(1,m),J_sc(1,m),J_0_r(1,m),J_0_nr(1,m),ERE(1,m));
endfor
# Close file, save data
fclose(fid);
