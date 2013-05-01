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
load("am15short");
# Renormalise shortened spec to correct energy-density: 1kW/m^2
# Do I have to normalise? W/o normalising: J_sc=28mA/cm^2; with: J_sc=42mA/cm^2
## As is: normalised to 683.56W/m^2, i.e. int 350 to 920 in am15full
wave=am15short(:,1); #wavelength in m (useful for plotting, finding, integrating)
wavem=wave*10^-9;
energy=flipud(wavetoen(wavem));
start=find(am15full==(min(wave)));
stop=find(am15full==(max(wave)));
shouldpower=trapz(am15full(start:stop,1),am15full(start:stop,3));
specpower=am15short(:,3)*shouldpower/trapz(am15short(:,1),am15short(:,3));

clear am15full am15short
specphot= specpower.*wavem/(h*c);
#specphot=numofphot(wave,specpower)*10^9; # to convert to SI, i.e. photons per m^2 *per m* instead of ...per nm
############## Get cell data ###################################
names=cellstr('');
nofcells=input("How many cells?");
nofpoints=size(wave)(1);
eqearray=zeros(nofpoints,nofcells);
Csize=zeros(1,nofcells);
Voc=zeros(1,nofcells);
Jscmeas=zeros(1,nofcells);
Sco=0;
for m=1:nofcells
	names{1,m}=input("Path to cell-eqe Data? ","s");
	Voc(1,m)=input("Voc? (V)");
	Jscmeas(1,m)=input("Jsc? (mA/cm^2)");
	%Sco=input("Surface-Coverage? ");
	data=dlmread(names{1,m});
	eqearray(:,m)=data(:,2);%*Sco;
end
Csize=Csize*10^-4; #From cm^2 to m^2
clear Sco data
############# Plot Eqe-data ####################################
plot(wave,eqearray)
xlabel "Wavelength (nm)"
ylabel "External Quantum Efficiency, surface coverage cleared"
legend(names,"location","south")

print(gen_pname(names,'_eqe','.png'),'-dpng');
print(gen_pname(names,'_eqe','.eps'),'-deps');
clf;

########################## Now do the real calculations #################
# Preallocate vectors and matrices for calculated parameters
J_sc=zeros(1,nofcells);
J_rad=zeros(1,nofcells);
Q_led=zeros(1,nofcells);
Voc_lim=zeros(1,nofcells);
ERE=zeros(1,nofcells);
lum=zeros(nofpoints,nofcells);
# Open txt file for saving parameters
[fid,msg]=fopen(gen_pname(names,'','.txt'),'w');
fprintf(fid,'cell\tVoc(V)\tVlim(V)\tJ_sc\tJ_sclim(mA/cm^2)\tJ_rad(mA/cm^2)\tERE(%%)\n');

for m=1:nofcells
	J_sc(1,m)=q*trapz(wave,specphot.*eqearray(:,m))/10;
	J_rad(1,m)=q*trapz(energy,blackbodyE(energy,Troom).*flipud(eqearray(:,m)))/10;
	Voc_lim(1,m)=k*Troom/q * log(J_sc(1,m) / J_rad(1,m) +1);
	lum(:,m)=q*exp(q*Voc(m)/(k*Troom))*blackbodyL(wavem,Troom).*eqearray(:,m); #check units!!! WHAT THE Fuck?!?
	ERE(1,m)=J_rad(1,m)*exp(q*Voc(m)/(k*Troom))/J_sc(1,m)*100;
	#dVoc=-k*Troom/q*log(ERE(1,m)/100); can be used as sanity check: Voc+dVoc=Voc_lim -> works!
	fprintf(fid,'%s\t%.4g\t%g\t%g\t%g\t%g\t%g\n',names{1,m},Voc(m),Voc_lim(1,m),Jscmeas(1,m),J_sc(1,m),J_rad(1,m),ERE(1,m));
end
# Close file, save data
fclose(fid);
