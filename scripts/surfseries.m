############## Define some physical constants 	##################
global k=1.3806488e-23;
global q=1.60217646e-19;
global h=6.62606957e-34;
global c=299792458;
global Troom=300;
global Tsun=5777;
############## Defined some physical constants	#################

# utility functions
function [val,pos]=findclosest(array,target)
	[val,pos]=min(abs(array-target));
end


############## Load necessary data		##################
fspec=load("am15full");
fspec=fspec(:,[1,3]);

fid=fopen('VOC');
Voc=fscanf(fid,'%f\n');
fclose(fid);
fid=fopen('A');
A=fscanf(fid,'%f\n');
fclose(fid);
P=fopen('P');
Voc=fscanf(fid,'%f\n');
fclose(fid);
R=P./A;

nofcells=input("How many cells?");
names=cellstr('');
#create cell to contain {1,1}=Voc; {2,:}=wavelength (m); 
#{3,:}=eqe; {4,:}=energy; {5,:}=eqeflipped; {6,:}=numberofphotons
datacell=cell(6,nofcells);

detectionlimit=1.5360e-03;

for m=1:nofcells
	names{1,m}=input("Path to cell-eqe Data? ","s");
	datacell{1,m}=Voc(m);
	data=dlmread(names{1,m});
	data(find(data(:,2)<=detectionlimit),2)=0; #set all EQEs below detection limit to zero
	datacell{2,m}=data(:,1)*10^(-9);
	datacell{3,m}=data(:,2);
	datacell{4,m}=wavetoen(flipud(datacell{2,m}));		
	datacell{5,m}=flipud(datacell{3,m});
	datacell{6,m}=interp1(fspec(:,1),fspec(:,2),datacell{2,m}*10^9);
	[val,start]=findclosest(fspec(:,1),min(datacell{2,m}*10^9));
	[val,stop]=findclosest(fspec(:,1),max(datacell{2,m}*10^9));
	should=trapz(fspec(start:stop,1),fspec(start:stop,2));
	is=trapz(datacell{2,m}*10^9,datacell{6,m});
	datacell{6,m}=datacell{6,m}*should/is .* datacell{2,m} / (h*c) ;
	clear data
endfor
############## Loaded necessary data		##################

############## do the calculations		##################
#create cell to contain results: {1,1}=J_sc; 
#{2,1}=J_0_r; {3,1}=Voc_limit; {4,1}=ERE; {5,:}=lum
results=cell(5,nofcells);

for m=1:nofcells
	results{1,m}=q*trapz(datacell{2,m}*10^9,datacell{6,m} .* datacell{3,m}) / 10 ;	
	results{2,m}=q*trapz(datacell{4,m}, blackbodyE(datacell{4,m},Troom) .* datacell{5,m} ) / 10 ;
	results{3,m}=k*Troom/q * log( results{1,m}/results{2,m} +1 ) ;
	results{4,m}=exp( q * datacell{1,m} / (k*Troom) ) * results{2,m} / results{1,m} * 100;
	results{5,m}=q*exp(q*results{3,m}/(k*Troom)-1) * blackbodyE(datacell{4,m},Troom) .* flipud(datacell{3,m}) /10;
	#this one goes wrong because of energy vs nm!
endfor
############## done the calculations		##################


############## print human readable output	##################
[fid,msg]=fopen(gen_pname(names,'','.txt'),'w');
fprintf(fid,'cell\t\tArea\tPerimeter\tRatio\tVoc(V)\t\tVlim(V)\tJ_sc(mA/cm^2)\tJ_0_r\tERE(%%)\n\n');
for m=1:nofcells
	fprintf(fid,'%s\t\t%.3g\t%.3g\t%.3g\t%.4g\t\t%g\t%g\t%g\t%g\n',names{1,m},A(m),P(m),R(m),datacell{1,m},results{3,m},results{1,m},results{2,m},results{4,m});
endfor
# Close file, save data
fclose(fid);
