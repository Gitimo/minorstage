clear all; clf; format long;

# Define some physical constants
global k=1.3806488e-23;
global q=1.60217646e-19;
global h=6.62606957e-34;
global c=299792458;
global Troom=293;
global Tsun=6000;

# Black-body energy density emission spectrum as function of wavelength in m
function [y] = bb(lam,Fs,T)
global k q h c;
	y = 8*Fs*h*c./(lam.^5) ./(exp(h*c./(lam*k *T)-1) ); # not sure about the 8
end
# Fuction to compute the number of photons
function [y] = numofphot(lam,Fs,T)
global k q h c;
	y = bb(lam,Fs,T).*lam/(h*c);
end
# Example EQE for W1689
EQE_1689=[350,0.2106928336
360,0.2301457638
370,0.2528985274
380,0.3048512132
390,0.3566626976
400,0.4071348396
410,0.4465933858
420,0.4896145631
430,0.5348703435
440,0.5753443344
450,0.5985453333
460,0.6223220811
470,0.6454429784
480,0.667584027
490,0.6868473876
500,0.70255083
510,0.7062906004
520,0.7132337738
530,0.7159509645
540,0.719052945
550,0.7202192665
560,0.7217732254
570,0.7223337092
580,0.7240661553
590,0.7237737724
600,0.7251560248
610,0.7257512675
620,0.7255840815
630,0.7268774147
640,0.7269128676
650,0.7265505948
660,0.7256640326
670,0.7250425359
680,0.7233006819
690,0.7208847828
700,0.7177654633
710,0.7161511134
720,0.7090033738
730,0.7128596609
740,0.7049346516
750,0.6982204498
760,0.7123020726
770,0.6948112422
780,0.6787893458
790,0.7078854339
800,0.736981522
810,0.7104851513
820,0.724192147
830,0.7444049968
840,0.7042466067
850,0.6429558937
860,0.637958733
870,0.5812911134
880,0.3144970891
890,0.0714055785
900,0.0192755004
910,0.0149613806
920,0.0056980765];
#W1688
EQE_1688=[350,0.0200133328
360,0.0109828343
370,0.0186919176
380,0.0155512493
390,0.0110166065
400,0.0138951277
410,0.0088050184
420,0.0079768259
430,0.0046475948
440,0.0057629112
450,0.0087415122
460,0.008832778
470,0.0145056795
480,0.0218672626
490,0.0279607592
500,0.0372175636
510,0.0485117289
520,0.0585402983
530,0.0714390316
540,0.0833273472
550,0.0964207185
560,0.1087318976
570,0.121922677
580,0.135874046
590,0.1489412532
600,0.1644164281
610,0.180547777
620,0.1964294488
630,0.211892315
640,0.2288818259
650,0.2453568681
660,0.2668926
670,0.2839491808
680,0.297600672
690,0.3107953973
700,0.3237631487
710,0.339343406
720,0.3550334838
730,0.3683078815
740,0.3801417336
750,0.3915064491
760,0.3995952867
770,0.4062665313
780,0.4114035278
790,0.4182320876
800,0.4250606475
810,0.4276094576
820,0.4271175156
830,0.425488329
840,0.422768264
850,0.4172130983
860,0.3894985062
870,0.2880523564
880,0.1365825386
890,0.0386944316
900,0.0098643869
910,0.0022668211
920,0.0008279732];

# Utility functions
function [y] = wavetoen(lambda)
global k q h c;
	y = h*c ./ lambda ;
end
function [y] =entowave(E);
global k q h c;
	y = h*c ./ E ;
end
function [ax]=replotyy(x1,y1,x2,y2)
	facx=1;	
	facy=1.1;
	[ax]=plotyy(x1,y1,x2,y2);
	set(ax(:),"xlim",facx*[750,950]);
	set(ax(1),"ylim",facy*[min(y1),max(y1)]);
	set(ax(2),"ylim",facy*[min(y2),max(y2)]);
	set(ax(:),"dataaspectratiomode","auto");
end

function ylabelyy(str,ax)
	set(ax(:),"ylabel",str)
end

### actual code

plot(EQE_1689(:,1),EQE_1689(:,2)*100,"b",EQE_1688(:,1),EQE_1688(:,2)*100,"g");
title "Quantum efficies. Blue is W1689, green W1688"
xlabel "wavelength in nm"
ylabel "EQE in %"
print -dpng eqes.png
disp("press return for next plot");
pause();

dlambda=(10:10:920000)*10^(-9); #wavelength-vector in m	
EQElong_1689=zeros(size(dlambda)(2),2);
EQElong_1688=zeros(size(dlambda)(2),2);
EQElong_1689(:,1)=dlambda;
EQElong_1688(:,1)=dlambda;
# fill new matrix with EQE data, zero where no data is given
for m=1:size(EQE_1689)(1)
	EQElong_1689(EQE_1689(m,1)/10,2)=EQE_1689(m,2);	#first row of EQE: lambda =350 assigned to 35th row of EQElong
	EQElong_1688(EQE_1688(m,1)/10,2)=EQE_1688(m,2);
end

F_sun=2.16e-5*pi;			#geometric factors, angular dependence of bb-law
F_bb=pi;

dlam=EQElong_1689(:,1);
J_sc_1689=q*trapz(dlam,numofphot(dlam,F_sun,Tsun).*EQElong_1689(:,2))
J_sc_1688=q*trapz(dlam,numofphot(dlam,F_sun,Tsun).*EQElong_1688(:,2))
J_rad_1689=q*trapz(dlam,numofphot(dlam,F_bb,Troom).*EQElong_1689(:,2))
J_rad_1688=q*trapz(dlam,numofphot(dlam,F_bb,Troom).*EQElong_1688(:,2))


Voc_1689=0.6;
Voc_1688=1.01;
Q_led_1689=J_rad_1689*exp(q*Voc_1689/(k*Troom))/J_sc_1689;
lum_1689=exp(q*Voc_1689/(k*Troom))*bb(dlam,F_bb,Troom).*EQElong_1689(:,2);
Q_led_1688=J_rad_1688*exp(q*Voc_1688/(k*Troom))/J_sc_1688
#Voc_lim_1689=k*Troom/q * log(J_sc_1689 / J_rad_1689 +1)
Voc_lim_1688=k*Troom/q * log(J_sc_1688 / J_rad_1688 +1)

lum_1689=exp(q*Voc_1689/(k*Troom))*bb(dlam,F_bb,Troom).*EQElong_1689(:,2);
lum_1688=exp(q*Voc_1688/(k*Troom))*bb(dlam,F_bb,Troom).*EQElong_1688(:,2);
ax=replotyy(dlam*10^(9),lum_1689,dlam*10^(9),lum_1688); #plotting in nm
title "Blue is W1689, green W1688, V_{OC,1689}=0.6V,V_{OC,1688}=1.01V "
xlabel "wavelength (nm)";
ylabel "luminescence (A/m^2/nm)";
print -dpng compareW1689_W1688.png
#disp("press return for next plot");
#pause();
