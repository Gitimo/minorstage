%% Shockley-Queisser (using Timo's)
% Complicated matmathics to calculate the Shockley Queisser limit of
% efficiency over bandgap for photovoltaic single layer

clear all
close all
clc
%% constants
T=296;
Tc=273;%K
Ts=6000;%K
Planck=6.626E-34; %Js
PlanckeV=4.135E-15; %eV s
LS=2.99E8; %lightspeed m/s
q=1.6021E-19; %Coulomb
kb=1.380E-23; % J K

%Load Data variable
load AM15.mat


%% Convert to photon density using E*lamdba/(Planck*LS)
    for n=1:1:size(Data,1);
        SPhotonsPerTEA(n,1)=Data(n,1); %wavelength
        SPhotonsPerTEA(n,2)=Data(n,2)*(Data(n,1)*10^-9 )/(Planck*LS); %photonflux
        %SPhotonsPerTEA(n,2) = Data(n,2).*1/EPhotons(n,2)*Planck*LS/(EPhotons(n,2)^2);
    end

%% Amount of photons per bandgap meaning the sum of the amount of photons from that band gap
for k=1:size(Data,1);
    Amount(k,1)=Data(k,1); %wavelength
    Amount(k,2)=sum(SPhotonsPerTEA(1:k,2));% amount of photons from wavelength 280nm till specific bandgap
end
%Just some renaming and combination
SolarPhotonsAboveGap(:,1)=Data(:,1); %wavelength
SolarPhotonsAboveGap(:,2)=1240./Data(:,1); %Eg (eV)
SolarPhotonsAboveGap(:,3)=Amount(:,2); %amount of photons
%Plot amount of Photons
plot(SolarPhotonsAboveGap(:,1),SolarPhotonsAboveGap(:,3))
xlabel('Wavelength (nm)');
ylabel('Amount of Photons')

%Calculate Isc directly from amount of photons
Isc=q.*SolarPhotonsAboveGap(:,3); % (A)
%plot
figure
plot(SolarPhotonsAboveGap(:,1),Isc)
xlabel('Wavelength (nm)','fontsize',16);
ylabel('Isc (A)','fontsize',16)

        %possibility of calculating energy
        %Energy=(Amount(:,2)*(Planck*LS))./(Data(:,1)*10^-9);

        %plot(SolarPhotonsAboveGap(:,1),Energy);


%% IV characteristic
Fg=2.19E-5;% constant

Pin=sum(Data(:,2));%Sum move all energy from AM1.5

Eg(:,1)=Planck*LS./(Data(:,1)*10^-9);% Calculating Eg (J) from wavelength

%Calculating Voc, voc and eff for all wavelengths/bandgaps A lot of steps
%for a long equations only the last is interesting
for l=1:size(Data,1);
       
    Voc(l,1)=Data(l,1);
    Voc(l,2)=((kb*T)/q);
    Voc(l,3)=Fg*(Ts/Tc)*((Eg(l,1)+kb*Ts)^2/(Eg(l,1)+kb*Tc)^2);
    Voc(l,4)=((Eg(l,1)/kb)*((1/Tc)-(1/Ts)));
    Voc(l,5)=exp(Voc(l,4));
    Voc(l,6)=Voc(l,3)*Voc(l,5);
    Voc(l,7)=Voc(l,6)+1;
    Voc(l,8)=log(Voc(l,7));
    Voc(l,9)=Voc(l,2)*Voc(l,8);

    voc(l,1)=Voc(l,9)*q/(kb*200);

    eff(l,1)=(Voc(l,9)*Isc(l,1)*(voc(l,1)-log(voc(l,1)+0.72)/(voc(l,1)+1)))/Pin;% Eff for different wavelengths


end
%plot
figure
plot(SolarPhotonsAboveGap(:,1),eff);
xlabel('Wavelength (nm)','fontsize',16);
ylabel('Effiecieny (%)','fontsize',16)

% Now transfer to eV's
EVeff(:,1)=flipud(1240./Data(:,1));
EVeff(:,2)=flipud(eff);

%plot
figure
plot(EVeff(:,1),EVeff(:,2));
xlabel('Band Gap (eV)','fontsize',16);
ylabel('Effiecieny (%)','fontsize',16)

