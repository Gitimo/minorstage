Ag=dlmread('refindAg.txt');
Au=dlmread('refindAu.txt');
GaAs=dlmread('refindGaAs.txt');
AlGaAs=dlmread('refindAlGaAs.txt');

Ag(:,1)=Ag(:,1)*1000;
Au(:,1)=Au(:,1)*1000;
GaAs(:,1)=flipud(GaAs(:,1)*1000);
AlGaAs(:,1)=flipud(AlGaAs(:,1)*1000);

wave=350:10:820;
indAg(1,:)=interp1(Ag(:,1),Ag(:,2),wave);
indAg(2,:)=interp1(Ag(:,1),Ag(:,3),wave);

indAu(1,:)=interp1(Au(:,1),Au(:,2),wave);
indAu(2,:)=interp1(Au(:,1),Ag(:,3),wave);

indGaAs(1,:)=interp1(GaAs(:,1),GaAs(:,2),wave);
indGaAs(2,:)=interp1(GaAs(:,1),GaAs(:,3),wave);

indAlGaAs(1,:)=interp1(AlGaAs(:,1),AlGaAs(:,2),wave);
indAlGaAs(2,:)=interp1(AlGaAs(:,1),AlGaAs(:,3),wave);

for m=1:48
nAg(m)=indAg(1,m)+i*indAg(2,m);
nAu(m)=indAu(1,m)+i*indAu(2,m);
nGaAs(m)=indGaAs(1,m)+i*indGaAs(2,m);
nAlGaAs(m)=indGaAs(1,m)+i*indAlGaAs(2,m);
lAlGaAs(m)=1000*indAlGaAs(1,m)/wave(m);
lGaAs(m)=10*indGaAs(1,m)/wave(m);

rAg(m)=abs(multidiel([1,nAlGaAs(m),nGaAs(m),nAg(m)],[lAlGaAs(m),lGaAs(m)],1))^2;
rAu(m)=abs(multidiel([1,nAlGaAs(m),nGaAs(m),nAu(m)],[lAlGaAs(m),lGaAs(m)],1))^2;
end
