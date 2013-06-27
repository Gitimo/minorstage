# Define some physical constants
global k=1.3806488e-23;
global q=1.60217646e-19;
global T=293;

# Single Diode Model
function [y]=single_diode(x,par)
global k q T
  y=par(1)*(exp(q*x/(par(2)*k*T))-1);
end

# Double Diode Model
function [y]=double_diode(x,par)
global k q T
  y=par(1)*(exp(q*x/(k*T))-1) + par(2)*(exp(q*x/(2*k*T))-1);
end

ca=[14.025, 12.637, 9.489, 7.743, 4.925, 4.168];
CA=[.997, .973, .962, .952, .935, .928];
cb=[7.35, 6.56, 5.18, 4.24, 2.45, 2.13]*0.2858/0.1429;
CB=[.977, .972, .962, .954, .93, .923];
cc=[17.1,15.46,12.84, 9.81, 6.02 5.27];
CC=[.992,.989,.981,.969,.947,.941];
cD=[16.27,14.53,11.4,9.32,5.41,4.81];
CD=[.986,.981,.97,.962,.937,.931];
ce=[18.12,16.12,12.7,10.35,6.09,5.21];
CE=[1.002,.996,.985,.975,.95,.942];
cf=[19.06,17,13.3,10.89,6.28,5.47];
CF=[1.025,1.02,1.01,1.001,.975,.968];
cg=[19.14,17.19,13.34,10.86,6.4,5.6];
CG=[1.035,1.031,1.022,1.014,.991,.984];

%[f,p,cvg,iter,corp,covp,covr,stdresid,Z,r2]=leasqr(x,y,pin,F,{stol,niter,wt,dp,dFdp,options})

[fA,pA,cvg,iter,corp,covp,covr,stdresid,Z,r2A]=leasqr(CA,ca,[3.4006e-08,1.9],'single_diode');
disp('sd A')
[fB,pB,cvg,iter,corp,covp,covr,stdresid,Z,r2B]=leasqr(CB,cb,[1.1387e-09,1.7],'single_diode');
disp('sd B')
[fC,pC,cvg,iter,corp,covp,covr,stdresid,Z,r2C]=leasqr(CC,cc,[2.6576e-09,1.7],'single_diode');
disp('sd C')
[fD,pD,cvg,iter,corp,covp,covr,stdresid,Z,r2D]=leasqr(CD,cD,[4.6661e-09,1.7],'single_diode');
disp('sd D')
[fE,pE,cvg,iter,corp,covp,covr,stdresid,Z,r2E]=leasqr(CE,ce,[1.4502e-08,1.9],'single_diode');
disp('sd E')
[fF,pF,cvg,iter,corp,covp,covr,stdresid,Z,r2F]=leasqr(CF,cf,[3.3800e-09,1.8],'single_diode');
disp('sd F')
[fG,pG,cvg,iter,corp,covp,covr,stdresid,Z,r2G]=leasqr(CG,cg,[2.5101e-10,1.6],'single_diode');
disp('sd G')

[fA2,pA2,cvg,iter,corp,covp,covr,stdresid,Z,r2A2]=leasqr(CA,ca,[2e-17,6e-8],'double_diode');
disp('dd A')
[fB2,pB2,cvg,iter,corp,covp,covr,stdresid,Z,r2B2]=leasqr(CB,cb,[pB(1)*0.999,pB(1)*0.001],'double_diode');
disp('dd B')
[fC2,pC2,cvg,iter,corp,covp,covr,stdresid,Z,r2C2]=leasqr(CC,cc,[pC(1)*0.999,pC(1)*0.001],'double_diode');
disp('dd C')
[fD2,pD2,cvg,iter,corp,covp,covr,stdresid,Z,r2D2]=leasqr(CD,cD,[pD(1)*0.999,pD(1)*0.001],'double_diode');
disp('dd D')
[fE2,pE2,cvg,iter,corp,covp,covr,stdresid,Z,r2E2]=leasqr(CE,ce,[pE(1)*0.999,pE(1)*0.001],'double_diode');
disp('dd E')
[fF2,pF2,cvg,iter,corp,covp,covr,stdresid,Z,r2F2]=leasqr(CF,cf,[pF(1)*0.999,pF(1)*0.001],'double_diode');
disp('dd G')
[fG2,pG2,cvg,iter,corp,covp,covr,stdresid,Z,r2G2]=leasqr(CG,cg,[pG(1)*0.999,pG(1)*0.001],'double_diode');

