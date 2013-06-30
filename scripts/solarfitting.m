function [ddres_J]=solarfitting(Jscm,Vocm,fullJV)
############## Define some physical constants 	##################
	pkg load physicalconstants;
	global k T q;	
	T=298;
	k=physical_constant("Boltzmann constant");
	q=physical_constant("elementary charge");
	pkg unload physicalconstants;
############## Defined some physical constants	#################

	Jsc=max(Jscm);
	Voc=max(Vocm);
	printf("Measured Voc: %g V\n", Voc);

############## initial, linear fit of Jsc-Voc data ##############
	[lin_n,lin_J0,lin_r2]=fit_jscvoc(Jscm,Vocm);
	#check Voc
	lin_Voc=lin_n*k*T/q*log(Jsc/lin_J0+1);
	printf("Linear Voc: %g V\n", lin_Voc);

############## initial, non-linear sd-fit of Jsc-Voc data #######
	[nlin_n_1,nlin_J0_1,nlin_J_1,nlin_r_1]=fit_sd(Jscm,Vocm,lin_J0,lin_n);	

############## second , non-linear sd-fit of Jsc-Voc data #######
	#	using initial's parameters as guesses
	[nlin_n_2,nlin_J0_2,nlin_J_2,nlin_r_2]=fit_sd(Jscm,Vocm,nlin_J0_1,nlin_n_1);
	if nlin_r_2>=nlin_r_1
		sd_n=nlin_n_2; sd_J0=nlin_J0_2; sd_J=nlin_J_2; sd_r2=nlin_r_2;
		fitinfo{1,1}="Second sd-fit was better than first"; 
	else
		sd_n=nlin_n_1; sd_J0=nlin_J0_1; sd_J=nlin_J_1; sd_r2=nlin_r_1;
		fitinfo{1,1}="Second sd-fit was better than first"; 
	endif
	clear nlin_n_1 nlin_n_2 nlin_J0_1 nlin_J0_2; 
	#check Voc
	sd_Voc=sd_n*k*T/q*log(Jsc/sd_J0 +1);	
	printf("Single Diode-Voc: %g V\n", sd_Voc);
	
############## initial, non-linear dd-fit of Jsc-Voc data #######
	#	using sd-fit parameters as guesses
	[dd_J01,dd_J02,dd_J,dd_r]=fit_dd(Jscm,Vocm,sd_J0);
	#check Voc, 'analytically' and numerically
	dd_Voc=[analvoc(Jsc,dd_J01,dd_J02),...
		fzero(@ (x)numvoc(Jsc,dd_J01,dd_J02,[0,0],x),[0.5,2])];
	printf("Double Diode-Voc. Analytical: %g V Numerical: %g V\n", dd_Voc(1),dd_Voc(2));

############## initial, non-linear ddres-fit of full JV data ########
	#check for full JV-data
	if !isempty (fullJV)
		fullV=fullJV(:,1);
		fullJ=fullJV(:,2);
		guess=sqrt([dd_J01,dd_J02,1,10000,Jsc]);
		[par,retJ,solverinfo]=fsolve(@(x) fit_dd_res(fullJ,fullV,x),guess,...
					optimset("FunValCheck","on","TolX",1e-5));
		[par,retJ,solverinfo]=fsolve(@(x) fit_dd_res(fullJ,fullV,x),par,...
					optimset("FunValCheck","on","TolX",1e-5));

		ddres_J01=par(1)^2; ddres_J02=par(2)^2; ddres_Rs=par(3)^2; ddres_Rsh=par(4)^2; 
		ddres_Jsc=par(5)^2; ddres_J=retJ+fullJ;		
		ddres_Voc=[fzero(@ (x)numvoc(ddres_Jsc,ddres_J01,ddres_J02,[1,ddres_Rsh],x),[0.5,2])];
		printf("Double Diode-Res-Voc. Numerical: %g V\n", ddres_Voc);
	endif
endfunction
