function [p,Rsq]=linfit(X,Y)
#Compute linear fit and R² of the same
	p=polyfit(X,Y,1);
	yfit=polyval(p,X);
	yresid= Y-yfit;
	SSres=sum(yresid.^2);
	SStot=(length(Y)-1)*var(Y);
	Rsq=1-SSres/SStot;
end
