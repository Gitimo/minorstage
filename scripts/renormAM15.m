function [y]=renormAM15(am15short,am15full)
	start=find(am15full==min(am15short)(1));
	stop=find(am15full==max(am15short)(1));
	
	is=trapz(am15short(:,1),am15short(:,3));
	should=trapz(am15full(start:stop,1),am15full(start:stop,3));

	y=am15short(:,3)*should/is;
end
