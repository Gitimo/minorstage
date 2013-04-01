function [pname]=gen_pname(fnames,add,fform)
	#Utility-function to generate nice filenames
	pname=fnames{1,1};
	if size(fnames)(2)<4
		for m=2:size(fnames)(2)
			pname=strcat(pname,fnames{1,m});	
		end
		pname=strcat(pname,add,fform);
	else
	pname=strcat(pname,'_to_',fnames{1,size(fnames)(2)},add,fform);
	end
end
