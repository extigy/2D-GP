function [comx,comy] = gpege2dtcenterofmass_dt(dirarg,startno,endno)
for i=startno:1:endno
    [gridx,gridy,psi,potential] = gpeget2dPSI(dirarg,i);
	dens = psi.*conj(psi);
    fprintf('read %d\n',i);
    [comx(i/1),comy(i/1)] = gpegetcenterofmass(dens,gridx,gridy);
end

end
