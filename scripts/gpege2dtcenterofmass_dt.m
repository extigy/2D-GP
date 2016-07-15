function [comx,comy] = gpege2dtcenterofmass_dt(dirarg,startno,endno,speed,nx,ny)
for i=startno:1:endno
    [gridx,gridy,dens,phase,potential] = gpeget2dWF(dirarg,i,speed,nx,ny);
    fprintf('read %d\n',i);
    [comx(i/1),comy(i/1)] = gpegetcenterofmass(dens,gridx,gridy);
end

end
