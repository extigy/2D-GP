function [total,pos,neg] = gpeget2dvort_homg_dt(dirarg,startno,stride,endno,speed,nx,ny)
total = [];
pos = [];
neg = [];
for i=startno:stride:endno
    [gridx,gridy,dens,phase,potential] = gpeget2dWF(dirarg,i,speed,nx,ny);
    fprintf('read %d\n',i);
    j = (i+(stride-startno))/stride;
    [xlocs,ylocs,pol] = gpeget2dvort_homg(dens,phase,gridx,gridy,potential);
    total(j) = length(pol);
    neg(j) = -sum(pol-1)/2;
    pos(j) = length(pol) - neg(j);
end
%imagesc(dens)

end