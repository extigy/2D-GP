function [FX,FY] = gpeget2dforce_dt(dirarg,startno,stride,endno,speed,nx,ny)
total = [];
pos = [];
neg = [];
for i=startno:stride:endno
    [gridx,gridy,dens,phase,potential] = gpeget2dWF(dirarg,i,speed,nx,ny);
    fprintf('read %d\n',i);
    j = (i+(stride-startno))/stride;
    [fxt,fyt] = gpeget2dforce(dens,phase,gridx,gridy,potential);
    FX(j) = fxt;
    FY(j) = fyt;
end
end