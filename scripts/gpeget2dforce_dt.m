function [FX,FY] = gpeget2dforce_dt(dirarg,startno,stride,endno,speed,nx,ny)
total = [];
pos = [];
neg = [];
for i=startno:stride:endno
    try 
        [gridx,gridy,psi,potential] = gpeget2dPSI(dirarg,i,speed,nx,ny);
        fprintf('read %d\n',i);
        j = (i+(stride-startno))/stride;
        [fxt,fyt] = gpeget2dforce(psi,gridx,gridy,potential);
        FX(j) = fxt;
        FY(j) = fyt;
    catch E
      fprintf('%d not found\n',i);  
    end
end
end