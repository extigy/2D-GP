function [total,pos,neg] = gpeget2dvort_dt(dirarg,startno,stride,endno)
total = [];
pos = [];
neg = [];
    for i=startno:stride:endno
        [gridx,gridy,psi,pot] = gpeget2dPSI(dirarg,i);
        fprintf('read %d\n',i);
        [xlocs,ylocs,pol] = gpeget2dvort(psi,gridx,gridy,'potential',pot);
        j = (i+(stride-startno))/stride;
        total(j) = length(pol);
        neg(j) = -sum(pol-1)/2;
        pos(j) = length(pol) - neg(j);
    end
end