function [totalE,kinE,potE,interE] = gpeget2denergy_dt(dirarg,startno,stride,endno)
totalE = [];
kinE = [];
potE = [];
for i=startno:stride:endno
    [gridx,gridy,psi,potential] = gpeget2dPSI(dirarg,i);
    fprintf('read %d\n',i);
    [tE,kE,pE,lE] = gpeget2denergy(gridx,gridy,psi,potential);
    j = (i+(stride-startno))/stride;
    totalE(j) = tE;
    kinE(j) = kE;
    potE(j) = pE;
    interE(j) = lE;
end
fclose('all');
end