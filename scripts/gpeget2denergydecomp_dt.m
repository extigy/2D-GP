function [Et,Ec, Ei] = gpeget2denergydecomp_dt(dirarg,startno,stride,endno)
Et = [];
Ec = [];
Ei = [];
for i=startno:stride:endno
    [gridx,gridy,psi] = gpeget2dPSI(dirarg,i);
    fprintf('read %d\n',i);
    [Ekinsq_t,Ekinsq_c, Ekinsq_i] = gpeget2denergydecomp(gridx,psi);
    j = (i+(stride-startno))/stride;
    Et(j) = trapz(gridy,trapz(gridx,Ekinsq_t));
    Ec(j) = trapz(gridy,trapz(gridx,Ekinsq_c));
    Ei(j) = trapz(gridy,trapz(gridx,Ekinsq_i));
end
fclose('all');
end