function norm = gpeget2dnorm_dt(dirarg,startno,stride,endno)
norm = [];
for i=startno:stride:endno
    [gridx,gridy,psi,~] = gpeget2dPSI(dirarg,i);
	dens = psi.*conj(psi)
    fprintf('read %d\n',i);
    normt = gpeget2dnorm(gridx,gridy,dens);
    j = (i+(stride-startno))/stride;
    norm(j) = normt;
end
fclose('all');
end