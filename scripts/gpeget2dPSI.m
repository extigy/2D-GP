function [gridx,gridy,psi,pot] = gpeget2dPSI(dirarg,startno,speed,nx,ny)
gridx=0;
gridy=0;
FX = 0;
FY = 0;
dirarg = regexprep(dirarg, '/$', '');

%datalocation = strcat(dirarg, '/%04d.dumpwf.%04d');
datalocation = strcat(dirarg, '/%dplotwf.%04d');
fname = sprintf(datalocation,speed,startno);
densn = fopen(fname);
A = fscanf(densn, '%g %g %g %g %g\n', [5 inf]);
fclose(densn);
gridy = A(2,1:ny+1);
gridx = A(1,1:ny+1:((nx+1)*(ny+1)));

psi = A(3,:) + 1i*A(4,:);
psi = reshape(psi,(ny+1),(nx+1));
pot = reshape(A(5,:),(ny+1),(nx+1));
fclose('all');
end