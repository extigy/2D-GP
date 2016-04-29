function [gridx,gridy,dens,phase,potential] = gpeget2dWF(dirarg,startno,speed,nx,ny)
gridx=0;
gridy=0;
FX = 0;
FY = 0;
dirarg = regexprep(dirarg, '/$', '');

datalocation = strcat(dirarg, '/%dplotwf.%04d');
fname = sprintf(datalocation,speed,startno);
densn = fopen(fname);
A = fscanf(densn, '%f %f %f %f %f\n', [5 inf]);
fclose(densn);
gridy = A(2,1:ny+1);
gridx = A(1,1:ny+1:((nx+1)*(ny+1)));

psi = A(3,:) + 1i*A(4,:);
psi = reshape(psi,(ny+1),(nx+1));
dens = psi.*conj(psi);
phase = atan2(imag(psi),real(psi));
%phase = 2*atan(imag(psi)./sqrt(real(psi).*real(psi) + imag(psi).*imag(psi))+real(psi));
potential = reshape(A(5,:),(ny+1),(nx+1));

fclose('all');
end
