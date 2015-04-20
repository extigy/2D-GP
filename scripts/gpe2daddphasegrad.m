function [psi] = gpe2daddphasegrad(psi,gridx,gridy,dtheta)
dims = size(real(psi));
DSPACE = gridx(2)-gridx(1);

for k = 1:dims(1)
for j = 1:dims(2)
    psi(k,j) = psi(k,j)*exp(1i*(j*DSPACE*dtheta));
end
end
end