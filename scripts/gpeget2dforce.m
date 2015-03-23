function [FX,FY] = gpeget2dforce(dens,phase,gridx,gridy,potential)
psi = dens + exp(1i*phase);

dx = gridx(5)-gridx(4);
dy = gridy(5)-gridy(4);
[FX,FY] = gradient(potential,dx);
FX = psi.*conj(psi).*FX;
FY = psi.*conj(psi).*FY;
FX = -trapz(trapz(FX))*dx*dx;
FY = -trapz(trapz(FY))*dy*dy;
end