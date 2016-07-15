function [E,P] = gpeget2denergymomentum(gridx,gridy,inpsi)
dx=gridx(2)-gridx(1);
[FX, FY] = gradient(inpsi,dx,dx);
[FXS, FYS] = gradient(conj(inpsi),dx,dx);


modDpsisq = FX.*conj(FX)+FY.*conj(FY);
modpsisq = inpsi.*conj(inpsi);

kin = 0.5*modDpsisq;
inter = 0.5*(modpsisq-1).^2;
kinE=dx*dx*simpson(simpson(kin));
interE=dx*dx*simpson(simpson(inter));
E=kinE+interE;

momx = ((conj(inpsi)-1).*FX - (inpsi-1).*FXS);
momy = ((conj(inpsi)-1).*FY - (inpsi-1).*FYS);
momxP = 0.5./1i.*dx*dx*simpson(simpson(momx));
momyP = 0.5./1i.*dx*dx*simpson(simpson(momy));

P = sqrt(momxP.*conj(momxP)+momyP.*conj(momyP));
end 