function [totalE,kinE,potE,interE] = gpeget2denergy(gridx,gridy,inpsi,inpot)
dx=gridx(2)-gridx(1);
[FX, FY] = derivative5(inpsi, 'x', 'y');
FX = FX/dx;
FY = FY/dx;
%[FX, FY] = gradient(inpsi,dx,dx);
dpsi = FX.*conj(FX)+FY.*conj(FY);
kin=0.5*dpsi;

pot=inpot.*inpsi.*conj(inpsi);
inter=0.5.*inpsi.*conj(inpsi).*inpsi.*conj(inpsi);

kinE=abs(dx*dx*simpson(simpson(kin)));
potE=abs(dx*dx*simpson(simpson(pot)));
interE=abs(dx*dx*simpson(simpson(inter)));
totalE=abs(kinE+potE+interE);

end 