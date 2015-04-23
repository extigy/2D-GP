[gridx,gridy,psi,potential] = gpeget2dPSI('/data/redrocks/osc/osc10',0,0,512,512);
for i=0:1000,
[npsi] = gpe2daddphasegrad(psi,gridx,gridy,-0.5);
randy = rand(1)*30 - 15;
randx = rand(1)*30 - 15;
[npsi] = gpe2daddvortex(npsi,gridx,gridy,randx,randy,-1);
[FX0,FY0] = gpeget2dforce(npsi,gridx,gridy,potential);
X0(i+1) = randx;
Y0(i+1) = randy;
FFX0(i+1) = FX0;
FFY0(i+1) = FY0;
end
F = TriScatteredInterp(X0',Y0',FFX0');
ti = -15:.1:15;
[qx,qy] = meshgrid(ti,ti);
qz = F(qx,qy);
mesh(qx,qy,qz);
h=pcolor(qx,qy,qz);
set(h,'EdgeColor','none')