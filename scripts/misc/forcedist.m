[gridx,gridy,psi,potential] = gpeget2dPSI('/data/redrocks/osc/osc10',0,0,512,512);
for i=0:160,
[npsi] = gpe2daddphasegrad(psi,gridx,gridy,-0.5);
[npsi] = gpe2daddvortex(npsi,gridx,gridy,-i/4,i/4,1);
[FX0,FX1,FX2,FX3,FY0,FY1,FY2,FY3] = gpeget2dforcesplit(npsi,gridx,gridy,potential);
FFX0(i+1) = FX0;
FFX1(i+1) = FX1;
FFX2(i+1) = FX2;
FFX3(i+1) = FX3;
FFY0(i+1) = FY0;
FFY1(i+1) = FY1;
FFY2(i+1) = FY2;
FFY3(i+1) = FY3;
% [FX0,FY0] = gpeget2dforce(npsi,gridx,gridy,potential);
% FFX0(i+1) = FX0;
% FFY0(i+1) = FY0;
end