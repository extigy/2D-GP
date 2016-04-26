function [Lpos,Lneg] = gpegetL(dens,phase,gridx,gridy,potential)
[xlocs,ylocs,pol] = gpeget2dvort(dens,phase,gridx,gridy,potential);

stride = 1;
xlocsp = xlocs(pol > 0);
ylocsp = ylocs(pol > 0);

xlocsn = xlocs(pol < 0);
ylocsn = ylocs(pol < 0);

Kpos = RipleysK([xlocsp',ylocsp'],0:stride:200,[min(gridx) max(gridx) min(gridy) max(gridy)],0);
Kneg = RipleysK([xlocsn',ylocsn'],0:stride:200,[min(gridx) max(gridx) min(gridy) max(gridy)],0);
Lpos = sqrt(Kpos'/pi)-(0:stride:200);
Lneg = sqrt(Kneg'/pi)-(0:stride:200);
%plot((0:stride:200),Lpos,'r')
%hold on
%plot((0:stride:200),Lneg,'b--')
end
