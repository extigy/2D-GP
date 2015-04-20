function [psi] = gpe2daddvortex(psi,gridx,gridy,xloc,yloc,circ)
dims = size(real(psi));
DSPACE = gridx(2)-gridx(1);

for i = 1:dims(1)
for j = 1:dims(2)
    xx = ((j-((dims(1)+1)/2))*DSPACE);
	yy = ((i-((dims(2)+1)/2))*DSPACE);
	rs=(xx-xloc)^2.0d0 + (yy-yloc)^2.0d0;
	phse(i,j) = exp(circ*1i*(atan2(yy-yloc,xx-xloc)));
	R(i,j) = sqrt(rs*(0.3437+0.0286*rs)/(1+0.3333*rs+0.0286*rs*rs)); 
end
end
	psi = psi.*R.*phse;
end