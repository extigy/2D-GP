function [wx,wy] = gpeget2dJRwh(gridx,gridy,dens)
dx=gridx(2)-gridx(1);
dy=gridy(2)-gridy(1);
th = (max(dens(:))+min(dens(:)))/2;
bw = dens<th;
%imagesc(gridx,gridy,bw);
L = bwlabel(bw);
maxl = max(L(:));
for i=1:maxl
    [r,c] = find(L == i);
    wxi(i) = (max(c)-min(c))*dx;
    wyi(i) = (max(r)-min(r))*dy;
end
wx = max(wxi);
wy = max(wyi);

end 