function imagesc_comb(gridx,gridy,dens,G)
%C = flipud(hsv(256));
C = hsv;
L = size(C,1);
Gs = round(interp1(linspace(min(G(:)),max(G(:)),L),1:L,G));
H = reshape(C(Gs,:),[size(Gs) 3]);
mask = mat2gray(dens);
dims = size(Gs);
for i = 1:dims(1)
for j = 1:dims(2)
   H(i,j,:) = H(i,j,:).*mask(i,j);
end
end

image(gridx,gridy,H)
end