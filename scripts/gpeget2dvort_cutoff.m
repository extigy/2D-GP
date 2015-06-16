function [xlocs,ylocs,pol] = gpeget2dvort_cutoff(filename,Nx,Ny)
dirarg = regexprep(filename, '/$', '');
denslocation = strcat(dirarg, '/densities/iBC_2d050.dat');
phaselocation = strcat(dirarg, '/phase/phase001_001_050.dat');

density = importdata(denslocation);
ophase = importdata(phaselocation);
gridx = density(1:Nx+1,1);
gridy = density(1:Nx+1:(Nx+1)*(Ny+1),2);
ophase = reshape(ophase(:,3),Nx+1,Ny+1);
dens = reshape(density(:,3),Nx+1,Ny+1);
xlocs=[];
ylocs=[];
pol=[];

disp(['Grid spacing is ',num2str(gridx(2)-gridx(1)),'.']);
gf = 0.7; %Gaussian filter width. Should really be about the size of a vortex core in your units

disp(['Using gaussian filter of width ',num2str(gf),'.']);
li = 4; %Line integral size. Should be around the number of grid points in a vortex core

disp(['Doing line integrals over ',num2str(li),' points.']);
th = 2.2; %threshold value. I have set it to just over pi/2.
          %This can be freely tweaked as a way to control "sensitivity".
disp(['Using a threshold value of ',num2str(th),'.']);


dims = size(dens);
dspace=(gridx(2)-gridx(1));

psi = sqrt(dens).*exp(1i*ophase);
v1=fftn(psi);
v1=v1/(sqrt(Nx*Ny));
v1=ifftshift(v1);
    
kx=(-(Nx-1)/2:(Nx-1)/2);
ky=(-(Ny-1)/2:(Ny-1)/2);
for j=1:Nx
    for k=1:Ny
        k2=kx(j)^2+ky(k)^2;
        v1(j,k)=v1(j,k)*max(1-k2/30^2,0);
    end  
end

v1=ifftshift(v1);
v1=v1*(sqrt(Nx*Ny));
u1=ifftn(v1);
nphase = atan2(imag(u1),real(u1));
ndens=abs(u1).^2;
spatial_av1=mean(mean(mean(ndens)));

phase = unwrap(nphase);
dims = size(nphase);
for i = 2:dims(1)-1
for j = 2:dims(2)-1
    if (phase(i+1,j)-phase(i-1,j)<-(pi/2.0d0))
        temp1 = phase(i+1,j)-(phase(i-1,j) - pi);
    elseif (phase(i+1,j)-phase(i-1,j)>(pi/2.0d0))
        temp1 = phase(i+1,j)-(phase(i-1,j) + pi);
    else
        temp1 = phase(i+1,j)-phase(i-1,j);
    end
    mgx(i,j) = gridx(j);
    mgy(i,j) = gridy(i);
    vely(i,j) = real(temp1)/(2*dspace);
end
    vely(i,dims(2)) = vely(i,dims(2)-1);
    vely(i,1) = vely(i,2);
end
for j = 1:dims(2)
    vely(dims(1),j) = vely(dims(1)-1,j);
    vely(1,j) = vely(2,j);
end


phase = unwrap(nphase,[],2);
for i = 2:dims(1)-1
for j = 2:dims(2)-1
    if (phase(i,j+1)-phase(i,j-1)<-(pi/2.0d0))
        temp1 = phase(i,j+1)-(phase(i,j-1) - pi);
    elseif (phase(i,j+1)-phase(i,j-1)>(pi/2.0d0))
        temp1 = phase(i,j+1)-(phase(i,j-1) + pi);
    else
        temp1 = phase(i,j+1)-phase(i,j-1);
    end
    velx(i,j) = real(temp1)/(2*dspace);
end
    velx(i,dims(2)) = velx(i,dims(2)-1);
    velx(i,1) = velx(i,2);
end
for j = 1:dims(2)
    velx(dims(1),j) = velx(dims(1)-1,j);
    velx(1,j) = velx(2,j);
end

for i = li:1:dims(1)-li
for j = li:1:dims(2)-li
      presort(i,j)=LINEINTVF(velx,vely,i-li/2,i+li/2,j-li/2,j+li/2);
end
end
imagesc(ophase)
figure
h = fspecial('gaussian', size(nphase), 0.6);
nphase = imfilter(nphase, h);
imagesc(nphase)

h = fspecial('gaussian', size(presort), gf);
presort = imfilter(presort, h);
negareas = bwlabel(presort<-th);
posareas = bwlabel(presort>th);

for i = 1:max(max(posareas))
    [r,c] = find(posareas== i);
    if(length(r) > 1)
        xlocs = [xlocs,mean(gridx(c))];
        ylocs = [ylocs,mean(gridy(r))];
        pol = [pol,1];
    end
end

for i = 1:max(max(negareas))
    [r,c] = find(negareas== i);
    if(length(r) > 1)
        xlocs = [xlocs,mean(gridx(c))];
        ylocs = [ylocs,mean(gridy(r))];
        pol = [pol,-1];
    end
end
h=figure();
imagesc(gridx,gridy,nphase);
colormap(gray);
axis image;
axis xy;
hold on;
g = gscatter(xlocs,ylocs,pol,['b','r'],['^','o'],5,'off');
        if(length(g)==1 && pol(1)==1)
            set(g(1), 'MarkerFaceColor', 'r')
            set(g(1),'Marker','o');
            set(g(1),'MarkerEdgeColor','none');
        end
        if(length(g)==1 && pol(1)==-1)
            set(g(1), 'MarkerFaceColor', 'b')
            set(g(1),'Marker','^');
            set(g(1),'MarkerEdgeColor','none');
        end
        if(length(g)>1)
            set(g(1),'MarkerEdgeColor','none');
            set(g(1), 'MarkerFaceColor', 'b')
            set(g(2),'MarkerEdgeColor','none');
            set(g(2), 'MarkerFaceColor', 'r')
        end
xlabel('x', 'FontSize',16);
ylabel('y', 'FontSize',16);
ax = findobj(h,'type','axes','Tag','');
set(ax,'FontSize',16)


function ret = LINEINTVF(fieldx,fieldy,x,ex,y,ey)
	l1=0.0d0;
	l2=0.0d0;
	l3=0.0d0;
	l4=0.0d0;
	for t = y:ey
		l1 = l1 + dspace*fieldy(x,t);
    end
	for t = x:ex
		l2 = l2 + dspace*fieldx(t,y);
    end
	for t = y:ey
		l3 = l3 + dspace*fieldy(ex,t);
    end
	for t = x:ex
		l4 = l4 + dspace*fieldx(t,ey);
    end
	ret = l2+l3-l4-l1;
end


end

