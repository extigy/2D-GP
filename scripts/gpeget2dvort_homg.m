function [xlocs,ylocs,pol] = gpeget2dvort_homg(dens,ophase,gridx,gridy,potential)
xlocs=[];
ylocs=[];
pol=[];
dims = size(dens);
dspace=(gridx(2)-gridx(1));
velx(dims(1),dims(2)) = 0;
vely(dims(1),dims(2)) = 0;
presort(dims(1),dims(2)) = 0;
postsort(dims(1),dims(2)) = 0;

phase = unwrap(ophase);
for i = 2:dims(1)-1
for j = 2:dims(2)-1
	if (phase(i+1,j)-phase(i-1,j)<-(pi/2.0d0))
		temp1 = phase(i+1,j)-(phase(i-1,j) - pi);
    elseif (phase(i+1,j)-phase(i-1,j)>(pi/2.0d0))
		temp1 = phase(i+1,j)-(phase(i-1,j) + pi);
	else
		temp1 = phase(i+1,j)-phase(i-1,j);
    end
	velx(i,j) = real(temp1)/dspace;
end
end

phase = unwrap(ophase,[],2);

for i = 2:dims(1)-1
for j = 2:dims(2)-1
	if (phase(i,j+1)-phase(i,j-1)<-(pi/2.0d0))
		temp1 = phase(i,j+1)-(phase(i,j-1) - pi);
    elseif (phase(i,j+1)-phase(i,j-1)>(pi/2.0d0))
		temp1 = phase(i,j+1)-(phase(i,j-1) + pi);
	else
		temp1 = phase(i,j+1)-phase(i,j-1);
    end
	vely(i,j) = real(temp1)/dspace;
end
end

for i = 6:1:dims(1)-6
for j = 6:1:dims(2)-6
        presort(i,j)=LINEINTVF(velx,vely,i-3,i+3,j-3,j+3);
        if(potential(i,j)>0.5)
            presort(i,j) = 0;
        end
end
end

%densblob=dens>0.002;
  
    
%phase(dens<0.0005)=0;
%imagesc(sqrt(velx.*velx+vely.*vely))
%imagesc(gridx,gridy,dens)
%imagesc(dens)
%imagesc(presort);
h = fspecial('gaussian', size(presort), 0.2);
presort = imfilter(presort, h);
%imagesc(presort);

negareas = bwlabel(presort>2);
posareas = bwlabel(presort<-2);


for i = 1:max(max(posareas))
    [r,c] = find(posareas== i);
    if(length(r) > 2)
        xlocs = [xlocs,mean(gridx(c))];
        ylocs = [ylocs,mean(gridy(r))];
        pol = [pol,1];
    end
end

for i = 1:max(max(negareas))
    [r,c] = find(negareas== i);
    if(length(r) > 2)
        xlocs = [xlocs,mean(gridx(c))];
        ylocs = [ylocs,mean(gridy(r))];
        pol = [pol,-1];
    end
end

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

