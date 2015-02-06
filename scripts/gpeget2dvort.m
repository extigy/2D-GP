function [xlocs,ylocs,pol] = gpeget2dvort(dens,ophase,gridx,gridy,potential)
xlocs=[];
ylocs=[];
pol=[];
TFradius = 10;
TFpot = (TFradius*TFradius)/2;
[comx,comy] = gpegetcenterofmass(dens,gridx,gridy);

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

for i = 6:3:dims(1)-6
for j = 6:3:dims(2)-6
        if(sqrt((gridx(j)-comx).^2+(gridy(i)-comy).^2)<TFradius)
            presort(i,j)=LINEINTVF(velx,vely,i,i+5,j,j+5);
            if(sqrt(gridx(j).^2+gridy(i).^2)<5 && potential(i,j)>40)
                presort(i,j) = 0;
            end
        end
end
end



h = fspecial('gaussian', size(presort), 0.5);
presort = imfilter(presort, h);
%imagesc(presort);

negareas = bwlabel(presort>0.5);
posareas = bwlabel(presort<-0.5);

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

