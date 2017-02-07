function [xlocs,ylocs,pol] = gpeget2dvort(psi,gridx,gridy,pot,shouldplot)
    disp(['Grid spacing is ',num2str(gridx(2)-gridx(1)),'.']);
    dims = size(psi);
    dy = gridy(2)-gridy(1);
    gf = 2;
    disp(['Using gaussian filter of width ',num2str(gf),'.']);
    
    h = fspecial('gaussian', dims, gf);
    psiflt = imfilter(psi, h,'circular');
    phaseflt = angle(psiflt);

    xlocs=[];
    ylocs=[];
    pol=[];

    for i = 1:dims(2)
        for j = 1:dims(1)
            t1 = anglediff(phaseflt(mod(i+1-1,dims(2))+1,j),phaseflt(mod(i-1-1,dims(2))+1,j));
            vely(i,j) = real(t1)/(2*dy);
        end
    end

    for i = 1:dims(2)
        for j = 1:dims(1)
            t1 = anglediff(phaseflt(i,mod(j+1-1,dims(1))+1),phaseflt(i,mod(j-1-1,dims(1))+1));
            velx(i,j) = real(t1)/(2*dy);
        end
    end

    for i = 1:dims(2)
        for j = 1:dims(1)
            presort(i,j)=LINEINTVF(vely,velx,dy,i-1,i+1,j-1,j+1,dims(2),dims(1));
            if(pot(i,j)>80.0)
                presort(i,j) = 0;
            end
        end
    end

    negareas = bwlabel(presort<-6.2); %just under 2pi
    posareas = bwlabel(presort>6.2);

    %periodicity. find areas where connected component is over a
    %boundary and choose one
    r = find(posareas(1,1:dims(1))~=0 & posareas(dims(2),1:dims(1))~=0);
    if(~isempty(r))
        posareas(posareas==posareas(1,r(1))) = 0;
    end

    r = find(posareas(1:dims(2),1)~=0 & posareas(1:dims(2),dims(1))~=0);
    if(~isempty(r))
        posareas(posareas==posareas(r(1),1)) = 0;
    end

    r = find(negareas(1,1:dims(1))~=0 & negareas(dims(2),1:dims(1))~=0);
    if(~isempty(r))
        negareas(negareas==negareas(1,r(1))) = 0;
    end

    r = find(negareas(1:dims(2),1)~=0 & negareas(1:dims(2),dims(1))~=0);
    if(~isempty(r))
        negareas(negareas==negareas(r(1),1)) = 0;
    end

    for i = 1:max(max(posareas))
        [r,c] = find(posareas== i);
        if(~isempty(r))
            xlocs = [xlocs,mean(gridx(c))];
            ylocs = [ylocs,mean(gridy(r))];
            pol = [pol,1];
        end
    end
    for i = 1:max(max(negareas))
        [r,c] = find(negareas== i);
        if(~isempty(r))
            xlocs = [xlocs,mean(gridx(c))];
            ylocs = [ylocs,mean(gridy(r))];
            pol = [pol,-1];
        end
    end 
            
            
    if(shouldplot == 1)
        disp('Plotting...');
        figure();
        imagesc(gridx,gridy,presort);
        h=figure();
        imagesc(gridx,gridy,abs(psiflt).^2);
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
    end
    disp('Done!');
end

function ret = LINEINTVF(fieldx,fieldy,dy,x,ex,y,ey,ii,jj)
    l1=0;
    l2=0;
    l3=0;
    l4=0;
    for t = y:ey
        l1 = l1 + dy*fieldy(mod(x-1,ii)+1,mod(t-1,jj)+1);
    end
    for t = x:ex
        l2 = l2 + dy*fieldx(mod(t-1,ii)+1,mod(y-1,jj)+1);
    end
    for t = y:ey
        l3 = l3 + dy*fieldy(mod(ex-1,ii)+1,mod(t-1,jj)+1);
    end
    for t = x:ex
        l4 = l4 + dy*fieldx(mod(t-1,ii)+1,mod(ey-1,jj)+1);
    end
    ret = -l2-l3+l4+l1;
end

function d = anglediff(th1, th2)
    if nargin < 2
        d = th1;
    else
        d = th1 - th2;
    end
    d = mod(d+pi, 2*pi) - pi;
end