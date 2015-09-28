function [xlocs,ylocs,pol] = gpeget2dvort_homg(dens,ophase,gridx,gridy,potential)
    disp(['Grid spacing is ',num2str(gridx(2)-gridx(1)),'.']);

    gf = 1; %Gaussian filter width. Smooths out noise.
    disp(['Using gaussian filter of width ',num2str(gf),'.']);

    li = 4; %Line integral size. Set around the number of grid points in a vortex core radius
    disp(['Doing line integrals over ',num2str(li),' points.']);

    th = 5.0; %threshold value. This can be freely tweaked as a way to control "sensitivity".
    disp(['Using a threshold value of ',num2str(th),'.']);

    xlocs=[];
    ylocs=[];
    pol=[];

    dims = size(dens);
    dspace=(gridx(2)-gridx(1));
    velx(dims(1),dims(2)) = 0;
    vely(dims(1),dims(2)) = 0;
    presort(dims(1),dims(2)) = 0;

    phase = ophase;
    for i = 2:dims(1)-1
    for j = 2:dims(2)-1
        temp1 = anglediff(phase(i+1,j),phase(i-1,j));
        velx(i,j) = real(temp1)/dspace;
     end
    end

    for i = 2:dims(1)-1
    for j = 2:dims(2)-1
        temp1 = anglediff(phase(i,j+1),phase(i,j-1));
        vely(i,j) = real(temp1)/dspace;
    end
    end

    for i = li:1:dims(1)-li
    for j = li:1:dims(2)-li
          presort(i,j)=LINEINTVF(velx,vely,i-li/2,i+li/2,j-li/2,j+li/2);
          if(potential(i,j)>1.0)
            presort(i,j) = 0;
          end
    end
    end

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
    imagesc(gridx,gridy,dens);
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

end