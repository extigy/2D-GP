function [xlocs,ylocs,pol] = gpeget2dvort(psi,gridx,gridy,pot,shouldplot)
    disp(['Grid spacing is ',num2str(gridx(2)-gridx(1)),'.']);

    th = 0.2; %threshold value. This can be freely tweaked as a way to control "sensitivity".
    disp(['Using a threshold value of ',num2str(th),'.']);
    
    xlocs=[];
    ylocs=[];
    pol=[];
    
    disp('Getting pseudo-vorticity.');
    pv = gpe2dpseudovorticity(gridx,gridy,psi);
    pv(pot>0.1) = 0.0; %ignore areas where potential is high[x
    
    dspace=(gridx(2)-gridx(1));
    disp('Splitting by polarity.');
    negareas = bwlabel(pv<-th);
    posareas = bwlabel(pv>th);
    
    disp('Finding +ve vortices...');
    for i = 1:max(max(posareas))
        [r,c] = find(posareas== i);
        if(length(r) > 1)
            xlocs = [xlocs,mean(gridx(c))];
            ylocs = [ylocs,mean(gridy(r))];
            pol = [pol,1];
        end
    end
    
    disp('Finding -ve vortices...');
    for i = 1:max(max(negareas))
        [r,c] = find(negareas== i);
        if(length(r) > 1)
            xlocs = [xlocs,mean(gridx(c))];
            ylocs = [ylocs,mean(gridy(r))];
            pol = [pol,-1];
        end
    end
    if(shouldplot == 1)
        disp('Plotting...');
        figure();
        imagesc(gridx,gridy,pv);
        h=figure();
        imagesc(gridx,gridy,abs(psi).^2);
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