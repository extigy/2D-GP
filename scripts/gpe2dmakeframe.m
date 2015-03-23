function gpe2dmakeframe(dirarg,startno,speed,nx,ny)
    dirarg = regexprep(dirarg, '/$', '');
    pngfolder = strcat(dirarg, '/png');
    mkdir(pngfolder);
    i=startno;
    [gridx,gridy,dens,phase,potential] = gpeget2dWF(dirarg,i,speed,nx,ny);
    [xlocs,ylocs,pol] = gpeget2dvort(dens,phase,gridx,gridy,potential);
    fprintf('read %d\n',i);
    imagesc(gridx,gridy,dens)
    colormap 'gray'
    axis image
    axis xy
    hold on
    g = gscatter(xlocs,ylocs,pol,['b','r'],['^','o'],5,'off');
    if(length(g)==1)
        set(g(1), 'MarkerFaceColor', 'r')
        set(g(1),'Marker','o');
        set(g(1),'MarkerEdgeColor','none');
    end
    if(length(g)>1)
        set(g(1),'MarkerEdgeColor','none');
        set(g(1), 'MarkerFaceColor', 'b')
        set(g(2),'MarkerEdgeColor','none');
        set(g(2), 'MarkerFaceColor', 'r')
    end
    axis([-20 20 -20 20]);
    axis off
    %curtime = num2str(-0.4+roundn((i*100*0.0005)/(15*2*pi),-1));
    %text(0,-18,strcat('\color{white}',curtime,' s'),'FontSize',16,'HorizontalAlignment','center')
    
    curtime = num2str(-430+roundn((i*1000*0.0005)/(15*2*pi)*1000,0));
    text(0,-18,strcat('\color{white}',curtime,' ms'),'FontSize',16,'HorizontalAlignment','center')
    
    filename = strcat(pngfolder, '/p%04d.png');
    finalfname = sprintf(filename,round(i));
    %print (h,'-dpng','-r150',finalfname);
    %close(h);
    hold off
end
