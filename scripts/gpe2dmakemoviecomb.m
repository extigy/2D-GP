function gpe2dmakemoviecomb(dirarg,startno,stride,endno,speed,nx,ny)
    dirarg = regexprep(dirarg, '/$', ''); 
    pngfolder = strcat(dirarg, '/png');
    mkdir(pngfolder);
    for i=startno:stride:endno
        [gridx,gridy,dens,phase,potential] = gpeget2dWF(dirarg,i,speed,nx,ny);
        fprintf('read %d\n',i);
        j = i/stride -startno/stride;
        h=figure('visible','off');
        imagesc_comb(gridx,gridy,dens,phase);
        cmapr = cmap('C1');
        colormap(cmapr);
        hcb = colorbar(); 
        set(hcb,'YTick',[5.77,46.51,87.26,128,168.74,209.5,250.23]);
        set(hcb,'YTickLabel',[-3,-2,-1,0,1,2,3]);
        axis image;
        axis xy;
        hold on;
        xlabel('x/l', 'FontSize',16);
        ylabel('y/l', 'FontSize',16);
        %axis([190 410 -50 50]);
        ax = findobj(h,'type','axes','Tag','');
        %set(ax,'YTick',[-40,-20,0,20,40]);
       % set(ax,'FontSize',16)
        filename = strcat(pngfolder, '/p%04d.png');
        finalfname = sprintf(filename,j);
        print (h,'-dpng','-r150',finalfname);
        close(h);
    end
end
