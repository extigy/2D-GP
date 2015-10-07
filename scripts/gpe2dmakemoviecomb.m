function gpe2dmakemoviehomg(dirarg,startno,stride,endno,speed,nx,ny)
    dirarg = regexprep(dirarg, '/$', ''); 
    pngfolder = strcat(dirarg, '/png');
    mkdir(pngfolder);
    for i=startno:stride:endno
        [gridx,gridy,dens,phase,potential] = gpeget2dWF(dirarg,i,speed,nx,ny);
        fprintf('read %d\n',i);
        j = i/stride -startno/stride;
        h=figure('visible','off');
        imagesc_comb(gridx,gridy,dens,phase);
        axis image;
        axis xy;
        hold on;
        xlabel('x/\xi', 'FontSize',16);
        ylabel('y/\xi', 'FontSize',16);
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
