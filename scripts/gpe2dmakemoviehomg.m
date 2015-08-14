function gpe2dmakemoviehomg(dirarg,startno,stride,endno,speed,nx,ny)
    dirarg = regexprep(dirarg, '/$', ''); 
    pngfolder = strcat(dirarg, '/png');
    mkdir(pngfolder);
    for i=startno:stride:endno
        [gridx,gridy,dens,phase,potential] = gpeget2dWF(dirarg,i,speed,nx,ny);
        %[xlocs,ylocs,pol] = gpeget2dvort_homg(dens,phase,gridx,gridy,potential);
        fprintf('read %d\n',i);
        j = i/stride -startno/stride;
        h=figure('visible','off');
        imagesc(gridx(1:nx+1),gridy(1:ny+1),dens);
        %colormap(gray);
        x=linspace(0,1,128);
        r = sqrt(x);
        g = x.^3.0;
        b=sin(x*2*pi);
        b(b<0)=0;
        pm3d7515_bbry=[r;g;b]';
        colormap(pm3d7515_bbry);
        caxis([0.0,1.1])
        axis image;
        axis xy;
        hold on;
        % g = gscatter(xlocs,ylocs,pol,['b','r'],['^','o'],5,'off');
        % if(length(g) > 1)
        %     set(g(1), 'MarkerFaceColor', 'b');
        %     set(g(2), 'MarkerFaceColor', 'r');
        % end
        xlabel('x/\xi', 'FontSize',16);
        ylabel('y/\xi', 'FontSize',16);
        %axis([190 410 -50 50]);
        ax = findobj(h,'type','axes','Tag','');
        %set(ax,'YTick',[-40,-20,0,20,40]);
        set(ax,'FontSize',16)
        filename = strcat(pngfolder, '/p%04d.png');
        finalfname = sprintf(filename,j);
        print (h,'-dpng','-r150',finalfname);
        %close(h);
    end
end
