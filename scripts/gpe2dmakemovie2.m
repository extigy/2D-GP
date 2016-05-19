function gpe2dmakemovie2(dirarg,startno,stride,endno,speed,nx,ny)
    dirarg = regexprep(dirarg, '/$', '');
    pngfolder = strcat(dirarg, '/png');
    mkdir(pngfolder);
    for i=startno:stride:endno
        [gridx,gridy,dens,phase,potential] = gpeget2dWF(dirarg,i,speed,nx,ny);
        fprintf('read %d\n',i);
        j = i/stride;
        h=figure('visible','off');
        imagesc(gridx,gridy,dens);
        caxis([0 1.1*max(dens(:))])
        axis image
        axis xy
        axis off
        colormap(gray)
        filename = strcat(pngfolder, '/d/%04d.png');
        finalfname = sprintf(filename,round(j));
        %if(i > 87)
        %if(i > 146)
            %text(-40,-40,strcat('\gamma = 0.0'),'FontSize',16,'HorizontalAlignment','center')
            %text(-80,-90,strcat('\gamma = 0.00'),'FontSize',16,'HorizontalAlignment','center')
        %else
            %text(-40,-40,strcat('\gamma = 0.1'),'FontSize',16,'HorizontalAlignment','center')
            %text(-80,-90,strcat('\gamma = 0.1'),'FontSize',16,'HorizontalAlignment','center')
        %end
        print (h,'-dpng','-r150',finalfname);
        close(h);
        
        
%         h=figure('visible','off');
%         imagesc(gridx,gridy,dens);
%         caxis([0.95*min(dens(:)) 1.05*max(dens(:))])
%         axis image
%         axis xy
%         axis off
%         colormap(gray)
%         filename = strcat(pngfolder, '/d2/%04d.png');
%         finalfname = sprintf(filename,round(j));
%         print (h,'-dpng','-r150',finalfname);
%         close(h);
        
        
        h=figure('visible','off');
        imagesc(gridx,gridy,phase);
        axis image
        axis xy
        axis off
        colormap(hsv)
        filename = strcat(pngfolder, '/p/%04d.png');
        finalfname = sprintf(filename,round(j));
        print (h,'-dpng','-r150',finalfname);
        close(h);
        
%         h=figure('visible','off');
%         subplot(2,2,[3 4])
%         mdens = dens;
%         mindens(j+1) = min(mdens(:));
%         plot((1:(j+1)).*stride,mindens,'LineWidth',1);
%         axis([0 10000 0 1]);
%         xlabel('Time');
%         ylabel('Density minimum');
%         filename = strcat(pngfolder, '/m/%04d.png');
%         finalfname = sprintf(filename,round(j));
%         print (h,'-dpng','-r150',finalfname);
%         close(h);
    end
end
