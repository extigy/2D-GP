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
        caxis([0 1.1])
        axis image
        axis xy
        axis off
        colormap(gray)
        filename = strcat(pngfolder, '/d%04d.png');
        finalfname = sprintf(filename,round(j));
        print (h,'-dpng','-r150',finalfname);
        close(h);
        h=figure('visible','off');
        imagesc(gridx,gridy,phase);
        axis image
        axis xy
        axis off
        colormap(hsv)
        filename = strcat(pngfolder, '/p%04d.png');
        finalfname = sprintf(filename,round(j));
        print (h,'-dpng','-r150',finalfname);
        close(h);
    end
end
