function gpe2dmakemoviecomb(dirarg,startno,stride,endno,speed,nx,ny)
    dirarg = regexprep(dirarg, '/$', ''); 
    pngfolder = strcat(dirarg, '/png');
    mkdir(pngfolder);
    pngfolder = strcat(dirarg, '/png/comb');
    mkdir(pngfolder);
    for i=startno:stride:endno
        [gridx,gridy,dens,phase,potential] = gpeget2dWF(dirarg,i,speed,nx,ny);
        fprintf('read %d\n',i);
        j = i/stride -startno/stride;
        h=figure('visible','off');
        imagesc_comb(gridx,gridy,dens,phase);
        cmapr = hsv;
        colormap(cmapr);
        axis image;
        axis xy;
        axis off;
        hold on;
        filename = strcat(pngfolder, '/%04d.png');
        finalfname = sprintf(filename,j);
        print (h,'-dpng','-r150',finalfname);
        close(h);
    end
end
