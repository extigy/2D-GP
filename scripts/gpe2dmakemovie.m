function gpe2dmakemovie(dirarg,startno,stride,endno)
    dirarg = regexprep(dirarg, '/$', '');
    pngfolder = strcat(dirarg, '/png');
    mkdir(pngfolder);
    for i=startno:stride:endno
        [gridx,gridy,psi,potential] = gpeget2dPSI(dirarg,i);
        [xlocs,ylocs,pol] = gpeget2dvort(psi,gridx,gridy,'potential',potential);
        dens = abs(psi).^2;
        fprintf('read %d\n',i);
        j = i/stride;
        h=figure('visible','off');
        imagesc(gridx,gridy,dens)
        colormap 'gray'
        axis image
        axis xy
        hold on
        g = gscatter(xlocs,ylocs,pol,['b','r'],['^','o'],5,'off');
        if(length(g)==1 && pol(1)==1)
            set(g(1), 'MarkerFaceColor', 'r')
            set(g(1),'Marker','o');
            set(g(1),'MarkerEdgeColor','none');
        elseif(length(g)==1 && pol(1)==-1)
            set(g(1), 'MarkerFaceColor', 'b')
            set(g(1),'Marker','^');
            set(g(1),'MarkerEdgeColor','none');
        else
            set(g(1),'MarkerEdgeColor','none');
            set(g(1), 'MarkerFaceColor', 'b')
            set(g(2),'MarkerEdgeColor','none');
            set(g(2), 'MarkerFaceColor', 'r')
        end
        axis off
        filename = strcat(pngfolder, '/d%06d.png');
        finalfname = sprintf(filename,round(j));
        print (h,'-dpng','-r150',finalfname);
        close(h);
    end
end
