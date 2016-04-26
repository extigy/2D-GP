function clus_dt = gpeget2dcluster_dt(dirarg,startno,endno,speed,nx,ny)
    stride = 1;
    for i=startno:stride:endno
        [gridx,gridy,dens,phase,potential] = gpeget2dWF(dirarg,i,speed,nx,ny);
        fprintf('read %d\n',i);
        k = (i-startno)/stride + 1;
        [xlocs,ylocs,pol] = gpeget2dvort(dens,phase,gridx,gridy,potential);
        [clusters,~,~] = gpeget2dcluster(xlocs,ylocs,pol);
        [clussize,~] = cellfun(@size,clusters);
        clus_dt(k) = mean(clussize);
    end

end