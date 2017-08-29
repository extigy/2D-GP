function clus_dt = gpeget2dcluster_dt(dirarg,startno,endno)
    stride = 1;
    for i=startno:stride:endno
        [gridx,gridy,psi,potential] = gpeget2dPSI(dirarg,i);
        fprintf('read %d\n',i);
        k = (i-startno)/stride + 1;
        [xlocs,ylocs,pol] = gpeget2dvort(psi,gridx,gridy,'potential',potential);
        [clusters,~,~] = gpeget2dcluster(xlocs,ylocs,pol);
        [clussize,~] = cellfun(@size,clusters);
        clus_dt(k) = mean(clussize);
    end

end