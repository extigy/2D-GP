function [pos_dv,neg_dv,total_dv] = gpeget2dvort_homg_nearobj_dv(dirarg,startno,stride,endno,nx,ny)
    datalocation = strcat(dirarg, '%02d');
    for run=1:40
        fname = sprintf(datalocation,run);
        speed = 0;

        for i=startno:stride:endno
            [gridx,gridy,dens,phase,potential] = gpeget2dWF(fname,i,speed,nx,ny);
            fprintf('read %d\n',i);
            j = (i+(stride-startno))/stride;
            [xlocs,ylocs,pol] = gpeget2dvort_homg_nearobj(dens,phase,gridx,gridy,potential);
            total_dv(j,run) = length(pol);
            neg_dv(j,run) = -sum(pol-1)/2;
            pos_dv(j,run) = length(pol) - neg_dv(j,run);
        end
    end
end
