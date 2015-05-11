function [FX,FY] = gpeget2dforce_rms_dv(dirarg,startno,stride,endno,speed,nx,ny)
    datalocation = strcat(dirarg, '%02d');
    for j=1:40
        fname = sprintf(datalocation,j);
        sp =  300 + 17*(j-1);
        [ffx,ffy] = gpeget2dforce_dt(fname,startno,stride,endno,sp,nx,ny);
        FX(j) = mean(ffx)
        FY(j) = mean(ffy)
    end
end
