module output
    use params
    use netcdf
    implicit none
    include 'netcdf.inc'
    integer :: ncdf_id,x_id,y_id,re_id,im_id,pot_id
    integer       :: icount(3)
    integer       :: istarting(3)

Contains

    subroutine make_file(fname)
        implicit none
        integer :: dims(2)
        integer :: r,x_dim_id,y_dim_id
        character(len=80) fname
        r=NF90_create(fname , IOR(NF90_CLOBBER,NF90_64BIT_OFFSET), ncdf_id)
        r=NF90_def_dim(ncdf_id, 'x_dim', NX+1, x_dim_id)
        r=NF90_def_dim(ncdf_id, 'y_dim', NY+1, y_dim_id)

        r=NF90_def_var(ncdf_id, 'x', NF90_DOUBLE,x_dim_id, x_id)
        r=NF90_def_var(ncdf_id, 'y', NF90_DOUBLE,y_dim_id, y_id)

        dims(1) = x_dim_id
        dims(2) = y_dim_id

        r=NF90_def_var(ncdf_id, 'real', NF90_DOUBLE, dims, re_id)
        r=NF90_def_var(ncdf_id, 'imag', NF90_DOUBLE, dims, im_id)
        r=NF90_def_var(ncdf_id, 'pot' , NF90_DOUBLE, dims, pot_id)

        r=NF90_enddef(ncdf_id)
        r=NF90_sync(ncdf_id)

    end subroutine


    subroutine write_wf_file()
        implicit none
        integer :: r,i
        istarting(1) = 1
        istarting(2) = 1

        icount(1) = NX+1
        r=NF90_put_var(ncdf_id, x_id, (/(i*DSPACE,i=-NX/2,NX/2)/),istarting,icount)
        icount(1) = NY+1
        r=NF90_put_var(ncdf_id, y_id, (/(i*DSPACE,i=-NY/2,NY/2)/),istarting,icount)

        icount(1) = NX+1
        icount(2) = NY+1

        r=NF90_put_var(ncdf_id, re_id,    dble(GRID),istarting,icount)
        r=NF90_put_var(ncdf_id, im_id,   aimag(GRID),istarting,icount)
        r=NF90_put_var(ncdf_id, pot_id, OBJPOT,istarting,icount)
    end subroutine

    subroutine close_file()
        implicit none
        integer :: r
        r=NF90_close(ncdf_id)
    end subroutine

    subroutine read_wf_file(fname)
        implicit none
        integer :: r,rwf_ncid,rwf_re_id,rwf_im_id,rwf_pot_id
        double precision, dimension(-NX/2:NX/2,-NY/2:NY/2) :: realgrid,imaggrid,potgrid
        character(len=2048) fname
        r = NF90_open(path = fname, mode = nf90_nowrite, ncid = rwf_ncid)
        r = NF90_inq_varid(rwf_ncid, "real",  rwf_re_id)
        r = NF90_inq_varid(rwf_ncid, "imag",  rwf_im_id)
        r = NF90_inq_varid(rwf_ncid,  "pot", rwf_pot_id)
        r = NF90_get_var(rwf_ncid, rwf_re_id, realgrid)
        r = NF90_get_var(rwf_ncid, rwf_im_id, imaggrid)
        r = NF90_get_var(rwf_ncid, rwf_pot_id, potgrid)
        GRID = realgrid + EYE*imaggrid
        OBJPOT = potgrid
        r = NF90_close(rwf_ncid)
        write(6,*) "Loaded saved grid. Center point is: ", GRID(0,0)
    end subroutine
    
    subroutine dump_wavefunction_ncdf (II,steptype)
        use params
        implicit none
        integer :: II,steptype
        character(len=80) fname
        if(steptype == 0) then
            write(fname, '(a,i0.6,a)') 'imag.',II/dumpwf,'.nc'
        else if (steptype == 1) then
            write(fname, '(a,i0.6,a)') 'smth.',II/dumpwf,'.nc'
        else if (steptype == 2) then
            write(fname, '(a,i0.6,a)') 'psi.',II/dumpwf,'.nc'
        end if
        call make_file(fname)
        call write_wf_file
        call close_file
    end subroutine
end module
