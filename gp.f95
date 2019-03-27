program gp
    use params
    use output
    CALL init_params
    DT = -EYE*DTSIZE
    call calc_OBJPOT
    call approx
    open (8, FILE = 'utils.dat')
    call runit(ISTEPS,0)
    include 'ic.in'
    DT = DTSIZE
    call runit(DSTEPS,1)
    call add_noise
    call runit(NSTEPS,2)
    close(8)
end PROGRAM gp
