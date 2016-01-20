program gp
	use params
	implicit none
	character(len=80) fname
	integer :: ji,ii
	double precision :: ret,tmptime
	CALL init_params

	do LOOPNO = VOBS,VOBE,VOBST
		!Initialise
		GRID = 0
		TIME = 0.0d0
		IMAGTIME = 0.0d0
		VOB = DBLE(LOOPNO)/VOBSCALE
		DT = -EYE*DTSIZE
		call calc_OBJPOT
		call approx
		write(fname, '(a,i0)') 'utils.',LOOPNO
		open (8, FILE = fname)
		!!!!!!!!!!!!
		call runit(ISTEPS,0,PLOTIT)
		include 'ic.in'
		call runit(VSTEPS,2,0)
		DT = DTSIZE
		call add_noise
		call runit(NSTEPS,1,1)
		close(8)
	end do
end PROGRAM gp

