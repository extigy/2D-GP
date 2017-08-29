program gp
	use params
	use output
	implicit none
	character(len=80) fname
	integer :: ji,ii
	double precision :: ret,tmptime
	CALL init_params
	!Initialise
	GRID = 0
	TIME = 0.0d0
	IMAGTIME = 0.0d0
	DT = -EYE*DTSIZE
	call calc_OBJPOT
	call approx
	open (8, FILE = 'utils.dat')
	call runit(ISTEPS,0,PLOTIT)
	include 'ic.in'
	call runit(VSTEPS,2,0)
	DT = DTSIZE
	call add_noise
	call runit(NSTEPS,1,1)
	close(8)
end PROGRAM gp
