program mu
	use params
	implicit none
	character(len=80) fname
	double precision :: ret
	CALL init_params
	GRID = 0
	TIME = 0.0d0
	LOOPNO = 0
	GAMMAC = 0.1d0
	VOB = 0.0
	call calc_OBJPOT
	call approx
	open (9, FILE = 'mu.dat')
	call findmu()
	close(9)
end PROGRAM mu

subroutine findmu()
	use params
	implicit none
	double precision :: normgs,muold,muold2
	complex*16, dimension(-NX/2:NX/2,-NY/2:NY/2) :: GRIDGS
	do
		DT = -EYE*DTSIZE
		call runit(2000,0,0)
		DT = DTSIZE	
		call calc_norm
		normgs = NORM
		GRIDGS = GRID
		muold = 0
		muold2 = 0
		muold2 = harm_osc_mu
		do
			call iterate(1)
			call calc_norm
			muold = harm_osc_mu
			if(NORM > normgs) then
				harm_osc_mu = harm_osc_mu - 5000.0*abs(NORM - normgs)
			end if
			if(NORM < normgs) then
				harm_osc_mu = harm_osc_mu + 5000.0*abs(NORM - normgs)
			end if		
			write (unit=6,fmt="(f15.8)") harm_osc_mu
			GRID = GRIDGS
			call runit(1,0,0)
			if(abs(harm_osc_mu - muold) < 1d-10) then
				write (unit=6,fmt="(a,f15.8)") '-->', harm_osc_mu
				EXIT
			end if
		end do
		if(abs(harm_osc_mu - muold2) < 1d-9) then
			EXIT
		end if
	end do
end subroutine
