subroutine runit(steps,rt,plot)
	use params
	use output
	implicit none
	integer :: steps,rt,i,ii,ji,plot,isvorts,j,k,l,m,n
	double precision :: a,b,c,d,f
	if (plot == 1) then
		call dump_wavefunction_ncdf(0) !dump initial condition
	end if
	do i = STARTI, steps	
		call iterate(rt)

		if(rt == 2) then !it's vortex imprinting
			if(HOLDICV == 1 .and. modulo(i,10) == 9) then
				GRID = sqrt(GRID*conjg(GRID))*IMPOSEDPHASE
			end if
		end if
		TIME = TIME + dble(DT)
		IMAGTIME = IMAGTIME + imag(DT)
		if (modulo(i,dumputil) == 0) then
				call calc_misc
		end if
		if (modulo(i,dumputil) == 0 .and. i > 0) then
			open (10, FILE = "STATUS")
			if (rt == 1) then
					write (unit=10,fmt="(a,f6.2,a)")&
						"Simulating:     ",(dble(i)/dble(steps))*100.0d0,"%"
				else if (rt == 2) then
					write (unit=10,fmt="(a,f6.2,a)")&
						"Imaginary Time Propogation (after ic.in):     ",(dble(i)/dble(steps))*100.0d0,"%"
				else
					write (unit=10,fmt="(a,f6.2,a)")&
						"Imaginary Time Propogation:   " ,(dble(i)/dble(steps))*100.0d0,"%"
			end if
			close(10)
		end if

		if (modulo(i,vortexKillFreq) == 0) then
			if (plot == 1) then
				call kill_vortex_far()
			end if
		end if

		if (modulo(i,dumpwf) == 0) then
			if (plot == 1) then
				call dump_wavefunction_ncdf(i)
			end if
		end if

		if(potRep .eq. 1 .and. rt == 1) then
			call calc_OBJPOT
		end if

		if(i .eq. killgamma) then
			write(6,*) "Killing gamma on iteration", i
			GAMMAC = 0.0d0
		end if
	end do
end subroutine

!Runge-Kutta 4th Order!
subroutine iterate (rt)
	use params

	implicit none
	integer :: rt,BC
	complex*16, dimension(-NX/2:NX/2,-NY/2:NY/2) :: k1,k2,k3,k4
	double precision :: energy

	call rhs(GRID, k1,rt)
	call rhs(GRID + 0.5d0*DT*k1,k2,rt)
	call rhs(GRID + 0.5d0*DT*k2,k3,rt)
	call rhs(GRID + DT*k3,k4,rt)
	GRID = GRID + DT*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4)/6.0d0

	if ((rt .eq. 0) .or. (rt .eq. 2) .or. (rtNorm .and. (rt .eq. 1))) then
		if (RHSType .eq. 0) then
			call calc_norm
			GRID = GRID/sqrt(NORM)
 			!GRID = GRID*sqrt(259.1797d0)
			GRID = GRID*sqrt(DSPACE*DSPACE*NX*NY)
		end if
		if (RHSType .eq. 1) then
			call calc_norm
			GRID = GRID/sqrt(NORM)
		end if
	end if
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!

!Homogeneous Dimensionless - With Ramping!

subroutine rhs (gt, kk,rt)
	use params
	implicit none
	integer :: i,j,BC,rt
	double precision :: r
	complex*16, dimension(-NX/2:NX/2,-NY/2:NY/2) :: gt, kk
	kk=0
	if(RHSType .eq. 0) then
		!$OMP PARALLEL DO
		do i = -NX/2,NX/2
			do j = -NY/2,NY/2
				if(abs(i).eq.NX/2 .and. BCX.eq.2) gt(i,j) = 0.0d0
				if(abs(j).eq.NY/2 .and. BCY.eq.2) gt(i,j) = 0.0d0
				kk(i,j) = 0.5d0*(4.0d0*gt(BC(i,0),BC(j,1))- gt(BC(i,0),BC(j+1,1))&
						-gt(BC(i,0),BC(j-1,1))-gt(BC(i+1,0),BC(j,1))&
						-gt(BC(i-1,0),BC(j,1)))/(DSPACE**2.0d0)&	!laplacian
						+gt(i,j)*gt(i,j)*CONJG(gt(i,j))&	!Nonlinear
			 			-gt(i,j)&	!Chemical Potential
						+OBJPOT(i,j)*gt(i,j)&	!Obstacle potential
						+tanh(TIME/200.0d0)*VOB*EYE*(gt(BC(i+1,0),BC(j,1))&
						-gt(BC(i-1,0),BC(j,1)))/(2.0d0*DSPACE)	!Moving frame
			end do
		end do
		!$OMP END PARALLEL DO
	end if
!Harmonic Dimensionless
	if(RHSType .eq. 1) then
		!$OMP PARALLEL DO
		do i = -NX/2,NX/2
			do j = -NY/2,NY/2
				if(abs(i).eq.NX/2 .and. BCX.eq.2) gt(i,j) = 0.0d0
				if(abs(j).eq.NY/2 .and. BCY.eq.2) gt(i,j) = 0.0d0
				kk(i,j) = 0.5d0*(4.0d0*gt(BC(i,0),BC(j,1))- gt(BC(i,0),BC(j+1,1))&
						-gt(BC(i,0),BC(j-1,1))-gt(BC(i+1,0),BC(j,1))&
						-gt(BC(i-1,0),BC(j,1)))/(DSPACE**2.0d0)&	!laplacian
						+harm_osc_C*gt(i,j)*gt(i,j)*CONJG(gt(i,j))&	!Nonlinear
			 			+OBJPOT(i,j)*gt(i,j)&	!potential
						- harm_osc_mu*gt(i,j)&	!Chemical Potential
						-ROM*EYE*(i*DSPACE)*(gt(BC(i,0),BC(j+1,1))&
						-gt(BC(i,0),BC(j-1,1)))/(2.0d0*DSPACE)&
						+ROM*EYE*(j*DSPACE)*(gt(BC(i+1,0),BC(j,1))&
						-gt(BC(i-1,0),BC(j,1)))/(2.0d0*DSPACE) !Rotating frame
			end do
		end do
		!$OMP END PARALLEL DO		
	end if
	
	if (rt == 1) then ! do we need to do damping?
		if ((dampedX .eqv. .false.) .and. (dampedY .eqv. .false.) .and. (dampedR .eqv. .false.)) then
			kk = kk/(EYE-GAMMAC)
		else if (dampedR .eqv. .false.) then
			!regional damping
			!$OMP PARALLEL DO
			do i = -NX/2,NX/2
				do j = -NY/2,NY/2
					if((dampedX .and. abs((i*DSPACE)) > dampedXDist) .and.&
			 		   (dampedY .and. abs((j*DSPACE)) > dampedYDist)) then
						kk(i,j) = kk(i,j)/(EYE-GAMMAC-&
								  max(0.0d0,min(dampedgamma,dampedgamma*0.5*(tanh((abs((i*DSPACE))-dampedXDist-3.0)/2.0)+1)+&
								  dampedgamma*0.5*(tanh((abs((j*DSPACE))-dampedYDist-3.0)/2.0)+1))))
					else if(dampedX .and. abs((i*DSPACE)) > dampedXDist) then
						kk(i,j) = kk(i,j)/(EYE-GAMMAC-&
								  dampedgamma*0.5*(tanh((abs((i*DSPACE))-dampedXDist-3.0)/2.0)+1))
					else if(dampedY .and. abs((j*DSPACE)) > dampedYDist) then
						kk(i,j) = kk(i,j)/(EYE-GAMMAC-&
								  dampedgamma*0.5*(tanh((abs((j*DSPACE))-dampedYDist-3.0)/2.0)+1))
					else 
						kk(i,j) = kk(i,j)/(EYE-GAMMAC)
					end if
				end do
			end do
			!$OMP END PARALLEL DO
		else
			!regional damping
			do i = -NX/2,NX/2
				do j = -NY/2,NY/2
					r = SQRT((i*DSPACE)**2.0d0+(j*DSPACE)**2.0d0)
					if (r > dampedRDist) then
						kk(i,j) = kk(i,j)/(EYE-GAMMAC-0.5d0*dampedGamma*(1.0d0+tanh(pi*(-1.0d0+2.0d0*(r-dampedRDist)/dampedRWidth))))
					else 
						kk(i,j) = kk(i,j)/(EYE-GAMMAC)
					end if
				end do
			end do

		end if
	else
		kk = kk/EYE
	end if
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function BC(s,n)
	use params
	implicit none
	integer :: s,n
	BC=s
	select case (n)
    	case (0)
    		if(s.eq.NX/2+1  .and. BCX.eq.0)BC=NX/2
		if(s.eq.NX/2+1  .and. BCX.eq.1)BC=-NX/2
		if(s.eq.NX/2+1  .and. BCX.eq.2)BC=NX/2
    		if(s.eq.-NX/2-1 .and. BCX.eq.0)BC=-NX/2
		if(s.eq.-NX/2-1 .and. BCX.eq.1)BC=NX/2
		if(s.eq.-NX/2-1 .and. BCX.eq.2)BC=-NX/2
    	case (1)
    		if(s.eq.NY/2+1  .and. BCY.eq.0)BC=NY/2
		if(s.eq.NY/2+1  .and. BCY.eq.1)BC=-NY/2
		if(s.eq.NY/2+1  .and. BCY.eq.2)BC=NY/2
    		if(s.eq.-NY/2-1 .and. BCY.eq.0)BC=-NY/2
		if(s.eq.-NY/2-1 .and. BCY.eq.1)BC=NY/2
		if(s.eq.-NY/2-1 .and. BCY.eq.2)BC=-NY/2	
	end select
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!
