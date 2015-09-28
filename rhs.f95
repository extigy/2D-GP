!Runge-Kutta 4th Order!

subroutine iterate (rt)
	use params

	implicit none
	integer :: rt,BC
	complex*16, dimension(-NX/2:NX/2,-NY/2:NY/2) :: k1,k2,k3,k4
	double precision :: energy

	if (rt == 0) then
		!call calc_norm
		!write(6,*) NORM
	end if

	call rhs(GRID, k1,rt)
	call rhs(GRID + 0.5d0*DT*k1,k2,rt)
	call rhs(GRID + 0.5d0*DT*k2,k3,rt)
	call rhs(GRID + DT*k3,k4,rt)
	GRID = GRID + DT*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4)/6.0d0

	if (rt == 0 .or. rtNorm) then
		if (RHSType .eq. 0) then
			call calc_norm
			GRID = GRID/sqrt(NORM)
			GRID = GRID*sqrt(DSPACE*DSPACE*NX*NY)
		end if
		if (RHSType .eq. 1) then
			call calc_norm
			GRID = GRID/sqrt(NORM)
		end if
	end if
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!

!Homogeneous Dimentionless - With Ramping!

subroutine rhs (gt, kk,rt)
	use params

	implicit none
	integer :: i,j,BC,rt
	complex*16, dimension(-NX/2:NX/2,-NY/2:NY/2) :: gt, kk
	kk=0
	if(RHSType .eq. 0) then
		!$OMP PARALLEL DO
		do i = -NX/2,NX/2
			do j = -NY/2,NY/2
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
!Harmonic Dimentionless
	if(RHSType .eq. 1) then
		!$OMP PARALLEL DO
		do i = -NX/2,NX/2
			do j = -NY/2,NY/2
				kk(i,j) = 0.5d0*(4.0d0*gt(BC(i,0),BC(j,1))- gt(BC(i,0),BC(j+1,1))&
						-gt(BC(i,0),BC(j-1,1))-gt(BC(i+1,0),BC(j,1))&
						-gt(BC(i-1,0),BC(j,1)))/(DSPACE**2.0d0)&	!laplacian
						+harm_osc_C*gt(i,j)*gt(i,j)*CONJG(gt(i,j))&	!Nonlinear
			 			+OBJPOT(i,j)*gt(i,j)&	!potential
						- harm_osc_mu*gt(i,j)	!Chemical Potential
			end do
		end do
		!$OMP END PARALLEL DO		
	end if

	if (rt == 1) then ! do we need to do damping?
		if ((dampedX .eqv. .false.) .and. (dampedY .eqv. .false.)) then
			kk = kk/(EYE-GAMMAC)
		else 
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
						kk(i,j) = kk(i,j)/EYE 
					end if
				end do
			end do
			!$OMP END PARALLEL DO

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
    		if(s.eq.-NX/2-1 .and. BCX.eq.0)BC=-NX/2
				if(s.eq.-NX/2-1 .and. BCX.eq.1)BC=NX/2
    	case (1)
    		if(s.eq.NY/2+1  .and. BCY.eq.0)BC=NY/2
				if(s.eq.NY/2+1  .and. BCY.eq.1)BC=-NY/2
    		if(s.eq.-NY/2-1 .and. BCY.eq.0)BC=-NY/2
				if(s.eq.-NY/2-1 .and. BCY.eq.1)BC=NY/2
	end select
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!
