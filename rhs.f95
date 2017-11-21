subroutine runit(steps,steptype)
    use params
    use output
    implicit none
    integer :: steps,steptype,i
    open (10, FILE = "STATUS")
    if (steptype == 2 .OR. PLOTALL) then
        call dump_wavefunction_ncdf(0,steptype) !dump initial condition
    end if

    do i = STARTI, steps

        if(steptype == 0) then
            call iterate(.true.,0.0d0)
            IMAGTIME = IMAGTIME + imag(DT)
        else if (steptype == 1) then
            call iterate(NORMALL,1.0d0)
        else if (steptype == 2) then
            call iterate(NORMALL,GAMMAC)
            TIME = TIME + dble(DT)
        end if

        if(steptype == 1 .and. HOLDIC .and. modulo(i,10) == 0) then !vortex imprinting
            GRID = sqrt(GRID*conjg(GRID))*IMPOSEDPHASE
        end if

        if (modulo(i,dumputil) == 0) then
                call calc_misc
        end if

        if (modulo(i,dumputil) == 0 .and. i > 0) then
            if (steptype == 2) then
                write (unit=10,fmt="(a,f6.2,a)") "Simulating:  ",(dble(i)/dble(steps))*100.0d0,"%"
            end if
            if (steptype == 1) then
                write (unit=10,fmt="(a,f6.2,a)") "Smoothing (after ic.in):  ",(dble(i)/dble(steps))*100.0d0,"%"
            end if
            if (steptype == 0) then
                write (unit=10,fmt="(a,f6.2,a)") "Imaginary Time Propogation:  " ,(dble(i)/dble(steps))*100.0d0,"%"
            end if
        end if

        if (steptype == 2 .and. modulo(i,vortexKillFreq) == 0) then
            call kill_vortex_far()
        end if

        if (modulo(i,dumpwf) == 0) then
            if (steptype == 2 .OR. PLOTALL) then
                call dump_wavefunction_ncdf(i,steptype)
            end if
        end if

        if(steptype == 2 .and. potRep) then
            call calc_OBJPOT
        end if

        if(i .eq. killgamma) then
            write(10,*) "Killing gamma on iteration", i
            GAMMAC = 0.0d0
        end if
    end do
    close(10)
end subroutine

!Runge-Kutta 4th Order!
subroutine iterate (renorm,gammaVal)
    use params

    implicit none
    logical :: renorm
    complex*16, dimension(-NX/2:NX/2,-NY/2:NY/2) :: k1,k2,k3,k4
    double precision :: gammaVal

    call rhs(GRID, k1,gammaVal)
    call rhs(GRID + 0.5d0*DT*k1,k2,gammaVal)
    call rhs(GRID + 0.5d0*DT*k2,k3,gammaVal)
    call rhs(GRID + DT*k3,k4,gammaVal)
    GRID = GRID + DT*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4)/6.0d0

    if (renorm) then
        if (RHSType .eq. 0) then
            call calc_norm
            !Approximate what norm should be to get psi_inf=1
            GRID = GRID/sqrt(NORM)
            GRID = GRID*sqrt(DSPACE*DSPACE*COUNT(OBJPOT>1.0d0))
        end if
        if (RHSType .eq. 1) then
            call calc_norm
            GRID = GRID/sqrt(NORM)
        end if
    end if
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!

!Homogeneous Dimensionless - With Ramping!

subroutine rhs (gt,kk,gammaVal)
    use params
    implicit none
    integer :: i,j,BC
    double precision :: r,gammaVal
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
                        -gt(BC(i-1,0),BC(j,1)))/(DSPACE**2.0d0)&    !laplacian
                        +gt(i,j)*gt(i,j)*CONJG(gt(i,j))&    !Nonlinear
                        -gt(i,j)&   !Chemical Potential
                        +OBJPOT(i,j)*gt(i,j)&   !Obstacle potential
                        +tanh(TIME/200.0d0)*VOB*EYE*(gt(BC(i+1,0),BC(j,1))&
                        -gt(BC(i-1,0),BC(j,1)))/(2.0d0*DSPACE)  !Moving frame
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
                        -gt(BC(i-1,0),BC(j,1)))/(DSPACE**2.0d0)&    !laplacian
                        +harm_osc_C*gt(i,j)*gt(i,j)*CONJG(gt(i,j))& !Nonlinear
                        +OBJPOT(i,j)*gt(i,j)&   !potential
                        - harm_osc_mu*gt(i,j)&  !Chemical Potential
                        -ROM*EYE*(i*DSPACE)*(gt(BC(i,0),BC(j+1,1))&
                        -gt(BC(i,0),BC(j-1,1)))/(2.0d0*DSPACE)&
                        +ROM*EYE*(j*DSPACE)*(gt(BC(i+1,0),BC(j,1))&
                        -gt(BC(i-1,0),BC(j,1)))/(2.0d0*DSPACE) !Rotating frame
            end do
        end do
        !$OMP END PARALLEL DO       
    end if
    
    if ((dampedX .eqv. .false.) .and. (dampedY .eqv. .false.) .and. (dampedR .eqv. .false.)) then
        kk = kk/(EYE-gammaVal)
    else if (dampedR .eqv. .false.) then
        !sqaure regional damping
        !$OMP PARALLEL DO
        do i = -NX/2,NX/2
            do j = -NY/2,NY/2
                if((dampedX .and. abs((i*DSPACE)) > dampedXDist) .and.&
                   (dampedY .and. abs((j*DSPACE)) > dampedYDist)) then
                    kk(i,j) = kk(i,j)/(EYE-gammaVal-&
                              max(0.0d0,min(dampedgamma,dampedgamma*0.5*(tanh((abs((i*DSPACE))-dampedXDist-3.0)/2.0)+1)+&
                              dampedgamma*0.5*(tanh((abs((j*DSPACE))-dampedYDist-3.0)/2.0)+1))))
                else if(dampedX .and. abs((i*DSPACE)) > dampedXDist) then
                    kk(i,j) = kk(i,j)/(EYE-gammaVal-&
                              dampedgamma*0.5*(tanh((abs((i*DSPACE))-dampedXDist-3.0)/2.0)+1))
                else if(dampedY .and. abs((j*DSPACE)) > dampedYDist) then
                    kk(i,j) = kk(i,j)/(EYE-gammaVal-&
                              dampedgamma*0.5*(tanh((abs((j*DSPACE))-dampedYDist-3.0)/2.0)+1))
                else 
                    kk(i,j) = kk(i,j)/(EYE-gammaVal)
                end if
            end do
        end do
        !$OMP END PARALLEL DO
    else
        !circle regional damping
        do i = -NX/2,NX/2
            do j = -NY/2,NY/2
                r = SQRT((i*DSPACE)**2.0d0+(j*DSPACE)**2.0d0)
                if (r > dampedRDist) then
                    kk(i,j) = kk(i,j)/(EYE-gammaVal-0.5d0*dampedGamma*(1.0d0+tanh(pi*&
                                (-1.0d0+2.0d0*(r-dampedRDist)/dampedRWidth))))
                else 
                    kk(i,j) = kk(i,j)/(EYE-gammaVal)
                end if
            end do
        end do
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
