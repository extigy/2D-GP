subroutine calc_OBJPOT
    use params
    implicit none
    OBJPOT = 0.0d0
    !Shin Experiment
    if (doShin) then
        call calc_OBJPOT_shin
    else
        if (enablePot) then
            if(potType .eq. 0) then
                call calc_OBJPOT_obj
            end if
            if(potType .eq. 1) then
                call calc_OBJPOT_rot
            end if
            if(potType .eq. 2) then
                call calc_OBJPOT_pins_square
            end if
            if(potType .eq. 4) then
                call calc_OBJPOT_img
            end if
        end if
        if (enableTrap) then
            if (trapType .eq. 0) then
                call add_harmonic_trap
            end if
            if (trapType .eq. 1) then
                call add_hard_circle_trap
            end if
            if (trapType .eq. 2) then
                call add_hard_box_trap 
            end if
            if (trapType .eq. 3) then
                call add_soft_circle_trap
            end if
            if (trapType .eq. 4) then
                call add_soft_box_trap
            end if

        end if
    end if

end subroutine

subroutine add_harmonic_trap
    use params
    implicit none
    integer :: i,j
    !$OMP PARALLEL DO
    do i = -NX/2,NX/2
    do j = -NY/2,NY/2
        OBJPOT(i,j) = OBJPOT(i,j) + 0.5d0*TXSCALE*((dble(i)*DSPACE+TXDASH)*(dble(i)*DSPACE+TXDASH))&
                        +0.5d0*TYSCALE*((dble(j)*DSPACE+TYDASH)*(dble(j)*DSPACE+TYDASH))
    end do
    end do
    !$OMP END PARALLEL DO
end subroutine

subroutine add_soft_circle_trap
    use params
    implicit none
    integer :: i,j
    double precision :: rx,ry,r
    !$OMP PARALLEL DO
    do i = -NX/2,NX/2
    do j = -NY/2,NY/2
        rx = (dble(i)*DSPACE)-TXDASH
        ry = (dble(j)*DSPACE)-TYDASH
        r = SQRT(rx**2.0+ry**2.0)
        OBJPOT(i,j) = OBJPOT(i,j)+TRAPHEIGHT/(1.0d0+EXP(-TRAPBETA*(ABS(rx)-TRAPR)))
    end do
    end do
    !$OMP END PARALLEL DO
end subroutine

subroutine add_hard_circle_trap
    use params
    implicit none
    integer :: i,j
    double precision :: rx,ry,r
    !$OMP PARALLEL DO
    do i = -NX/2,NX/2
    do j = -NY/2,NY/2
        rx = (dble(i)*DSPACE)-TXDASH
        ry = (dble(j)*DSPACE)-TYDASH
        r = SQRT(rx**2.0+ry**2.0)
        if (r > TRAPR) then
            OBJPOT(i,j) = OBJPOT(i,j)+TRAPHEIGHT
        end if
    end do
    end do
    !$OMP END PARALLEL DO
end subroutine

subroutine add_soft_box_trap
    use params
    implicit none
    integer :: i,j
    double precision :: rx,ry
    !$OMP PARALLEL DO
    do i = -NX/2,NX/2
        do j = -NY/2,NY/2
            rx = (dble(i)*DSPACE)-TXDASH
            ry = (dble(j)*DSPACE)-TYDASH
            OBJPOT(i,j) = OBJPOT(i,j)+MIN(TRAPHEIGHT,TRAPHEIGHT/(1.0d0+EXP(-TRAPBETA*(ABS(rx)-TRAPR)))+&
                                          TRAPHEIGHT/(1.0d0+EXP(-TRAPBETA*(ABS(ry)-TRAPR))))
        end do
    end do
    !$OMP END PARALLEL DO
end subroutine

subroutine add_hard_box_trap
    use params
    implicit none
    integer :: i,j
    double precision :: rx,ry
    !$OMP PARALLEL DO
    do i = -NX/2,NX/2
        do j = -NY/2,NY/2
            rx = (dble(i)*DSPACE)-TXDASH
            ry = (dble(j)*DSPACE)-TYDASH
            if (abs(rx) > TRAPR .OR. abs(ry) > TRAPR) then
                OBJPOT(i,j) = OBJPOT(i,j)+TRAPHEIGHT
            end if
        end do
    end do
    !$OMP END PARALLEL DO
end subroutine

subroutine calc_OBJPOT_obj
    use params
    implicit none
    integer :: i,j
    double precision :: rx,ry
    !$OMP PARALLEL DO
    do i = -NX/2,NX/2
        do j = -NY/2,NY/2
            rx = (dble(i)*DSPACE)-OBJXDASH-(OBJAMP*sin(OBJW*TIME))
            ry = (dble(j)*DSPACE)-OBJYDASH
            OBJPOT(i,j) = OBJHEIGHT*EXP(-(1.0d0/RRX**2.0d0)*(rx**2.0d0) - (1.0d0/RRY**2.0d0)*(ry**2.0d0))
        end do
    end do
    !$OMP END PARALLEL DO
end subroutine

subroutine calc_OBJPOT_pins_square
    use params
    implicit none
    integer :: i,j
    double precision :: rx,ry,rxy
    !$OMP PARALLEL DO
    do i = -NX/2,NX/2
        do j = -NY/2,NY/2
            rx = (dble(i)*DSPACE)-OBJXDASH
            ry = (dble(j)*DSPACE)-OBJYDASH
            rxy = (MODULO(rx,OBJPINSDIST)-(OBJPINSDIST/2.0d0))**2.0d0 + (MODULO(ry,OBJPINSDIST)-(OBJPINSDIST/2.0d0))**2.0d0
            OBJPOT(i,j) = OBJHEIGHT*EXP(-(1.0d0/OBJPINSSIG**2.0d0)*rxy)
        end do
    end do
    !$OMP END PARALLEL DO
end subroutine

subroutine calc_OBJPOT_rot
   use params
   implicit none
   integer :: i,j
   double precision :: rtheta
   double precision, dimension(2,2) :: rotmat
   double precision, dimension(2) :: point,rotpoint
   call calc_new_obj_angle
   rtheta = OBJANGLE
   rotmat(1,1) = cos(rtheta)
   rotmat(2,1) = -sin(rtheta)
   rotmat(1,2) = sin(rtheta)
   rotmat(2,2) = cos(rtheta)
    !$OMP PARALLEL DO
    do i = -NX/2,NX/2
        do j = -NY/2,NY/2
            point(1) = dble(i)*DSPACE
            point(2) = dble(j)*DSPACE
            rotpoint = matmul(rotmat,point)
            OBJPOT(i,j) = OBJHEIGHT*EXP( -(1.0d0/RRX**2.0d0)*(rotpoint(1)-OBJXDASH)**2.0d0&
                            -(1.0d0/RRY**2.0d0)*(rotpoint(2)-OBJYDASH)**2.0d0 )
       end do
   end do
   !$OMP END PARALLEL DO
end subroutine

subroutine calc_OBJPOT_img
    use params
    use bitmap
    implicit none
    integer :: i,j
    
    CALL load_bmp
    do i = -NX/2,NX/2
        do j = -NY/2,NY/2
            OBJPOT(i,j) = OBJHEIGHT*(Img(i,j)/255.0d0)
        end do
    end do
end subroutine

subroutine calc_OBJPOT_shin
!-Note:This is the shin experiment
    use params
    implicit none
    integer :: i,j
    double precision :: trVelx,trVely
    !Section 1 - Ramping up trap movement
    if(TIME .lt. (TTM/8.0d0)) then
        if (enablePot) then
            call calc_OBJPOT_obj
        end if
        if (enableTrap) then
            trVelx = (TVXDASH/2.0d0)*(tanh((6.0d0*TIME/(TTM*8.0d0))-3.0d0)+1.0d0)
            trVely = (TVYDASH/2.0d0)*(tanh((6.0d0*TIME/(TTM*8.0d0))-3.0d0)+1.0d0)
            TXDASH  = TXDASH  + (trVelx*DBLE(DT))
            TYDASH  = TYDASH  + (trVely*DBLE(DT))
            do i = -NX/2,NX/2
            do j = -NY/2,NY/2
                OBJPOT(i,j) = OBJPOT(i,j)&
                    +0.5d0*((dble(i)*DSPACE+TXDASH)*(dble(i)*DSPACE+TXDASH)&
                    +(dble(j)*DSPACE+TYDASH)*(dble(j)*DSPACE+TYDASH))
            end do
            end do
        end if
    end if
    !Section 2 - trap moving at terminal velocity
    if(TIME .gt. (TTM/8.0d0) .and. TIME .lt.  (7.0d0*TTM/8.0d0)) then
        if (enablePot) then
            call calc_OBJPOT_obj
        end if
        if (enableTrap) then
            TXDASH  = TXDASH  + (TVXDASH*DBLE(DT))
            TYDASH  = TYDASH  + (TVYDASH*DBLE(DT))
            do i = -NX/2,NX/2
            do j = -NY/2,NY/2
                OBJPOT(i,j) = OBJPOT(i,j)&
                    +0.5d0*((dble(i)*DSPACE+TXDASH)*(dble(i)*DSPACE+TXDASH)&
                    +(dble(j)*DSPACE+TYDASH)*(dble(j)*DSPACE+TYDASH))
            end do
            end do
        end if
    end if
    !Section 3 - Ramping down trap movement
    if(TIME .gt. (7.0d0*TTM/8.0d0) .and. TIME .lt. TTM) then
        if (enablePot) then
            call calc_OBJPOT_obj
        end if
        if (enableTrap) then
            trVelx = (TVXDASH/2.0d0)&
                *(tanh((-6.0d0*(TIME-(7.0d0*TTM/8.0d0))/(TTM/8.0d0))+3.0d0)+1.0d0)
            trVely = (TVYDASH/2.0d0)&
                *(tanh((-6.0d0*(TIME-(7.0d0*TTM/8.0d0))/(TTM/8.0d0))+3.0d0)+1.0d0)
            TXDASH  = TXDASH  + (trVelx*DBLE(DT))
            TYDASH  = TYDASH  + (trVely*DBLE(DT))
            do i = -NX/2,NX/2
            do j = -NY/2,NY/2
                OBJPOT(i,j) = OBJPOT(i,j)&
                    +0.5d0*((dble(i)*DSPACE+TXDASH)*(dble(i)*DSPACE+TXDASH)&
                    +(dble(j)*DSPACE+TYDASH)*(dble(j)*DSPACE+TYDASH))
            end do
            end do
        end if
    end if

    !Section 4 - Ramping down laser beam
    if(TIME .gt. TTM) then
        if (enablePot) then
            OBJHEIGHT  = OBJHEIGHT  - ((DBLE(DT)*15*harm_osc_mu)/(2.0*PI*15.0*0.4))
            if(OBJHEIGHT .gt. 0.0d0) then
                call calc_OBJPOT_obj
            end if
            if(OBJHEIGHT .lt. 0.0d0) then
                    OBJPOT = 0.0d0
                    OBJHEIGHT = -1.0d0
                    potRep = .false. !Should stop recalculating potential now
            end if
        end if
        if (enableTrap) then
            do i = -NX/2,NX/2
            do j = -NY/2,NY/2
                OBJPOT(i,j) = OBJPOT(i,j)&
                    +0.5d0*((dble(i)*DSPACE+TXDASH)*(dble(i)*DSPACE+TXDASH)&
                    +(dble(j)*DSPACE+TYDASH)*(dble(j)*DSPACE+TYDASH))
            end do
            end do
        end if
    end if

end subroutine