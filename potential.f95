subroutine calc_OBJPOT
	use params
	implicit none
	integer :: i,j
	OBJPOT = 0.0d0

	!Normal potential and Trap
	if (.NOT. doShin) then
		if (enablePot) then
			if(potType .eq. 0) then
				call calc_OBJPOT_obj
			end if
			if(potType .eq. 1) then
				call calc_OBJPOT_rot
			end if
			if(potType .eq. 2) then
				call calc_OBJPOT_osc
			end if
			if(potType .eq. 3) then
				call calc_OBJPOT_afm
			end if
			if(potType .eq. 4) then
				call calc_OBJPOT_img
			end if
			if(potType .eq. 5) then
				call calc_OBJPOT_bump_dt
			end if
		end if
		if (enableTrap) then
			!$OMP PARALLEL DO
			do i = -NX/2,NX/2
			do j = -NY/2,NY/2
				OBJPOT(i,j) = OBJPOT(i,j)&
					+0.5d0*((dble(i)*DSPACE+TXDASH)*(dble(i)*DSPACE+TXDASH)&
					+(dble(j)*DSPACE+TYDASH)*(dble(j)*DSPACE+TYDASH))
			end do
			end do
			!$OMP END PARALLEL DO
		end if
	end if
	!Shin Experiment
	if (doShin) then
		call calc_OBJPOT_shin
	end if
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
			OBJHEIGHT  = OBJHEIGHT  - ((DT*15*harm_osc_mu)/(2.0*PI*15.0*0.4))
			if(OBJHEIGHT .gt. 0.0d0) then
				call calc_OBJPOT_obj
			end if
			if(OBJHEIGHT .lt. 0.0d0) then
					OBJPOT = 0.0d0
					OBJHEIGHT = -1.0d0
					potRep = 0 !Should stop recalculating potential now
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

subroutine calc_OBJPOT_bump_dt
	use params
	implicit none
	integer :: i,j,obj
	double precision :: rx,ry,r2

	if(TIME .gt. 100.and. TIME .lt. 101) then
		do i = -NX/2,NX/2
			do j = -NY/2,NY/2
				rx = (dble(i)*DSPACE)-OBJXDASH
				OBJPOT(i,j) = OBJHEIGHT*EXP(-(1.0d0/RRX**2.0d0)*(rx**2.0d0))
			end do
		end do
	end if

end subroutine

subroutine calc_OBJPOT_obj
	use params
	implicit none
	integer :: i,j,obj
	double precision :: rx,ry,r2
	!$OMP PARALLEL DO
	do i = -NX/2,NX/2
		do j = -NY/2,NY/2
			rx = (dble(i)*DSPACE)-OBJXDASH
			ry = (dble(j)*DSPACE)-OBJYDASH
			OBJPOT(i,j) = OBJHEIGHT*&
				EXP(-(1.0d0/RRX**2.0d0)*(rx**2.0d0) - (1.0d0/RRY**2.0d0)*(ry**2.0d0))
		end do
	end do
	!$OMP END PARALLEL DO
end subroutine


subroutine calc_OBJPOT_osc
	use params
	implicit none
	integer :: i,j,obj
	double precision :: rx,ry,r2
	!$OMP PARALLEL DO
	do i = -NX/2,NX/2
		do j = -NY/2,NY/2
			rx = (dble(i)*DSPACE)-OBJXDASH-(OBJAMP*sin(OBJW*TIME))
			ry = (dble(j)*DSPACE)-OBJYDASH
			OBJPOT(i,j) = OBJHEIGHT*&
				EXP(-(1.0d0/RRX**2.0d0)*(rx**2.0d0) - (1.0d0/RRY**2.0d0)*(ry**2.0d0))
		end do
	end do
	!$OMP END PARALLEL DO
end subroutine

subroutine calc_OBJPOT_rot
   use params
   implicit none
   integer :: i,j,obj
   double precision :: rx,ry,r2,rtheta
   double precision, dimension(2,2) :: rotmat
   double precision, dimension(2) :: point,rotpoint
   call calc_new_obj_angle
   rtheta = OBJANGLE
   rotmat(1,1) = cos(rtheta)
   rotmat(2,1) = -sin(rtheta)
   rotmat(1,2) = sin(rtheta)
   rotmat(2,2) = cos(rtheta)

	do i = -NX/2,NX/2
		do j = -NY/2,NY/2
			point(1) = dble(i)*DSPACE
      point(2) = dble(j)*DSPACE
			rotpoint = matmul(rotmat,point)
      OBJPOT(i,j) = OBJHEIGHT*&
				EXP( -(1.0d0/RRX**2.0d0)*(rotpoint(1)-OBJXDASH)**2.0d0&
				-(1.0d0/RRY**2.0d0)*(rotpoint(2)-OBJYDASH)**2.0d0 )
       end do
   end do
end subroutine

subroutine calc_OBJPOT_afm
	use params
	implicit none
	integer :: i,j,obj,slicei,slicej
	double precision :: di,dj,ity,xtt,ytt,ztt
	double precision, dimension(afmRES) :: xdat,ydat
	double precision, dimension(afmRES,afmRES) :: ydatb
	open (unit = 5, file = TRIM(afm_filename))
	do j = 1,afmRES
	do i = 1,afmRES
		read (5,*) xtt,ytt,ztt
			xdat(i) = ((xtt-0.5d0)*1000.0d0*afmXScale)/xi1
			ydatb(i,j) = (ztt/xi1)
	end do
	end do
	ydat = ydatb(afmSlice,:)
	ydat = ydat - MINVAL(ydat(3:afmRES-3)) - (NY/2)*DSPACE
	ydat = ydat*afmYscale
	OBJPOT = 0.0d0

	do j = -NY/2,NY/2
		dj = j*DSPACE
		do i = -NX/2,NX/2
			di = i*DSPACE
			if(di>maxval(xdat(1:afmRES-3))) then
				di = xdat(afmRES-3)
			end if
			CALL interp(di,xdat,ydat,ity)
			if(dj < (ity+OBJYDASH)) then
				OBJPOT(i,j) = OBJHEIGHT
			end if
		end do
	end do
	close(5)
	
	OBJPOT(:,-NY/2+NINT(NY/TRUNCPARAM):NY/2) = 0.0d0
	
end subroutine

subroutine calc_OBJPOT_img
	use params
	use bitmap
	implicit none
	integer :: i,j,obj
	
	CALL load_bmp
	do i = -NX/2,NX/2
		do j = -NY/2,NY/2
			OBJPOT(i,j) = OBJHEIGHT*(Img(i,j)/255)
		end do
	end do
end subroutine


subroutine interp(xloc,xdat,ydat,rety)
	use params
	implicit none
	integer :: i,si
	double precision :: xloc,rety,low,high,diff,alpha
	double precision, dimension(afmRES) :: xdat,ydat

	si = 1
	do i = 1,afmRES-1
		if(xdat(i)>xloc) then
			si = i
			exit
		end if
	end do

	low = xdat(si-1)
	high = xdat(si)
	diff = high-low
	alpha = (xloc-low)/diff
	low = ydat(si-1)
	high = ydat(si)
	rety = low+((high-low)*alpha)
end subroutine
