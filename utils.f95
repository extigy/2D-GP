subroutine add_noise
	use params
	implicit none
	integer :: i,j,seed
	call SYSTEM_CLOCK(COUNT=seed)
	call srand(seed)
	do i = -NX/2, NX/2
		do j = -NY/2, NY/2
			GRID(i,j) = GRID(i,j) + GRID(i,j)*CMPLX((RAND()*noiseamp)&
									-(noiseamp/2.0d0),(RAND()*noiseamp)-(noiseamp/2.0d0))
		end do
	end do
end subroutine

subroutine load_all(t)
	use params
	implicit none
	integer :: i,j,t
	double precision :: x, y, real, imag, pot
	character(2048) :: tmp,fname
	fname = './xxxx.dumpwf.xxxx'
	WRITE (fname(3:6), '(I4)') LOOPNO
	WRITE (fname(15:18), '(I4)') t
	open (999, FILE = fname)
	do i = -NX/2, NX/2
		do j = -NY/2, NY/2
			read (999,"(f10.2,f10.2,3F20.10)") x, y, real, imag, pot
			GRID(i,j) = real + EYE*imag
			OBJPOT(i,j) = pot
		end do
		read (999,"(a)") tmp
	end do
	close(999)
	TIME = t*dumpwf*DTSIZE
	STARTI = t*dumpwf + 1
	write (8, *) 'Loading file:',fname,' with time: ',TIME, ' and frame: ', STARTI
	call calc_OBJPOT
end subroutine


!Approx - homogeneous!
subroutine approx
	use params
	implicit none
	integer :: i,j
	if (RHSType .eq. 0) then
		GRID =  1.00d0
	end if
	if (RHSType .eq. 1) then
		do i = -NX/2, NX/2
			do j = -NY/2, NY/2
				GRID(i,j) = (harm_osc_mu - OBJPOT(i,j))/harm_osc_C
				if(DBLE(GRID(i,j)) < 0) GRID(i,j) = 0.0d0
			end do
		end do
		GRID = SQRT(GRID)
	end if
end subroutine
!!!!!!!!!!!!!!!!!!!!!!

subroutine calc_norm
	use params
	implicit none
	integer :: i,j
	double precision :: simpsons_int_grid
	NORM=simpsons_int_grid(DBLE(GRID*CONJG(GRID)))
end subroutine

function simpsons_int_grid(intThing)
	use params
	implicit none
	double precision :: simpsons_int_grid,a(16)
	double precision, dimension(-NX/2:NX/2,-NY/2:NY/2) :: intThing

	a(1) = intThing(-NX/2,-NY/2)
	a(2) = intThing(-NX/2,NY/2)
	a(3) = intThing(NX/2,-NY/2)
	a(4) = intThing(NX/2,NY/2)

	a(5) = 4.0d0*sum(intThing(-NX/2,-NY/2+1:NY/2-1:2))
	a(6) = 2.0d0*sum(intThing(-NX/2,-NY/2+2:NY/2-1:2))
	a(7) = 4.0d0*sum(intThing(NX/2,-NY/2+1:NY/2-1:2))
	a(8) = 2.0d0*sum(intThing(NX/2,-NY/2+2:NY/2-1:2))

	a(9) = 4.0d0*sum(intThing(-NX/2+1:NX/2-1:2,-NY/2))
  a(10) = 2.0d0*sum(intThing(-NX/2+2:NX/2-1:2,-NY/2))
  a(11) = 4.0d0*sum(intThing(-NX/2+1:NX/2-1:2,NY/2))
  a(12) = 2.0d0*sum(intThing(-NX/2+2:NX/2-1:2,NY/2))

  a(13) = 16.0d0*sum(intThing(-NX/2+1:NX/2-1:2,-NY/2+1:NY/2-1:2))
  a(14) = 8.0d0* sum(intThing(-NX/2+1:NX/2-1:2,-NY/2+2:NY/2-1:2))
  a(15) = 8.0d0* sum(intThing(-NX/2+2:NX/2-1:2,-NY/2+1:NY/2-1:2))
  a(16) = 4.0d0*sum(intThing(-NX/2+2:NX/2-1:2,-NY/2+2:NY/2-1:2))

   simpsons_int_grid = (1.0d0/9.0d0)*(DSPACE*DSPACE)*sum(a)
end function

subroutine calc_misc
	use params
	implicit none
	double precision :: energy
	double precision, dimension(2) :: force
	call calc_force(force)
	call calc_energy(energy)
    write (unit=8,fmt="(f7.2,f15.8,5f15.10)") time,energy,force(1),force(2)
    flush(8)
end subroutine

subroutine calc_force(force)
	use params
	implicit none
	integer :: i,j
	COMPLEX*16 :: uux,uuy,uu
	double precision, dimension(2) :: force
	double precision, dimension(2) :: fvec
	fvec = 0.0d0

	do i = -NX/2+1, NX/2-1
		do j =  -NY/2+1, NY/2-1
			uu=GRID(i,j)
			uux=(GRID(i+1,j)-GRID(i-1,j))/(2.0d0*DSPACE)
			uuy=(GRID(i,j+1)-GRID(i,j-1))/(2.0d0*DSPACE)

			fvec(1) = fvec(1) + CONJG(uu)*EYE*uux
			fvec(2) = fvec(2) + CONJG(uu)*EYE*uuy
		end do
	end do
	fvec(1) = fvec(1)*DSPACE*DSPACE
	fvec(2) = fvec(2)*DSPACE*DSPACE

	force(1) = (fvec(1) - FVECOLD(1))/(2.0d0*DT*10.0d0)
	force(2) = (fvec(2) - FVECOLD(2))/(2.0d0*DT*10.0d0)

	FVECOLD(1) = fvec(1)
	FVECOLD(2) = fvec(2)

end subroutine

subroutine calc_force_2D(force2d)
  use params

  implicit none
  integer :: i,j
  COMPLEX*16 :: uu,uux,uuy
  double precision, dimension(2,-NX/2:NX/2,-NY/2:NY/2) :: force2d
  do i = -NX/2+1, NX/2-1
      do j = -NY/2+1, NY/2-1
        uu=GRID(i,j)*CONJG(GRID(i,j))
        uux=(OBJPOT(i+1,j)-OBJPOT(i-1,j))/(2.0d0*DSPACE)
        uuy=(OBJPOT(i,j+1)-OBJPOT(i,j-1))/(2.0d0*DSPACE)

  		  force2d(1,i,j) = uu*uux
		  	force2d(2,i,j) = uu*uuy
      end do
  end do
end subroutine

subroutine cross_2d(a, b,ret)
	implicit none
	double precision :: ret
    double precision, dimension(2) :: a,b
    ret = a(1)*b(2) - a(2)*b(1)
END subroutine

subroutine calc_net_torque(torque)
	use params
	implicit none
	integer :: i,j
	double precision, dimension(2,-NX/2:NX/2,-NY/2:NY/2) :: force2d
	double precision, dimension(2) :: r_rel
	double precision :: torque,crossp

	torque = 0.0d0
	call calc_force_2D(force2d)
	do i = -NX/2+1, NX/2-1
		do j = -NY/2+1, NY/2-1
			r_rel(1) = (dble(i)*DSPACE)-OBJXDASH
			r_rel(2) = (dble(j)*DSPACE)-OBJYDASH
			call cross_2d(r_rel,force2d(:,i,j),crossp)
		   torque = torque + crossp
		end do
	end do
end subroutine

subroutine calc_new_obj_angle
  use params
	implicit none
	double precision :: nettorque
	call calc_net_torque(nettorque)
	OBJANGLEV = OBJANGLEV + MOMINERTIA*(nettorque*DT)
	OBJANGLE  = OBJANGLE  + (OBJANGLEV*DT)
end subroutine

subroutine calc_energy(energy)
	use params
	implicit none
	integer :: i,j
	COMPLEX*16 :: uux,uuy,uu
	double precision :: energy,gg
	energy = 0.0d0
	if(RHSType .eq. 0) then
		gg=1.0d0
	end if
	if(RHSType .eq. 1) then
		gg=harm_osc_C
	end if
	do i = -NX/2+1, NX/2-1
		do j = -NY/2+1, NY/2-1
			uu=GRID(i,j)
			uux=(GRID(i+1,j)-GRID(i-1,j))/(2.0d0*DSPACE)
			uuy=(GRID(i,j+1)-GRID(i,j-1))/(2.0d0*DSPACE)

			energy = energy + 0.5d0*((uux*CONJG(uux))+(uuy*CONJG(uuy))) &
				+ 0.5d0*gg*uu*conjg(uu)*uu*conjg(uu)&
				+ OBJPOT(i,j)*uu*conjg(uu)
		end do
	end do
	energy = energy*DSPACE*DSPACE
end subroutine

subroutine calc_phase(phase)
	use params
	implicit none
	integer :: i,j
	double precision, dimension(-NX/2:NX/2,-NY/2:NY/2) :: phase
	do i = -NX/2,NX/2
		do j = -NY/2,NY/2
				phase(i,j) = atan2(imagpart(GRID(i,j)),realpart(GRID(i,j)))
		end do
	end do
end subroutine

subroutine insert_vortex(xloc,yloc,circ)
	use params
	implicit none
	integer :: i,j,k
	double precision :: xloc,yloc,rs,xx,yy,circ
	double precision, dimension(-NX/2:NX/2,-NY/2:NY/2) :: R
	complex*16, dimension(-NX/2:NX/2,-NY/2:NY/2) :: phse

	do i = -NX/2, NX/2
		do j = -NY/2, NY/2
			xx = (i*DSPACE)
			yy = (j*DSPACE)
			rs=(xx-xloc)**2.0d0 + (yy-yloc)**2.0d0;
			phse(i,j) = exp(circ*EYE*(atan2(yy-yloc,xx-xloc)))
			R(i,j) = sqrt(rs*(0.3437+0.0286*rs)/(1+0.3333*rs+0.0286*rs*rs)); 
		end do
	end do
	GRID = GRID*R*phse;
end subroutine

subroutine velxy(phase,velx,vely)
	use params
	implicit none
	integer :: i,j
	double precision, dimension(-NX/2:NX/2,-NY/2:NY/2) :: phase,velx,vely
	double precision :: temp1,temp2
	velx = 0.0d0
	vely = 0.0d0
	do i = -NX/2+1,NX/2-1
	do j = -NY/2+1,NY/2-1
		velx(i,j) = phase(i+1,j)-phase(i-1,j)
		do while (velx(i,j) > PI)
			velx(i,j) = velx(i,j) - 2.0d0*PI
		end do
		do while (velx(i,j) < -PI)
			velx(i,j) = velx(i,j) + 2.0d0*PI
		end do
	end do
	end do

	do i = -NX/2+1,NX/2-1
	do j = -NY/2+1,NY/2-1
		vely(i,j) = phase(i,j+1)-phase(i,j-1)
		do while (vely(i,j) > PI)
			vely(i,j) = vely(i,j) - 2.0d0*PI
		end do
		do while (vely(i,j) < -PI)
			vely(i,j) = vely(i,j) + 2.0d0*PI
		end do
	end do
	end do
end subroutine

subroutine lineintvf(fieldx,fieldy,x,ex,y,ey,ret)
	use params
	IMPLICIT NONE
	INTEGER :: t,x,ex,y,ey
	DOUBLE PRECISION, DIMENSION(-NX/2:NX/2,-NY/2:NY/2) :: fieldx,fieldy
	DOUBLE PRECISION l1,l2,l3,l4,ret
	l1=0.0d0
	l2=0.0d0
	l3=0.0d0
	l4=0.0d0
	do t = y,ey
		l1 = l1 + DSPACE*fieldy(x,t)
	end do
	do t = x,ex
		l2 = l2 + DSPACE*fieldx(t,y)
	end do
	do t = y,ey
		l3 = l3 + DSPACE*fieldy(ex,t)
	end do
	do t = x,ex
		l4 = l4 + DSPACE*fieldx(t,ey)
	end do
	ret = l2+l3-l4-l1
end subroutine

subroutine box_blur(field,x,y,ret)
	use params
	implicit none
	integer :: i,j,x,y
	double precision, DIMENSION(-NX/2:NX/2,-NY/2:NY/2) :: field
	double precision :: ret
	ret = 0.0d0
	do i = x-1,x+1
		do j = y-1,y+1
			ret = ret + field(i,j)
		end do
	end do
	ret = ret/9.0d0
end subroutine


subroutine calc_vortex_field(velx, vely, prevort)
	use params
	implicit none
	integer :: i,j,n
	double precision, dimension(-NX/2:NX/2,-NY/2:NY/2) :: tempfield,prevort,velx,vely
	prevort = 0.0d0
	do i = -NX/2+3,NX/2-3
		do j = -NY/2+3,NY/2-3
				call lineintvf(velx,vely,i-3,i+3,j-3,j+3,prevort(i,j))
				if(DBLE(OBJPOT(i,j))>5.0d0) then
		        	prevort(i,j) = 0.0d0;
		        end if
		end do
	end do
	do n = 1,3
		tempfield = prevort
		do i = -NX/2+1,NX/2-1
			do j = -NY/2+1,NY/2-1
				call box_blur(tempfield,i,j,prevort(i,j))
			end do
		end do
	end do
end subroutine

subroutine pos_neg_binary(field, threshold, posfield, negfield)
	use params
	implicit none
	integer :: i,j
	double precision, dimension(-NX/2:NX/2,-NY/2:NY/2) :: field
	double precision :: threshold
	integer, dimension(-NX/2:NX/2,-NY/2:NY/2) :: posfield, negfield

	posfield = 0;
	negfield = 0;
	do i = -NX/2,NX/2
		do j = -NY/2,NY/2
			if(field(i,j) > threshold) then
				posfield(i,j) = 1
			end if
			if(field(i,j) < -threshold) then
				negfield(i,j) = 1
			end if
		end do
	end do
end subroutine

subroutine bw_label(bwfield,labels)
	use params
	implicit none
	integer, dimension(-NX/2:NX/2,-NY/2:NY/2) :: bwfield,labels
	integer, dimension(4) :: check
	integer, dimension(NX*NY,4) :: linked
	integer :: i,j,ii,jj,rc,lc,m
	
	labels = -1
	linked = -1
	lc = 1
	rc = 0

	do j = -NY/2+1,NY/2-1
	do i = -NX/2+1,NX/2-1
		if (bwfield(i,j) .ne. 1) cycle
		check(1) = labels(i+1,j-1)
		check(2) = labels(i,j-1)
		check(3) = labels(i-1,j-1)
		check(4) = labels(i-1,j)

		if(check(1) >= 0) then 
			labels(i,j) = check(1)
			linked(lc,1) = check(1)
		end if
		if(check(2) >= 0) then
			labels(i,j) = check(2)
			linked(lc,2) = check(2)
		end if
		if(check(3) >= 0) then
			labels(i,j) = check(3)
			linked(lc,3) = check(3)
		end if
		if(check(4) >= 0) then
			labels(i,j) = check(4)
			linked(lc,4) = check(4)
		end if

		if(labels(i,j) >= 0) then
			lc = lc + 1
		else
			labels(i,j) = rc
			rc = rc + 1
		end if
	end do
	end do


	do jj = 1,NX*NY
	if(maxval(linked(jj,:)) .eq. -1) cycle
	m = minval(linked(jj,:),linked(jj,:) >= 0)
	do ii = 1,4
		if(linked(jj,ii) .ne. m .and. linked(jj,ii) >= 0 ) then
			do j = -NY/2,NY/2
			do i = -NX/2,NX/2
				if (labels(i,j) .eq. linked(jj,ii)) then
					labels(i,j) = m
				end if
			end do
			end do
		end if
	end do
	end do	
end subroutine

subroutine set_homg_sq(i,j,sizexy)
	use params
	implicit none
	integer :: ii,jj,i,j,sizexy
	double precision, dimension(-NX/2:NX/2,-NY/2:NY/2) :: phase
	complex*16, dimension(-NX/2:NX/2,-NY/2:NY/2) :: phse
	
	call calc_phase(phase)

	do ii = i-sizexy,i+sizexy
		do jj = j-sizexy, j+sizexy
			GRID(ii,jj) = 1.0d0*exp(EYE*phase(ii,jj));
		end do
	end do
end subroutine


subroutine kill_vortex_far()
	use params
	implicit none
	double precision, dimension(0:NX*NY,0:1) :: vortarray
	double precision, dimension(-NX/2:NX/2,-NY/2:NY/2) :: phase,velx,vely,vortfield
	integer, dimension(-NX/2:NX/2,-NY/2:NY/2) :: posfield, negfield,poslabels,neglabels
	integer :: i,j,n,m,maxlabel, circ
	double precision :: xloc, yloc

	if(doVortexKilling .eqv. .false.) return

	call calc_phase(phase)
	call velxy(phase,velx,vely)
	call calc_vortex_field(velx, vely, vortfield)
	call pos_neg_binary(vortfield, 0.5d0, posfield, negfield)
	call bw_label(posfield,poslabels)
	call bw_label(negfield,neglabels)

	maxlabel = maxval(poslabels)
	do n = 0, maxlabel
		xloc = 0.0d0
		yloc = 0.0d0
		m = 0
		do j = -NY/2,NY/2
		do i = -NX/2,NX/2
			if (poslabels(i,j) .eq. n) then
				xloc = xloc + (i*DSPACE)
				yloc = yloc + (j*DSPACE)
				m = m + 1
			end if
		end do
		end do

		xloc = xloc/DBLE(m)
		yloc = yloc/DBLE(m)

		if((vortexKillX .and. abs(xloc) > vortexKillXDist) .or.&
		   (vortexKillY .and. abs(yloc) > vortexKillYDist)) then
			write(6,*) "Killing positive vortex at: ", xloc, ", ", yloc
			call insert_vortex(xloc,yloc,-1.0d0)
			call set_homg_sq(nint(xloc/DSPACE),nint(yloc/DSPACE),8)
		end if
	end do

	maxlabel = maxval(neglabels)
	do n = 0, maxlabel
		xloc = 0.0d0
		yloc = 0.0d0
		m = 0
		do j = -NY/2,NY/2
		do i = -NX/2,NX/2
			if (neglabels(i,j) .eq. n) then
				xloc = xloc + (i*DSPACE)
				yloc = yloc + (j*DSPACE)
				m = m + 1
			end if
		end do
		end do

		xloc = xloc/DBLE(m)
		yloc = yloc/DBLE(m)
		if((vortexKillX .and. abs(xloc) > vortexKillXDist) .or.&
		   (vortexKillY .and. abs(yloc) > vortexKillYDist)) then
			write(6,*) "Killing negative vortex at: ", xloc, ", ", yloc
			call insert_vortex(xloc,yloc,1.0d0)
			call set_homg_sq(nint(xloc/DSPACE),nint(yloc/DSPACE),8)
		end if
	end do

end subroutine
