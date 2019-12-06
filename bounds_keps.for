!#############################################################################
        subroutine boundksgs(cmu)
!			Bruño Fraga Bugallo
!			Cardiff 2016
!##############################################################################
        use vars
        use multidata
        implicit none
        integer :: i,j,k,ib,ly
        integer :: is,ie,js,je,ks,ke
	  double precision :: cmu

        if (PERIODIC) call exchange_bc(6,pl_ex)

        do ly=0,pl_ex
        do ib=1,nbp

           is=dom(ib)%isp; ie=dom(ib)%iep
           js=dom(ib)%jsp; je=dom(ib)%jep
           ks=dom(ib)%ksp; ke=dom(ib)%kep
	
! Boundary Conditions for ksgs
!..............................................................................
!=== West ===> 
!..............................................................................
        if (dom(ib)%iprev.lt.0) then
           if (dom(ib)%bc_west.eq.4) then
           do k=ks-1,ke+1; do j=js-1,je+1
              dom(ib)%ksgs(is-1-ly,j,k)= 0.d0					!k=0 at the boundary?
           end do; end do
           else if (dom(ib)%bc_west.ge.63) then					!Wall functions
		if (ly.eq.0) call wall_function(1,ib)
           do k=ks-1,ke+1; do j=js-1,je+1
              dom(ib)%ksgs(is,j,k)=dom(ib)%tauww(j,k)/sqrt(cmu)		!k=ustar**2/cmu**0.5
              dom(ib)%ksgs(is-1-ly,j,k)=0.d0			
           end do; end do
           else if (dom(ib)%bc_west.eq.61
     &				.or.dom(ib)%bc_west.eq.62) then
	     if (ly.eq.0) call log_law(1,ib)
           do k=ks-1,ke+1; do j=js-1,je+1
              dom(ib)%ksgs(is,j,k)=dom(ib)%tauww(j,k)/sqrt(cmu)
              dom(ib)%ksgs(is-1-ly,j,k)=0.d0
           end do; end do
           else if (dom(ib)%bc_west.eq.1) then	
           do k=ks-1,ke+1; do j=js-1,je+1
              dom(ib)%ksgs(is-1-ly,j,k)= 
     &	  (3.d0/2.d0)*(ubulk*0.1)**2.0
           end do; end do
           else if (dom(ib)%bc_west.ne.5) then					!if 5->exchange, if 2,3 -> dk/dn=0  
           do k=ks-1,ke+1; do j=js-1,je+1
              dom(ib)%ksgs(is-1-ly,j,k)= dom(ib)%ksgs(is,j,k)
           end do; end do
           end if
        end if
!...............................................................................
!=== East ===> 
!...............................................................................
        if (dom(ib)%inext.lt.0) then
           if (dom(ib)%bc_east.eq.4) then
           do k=ks-1,ke+1; do j=js-1,je+1
              dom(ib)%ksgs(ie+1+ly,j,k)= 0.d0
           end do; end do
           else if (dom(ib)%bc_east.ge.63) then					!Wall functions
		if (ly.eq.0) call wall_function(2,ib)
           do k=ks-1,ke+1; do j=js-1,je+1
              dom(ib)%ksgs(ie,j,k)=dom(ib)%tauwe(j,k)/sqrt(cmu)
              dom(ib)%ksgs(ie+1+ly,j,k)=0.d0
           end do; end do
           else if (dom(ib)%bc_east.eq.61
     &				.or.dom(ib)%bc_east.eq.62) then
		if (ly.eq.0) call log_law(2,ib)
           do k=ks-1,ke+1; do j=js-1,je+1
              dom(ib)%ksgs(ie,j,k)=dom(ib)%tauwe(j,k)/sqrt(cmu)
              dom(ib)%ksgs(ie+1+ly,j,k)=0.d0
           end do; end do
           else if (dom(ib)%bc_east.ne.5) then
           do k=ks-1,ke+1; do j=js-1,je+1
              dom(ib)%ksgs(ie+1+ly,j,k)= dom(ib)%ksgs(ie,j,k)
           end do; end do
           end if
        end if
!...............................................................................
!=== South ===>
!...............................................................................
        if (dom(ib)%jprev.lt.0) then
           if (dom(ib)%bc_south.eq.4) then
           do k=ks-1,ke+1; do i=is-1,ie+1
              dom(ib)%ksgs(i,js-1-ly,k)= 0.d0
           end do; end do
           else if (dom(ib)%bc_south.ge.63) then					!Wall functions
		if (ly.eq.0) call wall_function(3,ib)
           do k=ks-1,ke+1; do i=is-1,ie+1
              dom(ib)%ksgs(i,js,k)=dom(ib)%tauws(i,k)/sqrt(cmu)
              dom(ib)%ksgs(i,js-1-ly,k)= 0.d0
           end do; end do
           else if (dom(ib)%bc_south.eq.61
     &				.or.dom(ib)%bc_south.eq.62) then
		if (ly.eq.0) call log_law(3,ib)
           do k=ks-1,ke+1; do i=is-1,ie+1
              dom(ib)%ksgs(i,js,k)=dom(ib)%tauws(i,k)/sqrt(cmu)
              dom(ib)%ksgs(i,js-1-ly,k)= 0.d0
           end do; end do
           else if (dom(ib)%bc_south.ne.5) then
           do k=ks-1,ke+1; do i=is-1,ie+1
              dom(ib)%ksgs(i,js-1-ly,k)= dom(ib)%ksgs(i,js,k)
           end do; end do
           end if
        end if
!.............................................................................
!=== North ===>
!.............................................................................
        if (dom(ib)%jnext.lt.0) then
           if (dom(ib)%bc_north.eq.4) then
           do k=ks-1,ke+1; do i=is-1,ie+1
              dom(ib)%ksgs(i,je+1+ly,k) = 0.d0
           end do; end do
           else if (dom(ib)%bc_north.ge.63) then					!Wall functions
		if (ly.eq.0) call wall_function(4,ib)
           do k=ks-1,ke+1; do i=is-1,ie+1
              dom(ib)%ksgs(i,je,k)=dom(ib)%tauwn(i,k)/sqrt(cmu)
              dom(ib)%ksgs(i,je+1+ly,k) = 0.d0
           end do; end do
           else if (dom(ib)%bc_north.eq.61
     &				.or.dom(ib)%bc_north.eq.62) then
		if (ly.eq.0) call log_law(4,ib)
           do k=ks-1,ke+1; do i=is-1,ie+1
              dom(ib)%ksgs(i,je,k)=dom(ib)%tauwn(i,k)/sqrt(cmu)
              dom(ib)%ksgs(i,je+1+ly,k) = 0.d0
           end do; end do
           else if (dom(ib)%bc_north.ne.5) then
           do k=ks-1,ke+1; do i=is-1,ie+1
              dom(ib)%ksgs(i,je+1+ly,k) = dom(ib)%ksgs(i,je,k)	
           end do; end do
           end if
        end if
!...............................................................................
!=== Bottom ===>
!...............................................................................
        if (dom(ib)%kprev.lt.0) then
           if (dom(ib)%bc_bottom.eq.4) then
           do j=js-1,je+1; do i=is-1,ie+1
              dom(ib)%ksgs(i,j,ks-1-ly)= 0.d0
           end do; end do
           else if (dom(ib)%bc_bottom.ge.63) then					!Wall functions
		if (ly.eq.0) call wall_function(5,ib)
           do j=js-1,je+1; do i=is-1,ie+1
              dom(ib)%ksgs(i,j,ks)=dom(ib)%tauwb(i,j)/sqrt(cmu)
              dom(ib)%ksgs(i,j,ks-1-ly)= 0.d0
           end do; end do
           else if (dom(ib)%bc_bottom.eq.61
     &				.or.dom(ib)%bc_bottom.eq.62) then
		if (ly.eq.0) then
			call log_law(5,ib)
		endif
           do j=js-1,je+1; do i=is-1,ie+1
              dom(ib)%ksgs(i,j,ks)=dom(ib)%tauwb(i,j)/sqrt(cmu)
              dom(ib)%ksgs(i,j,ks-1-ly)= 0.d0
           end do; end do
           else if (dom(ib)%bc_bottom.ne.5) then
           do j=js-1,je+1; do i=is-1,ie+1
              dom(ib)%ksgs(i,j,ks-1-ly)= dom(ib)%ksgs(i,j,ks)
           end do; end do
           end if
        end if
!.............................................................................
!=== Top ===>
!.............................................................................
        if (dom(ib)%knext.lt.0) then
           if (dom(ib)%bc_top.eq.4) then
           do j=js-1,je+1; do i=is-1,ie+1
              dom(ib)%ksgs(i,j,ke+1+ly) = 0.d0
           end do; end do
           else if (dom(ib)%bc_top.ge.63) then					!Wall functions
		if (ly.eq.0) call wall_function(6,ib)
           do j=js-1,je+1; do i=is-1,ie+1
              dom(ib)%ksgs(i,j,ke)=dom(ib)%tauwt(i,j)/sqrt(cmu)
              dom(ib)%ksgs(i,j,ke+1+ly) =0.d0
           end do; end do
           else if (dom(ib)%bc_top.eq.61
     &				.or.dom(ib)%bc_top.eq.62) then
		if (ly.eq.0) call log_law(6,ib)
           do j=js-1,je+1; do i=is-1,ie+1
             dom(ib)%ksgs(i,j,ke)=dom(ib)%tauwt(i,j)/sqrt(cmu)
              dom(ib)%ksgs(i,j,ke+1+ly) =0.d0
           end do; end do
           else if (dom(ib)%bc_top.ne.5) then
           do j=js-1,je+1; do i=is-1,ie+1
              dom(ib)%ksgs(i,j,ke+1+ly) = dom(ib)%ksgs(i,j,ke)
           end do; end do
           end if
        end if

!==============================================================================
        end do
        end do
        end subroutine boundksgs
!#############################################################################
        subroutine boundeps
!##############################################################################
        use vars
        use multidata
        implicit none
        integer :: i,j,k,ib,ly
        integer :: is,ie,js,je,ks,ke
	  double precision :: rrey,rk1,drkdy,lz
	  double precision :: dx,dy,dz,kappa,delta

        rrey=1.0/Re
	  kappa=0.41
	  lz=zen-zst
	        if (L_LSM) lz=length
	        if (L_LSMbase) lz=length

        if (PERIODIC) call exchange_bc(9,pl_ex)

        do ly=0,pl_ex 
        do ib=1,nbp

           dx=dom(ib)%dx
           dy=dom(ib)%dy
           dz=dom(ib)%dz

           is=dom(ib)%isp; ie=dom(ib)%iep
           js=dom(ib)%jsp; je=dom(ib)%jep
           ks=dom(ib)%ksp; ke=dom(ib)%kep
	
! Boundary Conditions for epsylon
!..............................................................................
!=== West ===> 
!..............................................................................
        if (dom(ib)%iprev.lt.0) then
           if (dom(ib)%bc_west.eq.4) then
           do k=ks-1,ke+1; do j=js-1,je+1
              rk1 = sqrt(dom(ib)%ksgs(is,j,k))
              drkdy = rk1 /(0.5*dom(ib)%dx)
              dom(ib)%eps(is-1-ly,j,k)= 2.d0*rrey*drkdy**2.0
           end do; end do
           else if (dom(ib)%bc_west.ge.63) then					!Wall functions
		if (ly.eq.0) call wall_function(1,ib)
           do k=ks-1,ke+1; do j=js-1,je+1
              dom(ib)%eps(is,j,k)= 
     &		dom(ib)%tauww(j,k)**1.5/(kappa*dx)
              dom(ib)%eps(is-1-ly,j,k)= 0.d0
           end do; end do
           else if (dom(ib)%bc_west.eq.61
     &				.or.dom(ib)%bc_west.eq.62) then
		if (ly.eq.0) call log_law(1,ib)
           do k=ks-1,ke+1; do j=js-1,je+1
              dom(ib)%eps(is,j,k)= 
     &		dom(ib)%tauww(j,k)**1.5/(kappa*dx)
              dom(ib)%eps(is-1-ly,j,k)= 0.d0
           end do; end do
           else if (dom(ib)%bc_west.eq.1) then	
           do k=ks-1,ke+1; do j=js-1,je+1			 
              dom(ib)%eps(is-1-ly,j,k) =
     &	  0.09**0.75*dom(ib)%ksgs(is-1-ly,j,k)**1.5/(0.07*lz)
           end do; end do
           else if (dom(ib)%bc_west.ne.5) then
           do k=ks-1,ke+1; do j=js-1,je+1
              dom(ib)%eps(is-1-ly,j,k)= dom(ib)%eps(is,j,k)
           end do; end do
           end if
        end if
!...............................................................................
!=== East ===> 
!...............................................................................
        if (dom(ib)%inext.lt.0) then
           if (dom(ib)%bc_east.eq.4) then
           do k=ks-1,ke+1; do j=js-1,je+1
              rk1 = sqrt(dom(ib)%ksgs(ie,j,k))
              drkdy = rk1 /(0.5*dom(ib)%dx)
              dom(ib)%eps(ie+1+ly,j,k)= 2.d0*rrey*drkdy**2.0
           end do; end do
           else if (dom(ib)%bc_east.ge.63) then					!Wall functions
		if (ly.eq.0) call wall_function(2,ib)
           do k=ks-1,ke+1; do j=js-1,je+1
              dom(ib)%eps(ie,j,k)= 
     &		dom(ib)%tauwe(j,k)**1.5/(kappa*dx)
              dom(ib)%eps(ie+1+ly,j,k)= 0.d0
           end do; end do
           else if (dom(ib)%bc_east.eq.61
     &				.or.dom(ib)%bc_east.eq.62) then
		if (ly.eq.0) call log_law(2,ib)
           do k=ks-1,ke+1; do j=js-1,je+1
              dom(ib)%eps(ie,j,k)= 
     &		dom(ib)%tauwe(j,k)**1.5/(kappa*dx)
              dom(ib)%eps(ie+1+ly,j,k)= 0.d0
           end do; end do
           else if (dom(ib)%bc_east.ne.5) then
           do k=ks-1,ke+1; do j=js-1,je+1
              dom(ib)%eps(ie+1+ly,j,k)= dom(ib)%eps(ie,j,k)
           end do; end do
           end if
        end if
!...............................................................................
!=== South ===>
!...............................................................................
        if (dom(ib)%jprev.lt.0) then
           if (dom(ib)%bc_south.eq.4) then
           do k=ks-1,ke+1; do i=is-1,ie+1
              rk1 = sqrt(dom(ib)%ksgs(i,js,k))
              drkdy = rk1 /(0.5*dom(ib)%dy)
              dom(ib)%eps(i,js-1-ly,k)= 2.d0*rrey*drkdy**2.0	
           end do; end do
           else if (dom(ib)%bc_south.ge.63) then					!Wall functions
		if (ly.eq.0) call wall_function(3,ib)
           do k=ks-1,ke+1; do i=is-1,ie+1
              dom(ib)%eps(i,js,k)= 
     &		dom(ib)%tauws(i,k)**1.5/(kappa*dy)
              dom(ib)%eps(i,js-1-ly,k)= 0.d0	
           end do; end do
           else if (dom(ib)%bc_south.eq.61
     &				.or.dom(ib)%bc_south.eq.62) then
		if (ly.eq.0) call log_law(3,ib)
           do k=ks-1,ke+1; do i=is-1,ie+1
              dom(ib)%eps(i,js,k)= 
     &		dom(ib)%tauws(i,k)**1.5/(kappa*dy)
              dom(ib)%eps(i,js-1-ly,k)= 0.d0
           end do; end do
           else if (dom(ib)%bc_south.ne.5) then
           do k=ks-1,ke+1; do i=is-1,ie+1
              dom(ib)%eps(i,js-1-ly,k)= dom(ib)%eps(i,js,k)
           end do; end do
           end if
        end if
!.............................................................................
!=== North ===>
!.............................................................................
        if (dom(ib)%jnext.lt.0) then
           if (dom(ib)%bc_north.eq.4) then
           do k=ks-1,ke+1; do i=is-1,ie+1
              rk1 = sqrt(dom(ib)%ksgs(i,je,k))
              drkdy = rk1 /(0.5*dom(ib)%dy)
              dom(ib)%eps(i,je+1+ly,k)= 2.d0*rrey*drkdy**2.0
           end do; end do
           else if (dom(ib)%bc_north.ge.63) then					!Wall functions
		if (ly.eq.0) call wall_function(4,ib)
           do k=ks-1,ke+1; do i=is-1,ie+1
              dom(ib)%eps(i,je,k)= 
     &		dom(ib)%tauwn(i,k)**1.5/(kappa*dy)
              dom(ib)%eps(i,je+1+ly,k) = 0.d0	
           end do; end do
           else if (dom(ib)%bc_north.eq.61
     &				.or.dom(ib)%bc_north.eq.62) then
		if (ly.eq.0) call log_law(4,ib)
           do k=ks-1,ke+1; do i=is-1,ie+1
              dom(ib)%eps(i,je,k)= 
     &		dom(ib)%tauwn(i,k)**1.5/(kappa*dy)
              dom(ib)%eps(i,je+1+ly,k) = 0.d0
           end do; end do
           else if (dom(ib)%bc_north.ne.5) then
           do k=ks-1,ke+1; do i=is-1,ie+1
              dom(ib)%eps(i,je+1+ly,k) = dom(ib)%eps(i,je,k)	
           end do; end do
           end if
        end if
!...............................................................................
!=== Bottom ===>
!...............................................................................
        if (dom(ib)%kprev.lt.0) then
           if (dom(ib)%bc_bottom.eq.4) then
           do j=js-1,je+1; do i=is-1,ie+1
              rk1 = sqrt(dom(ib)%ksgs(i,j,ks))
              drkdy = rk1 /(0.5*dom(ib)%dz)
              dom(ib)%eps(i,j,ks-1-ly)= 2.d0*rrey*drkdy**2.0	
           end do; end do
           else if (dom(ib)%bc_bottom.ge.63) then					!Wall functions
		if (ly.eq.0) call wall_function(5,ib)
           do j=js-1,je+1; do i=is-1,ie+1
              dom(ib)%eps(i,j,ks)= 
     &		dom(ib)%tauwb(i,j)**1.5/(kappa*dz)
              dom(ib)%eps(i,j,ks-1-ly)= 0.d0	
           end do; end do
           else if (dom(ib)%bc_bottom.eq.61
     &				.or.dom(ib)%bc_bottom.eq.62) then
		if (ly.eq.0) then
			call log_law(5,ib)
		endif
           do j=js-1,je+1; do i=is-1,ie+1
              dom(ib)%eps(i,j,ks)= 
     &		dom(ib)%tauwb(i,j)**1.5/(kappa*dz)
              dom(ib)%eps(i,j,ks-1-ly)= 0.d0
           end do; end do
           else if (dom(ib)%bc_bottom.ne.5) then
           do j=js-1,je+1; do i=is-1,ie+1
              dom(ib)%eps(i,j,ks-1-ly)= dom(ib)%eps(i,j,ks)
           end do; end do
           end if
        end if
!.............................................................................
!=== Top ===>
!.............................................................................
        if (dom(ib)%knext.lt.0) then
           if (dom(ib)%bc_top.eq.4) then
           do j=js-1,je+1; do i=is-1,ie+1
              rk1 = sqrt(dom(ib)%ksgs(i,j,ke))
              drkdy = rk1 /(0.5*dom(ib)%dz)
              dom(ib)%eps(i,j,ke+1+ly)= 2.d0*rrey*drkdy**2.0				
           end do; end do
           else if (dom(ib)%bc_top.ge.63) then					!Wall functions
		if (ly.eq.0) call wall_function(6,ib)
           do j=js-1,je+1; do i=is-1,ie+1
              dom(ib)%eps(i,j,ke)= 
     &		dom(ib)%tauwt(i,j)**1.5/(kappa*dz)
              dom(ib)%eps(i,j,ke+1+ly) = 0.d0
           end do; end do
           else if (dom(ib)%bc_top.eq.61
     &				.or.dom(ib)%bc_top.eq.62) then
		if (ly.eq.0) call log_law(6,ib)
           do j=js-1,je+1; do i=is-1,ie+1
              dom(ib)%eps(i,j,ke)= 
     &		dom(ib)%tauwt(i,j)**1.5/(kappa*dz)
              dom(ib)%eps(i,j,ke+1+ly) = 0.d0
           end do; end do
           else if (dom(ib)%bc_top.ne.5) then
           do j=js-1,je+1; do i=is-1,ie+1
              dom(ib)%eps(i,j,ke+1+ly) = dom(ib)%eps(i,j,ke)
           end do; end do
           end if
        end if

!==============================================================================
        end do
        end do
        end subroutine boundeps
!#############################################################################
