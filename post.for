!##########################################################################
        subroutine tecgrid(ki)
!##########################################################################
        use multidata
        use vars
        implicit none
        integer :: sn,sn1,i,j,k,toti,totj,totk
        integer :: is,ie,js,je,ks,ke
        integer :: ki,ib
        double precision ujkl,vikl,wijl
        character*8 :: chb,chb1
        character*25 :: gf

        do ib=1,nbp

	  if (L_anim_grd) then
          write(chb,'(i8)') dom_id(ib)
	    write(chb1,'(i8)') ki
          sn=len(trim(adjustl(chb)))
          sn1=len(trim(adjustl(chb1)))
          chb=repeat('0',(4-sn))//trim(adjustl(chb))
          chb1=repeat('0',(6-sn1))//trim(adjustl(chb1))
          gf='tecgrid'//trim(adjustl(chb))//'_'//
     & trim(adjustl(chb1))//'.dat'
	  else
          write(chb,'(i8)') dom_id(ib)
          sn=len(trim(adjustl(chb)))
          chb=repeat('0',(4-sn))//trim(adjustl(chb))
          gf='tecgrid'//trim(adjustl(chb))//'.dat'
	  endif

        open (unit=88, file=gf)

        is=pl+1; ie=dom(ib)%ttc_i-pl
        js=pl+1; je=dom(ib)%ttc_j-pl
        ks=pl+1; ke=dom(ib)%ttc_k-pl
        toti=ie-(is-1)+1
        totj=je-(js-1)+1
        totk=ke-(ks-1)+1
 
        write (88,*) 'title = ', 'grid quantities'
	  if (L_anim_grd) then
	    if (L_LSM) then
            write (88,*)
     & 'variables="x","y","z","U-velocity","V-velocity","W-velocity",',
     & '"dens"'
	      write(88,*)'zone ', 'STRANDID=', 1, 'SOLUTIONTIME=', ctime,
     &  'i=',toti,', ',' j=',totj,', k= ',totk, 
     &  'zonetype=', 'ordered',', DATAPACKING=point'
	    else
            write (88,*)
     & 'variables="x","y","z","U-velocity","V-velocity","W-velocity"'
	      write(88,*)'zone ', 'STRANDID=', 1, 'SOLUTIONTIME=', ctime,
     &  'i=',toti,', ',' j=',totj,', k= ',totk, 
     &  'zonetype=', 'ordered',', DATAPACKING=point'
	    endif
	  else
	    if (L_LSM) then
            write (88,*)
     & 'variables="x","y","z","U-velocity","V-velocity","W-velocity",',
     & '"dens"'
            write (88,*)  'zone ', ' i=',toti,', ',
     &' j=',totj,' k=',totk,', f=point'
	    else
            write (88,*)
     & 'variables="x","y","z","U-velocity","V-velocity","W-velocity"'
            write (88,*)  'zone ', ' i=',toti,', ',
     &' j=',totj,' k=',totk,', f=point'
	    endif
	  end if

        do k=ks-1,ke
           do j=js-1,je
              do i=is-1,ie

                 ujkl=0.25*(dom(ib)%u(i,j,k) +dom(ib)%u(i,j+1,k) +
     &dom(ib)%u(i,j,k+1) +dom(ib)%u(i,j+1,k+1))
                 vikl=0.25*(dom(ib)%v(i,j,k) +dom(ib)%v(i+1,j,k) + 
     &dom(ib)%v(i,j,k+1) +dom(ib)%v(i+1,j,k+1))
                 wijl=0.25*(dom(ib)%w(i,j,k) +dom(ib)%w(i+1,j,k) +
     &dom(ib)%w(i,j+1,k) +dom(ib)%w(i+1,j+1,k)) 

	          if (L_LSM) then
                  write (88,88) dom(ib)%x(i),dom(ib)%y(j),
     & dom(ib)%z(k),ujkl,vikl,wijl,dom(ib)%dens(i,j,k)
	          else
                  write (88,88) dom(ib)%x(i),dom(ib)%y(j),
     & dom(ib)%z(k),ujkl,vikl,wijl
	          end if

              end do
           end do
        end do
        close (88)

        end do
88      format (7e25.8)

        end
!##########################################################################
        subroutine tecplot_p(ki)
!##########################################################################
        use multidata
        implicit none
        integer :: sn,i,j,k,toti,totj,totk
        integer :: is,ie,js,je,ks,ke
        integer :: ki,ib
        character*8 :: chb
        character*25 :: gf

        do ib=1,nbp

        write(chb,'(i8)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        gf='tecout_p'//trim(adjustl(chb))//'.plt'
        open (unit=88, file=gf)

        is=dom(ib)%isp; ie=dom(ib)%iep
        js=dom(ib)%jsp; je=dom(ib)%jep
        ks=dom(ib)%ksp; ke=dom(ib)%kep
        toti=(ie+1)-(is-1)+1
        totj=(je+1)-(js-1)+1
        totk=(ke+1)-(ks-1)+1

        write (88,*) 'title = ', 'pgrid'
        write (88,*)'variables="x","y","z","Pressure","PMean",',
     & '"ppm","uvM","uwM","vwM","vis","ksgs","eps","T","Tm"'
        write (88,*)'zone ', ' i=',toti,', ',
     &  ' j=',totj,', k= ',totk,' f=point'

        do k=ks-1,ke+1
           do j=js-1,je+1
              do i=is-1,ie+1
                 write (88,88) dom(ib)%xc(i),dom(ib)%yc(j),
     & dom(ib)%zc(k),dom(ib)%p(i,j,k),dom(ib)%pm(i,j,k),
     & dom(ib)%ppm(i,j,k),dom(ib)%uvm(i,j,k),dom(ib)%uwm(i,j,k),
     & dom(ib)%vwm(i,j,k),dom(ib)%vis(i,j,k),dom(ib)%ksgs(i,j,k),
     & dom(ib)%eps(i,j,k),dom(ib)%T(i,j,k),dom(ib)%Tm(i,j,k)
              end do
           end do
        end do
        close (88)

        end do
88      format (15e25.8)

        end
!##########################################################################
        subroutine tecplot_u(ki)
!##########################################################################
        use multidata
        use vars
        implicit none
        integer :: sn,i,j,k,toti,totj,totk
        integer :: is,ie,js,je,ks,ke
        integer :: ki,ib
        character*8 :: chb
        character*25 :: gf
	  double precision :: tau

        do ib=1,nbp

        write(chb,'(i8)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        gf='tecout_u'//trim(adjustl(chb))//'.plt'
        open (unit=88, file=gf)

        is=dom(ib)%isu; ie=dom(ib)%ieu
        js=dom(ib)%jsu; je=dom(ib)%jeu
        ks=dom(ib)%ksu; ke=dom(ib)%keu
        toti=(ie+1)-(is-1)+1
        totj=(je+1)-(js-1)+1
        totk=(ke+1)-(ks-1)+1

        write (88,*) 'title = ', 'u-grid'
        write (88,*)
     &  'variables="x","y","z","U-velocity","UMean","uuMean","tau"'
        write (88,*)'zone ', ' i=',toti,', ',
     &  ' j=',totj,', k= ',totk,' f=point'

        do k=ks-1,ke+1
           do j=js-1,je+1
              do i=is-1,ie+1

	  tau=0.0
	  if (i.eq.is-1) tau=dom(ib)%tauwe(j,k)
	  if (i.eq.ie) tau=dom(ib)%tauww(j,k)
	  if (j.eq.js-1) tau=dom(ib)%tauws(i,k)
	  if (j.eq.je) tau=dom(ib)%tauwn(i,k)
	  if (k.eq.ks-1) tau=dom(ib)%tauwb(i,j)
	  if (k.eq.ke) tau=dom(ib)%tauwt(i,j)

                 write (88,88) dom(ib)%x(i),dom(ib)%yc(j),
     & dom(ib)%zc(k),dom(ib)%u(i,j,k),
     & dom(ib)%um(i,j,k),dom(ib)%uum(i,j,k),tau
              end do
           end do
        end do
        close (88)


        end do

88      format (10e25.8)

        end
!##########################################################################
        subroutine tecplot_v(ki)
!##########################################################################
        use multidata
        implicit none
        integer :: sn,i,j,k,toti,totj,totk
        integer :: is,ie,js,je,ks,ke
        integer :: ki,ib
        character*8 :: chb
        character*25 :: gf

        do ib=1,nbp

        write(chb,'(i8)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        gf='tecout_v'//trim(adjustl(chb))//'.plt'
        open (unit=88, file=gf)

        is=dom(ib)%isv; ie=dom(ib)%iev
        js=dom(ib)%jsv; je=dom(ib)%jev
        ks=dom(ib)%ksv; ke=dom(ib)%kev
        toti=(ie+1)-(is-1)+1
        totj=(je+1)-(js-1)+1
        totk=(ke+1)-(ks-1)+1

        write (88,*) 'title = ', 'v-grid'
        write (88,*)
     &  'variables="x","y","z","V-velocity","VMean","vvMean"'
        write (88,*)'zone ', ' i=',toti,', ',
     &  ' j=',totj,', k= ',totk,' f=point'

        do k=ks-1,ke+1
           do j=js-1,je+1
              do i=is-1,ie+1
                 write (88,88) dom(ib)%xc(i),dom(ib)%y(j),
     & dom(ib)%zc(k),dom(ib)%v(i,j,k),
     & dom(ib)%vm(i,j,k),dom(ib)%vvm(i,j,k)
              end do
           end do
        end do
        close (88)

        end do

88      format (10e25.8)

        end
!##########################################################################
        subroutine tecplot_w(ki)
!##########################################################################
        use multidata
        implicit none
        integer :: sn,i,j,k,toti,totj,totk
        integer :: is,ie,js,je,ks,ke
        integer :: ki,ib
        character*8 :: chb
        character*25 :: gf

        do ib=1,nbp

        write(chb,'(i8)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        gf='tecout_w'//trim(adjustl(chb))//'.plt'
        open (unit=88, file=gf)

        is=dom(ib)%isw; ie=dom(ib)%iew
        js=dom(ib)%jsw; je=dom(ib)%jew
        ks=dom(ib)%ksw; ke=dom(ib)%kew
        toti=(ie+1)-(is-1)+1
        totj=(je+1)-(js-1)+1
        totk=(ke+1)-(ks-1)+1

        write (88,*) 'title = ', 'w-grid'
        write (88,*)
     &  'variables="x","y","z","W-velocity","WMean","wwMean"'
        write (88,*)'zone ', ' i=',toti,', ',
     &  ' j=',totj,', k= ',totk,' f=point'

        do k=ks-1,ke+1
           do j=js-1,je+1
              do i=is-1,ie+1
                 write (88,88) dom(ib)%xc(i),dom(ib)%yc(j),
     & dom(ib)%z(k),dom(ib)%w(i,j,k),
     & dom(ib)%wm(i,j,k),dom(ib)%wwm(i,j,k)
              end do
           end do
        end do
        close (88)

        end do

88      format (10e25.8)

        end
!##########################################################################
        subroutine tec_turb(ki)
!##########################################################################
        use multidata
        use vars
        implicit none
        integer :: sn,i,j,k,toti,totj,totk
        integer :: is,ie,js,je,ks,ke
        integer :: ki,ib
        double precision u_cn,v_cn,w_cn,um_cn,vm_cn,wm_cn,tau
        double precision uum_cn,vvm_cn,wwm_cn,uvml,uwml,vwml
        character*8 :: chb
        character*25 :: gf

        do ib=1,nbp

        write(chb,'(i8)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        gf='tecturb'//trim(adjustl(chb))//'.plt'
        open (unit=88, file=gf)

        is=pl+1; ie=dom(ib)%ttc_i-pl
        js=pl+1; je=dom(ib)%ttc_j-pl
        ks=pl+1; ke=dom(ib)%ttc_k-pl
        toti=ie-(is-1)+1
        totj=je-(js-1)+1
        totk=ke-(ks-1)+1

        write (88,*) 'title = ', 'turb'
        write (88,*)
     &  'variables="x","y","z","U","V","W"',
     &  ',"UM","VM","WM","uuM","vvM","wwM"',
     &  ',"uvM","uwM","vwM"'
        write (88,*)'zone ', ' i=',toti,', ',
     &  ' j=',totj,', k= ',totk,' f=point'

        do k=ks-1,ke
           do j=js-1,je
              do i=is-1,ie

	  tau=0.0
	  if (i.eq.is-1) tau=dom(ib)%tauwe(j,k)
	  if (i.eq.ie) tau=dom(ib)%tauww(j,k)
	  if (j.eq.js-1) tau=dom(ib)%tauws(i,k)
	  if (j.eq.je) tau=dom(ib)%tauwn(i,k)
	  if (k.eq.ks-1) tau=dom(ib)%tauwb(i,j)
	  if (k.eq.ke) tau=dom(ib)%tauwt(i,j)

                 u_cn  =0.25*(dom(ib)%u(i,j,k)+
     &dom(ib)%u(i,j+1,k)+dom(ib)%u(i,j,k+1)+
     &dom(ib)%u(i,j+1,k+1))
                 um_cn  =0.25*(dom(ib)%um(i,j,k)+
     &dom(ib)%um(i,j+1,k)+dom(ib)%um(i,j,k+1)+
     &dom(ib)%um(i,j+1,k+1))
                 uum_cn  =0.25*(dom(ib)%uum(i,j,k)+
     &dom(ib)%uum(i,j+1,k)+dom(ib)%uum(i,j,k+1)+
     &dom(ib)%uum(i,j+1,k+1))


                 v_cn  =0.25*(dom(ib)%v(i,j,k)+
     &dom(ib)%v(i+1,j,k)+dom(ib)%v(i,j,k+1)+
     &dom(ib)%v(i+1,j,k+1))
                 vm_cn  =0.25*(dom(ib)%vm(i,j,k)+
     &dom(ib)%vm(i+1,j,k)+dom(ib)%vm(i,j,k+1)+
     &dom(ib)%vm(i+1,j,k+1))
                 vvm_cn  =0.25*(dom(ib)%vvm(i,j,k)+
     &dom(ib)%vvm(i+1,j,k)+dom(ib)%vvm(i,j,k+1)+
     &dom(ib)%vvm(i+1,j,k+1))

                 w_cn  =0.25*(dom(ib)%w(i,j,k)+
     &dom(ib)%w(i+1,j,k)+dom(ib)%w(i,j+1,k)+
     &dom(ib)%w(i+1,j+1,k)) 
                 wm_cn  =0.25*(dom(ib)%wm(i,j,k)+
     &dom(ib)%wm(i+1,j,k)+dom(ib)%wm(i,j+1,k)+
     &dom(ib)%wm(i+1,j+1,k)) 
                 wwm_cn  =0.25*(dom(ib)%wwm(i,j,k)+
     &dom(ib)%wwm(i+1,j,k)+dom(ib)%wwm(i,j+1,k)+
     &dom(ib)%wwm(i+1,j+1,k)) 

                 uvml  =0.125*(dom(ib)%uvm(i,j,k)+
     &dom(ib)%uvm(i+1,j,k)    +dom(ib)%uvm(i,j+1,k)+
     &dom(ib)%uvm(i+1,j+1,k)  +dom(ib)%uvm(i,j,k+1)+
     &dom(ib)%uvm(i+1,j,k+1)  +dom(ib)%uvm(i,j+1,k+1)+
     &dom(ib)%uvm(i+1,j+1,k+1))
                 uwml  =0.125*(dom(ib)%uwm(i,j,k)+
     &dom(ib)%uwm(i+1,j,k)    +dom(ib)%uwm(i,j+1,k)+
     &dom(ib)%uwm(i+1,j+1,k)  +dom(ib)%uwm(i,j,k+1)+
     &dom(ib)%uwm(i+1,j,k+1)  +dom(ib)%uwm(i,j+1,k+1)+
     &dom(ib)%uwm(i+1,j+1,k+1))
                 vwml  =0.125*(dom(ib)%vwm(i,j,k)+
     &dom(ib)%vwm(i+1,j,k)    +dom(ib)%vwm(i,j+1,k)+
     &dom(ib)%vwm(i+1,j+1,k)  +dom(ib)%vwm(i,j,k+1)+
     &dom(ib)%vwm(i+1,j,k+1)  +dom(ib)%vwm(i,j+1,k+1)+
     &dom(ib)%vwm(i+1,j+1,k+1))

                 write (88,88) dom(ib)%x(i),dom(ib)%y(j),
     &  dom(ib)%z(k),u_cn,v_cn,w_cn,um_cn,vm_cn,wm_cn,
     &  uum_cn,vvm_cn,wwm_cn,uvml,uwml,vwml!,tau

              end do
           end do
        end do

        close (88)

        end do

88      format (20e25.8)

        end
!##########################################################################
        subroutine tecplot_phi(ki)
!##########################################################################
        use vars
        use multidata
        implicit none

        integer :: sn,sn1,i,j,k,toti,totj,totk
        integer :: is,ie,js,je,ks,ke
        integer :: ki,ib
        character*8 :: chb,chb1
        character*27 :: gf

        do ib=1,nbp

	  if (L_anim_phi) then
          write(chb,'(i8)') dom_id(ib)
          write(chb1,'(i8)') ki
          sn=len(trim(adjustl(chb)))
          sn1=len(trim(adjustl(chb1)))
          chb=repeat('0',(4-sn))//trim(adjustl(chb))
          chb1=repeat('0',(6-sn1))//trim(adjustl(chb1))
          gf='tecout_phi_'//trim(adjustl(chb))//'_'//
     & trim(adjustl(chb1))//'.dat'
	  else if (L_LSMinit) then
          write(chb,'(i8)') dom_id(ib)
          write(chb1,'(i8)') ki
          sn=len(trim(adjustl(chb)))
          chb=repeat('0',(4-sn))//trim(adjustl(chb))
          gf='tecout_phi_'//trim(adjustl(chb))//'_initial.dat'
	  else
          write(chb,'(i8)') dom_id(ib)
          sn=len(trim(adjustl(chb)))
          chb=repeat('0',(4-sn))//trim(adjustl(chb))
          gf='tecout_phi_'//trim(adjustl(chb))//'.dat'
	  endif      

	  open (unit=88, file=gf)

        is=dom(ib)%isp; ie=dom(ib)%iep
        js=dom(ib)%jsp; je=dom(ib)%jep
        ks=dom(ib)%ksp; ke=dom(ib)%kep
        toti=(ie)-(is-1)+1
        totj=(je)-(js-1)+1
        totk=(ke)-(ks-1)+1

        write (88,*) 'title = ', 'phi'
        write (88,*)'variables="x","y","z","phi","phi_reinit","phim",',
     & '"dens","mu"'
	  if (L_anim_phi) then
          write (88,*)'zone ','STRANDID=', 1, 'SOLUTIONTIME=', ctime,
     &    ' i=',toti,', ',' j=',totj,', k= ',totk,
     &    'zonetype=', 'ordered',', DATAPACKING=point'
	  else if (L_LSMinit) then
          write (88,*)'zone ','STRANDID=', 1, 'SOLUTIONTIME=', 0.00000, 
     &    ' i=',toti,', ',' j=',totj,', k= ',totk,
     &    'zonetype=', 'ordered',', DATAPACKING=point'
	  else
          write (88,*)'zone ', ' i=',toti,', ',
     &    ' j=',totj,', k= ',totk,' f=point'
	  endif

        do k=ks-1,ke
           do j=js-1,je
              do i=is-1,ie
                write (88,88) dom(ib)%xc(i),dom(ib)%yc(j),
     & dom(ib)%zc(k),dom(ib)%phi(i,j,k),dom(ib)%phi_reinit(i,j,k),
     & dom(ib)%phim(i,j,k),
     & dom(ib)%dens(i,j,k),dom(ib)%mu(i,j,k)
              end do
           end do
        end do
        close (88)

        end do

88      format (8e25.8)

        end
!##########################################################################
        subroutine tecbin(ki)
!##########################################################################
        use multidata
        use vars
        implicit none
        integer :: sn,i,j,k,toti,totj,totk
        integer :: ki,ib,inind,jnind,knind
        character*8 :: chb
        character*25 :: gf

        do ib=1,nbp

        write(chb,'(i8)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        gf='tecbin'//trim(adjustl(chb))//'.bin'
        open (unit=88, file=gf, form='unformatted')

        toti=dom(ib)%ttc_i
        totj=dom(ib)%ttc_j
        totk=dom(ib)%ttc_k

        write (88) toti,totj,totk
        write (88) pl
!====================================================================
        inind=0; jnind=0; knind=0
        if (dom(ib)%inext.lt.0 .and. dom(ib)%bc_east.ne.5) inind=-1
        if (dom(ib)%jnext.lt.0 .and. dom(ib)%bc_north.ne.5) jnind=-1
        if (dom(ib)%knext.lt.0 .and. dom(ib)%bc_top.ne.5) knind=-1
        write (88) inind,jnind,knind
!====================================================================

        do k=1,totk
           do j=1,totj
              do i=1,toti
        write (88) dom(ib)%x(i),dom(ib)%y(j),dom(ib)%z(k),
     & dom(ib)%p(i,j,k),dom(ib)%pm(i,j,k),dom(ib)%ppm(i,j,k),
     & dom(ib)%vis(i,j,k),
     & dom(ib)%u(i,j,k),dom(ib)%um(i,j,k),dom(ib)%uum(i,j,k),
     & dom(ib)%v(i,j,k),dom(ib)%vm(i,j,k),dom(ib)%vvm(i,j,k),
     & dom(ib)%w(i,j,k),dom(ib)%wm(i,j,k),dom(ib)%wwm(i,j,k),
     & dom(ib)%uvm(i,j,k),dom(ib)%uwm(i,j,k),dom(ib)%vwm(i,j,k),
     & dom(ib)%S(i,j,k),dom(ib)%Sm(i,j,k),
     & dom(ib)%ksgs(i,j,k),dom(ib)%eps(i,j,k),
     & dom(ib)%T(i,j,k),dom(ib)%Tm(i,j,k),dom(ib)%Ttm(i,j,k)
              end do
           end do
        end do
        close (88)

        end do
        end
!##########################################################################
