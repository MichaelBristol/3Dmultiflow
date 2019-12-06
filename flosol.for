!##########################################################################
        subroutine flosol
!##########################################################################
        use vars
        use mpi
        use multidata
	  use vars_pt
        implicit none
        real ( kind =8 )  :: wtimedum,wtime_total,wtime_solver,wtime_ib
	real (kind =8) :: wtime_cd,wtime_lpt
        integer ib,i,j,k,kutta,kuttacond,jtime,ii
        double precision :: alfark(3)

! Set some constants ---------------------
! ..... 3-STEP RUNGE KUTTA
        if (conv_sch.eq.4) then
           alfark(1)=1./3.
           alfark(2)=0.5
           alfark(3)=1.0
           kuttacond=3
! ..... 2-STEP RUNGE KUTTA
        else if (conv_sch.eq.3) then
           alfark(1)=0.5
           alfark(2)=1.0
           alfark(3)=1.0
           kuttacond=2
        end if

        alfabc = 1
        alfapr = 1.0

        if (LRESTART.eq..false.) cnt_pt = 1 

        numfile1=1002; numfile2=1003; numfile3=1004
        if(myrank.eq.0) then
           if (pressureforce) then
              open (unit=numfile1, file='forcn.dat')
              write (numfile1,*)
     & 'variables="time","forcn","qsttp","flwsum"'
           end if
           open (unit=numfile2, file='rms.dat')
           write (numfile2,*)
     & 'variables="iter","time","rmax","dt","Mdef"'
           if (L_LSM) then
              open (unit=numfile3, file='normV.dat')
              write (numfile3,*)
     & 'variables="ntime","normV","steps","ctime","dt"'
	     end if
        end if

        itime_start=ntime+1
	  jtime = 0
	  ireadinlet = 0
	  iaddinlet = 1

!	if(LSCALAR)   call sediment_init

	if (myrank.eq.0) then
	 open(unit=203, file='worktime.dat')
	 write(203,*)'Variables=it,C-D,PSolver,IBM,LPT,Total'
	endif

!================== Start Time Loop ====================================
        do itime=itime_start,itime_end

	   if (myrank.eq.1) wtime_total = MPI_WTIME ( )

           if (L_dt)	call checkdt

           ctime=ctime+dt
           ntime = ntime + 1

!----------------------reading inflow data------------------------------!brunho2014
	   if (read_inflow.eq..true.) then
           		if (ireadinlet.eq.ITMAX_PI.and.iaddinlet.eq.1) then
				iaddinlet=-1
          		elseif (ireadinlet.eq.1.and.iaddinlet.eq.-1) then
				iaddinlet=1
	     		endif
          		ireadinlet=ireadinlet+iaddinlet
!---------------------------reading SEM---------------------------------!Pablo2015
	   elseif ((bc_w.eq.8)) then 
           		if (ireadinlet.eq.ITMAX_SEM.and.iaddinlet.eq.1) then 
				iaddinlet=-1 
          		elseif (ireadinlet.eq.1.and.iaddinlet.eq.-1) then 
				iaddinlet=1 
	     		endif 
          		ireadinlet=ireadinlet+iaddinlet 
	   endif
!-----------------------------------------------------------------------

           if (LENERGY) call boundT
           if (LSCALAR) call boundS

           do ib=1,nbp
              do k=1,dom(ib)%ttc_k
              do i=1,dom(ib)%ttc_i
              do j=1,dom(ib)%ttc_j
                 dom(ib)%uoo(i,j,k)=dom(ib)%u(i,j,k)
                 dom(ib)%voo(i,j,k)=dom(ib)%v(i,j,k)
                 dom(ib)%woo(i,j,k)=dom(ib)%w(i,j,k)
                 dom(ib)%To(i,j,k)=dom(ib)%T(i,j,k)
	     	     if (LSCALAR) 	dom(ib)%So(i,j,k)=dom(ib)%S(i,j,k)
              end do
              end do
              end do
           end do

           if (LENERGY) call energy
           if (LSCALAR) call sediment_4thtest
           if (L_LSM)  then
             call LSM_3D
             call heaviside  
           end if

           if (LIMB)  call IB_previous						!Pablo2015    
            
	     if (LPT) then								!Brunho2013
		if (myrank.eq.0) then
			if (mod(itime,tsnr).eq.0)	call release_pt			
		 	call alloc_pt
		endif
		call MPI_pt					
	     endif

           if(SGS) then
              if(sgs_model.eq.1) then
                 call eddyv_smag
              else if(sgs_model.eq.2) then
                 call eddyv_wale
              else if(sgs_model.eq.3) then
                 call eddyv_1eqn
              else if(sgs_model.eq.4) then
                 call eddyv_keps
              end if
           end if

           if (pressureforce) call pressure_forcing

	   if (myrank.eq.0) wtime_cd = MPI_WTIME ( )

           if(conv_sch.eq.3 .or. conv_sch.eq.4) then
              do kutta=1,kuttacond
                 if (kutta.eq.1) then
                    alfabc = 1
                 else
                    alfabc = 0
                 endif
                 alfapr=alfark(kutta)

                 select case (differencing)
                    case (1) 
                       call rungek_conv2nd
                    case (2) 
                       call rungek_conv4th
                    case (3)
                       call rungek_convWENO
                 end select
                 if (diff_sch.eq.3) then
                    select case (differencing)
                       case (1) 
                          call rungek_diff2nd
                       case (2) 
                          call rungek_diff4th
                       case (3)
                          call rungek_diff2nd
                    end select
                 	else
                    call diffusion
                 	end if
                 	if (kutta.lt.kuttacond) call calvel
              end do
           else
              call convection
              call diffusion
           end if

	   CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	   if (myrank.eq.0) wtime_cd=MPI_WTIME( )-wtime_cd

	   if (myrank.eq.0) wtime_lpt = MPI_WTIME ( )					!Brunho2013
           if (LPT) then			
        		call exchange(11)
        		call exchange(22)
        		call exchange(33)
        		call exchange(10)
			if (np_loc.gt.0) call particle_tracking			!Procs without particles do not enter
			call final_LPT
	     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	   endif
	   if (myrank.eq.0) wtime_lpt = MPI_WTIME ( ) - wtime_lpt

	     if (myrank.eq.1) wtime_ib = MPI_WTIME ( )
           if (LIMB)  call IBM							!Pablo2015
	     if (myrank.eq.1) wtime_ib = MPI_WTIME ( ) - wtime_ib

           call correctoutflux

	     if (myrank.eq.1) wtime_solver = MPI_WTIME ( )
           if(solver.eq.1) then
              call pressure_1sweep
           else if(solver.eq.2) then
              call newsolv_mg
           end if
	     if (myrank.eq.1) wtime_solver = MPI_WTIME ( ) - wtime_solver

           if (time_averaging) call update_mean

           if (MOD(itime,10).eq.0 .and. myrank.eq.0)  then
              wtimedum = MPI_WTIME ( ) - wtime

              write (6,5000) itime,ctime,dt,dtavg
              write (6,*) ' '
              write (6,5500) wtimedum
              write (6,*) ' '

              write (numfile,5000) itime,ctime,dt,dtavg
              write (numfile,*) ' '
              write (numfile,5500) wtimedum
              write (numfile,*) ' '
           end if

!           call MPI_BARRIER (MPI_COMM_WORLD,ierr)

!-----------Saving unsteady variables-----------------------------------
	if (ctime.ge.t_start_averaging2) then					!Brunho 2014 
	      jtime = jtime + 1
		do ii=1,n_unstpt
		DO ib=1,nbp
		if (dom_id(ib).eq.id_unst(ii)) then
			dom(ib)%u_unst(ii,jtime)=
     &		dom(ib)%u(i_unst(ii),j_unst(ii),k_unst(ii))
			dom(ib)%v_unst(ii,jtime)=
     &		dom(ib)%v(i_unst(ii),j_unst(ii),k_unst(ii))
			dom(ib)%w_unst(ii,jtime)=
     &		dom(ib)%w(i_unst(ii),j_unst(ii),k_unst(ii))
			dom(ib)%um_unst(ii,jtime)=
     &		dom(ib)%um(i_unst(ii),j_unst(ii),k_unst(ii))
			dom(ib)%vm_unst(ii,jtime)=
     &		dom(ib)%vm(i_unst(ii),j_unst(ii),k_unst(ii))
			dom(ib)%wm_unst(ii,jtime)=
     &		dom(ib)%wm(i_unst(ii),j_unst(ii),k_unst(ii))
			dom(ib)%p_unst(ii,jtime)=
     &		dom(ib)%p(i_unst(ii),j_unst(ii),k_unst(ii))
			dom(ib)%pm_unst(ii,jtime)=
     &		dom(ib)%pm(i_unst(ii),j_unst(ii),k_unst(ii))
			dom(ib)%ksgs_unst(ii,jtime)=
     &		dom(ib)%ksgs(i_unst(ii),j_unst(ii),k_unst(ii))
			dom(ib)%eps_unst(ii,jtime)=
     &		dom(ib)%eps(i_unst(ii),j_unst(ii),k_unst(ii))
		endif
		ENDDO
		enddo	
!-----------Saving inflow-----------------------------------------------!Brunho2014
	DO ib=1,nbp				
	if ((save_inflow.eq..true.).and.(mod(dom_id(ib),idom).eq.0)) then	
			call write_inflow(ib)
	endif
	ENDDO
	endif


!------write solution in tecplot format------------------------------
           if ((mod(itime,n_out).eq.0).and.(itime.ge.itime_start)) then
!              call tecgrid(itime)
!              call tec_turb(itime)
!              call tecplot_p(itime)
!              call tecplot_u(itime)
!              call tecplot_v(itime)
!              call tecplot_w(itime)
      	  ! call TECPLOT(itime)
		  if (L_LSM) call tecplot_phi(itime)
              call tecbin(itime)
      	     	!call TECPLOT(itime)
		  if (ctime.ge.t_start_averaging2)  call timesig
              open (unit=101, file='final_ctime.dat')
              if(myrank.eq.0) write (101,'(i8,3F15.6)') 
     &   ntime,ctime,forcn,qstpn
              close(101)
		  if (LPT) then
        		if (myrank.eq.0) then          			
		    		open(30,file='final_particle.dat')
		   	 	write(30,*) np
		    		write(30,*) cnt_pt
		    		do k=1,np
		      		write(30,*) xp_pt(k),yp_pt(k),zp_pt(k)
		      		write(30,*) uop_pt(k),vop_pt(k),wop_pt(k)  
		 			write(30,*) dp_pt(k)
		    		end do
		    		close(30)
        		endif
		  endif
        end if
        if (LPT) then									!Brunho 2013 
		if ((mod(itime,tsteps_pt).eq.0).and.
     &		(itime.gt.itime_start)) then 			
 	           	if (myrank.eq.0) call TECPARTICLE(cnt_pt)
      	     	call TECPLOT(cnt_pt)
			cnt_pt = cnt_pt + 1
		endif
        end if
!-----------------------------------------------------------------------
	  if (myrank.eq.1) then
		wtime_total = MPI_WTIME ( ) - wtime_total
	    	write(203,*) 'solver',wtime_solver
		write(203,*) 'ibm',wtime_ib
		write(203,*) 'total',wtime_total
	  endif

        end do
! End Time Loop ---------------------------

        itime = itime - 1

	  close(203)

        if (mod(itime,n_out).ne.0) then
!              call tecgrid(itime)
!              call tec_turb(itime)
!              call tecplot_p(itime)
!              call tecplot_u(itime)
!              call tecplot_v(itime)
!              call tecplot_w(itime)
!              if (LSCALAR) call tecplot_S(itime)
              call tecbin(itime)
		  if (L_LSM) call tecplot_phi(itime)
		  if (ctime.ge.t_start_averaging2)  call timesig
           open (unit=101, file='final_ctime.dat')
           if(myrank.eq.0) write (101,'(i8,3F15.6)') 
     &   ntime,ctime,forcn,qstpn
           close(101)
        end if

	if (myrank.eq.0) write (6,*) 'ctime=' , ctime
	if (myrank.eq.0) write (numfile,*) 'ctime=' , ctime

5000  format(/1x,10(1h=),' nrtstp=',i8,2x,'ctime=',e14.6,2x,
     & 'dt=',e14.6,'  dtavg=',e14.6)
5500  format(/1x,'Work took ',e18.8,2x,' seconds')
        end subroutine flosol
!##########################################################################
