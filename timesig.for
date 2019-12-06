!##########################################################################
        subroutine timesig
!			Bruño Fraga Bugallo
!			Cardiff 2014
!##########################################################################
        use multidata
        use vars
	  use mpi
        implicit none
        integer :: i,j,idfile,ib,jtime
        character*25 :: unst
        character*8 :: numpt,x,y,z

		jtime = ntime-itime_start

	do i=1,n_unstpt
	do ib=1,nbp
		if (dom_id(ib).eq.id_unst(i)) then
			!crear nombre
			write(numpt,'(i2)') i
      		unst='unst_'//trim(adjustl(numpt))//'.dat'
			!abrir archivo
			idfile=499+i

			write(x,'(F6.3)')dom(ib)%x(i_unst(i))
			write(y,'(F6.3)')dom(ib)%y(j_unst(i))
			write(z,'(F6.3)')dom(ib)%z(k_unst(i))

      		open (unit=idfile, file=unst)

			!escribir

			write(idfile,*)'ZoneT = "Time series:',x,y,z,'"'
!		      write(idfile,*)'Zone= "Time step= ',itime,'"'
			write(idfile,*)'Variables = "itime","u","v","w"'
     &		,',"umean","vmean","wmean","p","pmean","k","eps"'
			do j=1,jtime
				write(idfile,*)j,dom(ib)%u_unst(i,j),
     &			dom(ib)%v_unst(i,j),dom(ib)%w_unst(i,j),
     &			dom(ib)%um_unst(i,j),dom(ib)%vm_unst(i,j),
     &			dom(ib)%wm_unst(i,j),dom(ib)%p_unst(i,j),
     &			dom(ib)%pm_unst(i,j),
     &			dom(ib)%ksgs_unst(i,j),dom(ib)%eps_unst(i,j)
			enddo
        		close (idfile)
		endif
	enddo
	enddo

!  500  format(A,F6.3,A)

	end subroutine
!##########################################################################
        subroutine write_inflow(i)
!			Bruño Fraga Bugallo
!			Cardiff 2014
!##########################################################################

        use vars
        use mpi
        use multidata
        implicit none

        integer i,j,k,jtime,kutta
        integer ib,l,strlen
	  character (LEN=19) :: filename
	  character (LEN=5) :: name_end
	  character (LEN=3) :: dominio

			if (ireadinlet.eq.10001) ireadinlet=0
            	ireadinlet=ireadinlet+1

        		write(name_end,'(I5)') ireadinlet
        		strlen=LEN(TRIM(ADJUSTL(name_end)))
       		name_end=REPEAT('0',(5-strlen))//
     &			TRIM(ADJUSTL(name_end)) ! e.g. "00001"

			write(dominio,'(I3)') dom_id(i)
        		strlen=LEN(TRIM(ADJUSTL(dominio)))
       		dominio=REPEAT('0',(3-strlen))//
     &			TRIM(ADJUSTL(dominio)) ! e.g. "001"

			filename='Inlet_'//dominio//'_'//name_end//'.dat'

            	open (unit=ireadinlet+500, file=filename)
			do k=dom(i)%ksu-1,dom(i)%keu+1
			do j=dom(i)%jsu-1,dom(i)%jeu+1
				write(ireadinlet+500,*)dom(i)%u(29,j,k)
			enddo
			enddo		
			do k=dom(i)%ksv-1,dom(i)%kev+1
			do j=dom(i)%jsv-1,dom(i)%jev+1
				write(ireadinlet+500,*)dom(i)%v(29,j,k)
			enddo
			enddo	
			do k=dom(i)%ksw-1,dom(i)%kew+1
			do j=dom(i)%jsw-1,dom(i)%jew+1
				write(ireadinlet+500,*)dom(i)%w(29,j,k)
			enddo
			enddo		
			close(ireadinlet+500)
	end subroutine
