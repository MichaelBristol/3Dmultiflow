!######################################################################
      module imb
!######################################################################
	  SAVE
      double precision  :: xt(5),yt(5),xdt(5),ydt(5),xddt(5),yddt(5),yto
	double precision  :: lambda,sigma,nxl
	
	INTEGER :: imp_proc_master,imb_block_master,forcefilej,a,b
        INTEGER :: bodynum,maxnode,master,maxnodeIBS,mdfsteps,yangcase	
	
	LOGICAL,allocatable,dimension (:):: rotating,LSELFST,intflow
	
      double precision,allocatable ::Cx(:),Cxor(:),Cy(:),Cyor(:)
      double precision,allocatable ::Cz(:),Czor(:),Czz(:),Cxx(:)	
      double precision,allocatable ::l2norm(:)				
	double precision,allocatable,dimension(:) :: xaero,yaero,zaero
	double precision,allocatable,dimension(:) :: iniT_selfST
	double precision,allocatable,dimension(:) :: SUMtorque_ST	
	double precision,allocatable,dimension(:) :: acc_selfST,massI	
	double precision,allocatable,dimension(:) :: radsin,rads,pitch
	double precision,allocatable,dimension(:) :: R,reddelta
	double precision,allocatable,dimension(:) :: FX1_MASTER,FX2_MASTER
	double precision,allocatable,dimension(:) :: FX3_MASTER,FX3_loc
        double precision,allocatable,dimension(:) :: nodex_loc,nodey_loc
        double precision,allocatable,dimension(:) :: nodez_loc
        double precision,allocatable,dimension(:) :: FX1_loc,FX2_loc
        double precision,allocatable,dimension(:) :: R0_loc,alpha0_loc
        double precision,allocatable,dimension(:) :: U_Beta1_loc
        double precision,allocatable,dimension(:) :: U_Beta2_loc,zend
        double precision,allocatable,dimension(:) :: U_Beta3_loc,zini
        double precision,allocatable,dimension(:,:) :: dh1_loc,dh2_loc
        double precision,allocatable,dimension(:,:) :: dh3_loc,delvol
        double precision,allocatable,dimension(:,:) :: nodexlocal
        double precision,allocatable,dimension(:,:) :: nodeylocal
        double precision,allocatable,dimension(:,:) :: nodezlocal
        double precision,allocatable,dimension(:,:) :: FX1NF,FX2NF,FX3NF
        double precision,allocatable,dimension(:,:) :: nodex,nodey,nodez
        double precision,allocatable,dimension(:,:) :: torque,U_Beta3 
        double precision,allocatable,dimension(:,:) :: U_Beta1,U_Beta2
        double precision,allocatable,dimension(:,:) :: FX1,FX2,FX3      
        double precision,allocatable,dimension(:,:) :: alpha0,R0        
        double precision,allocatable,dimension(:) :: dxm,dym,dzm
	INTEGER,allocatable,dimension(:,:) :: I_nr_V,J_nr_V,K_nr_V
	INTEGER,allocatable,dimension(:,:) :: I_nr_U,J_nr_U,K_nr_U
	INTEGER,allocatable,dimension(:,:) :: I_nr_W,J_nr_W,K_nr_W
	INTEGER,allocatable,dimension(:) :: imbnumber,imb_shape,IBMnum
	INTEGER,allocatable,dimension(:) :: kmaxU,kmaxV,kmaxW,nodes
	INTEGER,allocatable,dimension(:) :: imb_proc,imb_block,turax	
	INTEGER,allocatable,dimension(:) :: lag_bod_loc,cmax,linfin
        INTEGER,allocatable,dimension(:) :: domtemp,imb_block_loc,axis
        INTEGER,allocatable,dimension(:) :: imbinblock_loc,rott_loc
        double precision,allocatable,dimension(:) :: rdiv_imb
        CHARACTER*32, allocatable, dimension (:) :: filepoints
        
      end module imb
!#############################################################
      SUBROUTINE IMB_INITIAL
!#############################################################
      use vars
      use multidata
      use imb
      use mpi
      implicit none
      double precision         :: PI,revoltime
      INTEGER      :: L,I,strlen,maxn,K
      CHARACTER*8  :: char_block
      CHARACTER*31 :: gridfile
      character*80 :: dummyline

        PI = 4.D0*DATAN(1.D0)
        xt = 0.d0	; xdt = 0.d0	; xddt = 0.d0
	yt = 0.d0	; ydt = 0.d0	; yddt = 0.d0

!        imb_proc = -1  ; imb_block = -1; imb_mastercpu = -1

	master=0 ! 0 is going to be always the master processor

       open (unit=1, file='geom.cin')
       read (1,*) dummyline
       read (1,*) yangcase
       read (1,*) mdfsteps
       read (1,*) bodynum

!Allocate variables that all MPI need to know:
        allocate(nodes(bodynum),rotating(bodynum),rads(bodynum))
	allocate(imb_shape(bodynum),imbnumber(bodynum))
	allocate(radsin(bodynum),LSELFST(bodynum),filepoints(bodynum))
	allocate(turax(bodynum),reddelta(bodynum),rdiv_imb(bodynum))        	
	allocate(dxm(bodynum),dym(bodynum),dzm(bodynum)) 

	IF (myrank.ne.master) GOTO 545
!Allocate variables only needed by the master:
	allocate(Cx(bodynum),Cxor(bodynum),Cy(bodynum),Cyor(bodynum))
	allocate(Cz(bodynum),Czor(bodynum),pitch(bodynum))
	allocate(R(bodynum),l2norm(bodynum),IBMnum(bodynum))	
	allocate(xaero(bodynum),yaero(bodynum),zaero(bodynum))
        allocate(acc_selfST(bodynum),torque(itime_end,bodynum))
        allocate(SUMtorque_ST(bodynum),massI(bodynum))
        allocate(iniT_selfST(bodynum),cmax(bodynum),axis(bodynum))
	allocate(linfin(bodynum),zini(bodynum),zend(bodynum))        
	allocate(intflow(bodynum))        

	maxn=2999000	!maximum number of Lagrangian allowed

         allocate (nodex(bodynum,maxn),nodexlocal(bodynum,maxn))
	 allocate (nodey(bodynum,maxn),nodeylocal(bodynum,maxn))
	 allocate (nodez(bodynum,maxn),nodezlocal(bodynum,maxn))
         allocate (U_Beta1(bodynum,maxn),U_Beta2(bodynum,maxn))
         allocate (U_Beta3(bodynum,maxn))
         allocate (FX1NF(bodynum,maxn),FX1(bodynum,maxn))
         allocate (FX2NF(bodynum,maxn),FX2(bodynum,maxn))
         allocate (FX3NF(bodynum,maxn),FX3(bodynum,maxn))
         allocate (alpha0(bodynum,maxn),R0(bodynum,maxn))
         allocate (delvol(bodynum,maxn))

	 allocate (imb_proc(maxn),imb_block(maxn))
         nodex = 0.d0; nodey = 0.d0 ; nodez = 0.d0  ;maxnodeIBS=0

        i = 1
        DO WHILE (i.le.bodynum)
		read (1,*) dummyline 
		read (1,*) imb_shape(i)
		read (1,*) linfin(i),zini(i),zend(i)		
		read (1,*) Cx(i),Cy(i),Cz(i)
		read (1,*) R(i)
		read (1,*) cmax(i)
		read (1,*) axis(i)
		read (1,*) filepoints(i)
		read (1,*) rotating(i)
!		read (1,*) a,b
		read (1,*) intflow(i)	
		read (1,*) reddelta(i),rdiv_imb(i)		
!--- Turbine parameters: ---
		read (1,*) dummyline 
		read (1,*) turax(i)				
		read (1,*) xaero(i),yaero(i),zaero(i)		
		read (1,*) pitch(i)		
		read (1,*) imbnumber(i)
		read (1,*) radsin(i)		
		read (1,*) LSELFST(i)
		read (1,*) iniT_selfST(i)
		
	   if(imb_shape(i).ne.5) then
	    xaero(i)=0.d0 ; yaero(i)=0.d0 ; zaero(i)=0.d0 		
	    pitch(i)=0.d0 ; imbnumber(i)=1 ; radsin(i)=0.d0  
	   endif	
	   		
	 if(intflow(i) .eq. .true.) then
	 write(6,*)' This is not implemented in the current version'
	 STOP
	 endif
	 
           IBMnum(i) = i 

!Write possible combinations of movements/body types that doesn't work
        i = i + 1
      End do
      close (1)
	WRITE(6,*)' '
	WRITE(6,*)' '
	WRITE(6,*)'=================================================='
	WRITE(6,*)'===========  Immersed Boundary Details  =========='

!CALCULATE TO WHICH BLOCK IS THE CENTRE OF THE BODIES
545	CONTINUE

	call MPI_BARRIER(MPI_COMM_WORLD,ierr)


	IF (myrank.ne.master) RETURN
	
	Do i=1,bodynum	
       dxm(i)=g_dx/rdiv_imb(i)		!Minimum grid sizes
	 dym(i)=g_dy/rdiv_imb(i)	
	 dzm(i)=g_dz/rdiv_imb(i)		
	Enddo
	 
!	write(6,*)	     'Largest rdivmax  :',rdivmax	
!	write(6,'(a,3e12.4)')'Smallest gridsize: ',dxm,dym,dzm	 
		
	Do i=1,bodynum
          IF (imb_shape(i).eq.1) call imb_square(IBMnum(i))
          IF (imb_shape(i).eq.2) call imb_cylinder(IBMnum(i))
          IF (imb_shape(i).eq.3) call imb_cube(IBMnum(i))
	  IF (imb_shape(i).eq.4) call imb_sphere(IBMnum(i))
 	  IF (imb_shape(i).eq.5) call imb_file(IBMnum(i))
!          IF (imb_shape(i).eq.6) call imb_pipe(IBMnum(i))	  
   	  maxnodeIBS=maxnodeIBS+nodes(i)
	     IF (maxnodeIBS.gt.maxn) write(6,*)'Too many ib points'
	     IF (maxnodeIBS.gt.maxn) STOP
	Enddo

	call imb_alpha0	!-----> CHECK!!!!

!        call calpostn(0.d0,0.d0) !-----> CHECK!!!!
               
!      call imb_deternds	! Moved to flosol.for

!CREATE OUTPUT FILES FOR THE IMMERSED BOUNDARIES:

	L=0
	Do K=1,bodynum
	 if(imb_shape(K).ne.5) imbnumber(K)=1
	   Do i=1,imbnumber(K)
             	 L=L+1 ; forcefilej=399+L
	 IF (rotating(K) .AND. imb_shape(K).eq.5) then !Rotating VATT

	if(K.eq.1 .and. i.eq.1) then
	WRITE(6,*)' '
	WRITE(6,*)'=========== Rotating Parameters  ========= '
	endif

         write(char_block,'(I8)') L
         strlen=LEN(TRIM(ADJUSTL(char_block)))
         char_block=REPEAT('0',(3-strlen))//TRIM(ADJUSTL(char_block))
         gridfile='F_Blade_'//TRIM(ADJUSTL(char_block))//'.dat'
         open (unit=forcefilej, file=gridfile, status="unknown",
     &	action="write")
        write (forcefilej,*)'Variables="TIME","Deg","Fx","Fy","Fz"'
		IF(imb_shape(K).eq.5 .and. i.eq.1)then
			lambda=radsin(K)*R(K)/1.  	
			sigma=imbnumber(K)*1.d0/(R(K)*2*3.1416)
			revoltime=2.d0*PI/radsin(K)	
			write(6,'(a,i1)')   '        Turbine  ',K,''
			write(6,'(a,f12.3)')'        TSR     :',lambda
			write(6,'(a,f12.3)')'        Solidity:',sigma
			WRITE(6,*)' '
		ENDIF
	ELSE
         write(char_block,'(I8)') L
         strlen=LEN(TRIM(ADJUSTL(char_block)))
         char_block=REPEAT('0',(3-strlen))//TRIM(ADJUSTL(char_block))
         if(imb_shape(K).eq.1) then
         gridfile='F_Squ_'//TRIM(ADJUSTL(char_block))//'.dat'
     	 endif        
         if(imb_shape(K).eq.2) then
         gridfile='F_Cyl_'//TRIM(ADJUSTL(char_block))//'.dat'
     	 endif
         if(imb_shape(K).eq.3) then
         gridfile='F_Cub_'//TRIM(ADJUSTL(char_block))//'.dat'
     	 endif
         if(imb_shape(K).eq.4) then
         gridfile='F_Sph_'//TRIM(ADJUSTL(char_block))//'.dat'
     	 endif     	      
         if(imb_shape(K).eq.5) then
         gridfile='F_Body_'//TRIM(ADJUSTL(char_block))//'.dat'
     	 endif       
     
         open (unit=forcefilej, file=gridfile)
           write (forcefilej,*)'Variables="TIME","Fx","Fy","Fz"'
	 ENDIF
	Enddo !i

! IF SELF STARTING IS INTRODUCED IN THE CODE:
!	   if (LSELFST(K) .and. L.eq.1) then
!             open (unit=selfstarting, file=SelfStartingResults.dat')
!             write(selfstarting,*)'Variables="Accel","Veloc","Displ"'
!		radsin(K)=0.d0 ;  acc_selfST(K)=0.d0; rads(K)= 0.d0
!	   endif
	ENDDO !M

	open (unit=757, file= 'l2norm.dat')
	write(757,*)'Variables="Time","l2-norm","l1norm"'
	WRITE(6,*)' '
	WRITE(6,*)'Total # of IB POINTS.........',maxnodeIBS
	WRITE(6,*)' '
	WRITE(6,*)'=================================================='
!	write(6,*)'        Mesh sizes at the IB domain           '
!	write(6,'(a,f12.4,a,f12.4)')'    dx:',dx,'      dy:', dy
!	WRITE(6,*)' '


      RETURN
      end
!#############################################################
      SUBROUTINE imb_alpha0	!STILL NEED TO CHECK IT
!#############################################################
      use vars
      use multidata
      use imb
      use mpi
      implicit none
      INTEGER  :: M,L,iii,K
      double precision :: PI
       PI = 4.D0*DATAN(1.D0)
      
	Do M=1,bodynum
	 if(imb_shape(M).ne.5) then  !all but VATT

         do L=1,nodes(M) 
          alpha0(M,L)=atan((nodex(M,L)-Cx(M))/(nodey(M,L)-Cy(M)))
          R0(M,L)=sqrt((nodex(M,L)-Cx(M))**2+(nodey(M,L)-Cy(M))**2)
         enddo

	 else

	IF (turax(M).eq.1) then		! Vertical Axis Turbine
	   K=nodes(M)/imbnumber(M)
          do L=1,K!nodes(numIB) 
          alpha0(M,L)=atan(nodexlocal(M,L)/(nodeylocal(M,L)+R(M)))
          R0(M,L)=sqrt((nodexlocal(M,L))**2+(nodeylocal(M,L)+R(M))**2)
          enddo
          Do iii=1,imbnumber(M)-1
           do L=1,K!nodes(numIB)
                alpha0(M,L+K*iii)=alpha0(M,L)	
                R0(M,L+K*iii)= R0(M,L)
           enddo
          Enddo
        ENDIF
	IF (turax(M).eq.2) then	! Horizontal Axis Turbine
          do L=1,nodes(M) 
          alpha0(M,L)=atan(nodeylocal(M,L)/(nodezlocal(M,L)))
          
          if(nodeylocal(M,L).gt.0.d0 .and. nodezlocal(M,L).lt.0.d0) then
          alpha0(M,L)=PI+alpha0(M,L)
          endif
          if(nodeylocal(M,L).lt.0.d0 .and. nodezlocal(M,L).lt.0.d0) then
          alpha0(M,L)=PI+alpha0(M,L)
          endif
          if(nodeylocal(M,L).lt.0.d0 .and. nodezlocal(M,L).gt.0.d0) then
          alpha0(M,L)=2.D0*PI+alpha0(M,L)
          endif 
                                    
          R0(M,L)=sqrt((nodeylocal(M,L))**2+(nodezlocal(M,L))**2)
     
          enddo
        ENDIF
                 
	 endif
	Enddo

   88 FORMAT (i5)
   89 FORMAT (2e25.18)
        RETURN
        END SUBROUTINE
!######################################################################
      SUBROUTINE PartLocMPI
!######################################################################
      use vars
      use imb
      use mpi
      use multidata
      implicit none

       call MPI_BCAST(maxnodeIBS,1,MPI_INTEGER,
     &  master,MPI_COMM_WORLD,ierr)	!Total # IB points
       call MPI_BCAST(nodes,bodynum,MPI_INTEGER,
     &  master,MPI_COMM_WORLD,ierr)	!# IB points of each body
       call MPI_BCAST(reddelta,bodynum,MPI_DOUBLE_PRECISION,
     &  master,MPI_COMM_WORLD,ierr) !Reduction factor
       call MPI_BCAST(imb_shape,bodynum,MPI_INTEGER,
     &  master,MPI_COMM_WORLD,ierr)	!IB shape of each body
       call MPI_BCAST(imbnumber,bodynum,MPI_INTEGER,
     &  master,MPI_COMM_WORLD,ierr)	!# IB bodies of each body
       call MPI_BCAST(turax,bodynum,MPI_INTEGER,
     &  master,MPI_COMM_WORLD,ierr)	!# axis of rotation for turbines     
       call MPI_BCAST(radsin,bodynum,MPI_DOUBLE_PRECISION,
     &  master,MPI_COMM_WORLD,ierr)	!Rotational velocity of each IB body
       call MPI_BCAST(rotating,bodynum,MPI_LOGICAL,
     &  master,MPI_COMM_WORLD,ierr)	!If the body rotates
     
        if(yangcase.eq.1) nxl=1.4999d0
	if(yangcase.eq.2) nxl=2.4999d0
	if(yangcase.eq.3) nxl=1.9999d0	
	if(yangcase.eq.4) nxl=2.4999d0
	if(yangcase.eq.5) nxl=1.4999d0		         
	if(yangcase.eq.6) nxl=1.9999d0		!June 2015   
!nxl is the length of the kernel used for the delta functions.		
	  
         allocate(kmaxU(maxnodeIBS),kmaxV(maxnodeIBS),kmaxW(maxnodeIBS))
	  kmaxU=0 	; kmaxV=0 	; kmaxW=0	         		  
         allocate (nodex_loc(maxnodeIBS),nodey_loc(maxnodeIBS))
         allocate (nodez_loc(maxnodeIBS))
	 nodex_loc=0.d0 ; nodey_loc=0.d0 ; nodez_loc=0.d0 ;
	 allocate (U_Beta1_loc(maxnodeIBS),U_Beta2_loc(maxnodeIBS))
	 allocate (U_Beta3_loc(maxnodeIBS))
	  U_Beta1_loc =0.d0 	; U_Beta2_loc =0.d0 ; U_Beta3_loc =0.d0 
         allocate (FX1_loc(maxnodeIBS),FX1_MASTER(maxnodeIBS))
         allocate (FX2_loc(maxnodeIBS),FX2_MASTER(maxnodeIBS))
	 allocate (FX3_loc(maxnodeIBS),FX3_MASTER(maxnodeIBS))
	  FX1_loc =0.d0    ; FX2_loc =0.d0    ; FX3_loc =0.d0 
	  FX1_MASTER =0.d0 ; FX2_MASTER =0.d0 ; FX3_MASTER =0.d0 
         allocate (alpha0_loc(maxnodeIBS),R0_loc(maxnodeIBS))
	  alpha0_loc =0.d0 ; R0_loc =0.d0
	 allocate (imbinblock_loc(num_domains)) 
         allocate (rott_loc(maxnodeIBS))
	 allocate (imb_block_loc(maxnodeIBS),lag_bod_loc(maxnodeIBS))	
	  imbinblock_loc=0 ; imb_block_loc=0  ; lag_bod_loc = 0
	   rott_loc=0   
           
!        if (nbp.ge.1) then  !Only one block to procs where there are IBs!!
!        if (myrank.lt.40) then
        
           
	if (yangcase.eq.2 .or. yangcase.eq.4) then
         allocate (dh1_loc(maxnodeIBS,126),dh2_loc(maxnodeIBS,126))
         allocate (dh3_loc(maxnodeIBS,126))
	 allocate (I_nr_U(maxnodeIBS,126),J_nr_U(maxnodeIBS,126))
	 allocate (I_nr_V(maxnodeIBS,126),J_nr_V(maxnodeIBS,126))
	 allocate (I_nr_W(maxnodeIBS,126),J_nr_W(maxnodeIBS,126))
	 allocate (K_nr_U(maxnodeIBS,126),K_nr_V(maxnodeIBS,126))
	 allocate (K_nr_W(maxnodeIBS,126))
	endif
	if (yangcase.eq.1 .or. yangcase.eq.5) then
         allocate (dh1_loc(maxnodeIBS,28),dh2_loc(maxnodeIBS,28))
         allocate (dh3_loc(maxnodeIBS,28))
	 allocate (I_nr_U(maxnodeIBS,28),J_nr_U(maxnodeIBS,28))
	 allocate (I_nr_V(maxnodeIBS,28),J_nr_V(maxnodeIBS,28))
	 allocate (I_nr_W(maxnodeIBS,28),J_nr_W(maxnodeIBS,28))
	 allocate (K_nr_U(maxnodeIBS,28),K_nr_V(maxnodeIBS,28))
	 allocate (K_nr_W(maxnodeIBS,28))
	endif
	if (yangcase.eq.3 .or. yangcase.eq.6) then
         allocate (dh1_loc(maxnodeIBS,65),dh2_loc(maxnodeIBS,65))
         allocate (dh3_loc(maxnodeIBS,65))
	 allocate (I_nr_U(maxnodeIBS,65),J_nr_U(maxnodeIBS,65))
	 allocate (I_nr_V(maxnodeIBS,65),J_nr_V(maxnodeIBS,65))
	 allocate (I_nr_W(maxnodeIBS,65),J_nr_W(maxnodeIBS,65))
	 allocate (K_nr_U(maxnodeIBS,65),K_nr_V(maxnodeIBS,65))
	 allocate (K_nr_W(maxnodeIBS,65))
	endif		

	  dh1_loc=0.d0	; dh2_loc=0.d0	; dh3_loc=0.d0
	  I_nr_U=0	; J_nr_u=0 	; K_nr_U=0
	  I_nr_V=0 	; J_nr_V=0 	; K_nr_V=0
	  I_nr_W=0 	; J_nr_W=0 	; K_nr_W=0
!	endif           
                                         

	RETURN
	END
!######################################################################
      SUBROUTINE IB_previous
!######################################################################
      use vars
      use mpi
      use multidata
      use imb
      implicit none
      INTEGER :: K

	Do K=1,bodynum
	  IF (rotating(K).and.imb_shape(K).eq.5) call imb_moved(K)  !In shapes.for
	Enddo
	
		Call PartLoc  

	IF(itime.eq.itime_start) then
	 if(myrank.eq.master)write(6,*)'Delta functions initiating'
		Call Deltah
	 call MPI_BARRIER(MPI_COMM_WORLD,IERR)	
	 if(myrank.eq.master)write(6,*)'Delta functions generated'
	ENDIF

      END SUBROUTINE
!######################################################################
      SUBROUTINE PartLoc 
!######################################################################
      use vars
      use imb
      use mpi
      use multidata
      implicit none
      INTEGER :: M,L,ii,nxdom,nydom,nzdom,tnm,N,nx,ny,nz
      DOUBLE PRECISION :: lxdom(idom+1),lydom(jdom+1),lzdom(kdom+1)

       IF(myrank.eq.master)THEN
		lxdom=0 ; lydom=0 ; lzdom=0
		do N=2,idom+1
		 lxdom(N)=(xcor(N-2,2)-xcor(N-2,1))+lxdom(N-1)
		enddo
		do N=2,jdom+1
		 lydom(N)=(ycor((N-2)*idom,2)-ycor((N-2)*idom,1))+lydom(N-1)
		enddo
		do N=2,kdom+1
		lzdom(N)=
     & (zcor((N-2)*idom*jdom,2)-zcor((N-2)*idom*jdom,1))+lzdom(N-1)
		enddo
		imbinblk=0   !# Points in each block
		imb_block=0   !Block id to which every particle belongs
		ii=0  ; tnm=0
		do M=1,bodynum   !Perform this operation to all IB bodies
		   DO L=1,nodes(M) !Analyze all IB poins of the body.
		ii=ii+1 !; nxdom=0 ; nydom=0 ; nzdom=0

		 Do nx=1,idom 
		  if( (nodex(M,L)-1.d-11).gt.lxdom(nx) .and. 
     & (nodex(M,L)-1.d-11).le.lxdom(nx+1) )THEN
		nxdom=nx-1
		GOTO 490
		  endif
		 Enddo
490 	CONTINUE
		 Do ny=1,jdom 
		  if( (nodey(M,L)-1.d-11).gt.lydom(ny) .and. 
     & (nodey(M,L)-1.d-11).le.lydom(ny+1) )THEN
		nydom=ny-1
		GOTO 491
		  endif
		 Enddo
491 	CONTINUE
		 Do nz=1,kdom 
		  if( (nodez(M,L)-1.d-11).gt.lzdom(nz) .and. 
     & (nodez(M,L)-1.d-11).le.lzdom(nz+1) )THEN
		nzdom=nz-1
		GOTO 492
		  endif
		 Enddo
492 	CONTINUE
		imb_block(ii)=idom*jdom*nzdom+idom*nydom+nxdom
		imbinblk(imb_block(ii)+1)=imbinblk(imb_block(ii)+1)+1
		   ENDDO
		enddo 

		do L=1,num_domains !Check in all the domains
		tnm=tnm+imbinblk(L)
		IF (itime.eq.itime_start .AND. imbinblk(L).ne.0)  
     &   write(6,*)'Dom,#markrs',L-1,imbinblk(L),tnm  
		imbinblock_loc(L)=imbinblk(L) !New variable for all the other MPI
		enddo
		!Warning if some point is not assigned to some domain
		if(tnm.lt.maxnodeIBS) 
     &  write(6,*)'Some Lagrangian are not assigned to a domain!!!CHECK'
		 
		  ii=0
		DO M=1,bodynum
		   DO L=1,nodes(M)
		     ii=ii+1
		  nodex_loc(ii)=nodex(M,L) ; nodey_loc(ii)=nodey(M,L)
		  nodez_loc(ii)=nodez(M,L) ; imb_block_loc(ii)=imb_block(ii)

		    IF (itime.eq.itime_start) then !THIS IS DONE ONCE
		R0_loc(ii)=R0(M,L) 
		alpha0_loc(ii)=alpha0(M,L)
		lag_bod_loc(ii)=M 
		    ENDIF 
		rott_loc(ii)=1 !Moving Lagrangian
		IF(rotating(M).eq..false.)rott_loc(ii)=2   !Static Lagrangian
		 ENDDO
		Enddo
		ENDIF !master
		 
		IF(itime.eq.itime_start) then
		call MPI_BCAST(lag_bod_loc,maxnodeIBS,MPI_INTEGER,
     & master,MPI_COMM_WORLD,ierr) !# of the body to which the Lag is.
		call MPI_BCAST(alpha0_loc,maxnodeIBS,MPI_DOUBLE_PRECISION,
     & master,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(R0_loc,maxnodeIBS,MPI_DOUBLE_PRECISION,
     & master,MPI_COMM_WORLD,ierr)
		ENDIF

        call MPI_BCAST(rott_loc,maxnodeIBS,MPI_INTEGER,
     & master,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(imbinblock_loc,num_domains,MPI_INTEGER,
     & master,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(imb_block_loc,maxnodeIBS,MPI_INTEGER,
     & master,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(nodex_loc,maxnodeIBS,MPI_DOUBLE_PRECISION,
     & master,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(nodey_loc,maxnodeIBS,MPI_DOUBLE_PRECISION,
     & master,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(nodez_loc,maxnodeIBS,MPI_DOUBLE_PRECISION,
     & master,MPI_COMM_WORLD,ierr)

      RETURN
      END
!######################################################################
      SUBROUTINE Deltah
!######################################################################
      use vars
      use mpi
      use multidata
      use imb
      implicit none
      double precision :: dh,dhtotal
      INTEGER :: I,J,L,ib,nl,K

!	IF(nmls.eq.0) then

	Do ib=1,nbp  !Loop through all the blocks of one processor

       if (imbinblock_loc(dom_id(ib)+1).eq.0) GOTO 600 !IF THERE ARE NO POINTS IN THE BLOCK

      Do L = 1,maxnodeIBS !investigate all the IB points
      	nl=0 ;dhtotal=0.d0
	IF(imb_block_loc(L).ne.dom_id(ib)) GOTO 700 !If the IB point is not in the present block
	IF(rott_loc(L).ne.2) GOTO 700	!If the Lagrangian is dynamic:exit
!NEIGHBOURS FOR THE U-GRID
          DO I = 1, dom(ib)%ttc_i 
       IF (dom(ib)%x(i) .gt.(nodex_loc(L)+nxl*dom(ib)%dx) .or.
     &     dom(ib)%x(i) .lt.(nodex_loc(L)-nxl*dom(ib)%dx)) GOTO 210
           DO J = 1, dom(ib)%ttc_j 
       IF (dom(ib)%yc(j).gt.(nodey_loc(L)+nxl*dom(ib)%dy) .or.
     &     dom(ib)%yc(j).lt.(nodey_loc(L)-nxl*dom(ib)%dy)) GOTO 211
            DO K = 1, dom(ib)%ttc_k 
       IF (dom(ib)%zc(k).gt.(nodez_loc(L)+nxl*dom(ib)%dz) .or.
     &     dom(ib)%zc(k).lt.(nodez_loc(L)-nxl*dom(ib)%dz)) GOTO 212 
!nl indicates the number of the neighbour and dh1 the delta functions value.
	nl=nl+1
 	  dh1_loc(L,nl)=dh(dom(ib)%dx,dom(ib)%dy,dom(ib)%dz,
     &  dom(ib)%X(I),dom(ib)%YC(J),dom(ib)%ZC(K)
     & ,nodex_loc(L),nodey_loc(L),nodez_loc(L),yangcase) 
!The index of the neighbours number nl to the Lagrangian L are:    
	I_nr_U(L,nl)=I ;  J_nr_U(L,nl)=J ;  K_nr_U(L,nl)=K
	  dhtotal=dhtotal+dh1_loc(L,nl)
	  if(dhtotal.ge.0.9999) goto 876
212         CONTINUE
            END DO
211         CONTINUE
           END DO
210         CONTINUE
          END DO
!        dh1_loc(L,nl)=dh1_loc(L,nl)/dhtotal
876	continue       
	kmaxU(L)=nl !# of neighbours of the Lagrangian L
!NEIGHBOURS FOR THE V-GRID
      	nl=0 ;dhtotal=0.d0
          DO I = 1, dom(ib)%ttc_i 
       IF (dom(ib)%xc(i).gt.(nodex_loc(L)+nxl*dom(ib)%dx) .or.
     &     dom(ib)%xc(i).lt.(nodex_loc(L)-nxl*dom(ib)%dx))GOTO 220
           DO J = 1, dom(ib)%ttc_j 
       IF (dom(ib)%y(j) .gt.(nodey_loc(L)+nxl*dom(ib)%dy) .or.
     &     dom(ib)%y(j) .lt.(nodey_loc(L)-nxl*dom(ib)%dy))GOTO 221 
            DO K = 1, dom(ib)%ttc_k 
       IF (dom(ib)%zc(k).gt.(nodez_loc(L)+nxl*dom(ib)%dz) .or.
     &     dom(ib)%zc(k).lt.(nodez_loc(L)-nxl*dom(ib)%dz) )  GOTO 222 
	nl=nl+1
 	  dh2_loc(L,nl)=dh(dom(ib)%dx,dom(ib)%dy,dom(ib)%dz,
     &  dom(ib)%XC(I),dom(ib)%Y(J),dom(ib)%ZC(K)
     & ,nodex_loc(L),nodey_loc(L),nodez_loc(L),yangcase)  

	I_nr_V(L,nl)=I ;  J_nr_V(L,nl)=J ;  K_nr_V(L,nl)=K
	  dhtotal=dhtotal+dh2_loc(L,nl)
	  if(dhtotal.ge.0.9999) goto 877
222         CONTINUE
            END DO
221         CONTINUE
           END DO
220         CONTINUE
          END DO
!        dh2_loc(L,nl)=dh2_loc(L,nl)/dhtotal
877	continue              
	kmaxV(L)=nl 
!NEIGHBOURS FOR THE W-GRID
      	nl=0 ;dhtotal=0.d0
          DO I = 1, dom(ib)%ttc_i 
       IF (dom(ib)%xc(i).gt.(nodex_loc(L)+nxl*dom(ib)%dx) .or.
     &     dom(ib)%xc(i).lt.(nodex_loc(L)-nxl*dom(ib)%dx) )  GOTO 230 
           DO J = 1, dom(ib)%ttc_j 
       IF (dom(ib)%yc(j).gt.(nodey_loc(L)+nxl*dom(ib)%dy) .or.
     &     dom(ib)%yc(j).lt.(nodey_loc(L)-nxl*dom(ib)%dy) )  GOTO 231 
            DO K = 1, dom(ib)%ttc_k 
       IF (dom(ib)%z(k) .gt.(nodez_loc(L)+nxl*dom(ib)%dz) .or.
     &     dom(ib)%z(k) .lt.(nodez_loc(L)-nxl*dom(ib)%dz) )  GOTO 232 
	nl=nl+1
 	  dh3_loc(L,nl)=dh(dom(ib)%dx,dom(ib)%dy,dom(ib)%dz,
     &  dom(ib)%XC(I),dom(ib)%YC(J),dom(ib)%Z(K)
     & ,nodex_loc(L),nodey_loc(L),nodez_loc(L),yangcase)  

	I_nr_W(L,nl)=I ;  J_nr_W(L,nl)=J ;  K_nr_W(L,nl)=K
	  dhtotal=dhtotal+dh3_loc(L,nl)
	  if(dhtotal.ge.0.9999) goto 878
232         CONTINUE
            END DO
231         CONTINUE
           END DO
230         CONTINUE
          END DO
!        dh3_loc(L,nl)=dh3_loc(L,nl)/dhtotal
878	continue                     
	kmaxW(L)=nl 	
	
700	CONTINUE   
       Enddo
600	CONTINUE   

     	if(myrank.eq.master)	 write(6,*)'Ended',ib

      ENDDO

!	ELSE !MLS IS USED        !!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!	Do ib=1,nbp  
!          if (ptsinblock_loc(dom_id(ib)+1).eq.0) GOTO 601
!          Do L = 1,maxnodeIBS!*
!	   IF(imb_block_loc(L).ne.dom_id(ib)) GOTO 701!*
!	   IF(rott_loc(L).ne.2) GOTO 701!*
!	
!	   call ShapeFunction_MLS(1,L,ib)
!	   call ShapeFunction_MLS(2,L,ib)
!	 do M=1,neignum_mls
!	    write(6,*)L,M,Nu(L,M)
!	 enddo

!701	CONTINUE   
 !       Enddo
!       ENDDO
!601	 CONTINUE
!	ENDIF	

      RETURN
      END
!######################################################################
      SUBROUTINE IBM
!######################################################################
      use vars
      use mpi
      use multidata
      use imb
      implicit none
      INTEGER :: NF,M,L
      double precision :: sumvel,l1norm
   
!call exchange subroutines to fill the ghost cells with values
	!call exchange(1)
	!call exchange(2)
	!call exchange(3)
	call exchange(11)
	call exchange(22)
	call exchange(33)
	

        IF (Myrank.eq.master) THEN
	  FX1NF = 0.d0  ;  FX2NF=0.d0 	;  FX3NF=0.d0 
	  FX1 = 0.d0  	;  FX2=0.d0 	;  FX3=0.d0 
	ENDIF

	DO NF =1,mdfsteps+1	!MDF loops. +1 as the default loop for IB
     	 IF (Myrank.eq.master) THEN !Calculate the accumulated force
	  DO M=1,bodynum
	   Do L=1,nodes(M)
	    FX1NF(M,L) = FX1NF(M,L) + FX1(M,L)
	    FX2NF(M,L) = FX2NF(M,L) + FX2(M,L)
	    FX3NF(M,L) = FX3NF(M,L) + FX3(M,L)
	    Enddo
	   ENDDO
	  ENDIF
	  
	  
!	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!	if(myrank.eq.master)write(6,*)'entering UV'

	  call interpolate_UV 
	  
!	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!	if(myrank.eq.master)write(6,*)'entering CALFL'	  
 	  
	  call calfl

!	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!	if(myrank.eq.master)write(6,*)'entering DB'	   
 
	  call distfbeta
	  	   
	ENDDO

      IF (Myrank.eq.master) THEN !l2-norm is calculated in reference to the final velocitiy field
	sumvel=0.d0 ; l1norm=0.d0
	DO M=1,bodynum
	 Do L=1,nodes(M)
    	   sumvel=sumvel+
     &		((FX1(M,L)*dt)**2+(FX2(M,L)*dt)**2+(FX3(M,L)*dt)**2)
	   l1norm=l1norm+((FX1(M,L)+FX2(M,L)+FX3(M,L))*dt)**2
 	 End Do
	  l2norm(M)=SQRT(sumvel/nodes(M))
	  l1norm=SQRT(l1norm)/nodes(M)
	End do

! The final force is the sum of all
	DO M=1,bodynum
	 Do L=1,nodes(M)
	  FX1NF(M,L) = FX1NF(M,L) + FX1(M,L)
	  FX2NF(M,L) = FX2NF(M,L) + FX2(M,L)
	  FX3NF(M,L) = FX3NF(M,L) + FX3(M,L)
	 Enddo
	ENDDO

	DO M=1,bodynum
	 Do L=1,nodes(M)
	  FX1(M,L) = FX1NF(M,L) ;  FX2(M,L) = FX2NF(M,L)
	  FX3(M,L) = FX3NF(M,L)
	 Enddo
	enddo

	if (bodynum.eq.1) write(757,'(3f20.5)')CTIME,l2norm(1),l1norm
	if (bodynum.ge.2) write(757,'(3f20.5)')CTIME,l2norm(1),l2norm(2)
	
	ENDIF
	
	call caldrag


      END SUBROUTINE
!######################################################################
      SUBROUTINE interpolate_UV
!######################################################################
      use vars
      use mpi
      use multidata
      use imb
      implicit none
      double precision :: dh,dhtotal
      INTEGER :: I,J,K,L,M,ib,nl,nt

	U_Beta1_loc=0.d0 ; U_Beta2_loc=0.d0 ; U_Beta3_loc=0.d0	!Pablo
	Do ib=1,nbp

      if (imbinblock_loc(dom_id(ib)+1).eq.0) GOTO 600

	if (maxnodeIBS.le.100) nt = 1
	if (maxnodeIBS.gt.100) nt = OMP_threads

	call OMP_SET_NUM_THREADS(nt)

!$OMP parallel DEFAULT(SHARED)PRIVATE(I,J,M,L,K,nl,dhtotal)
!$OMP DO SCHEDULE(DYNAMIC,1)
      Do L = 1,maxnodeIBS
	 nl=0 ; dhtotal=0.d0 
	IF(imb_block_loc(L).ne.dom_id(ib)) GOTO 700 
	 IF( rott_loc(L).eq.1 )then
          DO I = 1, dom(ib)%ttc_i 
       IF (dom(ib)%x(i) .gt.(nodex_loc(L)+nxl*dom(ib)%dx) .or.
     &     dom(ib)%x(i) .lt.(nodex_loc(L)-nxl*dom(ib)%dx)) GOTO 210
           DO J = 1, dom(ib)%ttc_j 
       IF (dom(ib)%yc(j).gt.(nodey_loc(L)+nxl*dom(ib)%dy) .or.
     &     dom(ib)%yc(j).lt.(nodey_loc(L)-nxl*dom(ib)%dy)) GOTO 211
            DO K = 1, dom(ib)%ttc_k 
       IF (dom(ib)%zc(k).gt.(nodez_loc(L)+nxl*dom(ib)%dz) .or.
     &     dom(ib)%zc(k).lt.(nodez_loc(L)-nxl*dom(ib)%dz)) GOTO 212 

	 	IF (abs(dom(ib)%USTAR(I,J,K)).gt.1.d-3) then			!Brunho comp channel 2016

	 nl=nl+1
        dh1_loc(L,nl)= dh(dom(ib)%dx,dom(ib)%dy,dom(ib)%dz,
     &  dom(ib)%X(I),dom(ib)%YC(J),dom(ib)%ZC(K)
     & ,nodex_loc(L),nodey_loc(L),nodez_loc(L),yangcase)  			!June 2015	    

           U_Beta1_loc(L) = U_Beta1_loc(L) +
     &     dom(ib)%USTAR(I,J,K) * dh1_loc(L,nl)
        
		dhtotal=dhtotal+dh1_loc(L,nl)
      	KmaxU(L)=nl
        	I_nr_U(L,nl)=I ; J_nr_U(L,nl)=J ; K_nr_U(L,nl)=K
   	            
		ENDIF

	    IF (dhtotal.ge.0.9999) GOTO 700

212         CONTINUE
            END DO
211         CONTINUE
           END DO
210         CONTINUE
          END DO
          
          if (nl.eq.0) write(6,*)L,'nl is equal to 0!!'
                    
	U_Beta1_loc(L)=U_Beta1_loc(L)*1.0d0/dhtotal

	ELSE
	 Do nl=1,KmaxU(L)
	  I=I_nr_U(L,nl) ;  J=J_nr_U(L,nl) ;  K=K_nr_U(L,nl)
	 	IF (abs(dom(ib)%USTAR(I,J,K)).gt.1.d-3) then			!Brunho comp channel 2016
          U_Beta1_loc(L)=(U_Beta1_loc(L)+
     &    dom(ib)%USTAR(I,J,K)*dh1_loc(L,nl))     
		ENDIF
  	 Enddo  	 
	ENDIF
700	continue

      Enddo
!$OMP end DO NOWAIT
!$OMP DO SCHEDULE(DYNAMIC,1)
      Do L = 1,maxnodeIBS
	IF(imb_block_loc(L).ne.dom_id(ib)) GOTO 701 
	 nl=0 ; dhtotal=0.d0 
	IF( rott_loc(L).eq.1 )then
          DO I = 1, dom(ib)%ttc_i 
       IF (dom(ib)%xc(i).gt.(nodex_loc(L)+nxl*dom(ib)%dx) .or.
     &     dom(ib)%xc(i).lt.(nodex_loc(L)-nxl*dom(ib)%dx))GOTO 220
           DO J = 1, dom(ib)%ttc_j 
       IF (dom(ib)%y(j) .gt.(nodey_loc(L)+nxl*dom(ib)%dy) .or.
     &     dom(ib)%y(j) .lt.(nodey_loc(L)-nxl*dom(ib)%dy))GOTO 221 
            DO K = 1, dom(ib)%ttc_k 
       IF (dom(ib)%zc(k).gt.(nodez_loc(L)+nxl*dom(ib)%dz) .or.
     &     dom(ib)%zc(k).lt.(nodez_loc(L)-nxl*dom(ib)%dz) )  GOTO 222 
	 	IF (abs(dom(ib)%USTAR(I,J,K)).gt.1.d-3) then			!Brunho comp channel 2016
	 nl=nl+1
        dh2_loc(L,nl)= dh(dom(ib)%dx,dom(ib)%dy,dom(ib)%dz,
     &  dom(ib)%XC(I),dom(ib)%Y(J),dom(ib)%ZC(K)
     & ,nodex_loc(L),nodey_loc(L),nodez_loc(L),yangcase) 		!June 2015	    
           U_Beta2_loc(L) = U_Beta2_loc(L) +
     &     dom(ib)%VSTAR(I,J,K) * dh2_loc(L,nl)
     
	  	dhtotal=dhtotal+dh2_loc(L,nl)  
     
      	KmaxV(L)=nl
        I_nr_V(L,nl)=I ; J_nr_V(L,nl)=J ; K_nr_V(L,nl)=K

		ENDIF
           
	    IF (dhtotal.ge.0.9999) GOTO 701
222         CONTINUE
            END DO
221         CONTINUE
           END DO
220         CONTINUE
          END DO
          
	U_Beta2_loc(L)=U_Beta2_loc(L)*1.0d0/dhtotal          
	ELSE
	 Do nl=1,KmaxV(L)
	  I=I_nr_V(L,nl) ;  J=J_nr_V(L,nl) ;  K=K_nr_V(L,nl)
	 	IF (abs(dom(ib)%USTAR(I,J,K)).gt.1.d-3) then			!Brunho comp channel 2016
          U_Beta2_loc(L)=U_Beta2_loc(L)+
     &                 dom(ib)%VSTAR(I,J,K)*dh2_loc(L,nl)
		ENDIF
  	 Enddo
	ENDIF
701	continue
      Enddo
!$OMP end DO NOWAIT
!$OMP DO SCHEDULE(DYNAMIC,1)
      Do L = 1,maxnodeIBS
	IF(imb_block_loc(L).ne.dom_id(ib)) GOTO 702 
	 nl=0 ; dhtotal=0.d0 
	IF( rott_loc(L).eq.1 )then
          DO I = 1, dom(ib)%ttc_i 
       IF (dom(ib)%xc(i).gt.(nodex_loc(L)+nxl*dom(ib)%dx) .or.
     &     dom(ib)%xc(i).lt.(nodex_loc(L)-nxl*dom(ib)%dx) )  GOTO 230 
           DO J = 1, dom(ib)%ttc_j 
       IF (dom(ib)%yc(j).gt.(nodey_loc(L)+nxl*dom(ib)%dy) .or.
     &     dom(ib)%yc(j).lt.(nodey_loc(L)-nxl*dom(ib)%dy) )  GOTO 231 
            DO K = 1, dom(ib)%ttc_k 
       IF (dom(ib)%z(k) .gt.(nodez_loc(L)+nxl*dom(ib)%dz) .or.
     &     dom(ib)%z(k) .lt.(nodez_loc(L)-nxl*dom(ib)%dz) )  GOTO 232 
	 	IF (abs(dom(ib)%USTAR(I,J,K)).gt.1.d-3) then			!Brunho comp channel 2016
	 nl=nl+1
        dh3_loc(L,nl)= dh(dom(ib)%dx,dom(ib)%dy,dom(ib)%dz,
     &  dom(ib)%XC(I),dom(ib)%YC(J),dom(ib)%Z(K)
     & ,nodex_loc(L),nodey_loc(L),nodez_loc(L),yangcase)  		!June 2015	    
           U_Beta3_loc(L) = U_Beta3_loc(L) +
     &     dom(ib)%WSTAR(I,J,K) *  dh3_loc(L,nl)

	  	dhtotal=dhtotal+ dh3_loc(L,nl)
      	KmaxW(L)=nl
        I_nr_W(L,nl)=I ; J_nr_W(L,nl)=J ; K_nr_W(L,nl)=K
		ENDIF                
	    IF (dhtotal.ge.0.9999) GOTO 702
232         CONTINUE
            END DO
231         CONTINUE
           END DO
230         CONTINUE
          END DO

	U_Beta3_loc(L)=U_Beta3_loc(L)*1.0d0/dhtotal          
          
	ELSE
	 Do nl=1,KmaxW(L)
	  I=I_nr_W(L,nl) ;  J=J_nr_W(L,nl) ;  K=K_nr_W(L,nl)
	 	IF (abs(dom(ib)%USTAR(I,J,K)).gt.1.d-3) then			!Brunho comp channel 2016
          U_Beta3_loc(L)=U_Beta3_loc(L)+
     &    dom(ib)%WSTAR(I,J,K)*dh3_loc(L,nl)
		ENDIF
  	 Enddo
	ENDIF
702	continue
      Enddo
!$OMP end DO
!$OMP END PARALLEL

600	 CONTINUE
	Enddo !ib-loop
	
	
      RETURN
      END
!######################################################################
      SUBROUTINE calfl
!######################################################################
      use vars
      use multidata
      use mpi
      use imb
      implicit none
      INTEGER :: M,L,KK,ib,iii
      double precision :: PI,aplh,UIB_loc,VIB_loc,WIB_loc
       PI = 4.D0*DATAN(1.D0)
       
	 FX1_loc = 0.d0; FX2_loc = 0.d0; FX3_loc = 0.d0


	DO ib=1,nbp

	IF (imbinblock_loc(dom_id(ib)+1).eq.0) GOTO 333	!No points within the block

        Do L = 1,maxnodeIBS
	  IF(imb_block_loc(L).ne.dom_id(ib)) GOTO 800
	   UIB_loc = 0.d0; VIB_loc = 0.d0; WIB_loc = 0.d0
		M=lag_bod_loc(L) 
	    IF(imb_shape(M).eq.5.and.rott_loc(L).eq.1) then
	     IF (turax(M).eq.1) then	! Vertical Axis Turbine
		iii=INT((L-1)/(nodes(M)/imbnumber(M)))+1			
		aplh=rads(M)+(iii-1)*2.D0*PI/imbnumber(M)	
	      UIB_loc=-radsin(M)*R0_loc(L)*cos(aplh-alpha0_loc(L))
	      VIB_loc=-radsin(M)*R0_loc(L)*sin(aplh-alpha0_loc(L))
	      WIB_loc= 0.d0	
	     ENDIF	     
	     IF (turax(M).eq.2) then ! Horizontal Axis Turbine
!		iii=INT((L-1)/(nodes(M)/imbnumber(M)))+1			
!		aplh=rads(M)+(iii-1)*2.D0*PI/imbnumber(M)	
	      UIB_loc=0.d0			
	      VIB_loc= radsin(M)*R0_loc(L)*cos(rads(M)+alpha0_loc(L))
	      WIB_loc=-radsin(M)*R0_loc(L)*sin(rads(M)+alpha0_loc(L))
	     ENDIF	
	     
	    ENDIF
!Write here any other imposed movement in case.

            FX1_loc(L) = ( UIB_loc - U_Beta1_loc(L) )/dt
            FX2_loc(L) = ( VIB_loc - U_Beta2_loc(L) )/dt
            FX3_loc(L) = ( WIB_loc - U_Beta3_loc(L) )/dt  
	

800	CONTINUE
	ENDDO        
                     
333	CONTINUE    

       Enddo !ib-loop
              		

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!The force vectors are sum with the resultant located in the master
        call MPI_ALLREDUCE (FX1_loc,FX1_MASTER,maxnodeIBS,
     &            MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr )
        call MPI_ALLREDUCE (FX2_loc,FX2_MASTER,maxnodeIBS,
     &            MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr )
        call MPI_ALLREDUCE (FX3_loc,FX3_MASTER,maxnodeIBS,
     &            MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr )

	if (myrank.eq.master) then
	KK=0 
	Do M=1,bodynum
	 Do L=1,nodes(M)
	  KK=KK+1
	  FX1(M,L)=FX1_MASTER(KK)
	  FX2(M,L)=FX2_MASTER(KK)		
	  FX3(M,L)=FX3_MASTER(KK)	 
	 enddo
	enddo

	endif   
       
      RETURN
      END
!######################################################################
      SUBROUTINE distfbeta
!######################################################################
      use vars
      use multidata
      use mpi
      use imb
      implicit none
      INTEGER :: I,J,K,L,ib,nl
      double precision :: fbeta

	Do ib=1,nbp

       if(imbinblock_loc(dom_id(ib)+1).eq.0) GOTO 600

      Do L = 1,maxnodeIBS
	IF(imb_block_loc(L).ne.dom_id(ib)) GOTO 802
	 Do nl=1,KmaxU(L)
           I=I_nr_U(L,nl) ;  J=J_nr_U(L,nl) ;  K=K_nr_U(L,nl)
	 	IF (abs(dom(ib)%USTAR(I,J,K)).gt.1.d-3) then			!Brunho comp channel 2016
          fbeta = FX1_loc(L)*dh1_loc(L,nl)*reddelta(lag_bod_loc(L))  
          dom(ib)%USTAR(I,J,K) = dom(ib)%USTAR(I,J,K) + dt*alfapr*fbeta
		endif

  	 Enddo	
	 Do nl=1,KmaxV(L)
	  I=I_nr_V(L,nl) ;  J=J_nr_V(L,nl);  K=K_nr_V(L,nl)
	 	IF (abs(dom(ib)%USTAR(I,J,K)).gt.1.d-3) then			!Brunho comp channel 2016
          fbeta = FX2_loc(L)*dh2_loc(L,nl)*reddelta(lag_bod_loc(L))   
          dom(ib)%VSTAR(I,J,K) = dom(ib)%VSTAR(I,J,K) + dt*alfapr*fbeta
		endif
  	 Enddo
	 Do nl=1,KmaxW(L)
	  I=I_nr_W(L,nl) ;  J=J_nr_W(L,nl);  K=K_nr_W(L,nl)
	 	IF (abs(dom(ib)%USTAR(I,J,K)).gt.1.d-3) then			!Brunho comp channel 2016
          fbeta = FX3_loc(L)*dh3_loc(L,nl)*reddelta(lag_bod_loc(L))      
          dom(ib)%WSTAR(I,J,K) = dom(ib)%WSTAR(I,J,K) + dt*alfapr*fbeta
		endif
  	 Enddo
802	CONTINUE
      End do 

600     CONTINUE

	Enddo !ib-loop 
	
      RETURN
      END
!######################################################################
      SUBROUTINE caldrag
!######################################################################
        use vars
        use multidata
        use imb
        use mpi
        implicit none
        INTEGER :: M,L,j,iii,inipts,finpts,totalpoints
        double precision :: PI,fx_loc,fy_loc,fz_loc,alph,alpharads
        
        PI = 4.D0*DATAN(1.D0)

	IF (Myrank.ne.master) RETURN

	J=0
	Do M = 1,bodynum
	IF (imb_shape(M).eq.5 .and. rotating(M).eq..true.) then
	  Do iii=1,imbnumber(M)
		J=J+1 ; forcefilej=399+J    
  	        fx_loc = 0.d0   ; fy_loc = 0.d0 ; fz_loc = 0.d0 
		totalpoints=nodes(M)/imbnumber(M)
		inipts=(iii-1)*totalpoints+1 ;   finpts=iii*totalpoints

	  do L=inipts,finpts
	    fx_loc = fx_loc + FX1(M,L)*dxm(M)*dym(M)*dzm(M)
     &	*reddelta(lag_bod_loc(L)) 	
	    fy_loc = fy_loc + FX2(M,L)*dxm(M)*dym(M)*dzm(M)
     &	*reddelta(lag_bod_loc(L))  	
	    fz_loc = fz_loc + FX3(M,L)*dxm(M)*dym(M)*dzm(M)
     &	*reddelta(lag_bod_loc(L))    	    	
	  end do
       	
         alph  	    = rads(M)+(iii-1)*2.*PI/imbnumber(M)
	 alpharads=alph*180./PI

         write(forcefilej,88) CTIME,alpharads,fx_loc,fy_loc,fz_loc

	  Enddo !iii-loop

	ELSE
  	  fx_loc = 0.d0   ; fy_loc = 0.d0 ; fz_loc = 0.d0 
	  J=J+1 ;	 forcefilej=399+J    
	  do L=1,nodes(M)	
	    fx_loc = fx_loc + FX1(M,L)*dxm(M)*dym(M)*dzm(M)
     &	*reddelta(lag_bod_loc(L)) 	
	    fy_loc = fy_loc + FX2(M,L)*dxm(M)*dym(M)*dzm(M)
     &	*reddelta(lag_bod_loc(L))  	
	    fz_loc = fz_loc + FX3(M,L)*dxm(M)*dym(M)*dzm(M)
     &	*reddelta(lag_bod_loc(L))   
	  end do
         write(forcefilej,88) CTIME,fx_loc,fy_loc,fz_loc
	ENDIF

       End do !M loop

   88 FORMAT (10F13.5)
      RETURN
      END SUBROUTINE caldrag
