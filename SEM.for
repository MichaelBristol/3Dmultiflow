!#############################################################################
	subroutine SEM
!#############################################################################
!THIS PROGRAM IMPLEMENTS THE SEM METHOD DESCRIBED IN N. JARRIN THESIS, CHP. 4
	USE IFPORT
	USE module_SEM
	USE vars
	USE multidata
	USE mpi
	IMPLICIT NONE
	DOUBLE PRECISION :: VOL,Ly,Lz,ENNE,XMIN,XMAX,YMIN
	DOUBLE PRECISION :: YMAX,ZMIN,ZMAX,PI,U0,HU,RIZ,RDIVz,UAVE(3)
	DOUBLE PRECISION :: SIGMA_VALUE,MAXVSEM,MINVSEM
	INTEGER :: DIVY,DIVZ,IY,IZ,II,N,I,J,M,IT,IGLOBAL
	INTEGER,allocatable,dimension(:)::lsy,lsz,ley,lez
!THE BOX DIMENSIONS ARE DIFINED AS [XLENGHT] * [Ly] * [Lz]
![DIVX], [DIVY] & [DIVZ] ARE THE NUMBERS OF SPACIAL POINTS.
	PI = 2*ACOS(0.0D0) ; Ly  = yen-yst ; Lz= zen-zst ; U0 = ubulk 				

	ALLOCATE(iddom(jdom*kdom),ljdom(jdom*kdom),lkdom(jdom*kdom))
	ALLOCATE(elemyst(jdom*kdom),elemzst(kdom*jdom))
	ALLOCATE(elemyen(jdom*kdom),elemzen(kdom*jdom))
	ALLOCATE(ley(jdom),lez(kdom),lsy(jdom),lsz(kdom))

	 ljdom=0 	; lkdom=0 	;  DIVY = 0 	; DIVZ = 0
	 elemyst=0 	; elemyen=0	;  elemzst=0 	; elemzen=0
	 I=0 ;ley=0 ; lez=0	;  lsy=0 		; lsz=0

!ID OF THE BLOCK
	 DO N=1,kdom ; DO J=1,jdom
	 I=I+1 	
	  iddom(I)=idom*(J-1)+(jdom*idom*(N-1))
	 ENDDO ; ENDDO
!DIVISION ON THE Y AND Z DIRECTIONS 
	 Do J=1,jdom
	  ley(J)=ley(J-1)+
     &   NINT((ycor(iddom(J),2)-ycor(iddom(J),1))/g_dy) +1
		 IF(J.EQ.1) THEN	
	        lsy(J)=1
		 ELSE
	        lsy(J)=ley(J-1)+1
		 ENDIF
	   DIVY = DIVY + NINT((ycor(iddom(J),2)-ycor(iddom(J),1))/g_dy) +1
	 Enddo
	 DO N=1,kdom 
	  lez(N)=lez(N-1)+
     &   NINT((zcor(iddom(N),2)-zcor(iddom(N),1))/g_dz)+1
		 IF(N.EQ.1) THEN	
	        lsz(N)=1 
		 ELSE
	        lsz(N)=lez(N-1)+1 
		 ENDIF
	  DIVZ = DIVZ + NINT((zcor(iddom(N),2)-zcor(iddom(N),1))/g_dz) +1
	 ENDDO
!THE DIVISIONS FOR EACH OF THE DOMAINS IS DETERMINED:
	 I=0
	 DO N=1,kdom  ;  DO J=1,jdom
	 I=I+1 	
	  elemyst(I) = lsy(J)  ;  elemzst(I) = lsz(N)
	  elemyen(I) = ley(J)  ;  elemzen(I) = lez(N)
	  ljdom(I) = elemyen(I) -  elemyst(I) + 1
	  lkdom(I) = elemzen(I) -  elemzst(I) + 1
	 enddo ; enddo

	write(6,*)'Divisions :',DIVY,DIVZ	
	write(6,*)'# BLOCKS  :',jdom,kdom	
	write(6,*)lsy(1),ley(1),lsy(2),ley(2),ley(2)-lsy(2)+1
	write(6,*)lsz(1),lez(1),lsz(2),lez(2),lez(2)-lsz(2)+1

!THE AMPLITUDE OF THE VORTICES. This can be changed to have larger structures!!!!!!
	SIGMA_VALUE=MIN(16.D0*g_dy,Lz/2.D0,Ly/2.D0) !isotropic

!THE NUMBER OF EDDIES IS THE INLET SURFACE DIVIDED BY THE SURFACE OF EACH TURBULENT SPOT
	NE_SEM = (Ly*Lz)/SIGMA_VALUE**2	
![N](INTEGER) AND [ENNE](DOUBLE PRECISION) REPRESENT THE NUMBER OF EDDIES
	N = INT(NE_SEM)	;  ENNE = REAL(N)

	if(myrank.eq.0) write(6,*) 'The number of SEM eddies is :', N

!  [REYNOLDS(6)] IS A VECTOR WITH THE SIX ELEMENTS OF REYNOLDS STRESSES.
!  |REYNOLDS(1)  REYNOLDS(2)  REYNOLDS(4)|
!  |REYNOLDS(2)  REYNOLDS(3)  REYNOLDS(5)|
!  |REYNOLDS(4)  REYNOLDS(5)  REYNOLDS(6)|
      REYNOLDS=(/(TI_SEM*U0)**2, 0.0D0,(TI_SEM*U0)**2, 0.0D0, 
     & 0.0D0,(TI_SEM*U0)**2/)   ![M/S]

!ALLOCATION OF THE EDDIES VECTOR.
![VSEM(DIVX,DIVY,DIVZ,3)] IS THE INSTANTANEOUS VELOCITY VECTOR IN THE POINT WITH
!X,Y,Z COMPONENTS
![X_EDDY(3,N)] IS THE N-TH EDDY LOCATION X,Y,Z
![EPSILO(3,N)] IS THE N-TH EDDY INTENSITY IN X,Y,Z
![MOLT(3,N)]   IS THE MATRIX PRODUCT [R(3,3)] * [EPSILO(3,N)]
![Ksem(N)]     IS A VECTOR USED TO CHECK WHITCH EDDY IS OUTSIDE THE BOX AFTER THE CONVECTION

	ALLOCATE(Vsem(DIVY,DIVZ,3),Usem(DIVY,DIVZ,3),Ksem(N))
	ALLOCATE(X_EDDY(3,N),EPSILO(3,N),MOLT(3,N),SIGMA(DIVY,DIVZ))
	PRINT *, "TOTAL TIME INTERVAL = ", DT * ITMAX_SEM, " [S]"

	XMIN=1.D8 ; XMAX=1.D-8 ; YMIN=1.D8 ; YMAX=1.D-8
	ZMIN=1.D8 ; ZMAX=1.D-8

!DEFINITION OF EDDY LENGTH SCALE AND INITIAL VELOCITY FIELD
	DO IY=1,DIVY
	  DO IZ=1,DIVZ
	   SIGMA(IY,IZ)=SIGMA_VALUE    !!!  MIN(8.0*g_dy,0.20D0) !isotropic
!Set the inlet velocity prof.  
  	Usem(IY,IZ,:)=(/ U0 ,0.D0,0.D0/)			
!	if(UPROF_SEM.eq.1) Usem(IY,IZ,:)=(/ U0 ,0.D0,0.D0/)      !elli
!   	if(UPROF_SEM.eq.15) then
!        RIZ=IZ							!elli
!	   Usem(IY,IZ,:)=(/ U0-0.04324*exp(-11.1764*RIZ), 0.D0,0.D0/)
!	endif
	   UAVE(:) = UAVE(:) + Usem(IY,IZ,:)
!CALCULATION OF THE EDDY BOX PARAMETERS
	XMIN=MIN(0.0D0-SIGMA(IY,IZ),XMIN);XMAX=MAX(0.0D0+SIGMA(IY,IZ),XMAX)
	YMIN=MIN(0.0D0-SIGMA(IY,IZ),YMIN);YMAX=MAX(Ly   +SIGMA(IY,IZ),YMAX)
	ZMIN=MIN(0.0D0-SIGMA(IY,IZ),ZMIN);ZMAX=MAX(Lz   +SIGMA(IY,IZ),ZMAX)
	  END DO
	END DO

	UAVE = UAVE / (DIVY * DIVZ)
	VOL  = (XMAX - XMIN) * (YMAX - YMIN) * (ZMAX - ZMIN)
!GENERATION OF THE EDDY LOCATION INSIDE THE BOX AND INITIALIZATION OF THE [Ksem] VECTOR
	DO II=1,N
	  X_EDDY(1,II) = (XMAX - XMIN) * RAND() + XMIN
	  X_EDDY(2,II) = (YMAX - YMIN) * RAND() + YMIN
	  X_EDDY(3,II) = (ZMAX - ZMIN) * RAND() + ZMIN
	  Ksem(II) = 0
!INITIALIZATION OF THE INTENSITIES. FOR EVERY DIRECTION THE AVERAGE INTENSITY VALUE IS CALCULATED AND
!IT IS FORCED TO BE LOWER THAN THE [VLIM] VALUE
	  EPSILO(1,II) = (RAND()*2.0D0 - 1.0D0)
	  EPSILO(2,II) = (RAND()*2.0D0 - 1.0D0)
	  EPSILO(3,II) = (RAND()*2.0D0 - 1.0D0)
	END DO
!INITIALIZATION OF THE [R(3,3)] MATRIX WITH THE CHOLENSKY DECOMPOSITION
!OF THE REYNOLDS STRESS TENSOR
	R = 0
	R(1,1) = DSQRT(REYNOLDS(1))
	R(2,1) = REYNOLDS(2) / R(1,1)
	R(2,2) = DSQRT(REYNOLDS(3) - R(2,1)*R(2,1))
	R(3,1) = REYNOLDS(4) / R(1,1)
	R(3,2) = (REYNOLDS(5) - R(2,1)*R(3,1)) / R(2,2)
	R(3,3) = DSQRT(REYNOLDS(6) - R(3,1)*R(3,1) - R(3,2)*R(3,2))
!BEGINNING OF TIME ITERATIONS
	DO IT=1,ITMAX_SEM			!PARALLELIZE THIS LOOP
	 if(mod(IT,50).EQ.0) 
     &  WRITE(*,*)"ITERATION ",IT,"IN PROGRESS.TIME: ",(IT-1) * DT,'[S]'

!PRINTINGS OF GLOBAL VELOCITY AND OF CONVECTION VELOCITY
	 Do I=1,jdom*kdom
		  IGLOBAL=500+I
	   WRITE (FILEGLOBAL,'(A13,i4.4,a1,I6.6,A4)')
     &		'inflow/Inlet_',iddom(I),'_',IT,'.dat'
	   OPEN(UNIT=IGLOBAL, FILE=FILEGLOBAL,STATUS="UNKNOWN", ACTION="WRITE")
	   WRITE (IGLOBAL,*)'Variables=up,vp,wp'
	   WRITE (IGLOBAL,*)
     &  'zone ',' i=',ljdom(I),',',' j=',lkdom(I),', k= ',1,' f=point'
	 Enddo

	  MOLT = MATMUL(R,EPSILO) ! !aij*epsij   MATRIX MULTIPLICATION
	  MAXVSEM=0.D0 ; MINVSEM=1000000.D0
!------BEGINNING OF SPATIAL ITERATION
	  DO IY = 1,DIVY
	    DO IZ = 1,DIVZ
!X_POINT = GRID POINT COORDINATES
	   X_POINT=(/0.0D0,IY*Ly/DIVY+YMIN,IZ*Lz/DIVZ+ZMIN/)
	   Vsem(IY,IZ,:)=(/ 0.d0, 0.d0,  0.d0 /)
!------------BEGINNING OF EDDIES ITERATIONS
		DO II=1,N
		  TEMP(:) = DABS(X_POINT(:) - X_EDDY(:,II))
	  IF (TEMP(1).LT.SIGMA(IY,IZ) .AND. TEMP(2).LT.SIGMA(IY,IZ) .AND.
     &      TEMP(3).LT.SIGMA(IY,IZ)) THEN
	  Vsem(IY,IZ,:)=Vsem(IY,IZ,:)+MOLT(:,II)
     &    *(DSQRT(1.5D0)**3.0D0)*DSQRT(VOL)/DSQRT(SIGMA(IY,IZ)**3)
     &    *(1.0D0- DABS(X_POINT(:) - X_EDDY(:,II))/SIGMA(IY,IZ))
     &    *(1.0D0- DABS(X_POINT(:) - X_EDDY(:,II))/SIGMA(IY,IZ))
     &    *(1.0D0- DABS(X_POINT(:) - X_EDDY(:,II))/SIGMA(IY,IZ))
!f(x)=SQRT(3/2)*(1-ABS(X)) JARRIN ET AL. 2009
		 END IF
		END DO
! Instananeous velocity=mean velocity + SEM fluctuation velocity
	  Vsem(IY,IZ,:) = Vsem(IY,IZ,:) / DSQRT(ENNE) 
	  if(Vsem(IY,IZ,1).gt.MAXVSEM) MAXVSEM=Vsem(IY,IZ,1)
	  if(Vsem(IY,IZ,1).le.MINVSEM) MINVSEM=Vsem(IY,IZ,1)
	    END DO
	  END DO

!	 write(6,*)	'========='
	  write(6,'(i6,2f15.6)')IT,MAXVSEM,MINVSEM
!------END OF SPATIAL ITERATIONS
	      Do I=1,jdom*kdom
	        IGLOBAL=500+I
		   Do M=elemzst(I),elemzen(I)
		    Do J=elemyst(I),elemyen(I)
                  WRITE(500+I,'(3E15.6)')Vsem(J,M,:)				!turn down precision for large files
	  	  enddo ;enddo 
	    	 CLOSE (UNIT=IGLOBAL)
	  	Enddo
!--------BEGINNING OF EDDIES CONVECTION ITERATIONS
	  DO II=1,N
!RE-CALCULATION OF THE EDDIES POSITION. IF ANY EDDY GOES BEYOND THE BOX LIMITS
!IT IS RESTARTED AT THE SURFACE FACING THE EXIT.
	    X_EDDY(1,II) = X_EDDY(1,II) + UAVE(1) * dt
	    X_EDDY(2,II) = X_EDDY(2,II) + UAVE(2) * dt
	    X_EDDY(3,II) = X_EDDY(3,II) + UAVE(3) * dt
!AFTER THE EDDIES CONVECTIONS ARE NECESSARY SOME TESTS TO CHECK IF ANY EDDY IS NOW OUTSIDE
!THE SEM BOX DEFINED EARLIER
!WHEN A EDDY IS RE-GENERATE THE Ksem(I) FACTOR ASSUME THE 1 VALUE. THIS VALUE IS USED LATER TO
!GENERATE A NEW INTENSITY FOR THE NEW EDDY
	    IF (X_EDDY(1,II) > XMAX) THEN
		X_EDDY(1,II) = XMIN
		X_EDDY(2,II) = (YMAX - YMIN) * RAND() + YMIN
		X_EDDY(3,II) = (ZMAX - ZMIN) * RAND() + ZMIN
		Ksem(II) = 1
	    ELSE IF (X_EDDY(1,II) < XMIN) THEN
		X_EDDY(1,II) = XMAX
		X_EDDY(2,II) = (YMAX - YMIN) * RAND() + ZMIN
		X_EDDY(3,II) = (ZMAX - ZMIN) * RAND() + ZMIN
		Ksem(II) = 1
	    ELSE IF (X_EDDY(2,II) > YMAX) THEN
		X_EDDY(1,II) = (XMAX - XMIN) * RAND() + XMIN
		X_EDDY(2,II) = YMIN
		X_EDDY(3,II) = (ZMAX - ZMIN) * RAND() + ZMIN
		Ksem(II) = 1
	    ELSE IF (X_EDDY(2,II) < YMIN) THEN
		X_EDDY(1,II) = (XMAX - XMIN) * RAND() + XMIN
		X_EDDY(2,II) = YMAX
		X_EDDY(3,II) = (ZMAX - ZMIN) * RAND() + ZMIN
		Ksem(II) = 1
	    ELSE IF  (X_EDDY(3,II) > ZMAX) THEN
		X_EDDY(1,II) = (XMAX - XMIN) * RAND() + XMIN
		X_EDDY(2,II) = (YMAX - YMIN) * RAND() + ZMIN
		X_EDDY(3,II) = ZMIN 
		Ksem(II) = 1
	    ELSE IF (X_EDDY(3,II) < ZMIN) THEN
		X_EDDY(1,II) = (XMAX - XMIN) * RAND() + XMIN
		X_EDDY(2,II) = (YMAX - YMIN) * RAND() + ZMIN
		X_EDDY(3,II) = ZMAX
		Ksem(II) = 1
	    END IF
!INTENSITY GENERATION FOR THE RE-CREATED EDDIES. WE ARE USING THE Ksem FACTOR AS EXPLAINED FATOR.
	    IF (Ksem(II)== 1) THEN
		EPSILO(3,II) = (RAND()*2.0D0 - 1.0D0)
		EPSILO(2,II) = (RAND()*2.0D0 - 1.0D0)
		EPSILO(1,II) = (RAND()*2.0D0 - 1.0D0)
	    END IF
	    Ksem(II) = 0
	  END DO
!---------END OF EDDIES CONVECTIONS ITERATIONS
	END DO
!----END OF TIME ITERATIONS

	END SUBROUTINE SEM
