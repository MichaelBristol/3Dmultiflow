!#############################################################################
	module module_SEM
!#############################################################################
	DOUBLE PRECISION,DIMENSION(:,:,:), ALLOCATABLE :: Vsem ,Usem
	DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: X_EDDY,EPSILO,MOLT
	DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: SIGMA
	DOUBLE PRECISION,DIMENSION(:), ALLOCATABLE ::  Ksem
	DOUBLE PRECISION,DIMENSION(:) :: X_POINT(3),REYNOLDS(6)
	DOUBLE PRECISION,DIMENSION(:) :: TEMP(3),TEMP2(3)
	DOUBLE PRECISION, DIMENSION(:,:) ::  R(3,3)
	CHARACTER*44 :: FILEGLOBAL
	INTEGER,ALLOCATABLE:: elemyst(:),elemyen(:),elemzst(:),elemzen(:)
	INTEGER,ALLOCATABLE:: iddom(:),ljdom(:),lkdom(:)
	end
