!! COMBINE RANKS TO SINGLE OUTPUT

PROGRAM SBMES_3D_Cmbn_MPI_v1
IMPLICIT NONE


INTEGER :: nx, ny, nz, St, tfrm, fl, fb, fe, itv
INTEGER :: gy, gx, gz, rnb, cnb, ddR, ddC, lwy, upy, lwz, upz, np, rank
INTEGER :: rmvz, bivz, rmvy, bivy
REAL(KIND=4), ALLOCATABLE :: ReadIn(:,:,:)
REAL(KIND=4), ALLOCATABLE :: Und(:,:,:)
INTEGER, ALLOCATABLE :: IdxAry(:,:)
CHARACTER(LEN=2) :: pred
CHARACTER(LEN=4) :: frk, fsn
CHARACTER(LEN=18) :: flnm

gy = 186;	gx = 1460;	gz = 320


rnb = 6*2		!! rows of ranks
cnb = 8*2		!! columns of ranks

np = rnb*cnb

nx = gx

!! create index array for read-in data
ALLOCATE(IdxAry(0:np-1,6))
!! find the relationship between local index and global index
DO rank = 0, np-1

	ddC = INT(rank/rnb)+1
	ddR = MOD(rank,rnb)+1

	!! determine subdomain size along Z-axis
	bivz = INT(gz/rnb)
	rmvz = MOD(gz,rnb)
	IF ( (ddR-1) <= rmvz-1 ) THEN
		nz = bivz+1
	ELSE
		nz = bivz
	ENDIF
	!! global index of each subdomain
	lwz = (ddR-1)*nz+1
	IF ( (ddR-1) >= rmvz ) lwz = lwz+rmvz
	upz = lwz+(nz-1)

	!! determine subdomain size along X-axis
	bivy = INT(gy/cnb)
	rmvy = MOD(gy,cnb)
	IF ( (ddC-1) <= rmvy-1 ) THEN
		ny = bivy+1
	ELSE
		ny = bivy
	ENDIF
	!! global index of each subdomain
	lwy = (ddC-1)*ny+1
	IF ( (ddC-1) >= rmvy ) lwy = lwy+rmvy
	upy = lwy+(ny-1)


	!! global index of each subdomain
	IdxAry(rank,1) = lwy
	IdxAry(rank,2) = upy
	IdxAry(rank,3) = ny
	IdxAry(rank,4) = lwz
	IdxAry(rank,5) = upz
	IdxAry(rank,6) = nz

! 	print*, rank, IdxAry(rank,1:6)

ENDDO

fb = 111
fe = 131
fl = 135
itv = 10

!! read-in data
ALLOCATE(Und(1:gy,1:gx,1:gz))
DO St = 2, 5
! DO St = 8, 9

	IF ( St == 1 ) THEN
		pred = 'CG'
	ELSEIF ( St == 2 ) THEN
		pred = 'CP'
	ELSEIF ( St == 3 ) THEN
		pred = 'CE'
	ELSEIF ( St == 4 ) THEN
		pred = 'PE'
	ELSEIF ( St == 5 ) THEN
		pred = 'PP'
	ELSEIF ( St == 6 ) THEN
		pred = 'PD'
	ELSEIF ( St == 7 ) THEN
		pred = 'PG'	
	ELSEIF ( St == 8 ) THEN
		pred = 'PS'			
	ELSEIF ( St == 9 ) THEN		
		pred = 'PL'
	ENDIF

! 	DO tfrm = fb, fe, itv
	DO tfrm = fl, fl !, 10	 
! 	DO tfrm = 1, 1
	
		DO rank = 0, np-1

			WRITE(frk,"(I4)") 1000+rank
			WRITE(fsn,"(I4)") 1000+tfrm

			flnm = pred//frk//fsn//'.dat'

			ALLOCATE(ReadIn(1:IdxAry(rank,3),0:nx+1,1:IdxAry(rank,6)))

			OPEN(UNIT=20,FILE='/mnt/gs21/scratch/hcy/SBMES_2022/2022_0720_A/DataCH1/'// &
				flnm,FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')
				READ(20) ReadIn(1:IdxAry(rank,3),1:nx,1:IdxAry(rank,6))
			CLOSE(20)

			Und(IdxAry(rank,1):IdxAry(rank,2),1:nx,IdxAry(rank,4):IdxAry(rank,5)) = &
				ReadIn(1:IdxAry(rank,3),1:nx,1:IdxAry(rank,6))

			DEALLOCATE(ReadIn)
		ENDDO
		PRINT*,pred,tfrm,rank

		CALL OutPut3D_A(tfrm, gy, gx, gz, pred//'ASMB', 2, Und(1:gy,1:gx,1:gz))

	ENDDO

ENDDO



END PROGRAM SBMES_3D_Cmbn_MPI_v1
!! ==================== ==================== ==================== ==================== !!

!! SUBROUTINES
!! ============================================================================================== !!
!! Data IOs
!! ============================================================================================== !!

SUBROUTINE OutPut3D_A(fr,py,px,pz,prfn,AB,Cn1)
IMPLICIT NONE

INTEGER,INTENT(IN) :: fr, px, py, pz, AB
REAL(KIND=4),INTENT(IN),DIMENSION(1:py,1:px,1:pz) :: Cn1
CHARACTER(LEN=6),INTENT(IN) :: prfn
CHARACTER(LEN=4) :: fs
CHARACTER(LEN=14) :: flnm

WRITE(fs,"(I4)") 1000+fr

!! case 1 ==> ascii; case 2 ==> binary
SELECT CASE (AB)
	CASE (1)
		flnm = prfn//fs//'.txt'
		OPEN(UNIT=20,FILE=flnm,FORM='FORMATTED',STATUS='REPLACE',ACTION='WRITE')
			WRITE(20,*) SNGL(Cn1)
		CLOSE(20)
	CASE (2)
		flnm = prfn//fs//'.dat'
		OPEN(UNIT=20,FILE='/mnt/scratch/hcy/SBMES_2022/2022_0720_A/DataCHA/'// &	
			flnm,FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE')
! 			WRITE(20) SNGL(Cn1)
			WRITE(20) Cn1
		CLOSE(20)
END SELECT

END SUBROUTINE OutPut3D_A
