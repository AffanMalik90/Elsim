!! This program does distance fucntion and tanh function from 3D voxelated data for LSC.


PROGRAM LevelSet_3D_v03
IMPLICIT NONE
!! variables for general code
INTEGER,PARAMETER :: yg = 136, xg = 170*4+50, zg = 236
INTEGER,PARAMETER :: yo = 500, xo = 520, zo = 170
INTEGER :: px, py, pz, i, j, k, AB, iter, n
REAL(KIND=8) :: dh, dx, dy, dz, rad, tr
!! variables for level-set initialization
INTEGER :: fr

REAL(KIND=8), DIMENSION(0:yg+1,0:xg+1,0:zg+1) :: sgF, dsF
REAL(KIND=8), DIMENSION(yo,xo,zo) :: TmpA
REAL(KIND=8), DIMENSION(yo,zo,xo) :: TmpB
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: CdsF
REAL(KIND=8), DIMENSION(0:13) :: TuRd



CHARACTER(LEN=4) :: prex, frm

!! setups
px = xg
py = yg
pz = zg

dh = 1.0d0
AB = 2

DATA TuRd / 0.0, 9.384, 18.768, 28.151, 37.535, 46.919, 56.303,  &
	65.686, 75.070, 84.454, 93.838, 103.221, 112.605, 121.989 /

!! loading  sign function from smoothed voxelated LSC geometry ====================================
OPEN(UNIT=20,FILE='/mnt/scratch/hcy/SBMES_2025/2025_0818_GraMS/Geom_Vx/'// &
	'GraMicVox_II_1_500_520_170.txt',FORM='FORMATTED',STATUS='OLD',ACTION='READ')
	READ(20,*) TmpA(1:yo,1:xo,1:zo)
CLOSE(20)

prex = 'LHQP'

DO n = 0, 6, 6 
	fr = n
	tr = TuRd(n)
	
	!! load in the graphite microstructure
	dsF(1:yg,1:xg,1:zg) = 0.0

	WRITE(frm,"(I4)") 1000+fr
	print*, frm

	DO k = 1, zo
		TmpB(1:yo,k,1:xo) = TmpA(1:yo,1:xo,zo-k+1)
	ENDDO

	!! tunnel
	DO k = 1, pz; DO j = 1, px; DO i = 1, py  
		rad = SQRT((DBLE(i)-1.0)**2+(DBLE(k)-236.0)**2)	
		IF ( rad < tr ) THEN	!! the tunnel radius
			TmpB(i,j,k) = 0.0
		ENDIF
	
		rad = SQRT((DBLE(i)-136.0)**2+(DBLE(k)-1.0)**2)	
		IF ( rad < tr ) THEN	!! the tunnel radius
			TmpB(i,j,k) = 0.0
		ENDIF	
	ENDDO; ENDDO; ENDDO
    !! first thickness
	dsF(1:yg,50+1:220,1:zg) = TmpB(1:yg,1:zo,1:zg)
    !! second thickness 
	DO j = 221, 390
		dsF(1:yg,j,1:zg) = TmpB(1:yg, 390+1-j,1:zg)
	ENDDO
    ! third thickness
 	DO j = 391, 560
 		dsF(1:yg,j,1:zg) = TmpB(1:yg, j-390,1:zg)
 	ENDDO
    ! fourth thickness
 	DO j = 561, 730
 		dsF(1:yg,j,1:zg) = TmpB(1:yg, 390+1-j,1:zg)
 	ENDDO

	CALL NoGrdBc3D(py, px, pz, dsF(0:py+1,0:px+1,0:pz+1))

	!! pre-smooth with Allen-Cahn phase field method
	DO iter = 1, 50 !! 200
	
		sgF(1:yg,1:xg,1:zg) = 2*dsF(1:yg,1:xg,1:zg)*(1.0-dsF(1:yg,1:xg,1:zg))*(1.0-2*dsF(1:yg,1:xg,1:zg)) - &
			0.64*(dsF(2:yg+1,1:xg,1:zg) + dsF(0:yg-1,1:xg,1:zg) + &
				  dsF(1:yg,2:xg+1,1:zg) + dsF(1:yg,0:xg-1,1:zg) + &
				  dsF(1:yg,1:xg,2:zg+1) + dsF(1:yg,1:xg,0:zg-1) - 6*dsF(1:yg,1:xg,1:zg))
	
		dsF(1:yg,1:xg,1:zg) = dsF(1:yg,1:xg,1:zg) - 0.05*sgF(1:yg,1:xg,1:zg)
	
		CALL NoGrdBc3D(py, px, pz, dsF(0:py+1,0:px+1,0:pz+1))
		print*, 'AC', iter, n
	ENDDO			

	!! sign function
	WHERE ( dsF(0:py+1,0:px+1,0:pz+1) > 0.5d0 ) 
		sgF(0:py+1,0:px+1,0:pz+1) = 1.0d0
	ELSEWHERE 
		sgF(0:py+1,0:px+1,0:pz+1) = -1.0d0
	END WHERE


	!! generating distance function ===================================================================
	dsF(0:py+1,0:px+1,0:pz+1) = 0.0
	CALL ReDist3D_XYZ(175, py, px, pz, sgF(0:py+1,0:px+1,0:pz+1), dsF(0:py+1,0:px+1,0:pz+1))

	CALL OutPut3D1(fr, prex, py, px, pz, AB, dsF(1:py,1:px,1:pz))

ENDDO

END PROGRAM LevelSet_3D_v03
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!=======================================================================================
!! for general purpose
!!=======================================================================================
SUBROUTINE NoGrdBc3D(py, px, pz, f)
IMPLICIT NONE
INTEGER, INTENT(IN) :: py, px, pz
REAL(KIND=8), INTENT(INOUT), DIMENSION(0:py+1,0:px+1,0:pz+1) :: f

f(0,1:px,1:pz) = f(1,1:px,1:pz)
f(py+1,1:px,1:pz) = f(py,1:px,1:pz)
f(:,0,1:pz) = f(:,1,1:pz)
f(:,px+1,1:pz) = f(:,px,1:pz)
f(:,:,0) = f(:,:,1)
f(:,:,pz+1) = f(:,:,pz)

END SUBROUTINE NoGrdBc3D
!!=======================================================================================
!!=======================================================================================
SUBROUTINE OutPut3D1(fr, prfn, py, px, pz, AB, Cn1)
IMPLICIT NONE
INTEGER, INTENT(IN) :: py, px, pz
INTEGER,INTENT(INOUT) :: fr
REAL(KIND=8),INTENT(IN),DIMENSION(1:py,1:px,1:pz) :: Cn1
INTEGER,INTENT(IN) :: AB
CHARACTER(LEN=4),INTENT(IN) :: prfn
CHARACTER(LEN=4) :: fs
CHARACTER(LEN=12) :: flnm

WRITE(fs,"(I4)") 1000+fr

!! case 1 ==> ascii; case 2 ==> binary
SELECT CASE (AB)
	CASE (1)
		flnm = prfn//fs//'.txt'
		OPEN(UNIT=20,FILE=flnm,FORM='FORMATTED',STATUS='REPLACE',ACTION='WRITE')
			!WRITE(20,*) SNGL(Cn1)
			WRITE(20,*) Cn1
		CLOSE(20)
	CASE (2)
		flnm = prfn//fs//'.dat'
		OPEN(UNIT=30,FILE='/mnt/scratch/hcy/SBMES_2025/2025_0818_GraMS/DATs/'// &
			flnm,FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE')
			WRITE(30) Cn1
		CLOSE(30)
	CASE (3)
		flnm = prfn//fs//'.dat'
		OPEN(UNIT=30,FILE=flnm,FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE')
			WRITE(20) SNGL(Cn1)
		CLOSE(30)		
END SELECT

END SUBROUTINE OutPut3D1
!!=======================================================================================
!! for level set
!!=======================================================================================
SUBROUTINE ReDist3D_XYZ(lit, py, px, pz, sg, f)
IMPLICIT NONE

INTEGER, INTENT(IN) :: lit, py, px, pz
REAL(KIND=8),INTENT(INOUT),DIMENSION(0:py+1,0:px+1,0:pz+1) :: sg, f

INTEGER :: i, j, k, m
REAL(KIND=8) :: dt
REAL(KIND=8),DIMENSION(0:py+1,0:px+1,0:pz+1) :: Gf, nf

dt = 1.0d-1

DO m = 1, lit

	CALL AbsGrd3DUW_XYZ(py, px, pz, sg(0:py+1,0:px+1,0:pz+1), f(0:py+1,0:px+1,0:pz+1), &
		Gf(1:py,1:px,1:pz))
	CALL NoGrdBc3D(py, px, pz, Gf(0:py+1,0:px+1,0:pz+1))
	
	DO k = 0, pz+1
		DO j = 0, px+1
			DO i = 0, py+1
				nf(i,j,k) = f(i,j,k) + sg(i,j,k)*(1.0d0-Gf(i,j,k))*dt
				IF ( nf(i,j,k)*f(i,j,k) > 0.0 ) THEN
				ELSE
					nf(i,j,k) = f(i,j,k)+dt*1.0d-7*sg(i,j,k)
				ENDIF
			ENDDO
		ENDDO
	ENDDO	
	
	print*, m
	
	f(0:py+1,0:px+1,0:pz+1) = nf(0:py+1,0:px+1,0:pz+1)	
ENDDO

END SUBROUTINE ReDist3D_XYZ
!!=======================================================================================
!!=======================================================================================
SUBROUTINE AbsGrd3DUW_XYZ(py, px, pz, sg, f, Gf)
IMPLICIT NONE

INTEGER, INTENT(IN) :: py, px, pz
REAL(KIND=8),INTENT(IN),DIMENSION(0:py+1,0:px+1,0:pz+1) :: sg, f
REAL(KIND=8),INTENT(OUT),DIMENSION(py,px,pz) :: Gf
INTEGER :: i,j,k
REAL(KIND=8) :: ag,bg,cg,dg,eg,fg


DO k = 1, pz
	DO j = 1, px
		DO i = 1, py 
			
			ag = (f(i,j,k)-f(i-1,j,k))
			bg = (f(i+1,j,k)-f(i,j,k))
			cg = (f(i,j,k)-f(i,j-1,k))
			dg = (f(i,j+1,k)-f(i,j,k))
			eg = (f(i,j,k)-f(i,j,k-1))
			fg = (f(i,j,k+1)-f(i,j,k))
			
			IF ( sg(i,j,k) == 0.0 ) THEN
				Gf(i,j,k) = 0.0d0
			ELSEIF ( sg(i,j,k) > 0.0 ) THEN
				Gf(i,j,k) = SQRT(MAX(MAX(ag,0.0d0)**2,MAX(-bg,0.0d0)**2)+ &
					             MAX(MAX(cg,0.0d0)**2,MAX(-dg,0.0d0)**2)+ &
					             MAX(MAX(eg,0.0d0)**2,MAX(-fg,0.0d0)**2))
			ELSEIF ( sg(i,j,k) < 0.0 ) THEN
				Gf(i,j,k) = SQRT(MAX(MAX(-ag,0.0d0)**2,MAX(bg,0.0d0)**2)+ &
							     MAX(MAX(-cg,0.0d0)**2,MAX(dg,0.0d0)**2)+ &
							     MAX(MAX(-eg,0.0d0)**2,MAX(fg,0.0d0)**2))
			ENDIF
		ENDDO
	ENDDO
ENDDO	
			
END SUBROUTINE AbsGrd3DUW_XYZ
!!=======================================================================================




