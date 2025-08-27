PROGRAM RM_UCNCT_A1
IMPLICIT NONE

INTEGER :: py, px, pz, N, fr
REAL(KIND=4), ALLOCATABLE :: psA(:,:,:),psB(:,:,:),Cn(:,:,:)
REAL(KIND=4), ALLOCATABLE, DIMENSION(:,:,:) :: RdIn

REAL(KIND=4), ALLOCATABLE :: sgF(:,:,:), dsF(:,:,:)
CHARACTER(LEN=3) :: frm
CHARACTER(LEN=4) :: ser
CHARACTER(LEN=12) :: flnm

INTEGER,PARAMETER :: yg = 136, xg = 170*4+50, zg = 236

py = yg
px = xg
pz = zg

ALLOCATE(psA(py,px,pz),Cn(py,px,pz),psB(py,px,pz),RdIn(yg,xg,zg))
ALLOCATE(sgF(0:py+1,0:px+1,0:pz+1), dsF(0:py+1,0:px+1,0:pz+1))

DO N = 0, 6 , 6

	fr = N
	WRITE(ser,"(I4)") 1000+fr
! 	flnm = 'L24Z'//ser//'.dat'
	
	!! LFP phase removal of the disconnected particles
	OPEN(UNIT=20,FILE= & 
		'/mnt/scratch/hcy/SBMES_2025/2025_0818_GraMS/DATs/'// &
		'LSS_3D_OR_136x730x236_Z_II_'//ser//'.dat', &
		FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')
		READ(20) RdIn(1:yg,1:xg,1:zg)	
	CLOSE(20)
	psA(1:py,1:px,1:pz) = RdIn(1:py,1:px,1:pz)
	psA(1:py,1:px,1:pz) = 0.5*(1.0+TANH(psA(1:py,1:px,1:pz)/1.0))
	
	WRITE(frm,"(I3)") 100+fr	
	OPEN(UNIT=20,FILE= &
		'/mnt/scratch/hcy/SBMES_2025/2025_0818_GraMS/DATs/'// &
		'SOD_3D_CN_136x730x236_Z_II_'//frm//'.dat', &		
		FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')
		READ(20) Cn(1:py,1:px,1:pz)	
	CLOSE(20)

	psB(1:py,1:px,1:pz) = psA(1:py,1:px,1:pz)
	WHERE ( Cn(1:py,1:px,1:pz) < 1.0e-4 )
		psB(1:py,1:px,1:pz) = 0.0e0
	END WHERE

	!! distancing again
	dsF(1:py,1:px,1:pz) = psB(1:py,1:px,1:pz)
	CALL NoGrdBc3D(py, px, pz, dsF(0:py+1,0:px+1,0:pz+1))

	!! sign function
	WHERE ( dsF(0:py+1,0:px+1,0:pz+1) > 0.5 ) 
		sgF(0:py+1,0:px+1,0:pz+1) = 1.0
	ELSEWHERE 
		sgF(0:py+1,0:px+1,0:pz+1) = -1.0
	END WHERE

	!! generating distance function ========================================================
	dsF(0:py+1,0:px+1,0:pz+1) = 0.0
	CALL ReDist3D_XYZ(175, py, px, pz, sgF(0:py+1,0:px+1,0:pz+1), dsF(0:py+1,0:px+1,0:pz+1))

! 	OPEN(UNIT=30,FILE= &
! 		'/mnt/scratch/hcy/SBMES_2025/2025_0818_GraMS/DATs/'// &
! 		'GRP_LS_136x730x236_Z_II_'//frm//'.dat', &		
! 		FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE')
! 		WRITE(30) dsF(1:py,1:px,1:pz)	
! 	CLOSE(30)

	OPEN(UNIT=30,FILE= &
		'/mnt/scratch/hcy/SBMES_2025/2025_0818_GraMS/DATs/'// &
		'GRP_LS_136x730x236_Z_II_'//frm//'.txt', &		
		FORM='FORMATTED',STATUS='REPLACE',ACTION='WRITE')
		WRITE(30,*) dsF(1:py,1:px,1:pz)	
	CLOSE(30)

	PRINT*,'graphite done!', N


	!! electrolyte phase
	psA(1:py,1:px,1:pz) = 1.0e0 - psA(1:py,1:px,1:pz)

	psB(1:py,1:px,1:pz) = psA(1:py,1:px,1:pz) 
	OPEN(UNIT=20,FILE= &
		'/mnt/scratch/hcy/SBMES_2025/2025_0818_GraMS/DATs/'// &
		'LIQ_3D_CN_136x730x236_Z_II_'//frm//'.dat', &		
		FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')
		READ(20) Cn(1:py,1:px,1:pz)	
	CLOSE(20)

	WHERE ( Cn(1:py,1:px,1:pz) < 1.0e-4 )
		psB(1:py,1:px,1:pz) = 0.0e0
	END WHERE

	!! distancing again
	dsF(1:py,1:px,1:pz) = psB(1:py,1:px,1:pz)
	CALL NoGrdBc3D(py, px, pz, dsF(0:py+1,0:px+1,0:pz+1))

	!! sign function
	WHERE ( dsF(0:py+1,0:px+1,0:pz+1) > 0.5 ) 
		sgF(0:py+1,0:px+1,0:pz+1) = 1.0
	ELSEWHERE 
		sgF(0:py+1,0:px+1,0:pz+1) = -1.0
	END WHERE

	!! generating distance function ========================================================
	dsF(0:py+1,0:px+1,0:pz+1) = 0.0
	CALL ReDist3D_XYZ(175, py, px, pz, sgF(0:py+1,0:px+1,0:pz+1), dsF(0:py+1,0:px+1,0:pz+1))

! 	OPEN(UNIT=30,FILE= &
! 		'/mnt/scratch/hcy/SBMES_2022/2022_1225_A/Geom_T/'// &
! 		'POR_LS_114x220x198_Z_II_'//frm//'.dat', &		
! 		FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE')
! 		WRITE(30) dsF(1:py,1:px,1:pz)	
! 	CLOSE(30)

	OPEN(UNIT=40,FILE= &
		'/mnt/scratch/hcy/SBMES_2025/2025_0818_GraMS/DATs/'// &
		'POR_LS_136x730x236_Z_II_'//frm//'.txt', &		
		FORM='FORMATTED',STATUS='REPLACE',ACTION='WRITE')
		WRITE(40,*) dsF(1:py,1:px,1:pz)	
	CLOSE(40)

	PRINT*,'Liq done!', N

ENDDO


END PROGRAM RM_UCNCT_A1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!=======================================================================================
!! for general purpose
!!=======================================================================================
SUBROUTINE NoGrdBc3D(py, px, pz, f)
IMPLICIT NONE
INTEGER, INTENT(IN) :: py, px, pz
REAL(KIND=4), INTENT(INOUT), DIMENSION(0:py+1,0:px+1,0:pz+1) :: f

f(0,1:px,1:pz) = f(1,1:px,1:pz)
f(py+1,1:px,1:pz) = f(py,1:px,1:pz)
f(:,0,1:pz) = f(:,1,1:pz)
f(:,px+1,1:pz) = f(:,px,1:pz)
f(:,:,0) = f(:,:,1)
f(:,:,pz+1) = f(:,:,pz)

END SUBROUTINE NoGrdBc3D
!!=======================================================================================
!!=======================================================================================
!! for level set
!!=======================================================================================
SUBROUTINE ReDist3D_XYZ(lit, py, px, pz, sg, f)
IMPLICIT NONE

INTEGER, INTENT(IN) :: lit, py, px, pz
REAL(KIND=4),INTENT(INOUT),DIMENSION(0:py+1,0:px+1,0:pz+1) :: sg, f

INTEGER :: i, j, k, m
REAL(KIND=4) :: dt
REAL(KIND=4),DIMENSION(0:py+1,0:px+1,0:pz+1) :: Gf, nf

dt = 1.0e-1

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
REAL(KIND=4),INTENT(IN),DIMENSION(0:py+1,0:px+1,0:pz+1) :: sg, f
REAL(KIND=4),INTENT(OUT),DIMENSION(py,px,pz) :: Gf
INTEGER :: i,j,k
REAL(KIND=4) :: ag,bg,cg,dg,eg,fg


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
				Gf(i,j,k) = 0.0e0
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
