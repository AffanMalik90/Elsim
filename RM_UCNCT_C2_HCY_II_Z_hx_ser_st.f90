!! This program removes disconnected regions

PROGRAM RM_UCNCT_3D
   IMPLICIT NONE

   ! grid dimensions and frame variables
   INTEGER :: py, px, pz, N, fr
   ! domain parameters and
   REAL(KIND=8), ALLOCATABLE :: psA(:,:,:), psB(:,:,:), Cn(:,:,:)

   REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: RdIn

   REAL(KIND=8), ALLOCATABLE :: sgF(:,:,:), dsF(:,:,:)
   CHARACTER(LEN=3) :: frm
   CHARACTER(LEN=4) :: ser
   CHARACTER(LEN=12) :: flnm

   INTEGER,PARAMETER :: yg = 93, xg = 170*4+50, zg = 160

   py = yg
   px = xg
   pz = zg

   ALLOCATE(psA(py,px,pz), Cn(py,px,pz), psB(py,px,pz), RdIn(yg,xg,zg))
   ALLOCATE(sgF(0:py+1,0:px+1,0:pz+1), dsF(0:py+1,0:px+1,0:pz+1))

   DO N = 3, 13 , 1

      fr = N
      WRITE(ser,"(I4)") 1000+fr
      flnm = 'L24Z'//ser//'.dat'

      !! LFP phase removal of the disconnected particles
      OPEN(UNIT=20,FILE='/mnt/scratch/hcy/SBMES_2022/2022_1225_A/Geom_T/'//flnm, &
         FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')
      READ(20) RdIn(1:yg,1:xg,1:zg)
      CLOSE(20)
      psA(1:py,1:px,1:pz) = SNGL(RdIn(1:py,1:px,1:pz))
      psA(1:py,1:px,1:pz) = 0.5*(1.0+TANH(psA(1:py,1:px,1:pz)/1.0))

      WRITE(frm,"(I3)") 100+fr
      OPEN(UNIT=20,FILE= &
         '/mnt/scratch/hcy/SBMES_2022/2022_1225_A/Geom_T/'// &
         'SOD_3D_CN_93x730x160_Z_II_'//frm//'.dat', &
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

      OPEN(UNIT=30,FILE= &
         '/mnt/scratch/hcy/SBMES_2022/2022_1225_A/Geom_T/'// &
         'GRP_LS_93x730x160_Z_II_'//frm//'.dat', &
         FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE')
      WRITE(30) dsF(1:py,1:px,1:pz)
      CLOSE(30)

      PRINT*,'graphite done!', N


      !! electrolyte phase
      psA(1:py,1:px,1:pz) = 1.0e0 - psA(1:py,1:px,1:pz)

      psB(1:py,1:px,1:pz) = psA(1:py,1:px,1:pz)
      OPEN(UNIT=20,FILE= &
         '/mnt/scratch/hcy/SBMES_2022/2022_1225_A/Geom_T/'// &
         'LIQ_3D_CN_93x730x160_Z_II_'//frm//'.dat', &
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
         '/mnt/scratch/hcy/SBMES_2022/2022_1225_A/Geom_T/'// &
         'POR_LS_93x730x160_Z_II_'//frm//'.txt', &
         FORM='FORMATTED',STATUS='REPLACE',ACTION='WRITE')
      WRITE(40,*) dsF(1:py,1:px,1:pz)
      CLOSE(40)

      PRINT*,'Liq done!', N

   ENDDO


END PROGRAM RM_UCNCT_3D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINES                                                                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!=======================================================================================
!! Boundary conditions in 3D
!!=======================================================================================
SUBROUTINE NoGrdBc3D(py, px, pz, f)
   IMPLICIT NONE

   !-----------------------------------------------------------------
   ! Arguments
   !-----------------------------------------------------------------
   ! grid sizes in y, x, z directions
   INTEGER, INTENT(IN) :: py, px, pz
   ! scalar field
   REAL(KIND=8), INTENT(INOUT), DIMENSION(0:py+1,0:px+1,0:pz+1) :: f

   !-----------------------------------------------------------------
   ! No-gradient boundary conditions
   !-----------------------------------------------------------------
   ! y-direction boundaries at i = 0 and i = py + 1
   f(0,1:px,1:pz) = f(1,1:px,1:pz)
   f(py+1,1:px,1:pz) = f(py,1:px,1:pz)

   ! x-direction boundaries at j = 0 and j = px + 1
   f(:,0,1:pz) = f(:,1,1:pz)
   f(:,px+1,1:pz) = f(:,px,1:pz)
   ! z-direction boundaries at k = 0 and k = pz + 1
   f(:,:,0) = f(:,:,1)
   f(:,:,pz+1) = f(:,:,pz)

END SUBROUTINE NoGrdBc3D
!!=======================================================================================
!! Level set reinitialization in 3D
!!=======================================================================================
SUBROUTINE ReDist3D_XYZ(lit, py, px, pz, sg, f)
   IMPLICIT NONE

   !-----------------------------------------------------------------
   ! Arguments
   !-----------------------------------------------------------------
   ! number of iterations for reinitialization
   INTEGER, INTENT(IN) :: lit
   ! grid sizes in y, x, z directions
   INTEGER, INTENT(IN) :: py, px, pz
   ! sign field (guides upwinding) and scalar field
   REAL(KIND=8),INTENT(INOUT),DIMENSION(0:py+1,0:px+1,0:pz+1) :: sg, f

   !-----------------------------------------------------------------
   ! Local variables
   !-----------------------------------------------------------------
   ! iteration and loop variables
   INTEGER :: i, j, k, m
   ! time step
   REAL(KIND=8) :: dt
   ! absolute gradient magnitude and temporary update field
   REAL(KIND=8),DIMENSION(0:py+1,0:px+1,0:pz+1) :: Gf, nf

   !-----------------------------------------------------------------
   ! Initialization
   !-----------------------------------------------------------------
   ! specified time step
   dt = 1.0e-1

   !-----------------------------------------------------------------
   ! Iterative reinitialization loop
   !-----------------------------------------------------------------
   DO m = 1, lit

      !! compute upwind absolute gradient magnitude
      CALL AbsGrd3DUW_XYZ(py, px, pz, sg(0:py+1,0:px+1,0:pz+1), f(0:py+1,0:px+1,0:pz+1), &
         Gf(1:py,1:px,1:pz))

      !! impose no-gradient boundary conditions on Gf
      CALL NoGrdBc3D(py, px, pz, Gf(0:py+1,0:px+1,0:pz+1))

      !! update field with upwind values according to reinitialization equation:
      !! f_t = sign(f0) * (1.0 - |grad(f)|)
      !! field, f, evolves to a signed distance function
      DO k = 0, pz+1
         DO j = 0, px+1
            DO i = 0, py+1
               ! update
               nf(i,j,k) = f(i,j,k) + sg(i,j,k)*(1.0d0-Gf(i,j,k))*dt

               ! sign check
               IF ( nf(i,j,k)*f(i,j,k) > 0.0 ) THEN
               ELSE
                  ! small correction if nf changes sign
                  nf(i,j,k) = f(i,j,k)+dt*1.0d-7*sg(i,j,k)
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      ! print iteration count
      print*, m

      ! copy updated field, nf, back into f for the next iteration
      f(0:py+1,0:px+1,0:pz+1) = nf(0:py+1,0:px+1,0:pz+1)
   ENDDO

END SUBROUTINE ReDist3D_XYZ
!!=======================================================================================
!! Absolute gradient magnitude in 3D
!!=======================================================================================
SUBROUTINE AbsGrd3DUW_XYZ(py, px, pz, sg, f, Gf)
   IMPLICIT NONE

   !-----------------------------------------------------------------
   ! Arguments
   !-----------------------------------------------------------------
   ! grid sizes in y, x, z directions
   INTEGER, INTENT(IN) :: py, px, pz
   ! sign field (guides upwinding) and scalar field
   REAL(KIND=8),INTENT(IN),DIMENSION(0:py+1,0:px+1,0:pz+1) :: sg, f
   ! absolute gradient magnitude
   REAL(KIND=8),INTENT(OUT),DIMENSION(py,px,pz) :: Gf

   !-----------------------------------------------------------------
   ! Local variables
   !-----------------------------------------------------------------
   ! iteration variables
   INTEGER :: i, j, k
   ! directional finite differences
   REAL(KIND=8) :: ag, bg, cg, dg, eg, fg

   !-----------------------------------------------------------------
   ! Loop over all grid points
   !-----------------------------------------------------------------
   DO k = 1, pz
      DO j = 1, px
         DO i = 1, py

            !! compute finite differences in each direction at point (i,j,k)
            ag = (f(i,j,k) - f(i-1,j,k))			! backward difference in y
            bg = (f(i+1,j,k) - f(i,j,k))			! forward difference in y

            cg = (f(i,j,k) - f(i,j-1,k))			! backward difference in x
            dg = (f(i,j+1,k) - f(i,j,k))			! forward difference in x

            eg = (f(i,j,k) - f(i,j,k-1))			! backward difference in z
            fg = (f(i,j,k+1) - f(i,j,k))			! forward difference in z

            !! compute gradient magnitude using upwind differences
            IF ( sg(i,j,k) == 0.0 ) THEN
               ! No propagation: gradient is zero
               Gf(i,j,k) = 0.0e0
            ELSEIF ( sg(i,j,k) > 0.0 ) THEN
               ! Gudonov upwind scheme; Positive sign
               Gf(i,j,k) = SQRT(MAX(MAX(ag,0.0d0)**2,MAX(-bg,0.0d0)**2)+ &
                  MAX(MAX(cg,0.0d0)**2,MAX(-dg,0.0d0)**2)+ &
                  MAX(MAX(eg,0.0d0)**2,MAX(-fg,0.0d0)**2))
            ELSEIF ( sg(i,j,k) < 0.0 ) THEN
               ! Gudonov upwind scheme; Negative sign
               Gf(i,j,k) = SQRT(MAX(MAX(-ag,0.0d0)**2,MAX(bg,0.0d0)**2)+ &
                  MAX(MAX(-cg,0.0d0)**2,MAX(dg,0.0d0)**2)+ &
                  MAX(MAX(-eg,0.0d0)**2,MAX(fg,0.0d0)**2))
            ENDIF
         ENDDO
      ENDDO
   ENDDO

END SUBROUTINE AbsGrd3DUW_XYZ
!!=======================================================================================
