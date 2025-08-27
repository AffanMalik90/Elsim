!! This program calculates distance function
PROGRAM LevelSet_3D
   IMPLICIT NONE

   !-----------------------------------------------------------------
   ! Grid and geometry parameters
   !-----------------------------------------------------------------
   ! computational domain
   INTEGER,PARAMETER :: yg = 93, xg = 170*4+50, zg = 160
   ! input voxel data
   INTEGER,PARAMETER :: yo = 500, xo = 520, zo = 170
   ! working grid sizes
   INTEGER :: px, py, pz
   ! loop variables
   INTEGER :: i, j, k, iter, n
   ! flag for selecting output mode
   INTEGER :: AB
   ! grid spacing
   REAL(KIND=8) :: dh
   ! tunnel radius and number
   REAL(KIND=8) :: rad, tr

   !-----------------------------------------------------------------
   ! Variables for level-set initialization
   !-----------------------------------------------------------------
   ! frame index
   INTEGER :: fr

   ! sign function and level-set distance function
   REAL(KIND=8), DIMENSION(0:yg+1,0:xg+1,0:zg+1) :: sgF, dsF
   ! temporary arrays
   REAL(KIND=8), DIMENSION(yo,xo,zo) :: TmpA
   REAL(KIND=8), DIMENSION(yo,zo,xo) :: TmpB
   ! table of tunnel radii
   REAL(KIND=8), DIMENSION(0:13) :: TuRd

   !  file naming variables
   CHARACTER(LEN=4) :: prex, frm

   !-----------------------------------------------------------------
   ! Setup parameters
   !-----------------------------------------------------------------
   px = xg
   py = yg
   pz = zg

   ! grid spacing
   dh = 1.0d0
   ! output mode (2 = binary)
   AB = 2

   ! table of tunnel radii
   DATA TuRd /  0,   6.3826e+00,   1.2764e+01,   1.9171e+01,   2.5566e+01, &
      3.1941e+01,   3.8331e+01,   4.4729e+01,   5.1114e+01,   5.7506e+01, &
      6.3890e+01,   7.0280e+01,   7.6674e+01,   8.3061e+01 /

   !-----------------------------------------------------------------
   ! Load sign function from voxelated data
   !-----------------------------------------------------------------
   OPEN(UNIT=20,FILE='../Geom/'// &
      'GraMicVox_II_1_500_520_170.txt',FORM='FORMATTED',STATUS='OLD',ACTION='READ')
   READ(20,*) TmpA(1:yo,1:xo,1:zo)
   CLOSE(20)

   ! set prefix for output files
   prex = 'L24Z'

   !-----------------------------------------------------------------
   ! Loop over multiple tunnel radii
   !-----------------------------------------------------------------
   DO n = 3, 13
      ! frame index
      fr = n
      ! tunnel radius for this frame
      tr = TuRd(n)

      ! initialize distance function
      dsF(1:yg,1:xg,1:zg) = 0.0

      ! set frame name
      WRITE(frm,"(I4)") 1000+fr

      ! Reorient voxel data
      DO k = 1, zo
         TmpB(1:yo,k,1:xo) = TmpA(1:yo,1:xo,zo-k+1)
      ENDDO

      !! impose cylindrical tunnels by zeroing voxel data
      DO k = 1, pz
         DO j = 1, px
            DO i = 1, py
               ! compute radial distance
               rad = SQRT((DBLE(i)-1.0)**2+(DBLE(k)-160.0)**2)
               ! first tunnel
               IF ( rad < tr ) THEN
                  TmpB(i,j,k) = 0.0
               ENDIF

               ! compute radial distance
               rad = SQRT((DBLE(i)-93.0)**2+(DBLE(k)-1.0)**2)
               ! second tunnel
               IF ( rad < tr ) THEN
                  TmpB(i,j,k) = 0.0
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      !! stack repeated thickness in x-direction
      ! first thickness
      dsF(1:yg,50+1:220,1:zg) = TmpB(1:yg,1:zo,1:zg)
      ! second thickness
      DO j = 221, 390
         dsF(1:yg,j,1:zg) = TmpB(1:yg, j-((j-221)*2+51),1:zg)
      ENDDO
      ! third thickness
      dsF(1:yg,391:560,1:zg) = dsF(1:yg,51:220,1:zg)
      ! fourth thickness
      dsF(1:yg,561:730,1:zg) = dsF(1:yg,221:390,1:zg)

      ! impose no-gradient boundary conditions
      CALL NoGrdBc3D(py, px, pz, dsF(0:py+1,0:px+1,0:pz+1))

      !! pre-smooth with Allen-Cahn phase field method
      ! for denoising and smoother interface
      DO iter = 1, 50 !! 200

         sgF(1:yg,1:xg,1:zg) = 2*dsF(1:yg,1:xg,1:zg)*(1.0-dsF(1:yg,1:xg,1:zg))*(1.0-2*dsF(1:yg,1:xg,1:zg)) - &
            0.64*(dsF(2:yg+1,1:xg,1:zg) + dsF(0:yg-1,1:xg,1:zg) + &
            dsF(1:yg,2:xg+1,1:zg) + dsF(1:yg,0:xg-1,1:zg) + &
            dsF(1:yg,1:xg,2:zg+1) + dsF(1:yg,1:xg,0:zg-1) - 6*dsF(1:yg,1:xg,1:zg))

         dsF(1:yg,1:xg,1:zg) = dsF(1:yg,1:xg,1:zg) - 0.05 * sgF(1:yg,1:xg,1:zg)

         CALL NoGrdBc3D(py, px, pz, dsF(0:py+1,0:px+1,0:pz+1))
         print*, 'AC', iter, n
      ENDDO

      !! sign function: +1 outside, -1 inside
      WHERE ( dsF(0:py+1,0:px+1,0:pz+1) > 0.5d0 )
         sgF(0:py+1,0:px+1,0:pz+1) = 1.0d0
      ELSEWHERE
         sgF(0:py+1,0:px+1,0:pz+1) = -1.0d0
      END WHERE

      !! generating distance function by reinitialization
      dsF(0:py+1,0:px+1,0:pz+1) = 0.0
      CALL ReDist3D_XYZ(175, py, px, pz, sgF(0:py+1,0:px+1,0:pz+1), dsF(0:py+1,0:px+1,0:pz+1))

      !! output distance function to file
      CALL OutPut3D(fr, prex, py, px, pz, AB, dsF(1:py,1:px,1:pz))

   ENDDO

END PROGRAM LevelSet_3D


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
   dt = 1.0d-1

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
               nf(i,j,k) = f(i,j,k) + sg(i,j,k) * (1.0d0 - Gf(i,j,k))*dt

               ! sign check
               IF ( nf(i,j,k) * f(i,j,k) > 0.0 ) THEN
               ELSE
                  ! small correction if nf changes sign
                  nf(i,j,k) = f(i,j,k) + dt * 1.0d-7 * sg(i,j,k)
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
SUBROUTINE AbsGrd3DUW_XYZ(py, px, pz, sg, f, gr)
   IMPLICIT NONE

   !-----------------------------------------------------------------
   ! Arguments
   !-----------------------------------------------------------------
   ! grid sizes in y, x, z directions
   INTEGER, INTENT(IN) :: py, px, pz
   ! sign field (guides upwinding) and scalar field
   REAL(KIND=8),INTENT(IN),DIMENSION(0:py+1,0:px+1,0:pz+1) :: sg, f
   ! absolute gradient magnitude
   REAL(KIND=8),INTENT(OUT),DIMENSION(py,px,pz) :: gr

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
               gr(i,j,k) = 0.0d0
            ELSEIF ( sg(i,j,k) > 0.0 ) THEN
               ! Gudonov upwind scheme; Positive sign
               gr(i,j,k) = SQRT(MAX(MAX(ag,0.0d0)**2,MAX(-bg,0.0d0)**2)+ &
                  MAX(MAX(cg,0.0d0)**2,MAX(-dg,0.0d0)**2)+ &
                  MAX(MAX(eg,0.0d0)**2,MAX(-fg,0.0d0)**2))
            ELSEIF ( sg(i,j,k) < 0.0 ) THEN
               ! Gudonov upwind scheme; Negative sign
               gr(i,j,k) = SQRT(MAX(MAX(-ag,0.0d0)**2,MAX(bg,0.0d0)**2)+ &
                  MAX(MAX(-cg,0.0d0)**2,MAX(dg,0.0d0)**2)+ &
                  MAX(MAX(-eg,0.0d0)**2,MAX(fg,0.0d0)**2))
            ENDIF
         ENDDO
      ENDDO
   ENDDO

END SUBROUTINE AbsGrd3DUW_XYZ
!!=======================================================================================
!! Output 3D scalar field
!!=======================================================================================
SUBROUTINE OutPut3D(fr, prfn, py, px, pz, AB, Cn)
   IMPLICIT NONE

   !-----------------------------------------------------------------
   ! Arguments
   !-----------------------------------------------------------------
   ! grid sizes in y, x, z directions
   INTEGER, INTENT(IN) :: py, px, pz
   ! frame number
   INTEGER, INTENT(INOUT) :: fr
   ! 3D scalar field
   REAL(KIND=8), INTENT(IN), DIMENSION(1:py,1:px,1:pz) :: Cn
   ! flag for selecting output mode
   INTEGER, INTENT(IN) :: AB
   ! output file name prefix
   CHARACTER(LEN=4),INTENT(IN) :: prfn

   !-----------------------------------------------------------------
   ! Local variables
   !-----------------------------------------------------------------
   CHARACTER(LEN=4) :: fs
   CHARACTER(LEN=12) :: flnm

   ! create file name string suffix
   WRITE(fs,"(I4)") 1000+fr

   !! select output mode
   ! case 1 ==> ascii; case 2 & 3 ==> binary
   SELECT CASE (AB)
    CASE (1)         ! ASCII output
      flnm = prfn//fs//'.txt'
      OPEN(UNIT=20,FILE=flnm,FORM='FORMATTED',STATUS='REPLACE',ACTION='WRITE')
      !WRITE(20,*) SNGL(Cn)
      WRITE(20,*) Cn
      CLOSE(20)
    CASE (2)         ! binary output; double precision
      flnm = prfn//fs//'.dat'
      OPEN(UNIT=30,FILE='../Geom/'//flnm, &
         FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE')
      WRITE(30) Cn
      CLOSE(30)
    CASE (3)         ! binary output; single precision
      flnm = prfn//fs//'.dat'
      OPEN(UNIT=30,FILE=flnm,FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE')
      WRITE(30) SNGL(Cn)
      CLOSE(30)
   END SELECT

END SUBROUTINE OutPut3D
!!=======================================================================================




