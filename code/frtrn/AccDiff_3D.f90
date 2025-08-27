!! This program checks connectivity of both solid and electrolyte phases
!! uses accelerated diffusion
PROGRAM AccDiff_3D
   IMPLICIT NONE

   !-----------------------------------------------------------------
   ! Grid and geometry parameters
   !-----------------------------------------------------------------
   ! computational domain
   INTEGER,PARAMETER :: yg = 93, xg = 170*4+50, zg = 160
   ! working grid sizes
   INTEGER :: py, px, pz

   !-----------------------------------------------------------------
   ! Loop variables and numerical paramaters
   !-----------------------------------------------------------------
   INTEGER :: iter, N, term, fr
   REAL(KIND=8) :: dh, dt, idh2, Cnbv, res, abU, ttV, MxErr

   !-----------------------------------------------------------------
   ! Field variables
   !-----------------------------------------------------------------
   ! domain parameter
   REAL(KIND=8), DIMENSION(0:yg+1, 0:xg+1, 0:zg+1) :: psi
   ! diffusion field
   REAL(KIND=8), DIMENSION(0:yg+1, 0:xg+1, 0:zg+1) :: Cn
   ! geometry input
   REAL(KIND=8), DIMENSION(yg, xg, zg) :: RdIn
   ! diffusion divergence
   REAL(KIND=8), DIMENSION(yg, xg, zg) :: Div
   ! previous diffusion field
   REAL(KIND=8), DIMENSION(yg, xg, zg) :: Ub

   !-----------------------------------------------------------------
   ! File handling variables
   !-----------------------------------------------------------------
   ! frame identifier for output files
   CHARACTER(LEN=3) :: frm
   ! serial number in input filenames
   CHARACTER(LEN=4) :: ser
   ! geometry input filename
   CHARACTER(LEN=12) :: flnm

   !-----------------------------------------------------------------
   ! Initialization
   !-----------------------------------------------------------------
   ! grid spacing and inverse square
   dh = 1.0e0
   idh2 = 1.0e0/dh**2
   ! time step for diffusion
   dt = 1.25e-1 * 0.5

   ! setup working grid size
   py = yg
   px = xg
   pz = zg

   ! number of diffusion iterations
   term = 1501
   ! boundary value
   Cnbv = 1.0
   ! inverse total volume
   ttV = 1.0e0/(py*px*pz)

   !-----------------------------------------------------------------
   ! Loop over multiple frames (tunnel radii)
   !-----------------------------------------------------------------
   DO N = 3, 13, 1

      ! set frame name
      fr = N
      WRITE(ser,"(I4)") 1000+fr
      flnm = 'L24Z'//ser//'.dat'

      !! load geometry data
      OPEN(UNIT=20,FILE='../Geom/'//flnm, &
         FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')
      READ(20) RdIn(1:yg,1:xg,1:zg)
      CLOSE(20)

      ! convert distance field to domain parameter
      !   psi(1:py,1:px,1:pz) =  SNGL(RdIn(1:py,1:px,1:pz))
      psi(1:py,1:px,1:pz) =  RdIn(1:py,1:px,1:pz)
      psi(1:py,1:px,1:pz) = 0.5*(1.0+TANH(psi(1:py,1:px,1:pz)/0.75))

      !! Neumann boundary conditions
      ! y-direction boundaries at i = 0 and i = py + 1
      psi(0,1:px,1:pz) = psi(1,1:px,1:pz)
      psi(py+1,1:px,1:pz) = psi(py,1:px,1:pz)

      ! x-direction boundaries at j = 0 and j = px + 1
      psi(0:py+1,0,1:pz) = psi(0:py+1,1,1:pz)
      psi(0:py+1,px+1,1:pz) = psi(0:py+1,px,1:pz)

      ! z-direction boundaries at k = 0 and k = pz + 1
      psi(0:py+1,0:px+1,0) = psi(0:py+1,0:px+1,1)
      psi(0:py+1,0:px+1,pz+1) = psi(0:py+1,0:px+1,pz)

      !-----------------------------------------------------------------
      ! Solid phase connectivity using accelerated diffusion
      !-----------------------------------------------------------------
      ! initialize diffusion field and divergence
      Cn(0:py+1,0:px+1,0:pz+1) = 0.0e0
      Div = 0.0e0

      ! initial residual and max error
      res = 1.0
      MxErr = 1.0

      !! diffusion in solid phase
      DO iter = 2, term
         !DO WHILE ( res > 1.0d-6 .OR. MxErr > 1.0d-5 .OR. iter < 2001)

         ! save previous solution for residual calculation
         Ub(1:py,1:px,1:pz) = Cn(1:py,1:px,1:pz)

         !! diffusion update only within threshold of solid phase, psi > 0.4
         WHERE ( psi(1:py,1:px,1:pz) >= 4.0e-1 )
            Div(1:py,1:px,1:pz) = 0.5*idh2/psi(1:py,1:px,1:pz)* &
               ((psi(2:py+1,1:px,1:pz)+psi(1:py,1:px,1:pz))*(Cn(2:py+1,1:px,1:pz)-Cn(1:py,1:px,1:pz))- &
               (psi(1:py,1:px,1:pz)+psi(0:py-1,1:px,1:pz))*(Cn(1:py,1:px,1:pz)-Cn(0:py-1,1:px,1:pz))+ &
               (psi(1:py,2:px+1,1:pz)+psi(1:py,1:px,1:pz))*(Cn(1:py,2:px+1,1:pz)-Cn(1:py,1:px,1:pz))- &
               (psi(1:py,1:px,1:pz)+psi(1:py,0:px-1,1:pz))*(Cn(1:py,1:px,1:pz)-Cn(1:py,0:px-1,1:pz))+ &
               (psi(1:py,1:px,2:pz+1)+psi(1:py,1:px,1:pz))*(Cn(1:py,1:px,2:pz+1)-Cn(1:py,1:px,1:pz))- &
               (psi(1:py,1:px,1:pz)+psi(1:py,1:px,0:pz-1))*(Cn(1:py,1:px,1:pz)-Cn(1:py,1:px,0:pz-1)))

            Cn(1:py,1:px,1:pz) = Cn(1:py,1:px,1:pz) + dt * Div(1:py,1:px,1:pz)
         END WHERE

         !! accelerated diffusion every 5 iterations
         !! diffusion field set to boundary value (1.0) within threshold
         IF ( MOD(iter,5) == 1 ) THEN
            WHERE ( Cn(1:py,1:px,1:pz) > 1.0e-4 .AND. psi(1:py,1:px,1:pz) > 4.0e-1 ) Cn(1:py,1:px,1:pz) = 1.0e0
         ENDIF

         !! boundary conditions for solid diffusion
         ! Neumann boundary condition in y-direction
         Cn(0,1:px,1:pz) = Cn(1,1:px,1:pz)
         Cn(py+1,1:px,1:pz) = Cn(py,1:px,1:pz)

         ! Neumann boundary condition in x-direction at j = 0
         Cn(0:py+1,0,1:pz) = Cn(0:py+1,1,1:pz)
         ! Dirichlet boundary condition at j = px + 1
         Cn(0:py+1,px+1,1:pz) = Cnbv

         ! Neumann boundary condition in z-direction
         Cn(0:py+1,0:px+1,0) = Cn(0:py+1,0:px+1,1)
         Cn(0:py+1,0:px+1,pz+1) = Cn(0:py+1,0:px+1,pz)

         !! check convergence
         abU = SUM(ABS(Cn(1:py,1:px,1:pz)))*ttV
         IF ( abU == 0.0 ) abU = 1.0e-4
         ! compute residual and max error
         res = SQRT(SUM((Cn(1:py,1:px,1:pz) - Ub(1:py,1:px,1:pz))**2)*ttV)/abU
         MxErr = MAXVAL(ABS(Cn(1:py,1:px,1:pz) - Ub(1:py,1:px,1:pz)))/abU

         ! print residuals and max error for every frame for solid phase
         IF ( MOD(iter,10) == 1 ) PRINT*,'S', N, fr, iter, res, MxErr, abU

         ! save solid output at the end of iteration cycle
         IF ( MOD(iter,term-1) == 1 ) THEN
            WRITE(frm,"(I3)") 100+fr
            OPEN(UNIT=80,FILE='../Geom/'// &
               'SOD_3D_CN_93x730x160_Z_II_'//frm//'.dat', &
               FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE')
            WRITE(80) Cn(1:py,1:px,1:pz)
            CLOSE(80)
         ENDIF
      ENDDO

      !-----------------------------------------------------------------
      ! Electrolyte phase connectivity using accelerated diffusion
      !-----------------------------------------------------------------
      ! electrolyte domain parameter
      psi(1:py,1:px,1:pz) = 1.0e0 - psi(1:py,1:px,1:pz)
      ! reset diffusion field and divergence
      Cn(0:py+1,0:px+1,0:pz+1) = 0.0e0
      Div = 0.0e0

      ! reset residuals and max error
      res = 1.0
      MxErr = 1.0

      !! diffusion in electrolyte phase
      DO iter = 2, term
         !DO WHILE ( res > 1.0d-6 .OR. MxErr > 1.0d-5 .OR. iter < 2001)

         ! save previous solution for residual calculation
         Ub(1:py,1:px,1:pz) = Cn(1:py,1:px,1:pz)

         !! diffusion update only within threshold of electrolyte phase, psi > 0.4
         WHERE ( psi(1:py,1:px,1:pz) >= 4.0e-1 )
            Div(1:py,1:px,1:pz) = 0.5*idh2/psi(1:py,1:px,1:pz)* &
               ((psi(2:py+1,1:px,1:pz)+psi(1:py,1:px,1:pz))*(Cn(2:py+1,1:px,1:pz)-Cn(1:py,1:px,1:pz))- &
               (psi(1:py,1:px,1:pz)+psi(0:py-1,1:px,1:pz))*(Cn(1:py,1:px,1:pz)-Cn(0:py-1,1:px,1:pz))+ &
               (psi(1:py,2:px+1,1:pz)+psi(1:py,1:px,1:pz))*(Cn(1:py,2:px+1,1:pz)-Cn(1:py,1:px,1:pz))- &
               (psi(1:py,1:px,1:pz)+psi(1:py,0:px-1,1:pz))*(Cn(1:py,1:px,1:pz)-Cn(1:py,0:px-1,1:pz))+ &
               (psi(1:py,1:px,2:pz+1)+psi(1:py,1:px,1:pz))*(Cn(1:py,1:px,2:pz+1)-Cn(1:py,1:px,1:pz))- &
               (psi(1:py,1:px,1:pz)+psi(1:py,1:px,0:pz-1))*(Cn(1:py,1:px,1:pz)-Cn(1:py,1:px,0:pz-1)))

            Cn(1:py,1:px,1:pz) = Cn(1:py,1:px,1:pz)+dt*Div(1:py,1:px,1:pz)
         END WHERE

         !! accelerated diffusion every 5 iterations
         !! diffusion field set to boundary value (1.0) within threshold
         IF ( MOD(iter,5) == 1 ) THEN
            WHERE ( Cn(1:py,1:px,1:pz) > 1.0e-4 .AND. psi(1:py,1:px,1:pz) > 4.0e-1 ) Cn(1:py,1:px,1:pz) = 1.0e0
         ENDIF

         !! boundary conditions for electrolyte diffusion
         ! Neumann boundary condition in y-direction at i = 0 and i = py + 1
         Cn(0,1:px,1:pz) = Cn(1,1:px,1:pz)
         Cn(py+1,1:px,1:pz) = Cn(py,1:px,1:pz)

         ! Dirichlet boundary condition in x-direction at j = 0
         Cn(0:py+1,0,1:pz) =  Cnbv
         ! Neumann boundary condition in x-direction at j = px + 1
         Cn(0:py+1,px+1,1:pz) = Cn(0:py+1,px,1:pz)

         ! Neumann boundary condition in z-direction at k = 0 and k = pz + 1
         Cn(0:py+1,0:px+1,0) = Cn(0:py+1,0:px+1,1)
         Cn(0:py+1,0:px+1,pz+1) = Cn(0:py+1,0:px+1,pz)

         !! check convergence
         abU = SUM(ABS(Cn(1:py,1:px,1:pz)))*ttV
         IF ( abU == 0.0 ) abU = 1.0e-4
         ! compute residual and max error
         res = SQRT(SUM((Cn(1:py,1:px,1:pz) - Ub(1:py,1:px,1:pz))**2)*ttV)/abU
         MxErr = MAXVAL(ABS(Cn(1:py,1:px,1:pz) - Ub(1:py,1:px,1:pz)))/abU

         ! print residuals and max error for every frame for electrolyte phase
         IF ( MOD(iter,10) == 1 ) PRINT*,'L', N, fr, iter,res,MxErr,abU

         ! save electrolyte output at the end of iteration cycle
         IF ( MOD(iter,term-1) == 1 ) THEN
            WRITE(frm,"(I3)") 100+fr
            OPEN(UNIT=80,FILE='/mnt/scratch/hcy/SBMES_2022/2022_1225_A/Geom_T/'// &
               'LIQ_3D_CN_93x730x160_Z_II_'//frm//'.dat', &
               FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE')
            WRITE(80) Cn(1:py,1:px,1:pz)
            CLOSE(80)
         ENDIF
      ENDDO

   ENDDO

END PROGRAM AccDiff_3D
!!=======================================================================================
