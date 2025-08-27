!! ========== !! ========== !! ========== !! ========== !! ========== !! ========== !!
!! This program read the data file of LSC 3D microstructure. Check connenctivity.
!! ========== !! ========== !! ========== !! ========== !! ========== !! ========== !!
PROGRAM LCO3D_MicroStruc_DifSmth_C
IMPLICIT NONE

INTEGER,PARAMETER :: yg = 136, xg = 170*4+50, zg = 236

INTEGER :: py,px,pz,i,j,k,iter,BcCn,fr, N, term
REAL(KIND=4) :: dh,dt,idh2,zeta,Cnbv,res,abU,ttV,MxErr
REAL(KIND=4),DIMENSION(0:yg+1,0:xg+1,0:zg+1) :: psi, Cn
REAL(KIND=8),DIMENSION(yg,xg,zg) :: RdIn
REAL(KIND=4),DIMENSION(yg,xg,zg) :: Div,Ub
CHARACTER(LEN=3) :: frm
CHARACTER(LEN=4) :: ser
CHARACTER(LEN=12) :: flnm


dh = 1.0e0
dt = 1.25e-1 * 0.5
idh2 = 1.0e0/dh**2
zeta = 1.0

py = yg; px = xg; pz = zg
term = 1501

Cnbv = 1.0
ttV = 1.0e0/(py*px*pz)


DO N = 0, 6, 6

	fr = N
	WRITE(ser,"(I4)") 1000+fr
	flnm = 'LHQP'//ser//'.dat'

	!! loading data
	OPEN(UNIT=20,FILE='/mnt/scratch/hcy/SBMES_2025/2025_0818_GraMS/DATs/'//flnm, &
		FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')	
		READ(20) RdIn(1:yg,1:xg,1:zg)
	CLOSE(20)
	
	OPEN(UNIT=80,FILE='/mnt/scratch/hcy/SBMES_2025/2025_0818_GraMS/DATs/'// &
		'LSS_3D_OR_136x730x236_Z_II_'//ser//'.dat', &			
		FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE')
		WRITE(80) sngl(RdIn(1:py,1:px,1:pz))
	CLOSE(80)
				

	psi(1:py,1:px,1:pz) =  SNGL(RdIn(1:py,1:px,1:pz))
	psi(1:py,1:px,1:pz) = 0.5*(1.0+TANH(psi(1:py,1:px,1:pz)/0.75))

	!! BC
	psi(0,1:px,1:pz) = psi(1,1:px,1:pz)
	psi(py+1,1:px,1:pz) = psi(py,1:px,1:pz)
	psi(0:py+1,0,1:pz) = psi(0:py+1,1,1:pz)
	psi(0:py+1,px+1,1:pz) = psi(0:py+1,px,1:pz)
	psi(0:py+1,0:px+1,0) = psi(0:py+1,0:px+1,1)
	psi(0:py+1,0:px+1,pz+1) = psi(0:py+1,0:px+1,pz)

	!! solid phase
	!! check connectivity using accelerated diffusion
	Cn(0:py+1,0:px+1,0:pz+1) = 0.0e0
	Div = 0.0e0

	res = 1.0
	MxErr = 1.0
	!! diffusion in the domain
	DO iter = 2, term 
	!DO WHILE ( res > 1.0d-6 .OR. MxErr > 1.0d-5 .OR. iter < 2001)

		Ub(1:py,1:px,1:pz) = Cn(1:py,1:px,1:pz)

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
	

		!! accelerated diffusion
		IF ( MOD(iter,5) == 1 ) THEN
			WHERE ( Cn(1:py,1:px,1:pz) > 1.0e-4 .AND. psi(1:py,1:px,1:pz) > 4.0e-1 ) Cn(1:py,1:px,1:pz) = 1.0e0
		ENDIF
	
		!! BCs
		Cn(0,1:px,1:pz) = Cn(1,1:px,1:pz)
		Cn(py+1,1:px,1:pz) = Cn(py,1:px,1:pz)
		Cn(0:py+1,0,1:pz) = Cn(0:py+1,1,1:pz)
		Cn(0:py+1,px+1,1:pz) = Cnbv			
		Cn(0:py+1,0:px+1,0) = Cn(0:py+1,0:px+1,1)
		Cn(0:py+1,0:px+1,pz+1) = Cn(0:py+1,0:px+1,pz)
	
		!! check convergence
		abU = SUM(ABS(Cn(1:py,1:px,1:pz)))*ttV
		IF ( abU == 0.0 ) abU = 1.0e-4
		res = SQRT(SUM((Cn(1:py,1:px,1:pz)-Ub(1:py,1:px,1:pz))**2)*ttV)/abU
		MxErr = MAXVAL(ABS(Cn(1:py,1:px,1:pz)-Ub(1:py,1:px,1:pz)))/abU
	
		IF ( MOD(iter,10) == 1 ) PRINT*,'S', N, fr, iter, res, MxErr, abU


		IF ( MOD(iter,term-1) == 1 ) THEN
			WRITE(frm,"(I3)") 100+fr		
			OPEN(UNIT=80,FILE='/mnt/scratch/hcy/SBMES_2025/2025_0818_GraMS/DATs/'// &
				'SOD_3D_CN_136x730x236_Z_II_'//frm//'.dat', &			
				FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE')
				WRITE(80) Cn(1:py,1:px,1:pz)
			CLOSE(80)
		ENDIF
	ENDDO
	
	!! electrolyte phase
	psi(1:py,1:px,1:pz) = 1.0e0 - psi(1:py,1:px,1:pz)
	!! check connectivity using accelerated diffusion
	Cn(0:py+1,0:px+1,0:pz+1) = 0.0e0
	Div = 0.0e0

	res = 1.0
	MxErr = 1.0
	!! diffusion in the domain
	DO iter = 2, term
	!DO WHILE ( res > 1.0d-6 .OR. MxErr > 1.0d-5 .OR. iter < 2001)

		Ub(1:py,1:px,1:pz) = Cn(1:py,1:px,1:pz)

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
	

		!! accelerated diffusion
		IF ( MOD(iter,5) == 1 ) THEN
			WHERE ( Cn(1:py,1:px,1:pz) > 1.0e-4 .AND. psi(1:py,1:px,1:pz) > 4.0e-1 ) Cn(1:py,1:px,1:pz) = 1.0e0
		ENDIF
	
		!! BCs
		Cn(0,1:px,1:pz) = Cn(1,1:px,1:pz)
		Cn(py+1,1:px,1:pz) = Cn(py,1:px,1:pz)
		Cn(0:py+1,0,1:pz) =  Cnbv	
		Cn(0:py+1,px+1,1:pz) = Cn(0:py+1,px,1:pz)			
		Cn(0:py+1,0:px+1,0) = Cn(0:py+1,0:px+1,1)
		Cn(0:py+1,0:px+1,pz+1) = Cn(0:py+1,0:px+1,pz)
	
		!! check convergence
		abU = SUM(ABS(Cn(1:py,1:px,1:pz)))*ttV
		IF ( abU == 0.0 ) abU = 1.0e-4
		res = SQRT(SUM((Cn(1:py,1:px,1:pz)-Ub(1:py,1:px,1:pz))**2)*ttV)/abU
		MxErr = MAXVAL(ABS(Cn(1:py,1:px,1:pz)-Ub(1:py,1:px,1:pz)))/abU
	
		IF ( MOD(iter,10) == 1 ) PRINT*,'L', N, fr, iter,res,MxErr,abU

		IF ( MOD(iter,term-1) == 1 ) THEN
			WRITE(frm,"(I3)") 100+fr
			OPEN(UNIT=80,FILE='/mnt/scratch/hcy/SBMES_2025/2025_0818_GraMS/DATs/'// &
				'LIQ_3D_CN_136x730x236_Z_II_'//frm//'.dat', &			
				FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE')		
				WRITE(80) Cn(1:py,1:px,1:pz)
			CLOSE(80)
		ENDIF
	ENDDO

ENDDO

END PROGRAM LCO3D_MicroStruc_DifSmth_C
!!=======================================================================================




