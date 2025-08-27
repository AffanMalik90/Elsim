!! SMOOTHED BOUNDARY METHOD ELECTROCHEMICAL SIMULATION ANODE 3D

PROGRAM SBMES_3D_MPI_CarGr_ConjG_v01
IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER :: nx, ny, nz, i, j, k, iter, itlp, nlp, cnt, itle, nfr
REAL(KIND=8) :: dx, dt, rad, tm

REAL(KIND=8) :: rho, Xrf, eps
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: Tbl, iTbl, OTbl, DTbl
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: CnP, DvC, DfC, Lap, mu, Mb, pmN

REAL(KIND=8) :: zeta, tPP
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: psP, AvP, AvPx
!!REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: ppN
! ppN - 4th dim is average of psP in each degree of freedom

REAL(KIND=8) :: BvP, CrtP, CvgP, AvrgP
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: Kap, phP, ppO, tmp
!!REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: kpN

REAL(KIND=8) :: De, Dmp, Kpl, t_minus, infx, tPE, CrtE, avCnE
REAL(KIND=8) :: D0, tc1, tc2
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: psE
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: peN
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: CnE, CEo
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: Dab
!!REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: ceN

REAL(KIND=8) :: BvE, CrtL, CvgL, AvrgL, Vsr
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: phE, RHE, peO, pce

REAL(KIND=8) :: Cst1, RMS, Frd, TrgI, Cr, sCrnt
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: Rxn, RHS, RES, ioC, OCV, Kfw, RxC
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: Kbw, dPh

INTEGER :: gy, gx, gz, rnb, cnb, ddR, ddC, lwy, upy, lwz, upz
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: GLP, GLE

REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: AvE
REAL(KIND=8) :: Vsr0, dCrnt
REAL(KIND=4), ALLOCATABLE, DIMENSION(:,:,:) :: DmyT
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: DmyD

REAL(KIND=8) :: Lsm
REAL(KIND=8),ALLOCATABLE :: VaER(:)

!! input and output file names
CHARACTER(LEN=3) :: frm3 
CHARACTER(LEN=4) :: RkNb, frm4
CHARACTER(LEN=33) :: INPFS, INPFL, INPFO
CHARACTER(LEN=6) :: OUTPF

!! variables for material parameters
REAL(KIND=8),ALLOCATABLE :: Vrb(:)


!! for MPI functions
INTEGER :: errcode, rank, np
INTEGER,DIMENSION(MPI_STATUS_SIZE) :: nstatus

REAL(KIND=8) :: DMY, alp, bet, BTM
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: PRJ, WMX, TOR
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: pkN, pcN


!! -- initiate MPI
CALL MPI_INIT(errcode)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,errcode)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,np,errcode)




gy =  272 - 0
gx =  1460 - 0
gz =  472 - 0

!! ============================================ !!
!! -- domain decomposition, only in Y and Z directions
rnb = 16		!! rows of ranks
cnb = 16		!! columns of ranks
CALL MyDmnCmpn_YZ(rank, gy, gz, rnb, cnb, ddR, ddC, ny, lwy, upy, nz, lwz, upz)
!! y index along R; z index along C
nx = gx
!! ============================================ !!

! print*, rank, ny, lwy, upy, nz, lwz, upz

! dx = 3.25* 1.0d-5		!! cm	double precision
! dt = 0.01		!! s
dx = 3.25d-5 *0.5
dt = 1.05625e-02 *0.125	*0.75


zeta = 1.0 *0.5
!! optimize  convergence criteria for speed
CrtE = 1.0d-7	!! convergency criterion for CnE, concentration in electrolyte
CrtP = 1.0d-7	!! convergency criterion for phP, potential in particle
CrtL = 1.0d-6	!! convergency criterion for phE, potential in electrolyte

rho = 0.0312		!! mole/cm^3; lattice density of particle in graphite in Gallagher pg. 2
Cst1 = 1.6021766d-19/(1.3806488d-23*300.0)		!! constant F/(RT), e/(kT)
Frd = 96485.3365		!! Faraday constant

!!!!!!!!!!!!!!!!! tune these params !!!!!!!!!!!!!!!!!

BvP = -0.1					!! Boundary condition for phP
BvE = -0.4495 				!! Boundary condition for phE

Cr = 6.0 					!! charge rate (0.5 -> 2 hours to charge)
! Vsr = 1.0d-2 * 1.0d1     	!! voltage scanning rate
Vsr = 1.0d-2 * 1.0d1 /(3.25**2)    	!! voltage scanning rate
Vsr0 = 1.0d-2 * 1.0d1 /(3.25**2) *2.5


!!!!!!!!!!!!!!!!! tune these params !!!!!!!!!!!!!!!!!

!! D_PF6 = 4.0d-6 and D_Li = 1.25d-6 in Bernardo p 46
De = 2.0d0*(4.0d-6*1.25d-6)/(4.0d-6+1.25d-6)		!! ambipolar diffusivity
Dmp = 4.0d-6 - 1.25d-6								!! D_m - D_p
Kpl = Cst1 * (4.0d-6+1.25d-6)						!! conductivity of electrolyte
t_minus = 4.0d-6/(4.0d-6+1.25d-6)					!! transference number

!! normalize De, Dmp, and Kpl
D0 = (De/2.598606299407882d-06)*0.006667
tc1 = (2*t_minus-1.0)/(2*t_minus*(1.0-t_minus))
tc2 = 1.0/(2*t_minus*(1.0-t_minus))*Cst1

!! -- LiPF6 electrolyte used. D_PF6 = 1.5d-6 and D_Li = 0.73d-6 . Newman p 284
! De = 2.0d0*(0.73d-6*1.5d-6)/(0.73d-6+1.5d-6)		!! ambipolar diffusivity
! Dmp = 1.5d-6 - 0.73d-6								!! D_m - D_p
! Kpl = Cst1 * (0.73d-6+1.5d-6)						!! conductivity of electrolyte
! t_minus = 1.5d-6/(0.73d-6+1.5d-6)					!! transference number

!! Cahn-Hilliard
! eps = 0.64 * (dx**2) 
eps = 6.7600d-10

nfr = 6
WRITE(frm3,"(I3)") 100+nfr

INPFS = 'GrMic_II_DS_272x1460x472_H'//frm3//'.txt'
INPFL = 'GrMic_II_DL_272x1460x472_H'//frm3//'.txt'
INPFO = 'GrMic_II_LO_272x1460x472_H'//frm3//'.txt'
OUTPF = 'GHQ'//frm3

ALLOCATE(psP(0:ny+1,0:nx+1,0:nz+1), AvP(ny,nx,nz), &
	AvPx(ny,nx,nz), psE(0:ny+1,0:nx+1,0:nz+1), peN(ny,nx,nz,6))

!! domain parameter of particle (global)
ALLOCATE(GLP(0:gy+1,0:gx+1,0:gz+1) )
ALLOCATE(DmyT(1:gy,1:gx,1:gz) )
OPEN(UNIT=20,FILE='/mnt/scratch/hcy/SBMES_2025/2025_0818_GraMS/DATs/'//INPFS, &
	FORM='FORMATTED',STATUS='OLD',ACTION='READ')
	READ(20,*) DmyT(1:gy,1:gx,1:gz)
! 	FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')
! 	READ(20) DmyT(1:gy,1:gx,1:gz)
CLOSE(20)

GLP(1:gy,1:gx,1:gz) = DBLE(DmyT(1:gy,1:gx,1:gz))
GLP(1:gy,1:gx,1:gz) = 0.5*(1.0+TANH(GLP(1:gy,1:gx,1:gz)/zeta))

!! inverse phase
! GLP(1:gy,1:gx,1:gz)= 1.0d0 - GLP(1:gy,1:gx,1:gz)

WHERE ( GLP(1:gy,1:gx,1:gz) < 0.0 ) GLP(1:gy,1:gx,1:gz) = 0.0d0
WHERE ( GLP(1:gy,1:gx,1:gz) > 1.0 ) GLP(1:gy,1:gx,1:gz) = 1.0d0

GLP(0,1:gx,1:gz) = GLP(1,1:gx,1:gz)
GLP(gy+1,1:gx,1:gz) = GLP(gy,1:gx,1:gz)
GLP(0:gy+1,0,1:gz) = GLP(0:gy+1,1,1:gz)
GLP(0:gy+1,gx+1,1:gz) = GLP(0:gy+1,gx,1:gz)
GLP(0:gy+1,0:gx+1,0) = GLP(0:gy+1,0:gx+1,1)
GLP(0:gy+1,0:gx+1,gz+1) = GLP(0:gy+1,0:gx+1,gz)

!! local psi
psP(1:ny,1:nx,1:nz) = GLP(lwy:upy,1:gx,lwz:upz)

CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, psP(0:ny+1,0:nx+1,0:nz+1) )

psP(0:ny+1,0,0:nz+1) = psP(0:ny+1,1,0:nz+1)
psP(0:ny+1,nx+1,0:nz+1) = psP(0:ny+1,nx,0:nz+1)

!! domain parameter of electrolyte
! psE(0:ny+1,0:nx+1,0:nz+1) = 1.0 - psP(0:ny+1,0:nx+1,0:nz+1)


OPEN(UNIT=20,FILE='/mnt/scratch/hcy/SBMES_2025/2025_0818_GraMS/DATs/'//INPFL, &
	FORM='FORMATTED',STATUS='OLD',ACTION='READ')
	READ(20,*) DmyT(1:gy,1:gx,1:gz)
! 	FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')
! 	READ(20) DmyT(1:gy,1:gx,1:gz)	
CLOSE(20)

GLP(1:gy,1:gx,1:gz) = DBLE(DmyT(1:gy,1:gx,1:gz))
GLP(1:gy,1:gx,1:gz) = 0.5*(1.0+TANH(GLP(1:gy,1:gx,1:gz)/zeta))

WHERE ( GLP(1:gy,1:gx,1:gz) < 0.0 ) GLP(1:gy,1:gx,1:gz) = 0.0d0
WHERE ( GLP(1:gy,1:gx,1:gz) > 1.0 ) GLP(1:gy,1:gx,1:gz) = 1.0d0

GLP(0,1:gx,1:gz) = GLP(1,1:gx,1:gz)
GLP(gy+1,1:gx,1:gz) = GLP(gy,1:gx,1:gz)
GLP(0:gy+1,0,1:gz) = GLP(0:gy+1,1,1:gz)
GLP(0:gy+1,gx+1,1:gz) = GLP(0:gy+1,gx,1:gz)
GLP(0:gy+1,0:gx+1,0) = GLP(0:gy+1,0:gx+1,1)
GLP(0:gy+1,0:gx+1,gz+1) = GLP(0:gy+1,0:gx+1,gz)

psE(1:ny,1:nx,1:nz) = GLP(lwy:upy,1:gx,lwz:upz)

CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, psE(0:ny+1,0:nx+1,0:nz+1) )

psE(0:ny+1,0,0:nz+1) = psE(0:ny+1,1,0:nz+1)
psE(0:ny+1,nx+1,0:nz+1) = psE(0:ny+1,nx,0:nz+1)

!! particle surface  !! absolute value of gradient psi
AvP(1:ny,1:nx,1:nz) = SQRT(((psP(2:ny+1,1:nx,1:nz)-psP(0:ny-1,1:nx,1:nz))/2)**2 + &
	  				       ((psP(1:ny,2:nx+1,1:nz)-psP(1:ny,0:nx-1,1:nz))/2)**2 + &
						   ((psP(1:ny,1:nx,2:nz+1)-psP(1:ny,1:nx,0:nz-1))/2)**2)

ALLOCATE(AvE(ny,nx,nz))
AvE(1:ny,1:nx,1:nz) = SQRT(((psE(2:ny+1,1:nx,1:nz)-psE(0:ny-1,1:nx,1:nz))/2)**2 + &
	  				       ((psE(1:ny,2:nx+1,1:nz)-psE(1:ny,0:nx-1,1:nz))/2)**2 + &
						   ((psE(1:ny,1:nx,2:nz+1)-psE(1:ny,1:nx,0:nz-1))/2)**2)						   

AvP(1:ny,1:nx,1:nz) = SQRT(AvP(1:ny,1:nx,1:nz)*AvE(1:ny,1:nx,1:nz))

DEALLOCATE(AvE)

WHERE ( AvP(1:ny,1:nx,1:nz) < 1.0d-2 )
	AvP(1:ny,1:nx,1:nz) = 0.0d0
END WHERE

AvPx(1:ny,1:nx,1:nz) = AvP(1:ny,1:nx,1:nz)/dx


!! domain parameters without removing isolated regions
OPEN(UNIT=20,FILE='/mnt/scratch/hcy/SBMES_2025/2025_0818_GraMS/DATs/'//INPFO, &
	FORM='FORMATTED',STATUS='OLD',ACTION='READ')
	READ(20,*) DmyT(1:gy,1:gx,1:gz)
! 	FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')
! 	READ(20) DmyT(1:gy,1:gx,1:gz)	
CLOSE(20)

GLP(1:gy,1:gx,1:gz) = DBLE(DmyT(1:gy,1:gx,1:gz))
GLP(1:gy,1:gx,1:gz) = 0.5*(1.0+TANH(GLP(1:gy,1:gx,1:gz)/zeta))

!! inverse phase
! GLP(1:gy,1:gx,1:gz)= 1.0d0 - GLP(1:gy,1:gx,1:gz)

WHERE ( GLP(1:gy,1:gx,1:gz) < 0.0 ) GLP(1:gy,1:gx,1:gz) = 0.0d0
WHERE ( GLP(1:gy,1:gx,1:gz) > 1.0 ) GLP(1:gy,1:gx,1:gz) = 1.0d0

GLP(0,1:gx,1:gz) = GLP(1,1:gx,1:gz)
GLP(gy+1,1:gx,1:gz) = GLP(gy,1:gx,1:gz)
GLP(0:gy+1,0,1:gz) = GLP(0:gy+1,1,1:gz)
GLP(0:gy+1,gx+1,1:gz) = GLP(0:gy+1,gx,1:gz)
GLP(0:gy+1,0:gx+1,0) = GLP(0:gy+1,0:gx+1,1)
GLP(0:gy+1,0:gx+1,gz+1) = GLP(0:gy+1,0:gx+1,gz)

!! local psi
psP(1:ny,1:nx,1:nz) = GLP(lwy:upy,1:gx,lwz:upz)

CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, psP(0:ny+1,0:nx+1,0:nz+1) )

psP(0:ny+1,0,0:nz+1) = psP(0:ny+1,1,0:nz+1)
psP(0:ny+1,nx+1,0:nz+1) = psP(0:ny+1,nx,0:nz+1)

!! domain parameter of electrolyte
psE(0:ny+1,0:nx+1,0:nz+1) = 1.0 - psP(0:ny+1,0:nx+1,0:nz+1)

DEALLOCATE(DmyT)
DEALLOCATE(GLP)


ALLOCATE(Tbl(1:2,1:101))
OPEN(UNIT=10,FILE='/mnt/home/hcy/SBMES_3D_2023/2023_0727_Gr_Tun/prop/C_Li_X_101.txt', &
	FORM='FORMATTED',STATUS='OLD',ACTION='READ')
READ(10,*) Tbl(1,1:101)
CLOSE(10)

OPEN(UNIT=10,FILE='/mnt/home/hcy/SBMES_3D_2023/2023_0727_Gr_Tun/prop/C_Li_M6_101.txt', &
	FORM='FORMATTED',STATUS='OLD',ACTION='READ')
READ(10,*) Tbl(2,1:101)
CLOSE(10)

ALLOCATE(iTbl(1:2,1:101))
iTbl(1,1:101) = Tbl(1,1:101)
OPEN(UNIT=10,FILE='/mnt/home/hcy/SBMES_3D_2023/2023_0727_Gr_Tun/prop/C_Li_J2_101.txt', &
	FORM='FORMATTED',STATUS='OLD',ACTION='READ')
READ(10,*) iTbl(2,1:101)
CLOSE(10)

ALLOCATE(OTbl(1:2,1:101))
OTbl(1,1:101) = Tbl(1,1:101)
OPEN(UNIT=10,FILE='/mnt/home/hcy/SBMES_3D_2023/2023_0727_Gr_Tun/prop/C_Li_O3_101.txt', &
	FORM='FORMATTED',STATUS='OLD',ACTION='READ')
READ(10,*) OTbl(2,1:101)
CLOSE(10)

ALLOCATE(DTbl(1:2,1:101))
DTbl(1,1:101) = Tbl(1,1:101)
OPEN(UNIT=10,FILE='/mnt/home/hcy/SBMES_3D_2023/2023_0727_Gr_Tun/prop/C_Li_Mb5_101.txt', &
	FORM='FORMATTED',STATUS='OLD',ACTION='READ')
READ(10,*) DTbl(2,1:101)
CLOSE(10)

DTbl(2,1:101) = DTbl(2,1:101)*1.0d2*2/3		!! average
! DTbl(2,1:101) = DTbl(2,1:101)*1.0d2				!! in y and z directions


!! added 1.0d-7 to avoid numerical instability
psP(0:ny+1,0:nx+1,0:nz+1) = psP(0:ny+1,0:nx+1,0:nz+1) + 1.0d-7
Lsm = SUM(psP(1:ny,1:nx,1:nz))							!! total amount of local psi_P

ALLOCATE(VaER(0:np-1))
CALL MPI_ALLGATHER(Lsm,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)
tPP = SUM(VaER(0:np-1))

! print*, rank, Lsm, tPP


!! added 1.0d-7 to avoid numerical instability
psE(0:ny+1,0:nx+1,0:nz+1) = psE(0:ny+1,0:nx+1,0:nz+1) + 1.0d-6

Lsm = SUM(psE(1:ny,1:nx,1:nz))
CALL MPI_ALLGATHER(Lsm,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)
tPE = SUM(VaER(0:np-1))			!! total amount of psi_E
print*, rank, Lsm, tPE

! similar to above
! peN(1:ny,1:nx,1:nz,1) = (psE(0:ny-1,1:nx,1:nz)+psE(1:ny,1:nx,1:nz))*0.5d0
! peN(1:ny,1:nx,1:nz,2) = (psE(2:ny+1,1:nx,1:nz)+psE(1:ny,1:nx,1:nz))*0.5d0
! peN(1:ny,1:nx,1:nz,3) = (psE(1:ny,0:nx-1,1:nz)+psE(1:ny,1:nx,1:nz))*0.5d0
! peN(1:ny,1:nx,1:nz,4) = (psE(1:ny,2:nx+1,1:nz)+psE(1:ny,1:nx,1:nz))*0.5d0
! peN(1:ny,1:nx,1:nz,5) = (psE(1:ny,1:nx,0:nz-1)+psE(1:ny,1:nx,1:nz))*0.5d0
! peN(1:ny,1:nx,1:nz,6) = (psE(1:ny,1:nx,2:nz+1)+psE(1:ny,1:nx,1:nz))*0.5d0



!! Li fraction, Divergence of CnP, diffusivity of Li in particle
! ALLOCATE(CnP(0:ny+1,0:nx+1,0:nz+1), DvC(ny,nx,nz), mu(0:ny+1,0:nx+1,0:nz+1), Lap(nx,ny,nz))
ALLOCATE(CnP(0:ny+1,0:nx+1,0:nz+1), DvC(ny,nx,nz), mu(0:ny+1,0:nx+1,0:nz+1), Lap(ny,nx,nz), &
	pmN(0:ny+1,0:nx+1,0:nz+1), Mb(0:ny+1,0:nx+1,0:nz+1))
CnP(0:ny+1,0:nx+1,0:nz+1) = 2.02d-2	!! initial value of concentration
Lap(1:ny,1:nx,1:nz) = 0.0d0 

!! average concentration in the particle
Lsm = SUM(CnP(1:ny,1:nx,1:nz)*psP(1:ny,1:nx,1:nz))
CALL MPI_ALLGATHER(Lsm,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)
Xrf = SUM(VaER(0:np-1))/tPP

!! salt concentration, old value of CnE,
ALLOCATE(CnE(0:ny+1,0:nx+1,0:nz+1), CEo(ny,nx,nz), Dab(0:ny+1,0:nx+1,0:nz+1) )
CnE(0:ny+1,0:nx+1,0:nz+1) = 1.0d-3 		!! mol/cm^3

!! reaction rate
ALLOCATE(Rxn(ny,nx,nz), RHS(ny,nx,nz), RES(ny,nx,nz), ioC(ny,nx,nz), &
	OCV(ny,nx,nz), Kfw(ny,nx,nz), Kbw(ny,nx,nz), dPh(ny,nx,nz), RxC(ny,nx,nz))
! Rxn(1:ny,1:nx,1:nz) = 1.0d-9
Rxn(1:ny,1:nx,1:nz) = 0.0d0

!! potential in particle, Kap conductivity
ALLOCATE(Kap(0:ny+1,0:nx+1,0:nz+1), phP(0:ny+1,0:nx+1,0:nz+1), &
	ppO(ny,nx,nz), pkN(ny,nx,nz,6) )
phP(0:ny+1,0:nx+1,0:nz+1) = BvP


!! particle conductivity (variable Xi3)
! Kap(0:ny+1,0:nx+1,0:nz+1) = 3.3  ! perpendicular to the plane
! Kap(0:ny+1,0:nx+1,0:nz+1) = 3.0d3  ! parallel to the plane
Kap(0:ny+1,0:nx+1,0:nz+1) = 3.3*psP(0:ny+1,0:nx+1,0:nz+1) 

pkN(1:ny,1:nx,1:nz,1) = (Kap(0:ny-1,1:nx,1:nz)+Kap(1:ny,1:nx,1:nz))*0.5d0
pkN(1:ny,1:nx,1:nz,2) = (Kap(2:ny+1,1:nx,1:nz)+Kap(1:ny,1:nx,1:nz))*0.5d0
pkN(1:ny,1:nx,1:nz,3) = (Kap(1:ny,0:nx-1,1:nz)+Kap(1:ny,1:nx,1:nz))*0.5d0
pkN(1:ny,1:nx,1:nz,4) = (Kap(1:ny,2:nx+1,1:nz)+Kap(1:ny,1:nx,1:nz))*0.5d0
pkN(1:ny,1:nx,1:nz,5) = (Kap(1:ny,1:nx,0:nz-1)+Kap(1:ny,1:nx,1:nz))*0.5d0
pkN(1:ny,1:nx,1:nz,6) = (Kap(1:ny,1:nx,2:nz+1)+Kap(1:ny,1:nx,1:nz))*0.5d0



!! potential in electrolyte
ALLOCATE(RHE(ny,nx,nz), phE(0:ny+1,0:nx+1,0:nz+1), peO(ny,nx,nz), pcN(ny,nx,nz,6), pce(0:ny+1,0:nx+1,0:nz+1) )
phE(0:ny+1,0:nx+1,0:nz+1) = BvE

ALLOCATE(tmp(ny,nx,nz), PRJ(0:ny+1,0:nx+1,0:nz+1), WMX(ny,nx,nz), TOR(ny,nx,nz) )

trgI = 0.95*tPP*rho/(3600/Cr) !! target current

tm = 0.0d0

!! -- file names for output data
WRITE(RkNb,"(I4)") 1000+rank

CALL OutPut3D_A(1, ny, nx, nz, 'PS'//RkNb, 2, psP(1:ny,1:nx,1:nz))
CALL OutPut3D_A(1, ny, nx, nz, 'PL'//RkNb, 2, psE(1:ny,1:nx,1:nz))

sCrnt = trgI

!! time evolution
cnt = 1
iter = 1
! DO iter = 1, 1
DO WHILE ( Xrf < 0.971 )


	dCrnt = ABS(trgI-sCrnt)
	IF ( dCrnt < 0.07 * ABS(trgI) ) THEN
		Vsr = Vsr0*0.5
	ELSEIF ( dCrnt < 0.20 * ABS(trgI) ) THEN
		Vsr = Vsr0
	ELSEIF ( dCrnt < 0.40 * ABS(trgI) ) THEN
		Vsr = Vsr0*3.0
	ELSEIF ( dCrnt < 0.75 * ABS(trgI) ) THEN
		Vsr = Vsr0*5.0
	ELSE
		Vsr = Vsr0*8.0
	ENDIF


	IF ( Xrf < 0.12 ) THEN
		dt = 1.05625e-02 *0.125	*0.75
		Vsr = Vsr*2
	ELSEIF ( Xrf < 0.20 ) THEN
		dt = 1.05625e-02 *0.125	*0.5
	 	Vsr = Vsr*10
	ELSE
		dt = 1.05625e-02 *0.125	*0.5
		Vsr = Vsr*15
	ENDIF
	

	! !  ===============================================
	! !    _____      _   _               _      
	! !   / ____|    | | | |             | |     
	! !  | |     __ _| |_| |__   ___   __| | ___ 
	! !  | |    / _` | __| '_ \ / _ \ / _` |/ _ \
	! !  | |___| (_| | |_| | | | (_) | (_| |  __/
	! !   \_____\__,_|\__|_| |_|\___/ \__,_|\___|
	! ! ===============================================

	!! cathode particle

	!! Laplacian C
	Lap(1:ny,1:nx,1:nz) = ( CnP(1:ny,0:nx-1,1:nz) - 2*CnP(1:ny,1:nx,1:nz) + CnP(1:ny,2:nx+1,1:nz) + &
						    CnP(0:ny-1,1:nx,1:nz) - 2*CnP(1:ny,1:nx,1:nz) + CnP(2:ny+1,1:nx,1:nz) + &
						    CnP(1:ny,1:nx,0:nz-1) - 2*CnP(1:ny,1:nx,1:nz) + CnP(1:ny,1:nx,2:nz+1) )/dx**2

	!! tabulating chemical potential
	CALL Tbl_Chk_Pln_3DA(ny, nx, nz, 101, Tbl(1,1:101), Tbl(2,1:101), CnP(1:ny,1:nx,1:nz), mu(1:ny,1:nx,1:nz))

	CALL Tbl_Chk_Pln_3DA(ny, nx, nz, 101, DTbl(1,1:101), DTbl(2,1:101), CnP(1:ny,1:nx,1:nz), Mb(1:ny,1:nx,1:nz))

	CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, Mb(0:ny+1,0:nx+1,0:nz+1) )
	Mb(0:ny+1,0,0:nz+1) = Mb(0:ny+1,1,0:nz+1)
	Mb(0:ny+1,nx+1,0:nz+1) = Mb(0:ny+1,nx,0:nz+1)

	mu(1:ny,1:nx,1:nz) = mu(1:ny,1:nx,1:nz) -eps*Lap(1:ny,1:nx,1:nz)

	CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, mu(0:ny+1,0:nx+1,0:nz+1) )
	mu(0:ny+1,0,0:nz+1) = mu(0:ny+1,1,0:nz+1)
	mu(0:ny+1,nx+1,0:nz+1) = mu(0:ny+1,nx,0:nz+1)


	!! the term of psi*Mob
	pmN(0:ny+1,0:nx+1,0:nz+1) = psP(0:ny+1,0:nx+1,0:nz+1)*Mb(0:ny+1,0:nx+1,0:nz+1)

	!! calculate the divergency of grad (psi D grad C); Mb in x is reduced.
	DvC(1:ny,1:nx,1:nz) = &
		(0.5*(pmN(2:ny+1,1:nx,1:nz)+pmN(1:ny,1:nx,1:nz))*(mu(2:ny+1,1:nx,1:nz)-mu(1:ny,1:nx,1:nz)) - &
		 0.5*(pmN(1:ny,1:nx,1:nz)+pmN(0:ny-1,1:nx,1:nz))*(mu(1:ny,1:nx,1:nz)-mu(0:ny-1,1:nx,1:nz)) + &
   1.0d0*0.5*(pmN(1:ny,2:nx+1,1:nz)+pmN(1:ny,1:nx,1:nz))*(mu(1:ny,2:nx+1,1:nz)-mu(1:ny,1:nx,1:nz)) - &
   1.0d0*0.5*(pmN(1:ny,1:nx,1:nz)+pmN(1:ny,0:nx-1,1:nz))*(mu(1:ny,1:nx,1:nz)-mu(1:ny,0:nx-1,1:nz)) + &
		 0.5*(pmN(1:ny,1:nx,2:nz+1)+pmN(1:ny,1:nx,1:nz))*(mu(1:ny,1:nx,2:nz+1)-mu(1:ny,1:nx,1:nz)) - &
		 0.5*(pmN(1:ny,1:nx,1:nz)+pmN(1:ny,1:nx,0:nz-1))*(mu(1:ny,1:nx,1:nz)-mu(1:ny,1:nx,0:nz-1)))/dx**2

	WHERE (psP(1:ny,1:nx,1:nz) > 1.0d-5)
		CnP(1:ny,1:nx,1:nz) = CnP(1:ny,1:nx,1:nz) + dt*(DvC(1:ny,1:nx,1:nz) + &
			(Rxn(1:ny,1:nx,1:nz)/rho)*AvPx(1:ny,1:nx,1:nz))/psP(1:ny,1:nx,1:nz)
	END WHERE

	!! no flux BC on X sides
	CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, CnP(0:ny+1,0:nx+1,0:nz+1) )
	CnP(0:ny+1,0,0:nz+1) = CnP(0:ny+1,1,0:nz+1)
	CnP(0:ny+1,nx+1,0:nz+1) = CnP(0:ny+1,nx,0:nz+1)

	!! average concentration in the particle
	Lsm = SUM(CnP(1:ny,1:nx,1:nz)*psP(1:ny,1:nx,1:nz))
	CALL MPI_ALLGATHER(Lsm,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)
	Xrf = SUM(VaER(0:np-1))/tPP

	IF ( MOD(iter,2000) == 1 .AND. rank == 1 ) print*, 'Xrf = ', Xrf, iter

! ! print*, rank, iter, Xrf
 
	! ! ===============================================
	! !   ______ _           _             _       _       
	! !  |  ____| |         | |           | |     | |      
	! !  | |__  | | ___  ___| |_ _ __ ___ | |_   _| |_ ___ 
	! !  |  __| | |/ _ \/ __| __| '__/ _ \| | | | | __/ _ \
	! !  | |____| |  __/ (__| |_| | | (_) | | |_| | ||  __/
	! !  |______|_|\___|\___|\__|_|  \___/|_|\__, |\__\___|
	! !                                       __/ |        
	! !                                      |___/      
	! ! ===============================================  
 
! 	! electrolyte concentration
	CEo(1:ny,1:nx,1:nz) = CnE(1:ny,1:nx,1:nz)
	RHS(1:ny,1:nx,1:nz) = psE(1:ny,1:nx,1:nz)*CEo(1:ny,1:nx,1:nz) - dt*(Rxn(1:ny,1:nx,1:nz))*AvPx(1:ny,1:nx,1:nz)*t_minus

	!! influx of Li
	Lsm = SUM(Rxn(1:ny,1:nx,1:nz)*AvPx(1:ny,1:nx,1:nz))*(dx**3)   !!/(ny*nz*dx**2)
	CALL MPI_ALLGATHER(Lsm,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)
	infx = SUM(VaER(0:np-1))/(gy*gz*dx**2)

! ! 	IF ( MOD(iter,20) == 1 .AND. rank == 1 ) print*, 'infx = ', infx


	!! variable salt diffusivity
	! Dab(0:ny+1,0:nx+1,0:nz+1) = D0*EXP(-7.02 -0.83d3*CnE(0:ny+1,0:nx+1,0:nz+1) +0.05d3*CnE(0:ny+1,0:nx+1,0:nz+1)**2)
	Dab(0:ny+1,0:nx+1,0:nz+1) = D0*EXP(-7.02 -0.83d3*CnE(0:ny+1,0:nx+1,0:nz+1) +0.05*(1.0d3*CnE(0:ny+1,0:nx+1,0:nz+1))**2)
! 	Dab(0:ny+1,0:nx+1,0:nz+1) = De

	!! SBM average
	pce(0:ny+1,0:nx+1,0:nz+1) = psE(0:ny+1,0:nx+1,0:nz+1)*Dab(0:ny+1,0:nx+1,0:nz+1)	
	peN(1:ny,1:nx,1:nz,1) = (pce(0:ny-1,1:nx,1:nz)+pce(1:ny,1:nx,1:nz))*0.5d0
	peN(1:ny,1:nx,1:nz,2) = (pce(2:ny+1,1:nx,1:nz)+pce(1:ny,1:nx,1:nz))*0.5d0
	peN(1:ny,1:nx,1:nz,3) = (pce(1:ny,0:nx-1,1:nz)+pce(1:ny,1:nx,1:nz))*0.5d0
	peN(1:ny,1:nx,1:nz,4) = (pce(1:ny,2:nx+1,1:nz)+pce(1:ny,1:nx,1:nz))*0.5d0
	peN(1:ny,1:nx,1:nz,5) = (pce(1:ny,1:nx,0:nz-1)+pce(1:ny,1:nx,1:nz))*0.5d0
	peN(1:ny,1:nx,1:nz,6) = (pce(1:ny,1:nx,2:nz+1)+pce(1:ny,1:nx,1:nz))*0.5d0	
	
	
	!! back calculate boundary value
! 	CnE(0:ny+1,0,0:nz+1) = infx*t_minus*dx/(0.5*(Dab(0:ny+1,0,0:nz+1)+Dab(0:ny+1,1,0:nz+1))) + CnE(0:ny+1,1,0:nz+1)
	
	
	!! conjugate gradient method
	RES(1:ny,1:nx,1:nz) = RHS(1:ny,1:nx,1:nz) - &
		(psE(1:ny,1:nx,1:nz)*CnE(1:ny,1:nx,1:nz) - dt/dx**2 * &
		(peN(1:ny,1:nx,1:nz,2)*(CnE(2:ny+1,1:nx,1:nz)-CnE(1:ny,1:nx,1:nz)) - &
		 peN(1:ny,1:nx,1:nz,1)*(CnE(1:ny,1:nx,1:nz)-CnE(0:ny-1,1:nx,1:nz)) + &
		 peN(1:ny,1:nx,1:nz,4)*(CnE(1:ny,2:nx+1,1:nz)-CnE(1:ny,1:nx,1:nz)) - &
		 peN(1:ny,1:nx,1:nz,3)*(CnE(1:ny,1:nx,1:nz)-CnE(1:ny,0:nx-1,1:nz)) + &
		 peN(1:ny,1:nx,1:nz,6)*(CnE(1:ny,1:nx,2:nz+1)-CnE(1:ny,1:nx,1:nz)) - &
		 peN(1:ny,1:nx,1:nz,5)*(CnE(1:ny,1:nx,1:nz)-CnE(1:ny,1:nx,0:nz-1))))

	PRJ(1:ny,1:nx,1:nz) = RES(1:ny,1:nx,1:nz)
		
	RMS = 1.0
	itlp = 1
	IF ( iter > 1 ) THEN
	
		!! internal loop
! ! 	DO itlp = 1, 1001	!!0001
		DO WHILE ( RMS > CrtE )

			CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, PRJ(0:ny+1,0:nx+1,0:nz+1) )
			PRJ(0:ny+1,nx+1,0:nz+1) = PRJ(0:ny+1,nx,0:nz+1)
			PRJ(0:ny+1,0,0:nz+1) = PRJ(0:ny+1,1,0:nz+1)
		
			WMX(1:ny,1:nx,1:nz) = psE(1:ny,1:nx,1:nz)*PRJ(1:ny,1:nx,1:nz) - dt/dx**2 * &
				(peN(1:ny,1:nx,1:nz,2)*(PRJ(2:ny+1,1:nx,1:nz)-PRJ(1:ny,1:nx,1:nz)) - &
				 peN(1:ny,1:nx,1:nz,1)*(PRJ(1:ny,1:nx,1:nz)-PRJ(0:ny-1,1:nx,1:nz)) + &
				 peN(1:ny,1:nx,1:nz,4)*(PRJ(1:ny,2:nx+1,1:nz)-PRJ(1:ny,1:nx,1:nz)) - &
				 peN(1:ny,1:nx,1:nz,3)*(PRJ(1:ny,1:nx,1:nz)-PRJ(1:ny,0:nx-1,1:nz)) + &
				 peN(1:ny,1:nx,1:nz,6)*(PRJ(1:ny,1:nx,2:nz+1)-PRJ(1:ny,1:nx,1:nz)) - &
				 peN(1:ny,1:nx,1:nz,5)*(PRJ(1:ny,1:nx,1:nz)-PRJ(1:ny,1:nx,0:nz-1))) 

			RMS = SUM(RES(1:ny,1:nx,1:nz)**2)
			CALL MPI_ALLGATHER(RMS,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)			
			DMY = SUM(VaER(0:np-1))
	
			RMS = SUM(PRJ(1:ny,1:nx,1:nz)*WMX(1:ny,1:nx,1:nz))
			CALL MPI_ALLGATHER(RMS,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)			
			BTM = SUM(VaER(0:np-1))
								
			alp = DMY/BTM
		
			CnE(1:ny,1:nx,1:nz) = CnE(1:ny,1:nx,1:nz) + alp*PRJ(1:ny,1:nx,1:nz)
		
			RES(1:ny,1:nx,1:nz) = RES(1:ny,1:nx,1:nz) - alp*WMX(1:ny,1:nx,1:nz)

			RMS = SUM(RES(1:ny,1:nx,1:nz)**2)
			CALL MPI_ALLGATHER(RMS,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)			
			BTM = SUM(VaER(0:np-1))
		
			bet = BTM/DMY
			PRJ(1:ny,1:nx,1:nz) = RES(1:ny,1:nx,1:nz) + bet*PRJ(1:ny,1:nx,1:nz)
			
		
			Lsm = SUM(RES(1:ny,1:nx,1:nz)**2*psE(1:ny,1:nx,1:nz))
			CALL MPI_ALLGATHER(Lsm,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)
			RMS = SQRT(SUM(VaER(0:np-1))/tPE)			


			if ( mod(itlp,200) == 1 .AND. mod(iter,100) == 1 .AND. rank == 0 ) Print*, 'CnE', iter, itlp, RMS
			IF ( itlp > 50000 ) EXIT

			itlp = itlp + 1
		ENDDO
  	ENDIF
	!! boundary condition
	CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, CnE(0:ny+1,0:nx+1,0:nz+1) )
	CnE(0:ny+1,0,0:nz+1) = infx*t_minus*dx/(0.5*(Dab(0:ny+1,0,0:nz+1)+Dab(0:ny+1,1,0:nz+1))) + CnE(0:ny+1,1,0:nz+1)	
! 	CnE(0:ny+1,0,0:nz+1) = infx*t_minus*dx/De + CnE(0:ny+1,1,0:nz+1)
	CnE(0:ny+1,nx+1,0:nz+1) = CnE(0:ny+1,nx,0:nz+1)
	
	IF (MOD(iter,100) == 1) THEN
		Lsm = SUM(psE(1:ny,1:nx,1:nz)*CnE(1:ny,1:nx,1:nz))   !!/(ny*nz*dx**2)
		CALL MPI_ALLGATHER(Lsm,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)	
		avCnE = SUM(VaER(0:np-1))/tPE
		
		CnE(1:ny,1:nx,1:nz) = CnE(1:ny,1:nx,1:nz) - (avCnE - 1.0d-3)
		CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, CnE(0:ny+1,0:nx+1,0:nz+1) )
		CnE(0:ny+1,0,0:nz+1) = infx*t_minus*dx/(0.5*(Dab(0:ny+1,0,0:nz+1)+Dab(0:ny+1,1,0:nz+1))) + CnE(0:ny+1,1,0:nz+1)	
	! 	CnE(0:ny+1,0,0:nz+1) = infx*t_minus*dx/De + CnE(0:ny+1,1,0:nz+1)
		CnE(0:ny+1,nx+1,0:nz+1) = CnE(0:ny+1,nx,0:nz+1)		
	ENDIF
		

	!! electrolyte potential !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!! variable Dmp
	! Dab(0:ny+1,0:nx+1,0:nz+1) = tc1*D0*EXP(-7.02 -0.83d3*CnE(0:ny+1,0:nx+1,0:nz+1) +0.05d3*CnE(0:ny+1,0:nx+1,0:nz+1)**2)
	Dab(0:ny+1,0:nx+1,0:nz+1) = tc1*D0*EXP(-7.02 -0.83d3*CnE(0:ny+1,0:nx+1,0:nz+1) +0.05*(1.0d3*CnE(0:ny+1,0:nx+1,0:nz+1))**2)
! 	Dab(0:ny+1,0:nx+1,0:nz+1) = Dmp

	!! SBM average
	pce(0:ny+1,0:nx+1,0:nz+1) = psE(0:ny+1,0:nx+1,0:nz+1)*Dab(0:ny+1,0:nx+1,0:nz+1)	
	peN(1:ny,1:nx,1:nz,1) = (pce(0:ny-1,1:nx,1:nz)+pce(1:ny,1:nx,1:nz))*0.5d0
	peN(1:ny,1:nx,1:nz,2) = (pce(2:ny+1,1:nx,1:nz)+pce(1:ny,1:nx,1:nz))*0.5d0
	peN(1:ny,1:nx,1:nz,3) = (pce(1:ny,0:nx-1,1:nz)+pce(1:ny,1:nx,1:nz))*0.5d0
	peN(1:ny,1:nx,1:nz,4) = (pce(1:ny,2:nx+1,1:nz)+pce(1:ny,1:nx,1:nz))*0.5d0
	peN(1:ny,1:nx,1:nz,5) = (pce(1:ny,1:nx,0:nz-1)+pce(1:ny,1:nx,1:nz))*0.5d0
	peN(1:ny,1:nx,1:nz,6) = (pce(1:ny,1:nx,2:nz+1)+pce(1:ny,1:nx,1:nz))*0.5d0	
	
	RHE(1:ny,1:nx,1:nz) = 1.0/dx**2*(peN(1:ny,1:nx,1:nz,2)*(CnE(2:ny+1,1:nx,1:nz)-CnE(1:ny,1:nx,1:nz)) - &
									 peN(1:ny,1:nx,1:nz,1)*(CnE(1:ny,1:nx,1:nz)-CnE(0:ny-1,1:nx,1:nz)) + &
									 peN(1:ny,1:nx,1:nz,4)*(CnE(1:ny,2:nx+1,1:nz)-CnE(1:ny,1:nx,1:nz)) - &
									 peN(1:ny,1:nx,1:nz,3)*(CnE(1:ny,1:nx,1:nz)-CnE(1:ny,0:nx-1,1:nz)) + &
									 peN(1:ny,1:nx,1:nz,6)*(CnE(1:ny,1:nx,2:nz+1)-CnE(1:ny,1:nx,1:nz)) - &
									 peN(1:ny,1:nx,1:nz,5)*(CnE(1:ny,1:nx,1:nz)-CnE(1:ny,1:nx,0:nz-1)))


	!! variable Kpl
	! Dab(0:ny+1,0:nx+1,0:nz+1) = tc2*D0*EXP(-7.02 -0.83d3*CnE(0:ny+1,0:nx+1,0:nz+1) +0.05d3*CnE(0:ny+1,0:nx+1,0:nz+1)**2)
	Dab(0:ny+1,0:nx+1,0:nz+1) = tc2*D0*EXP(-7.02 -0.83d3*CnE(0:ny+1,0:nx+1,0:nz+1) +0.05*(1.0d3*CnE(0:ny+1,0:nx+1,0:nz+1))**2)
! 	Dab(0:ny+1,0:nx+1,0:nz+1) = Kpl

	pce(0:ny+1,0:nx+1,0:nz+1) = psE(0:ny+1,0:nx+1,0:nz+1)*CnE(0:ny+1,0:nx+1,0:nz+1)*Dab(0:ny+1,0:nx+1,0:nz+1)
	pcN(1:ny,1:nx,1:nz,1) = (pce(0:ny-1,1:nx,1:nz)+pce(1:ny,1:nx,1:nz))*0.5d0
	pcN(1:ny,1:nx,1:nz,2) = (pce(2:ny+1,1:nx,1:nz)+pce(1:ny,1:nx,1:nz))*0.5d0
	pcN(1:ny,1:nx,1:nz,3) = (pce(1:ny,0:nx-1,1:nz)+pce(1:ny,1:nx,1:nz))*0.5d0
	pcN(1:ny,1:nx,1:nz,4) = (pce(1:ny,2:nx+1,1:nz)+pce(1:ny,1:nx,1:nz))*0.5d0
	pcN(1:ny,1:nx,1:nz,5) = (pce(1:ny,1:nx,0:nz-1)+pce(1:ny,1:nx,1:nz))*0.5d0
	pcN(1:ny,1:nx,1:nz,6) = (pce(1:ny,1:nx,2:nz+1)+pce(1:ny,1:nx,1:nz))*0.5d0


	!! tabulating exchange current density
	CALL Tbl_Chk_Pln_3DA(ny, nx, nz, 101, iTbl(1,1:101), iTbl(2,1:101), CnP(1:ny,1:nx,1:nz), ioC(1:ny,1:nx,1:nz)) 
	! Change table function to accomodate 99 instead of 101

	ioC = ioC * 1.0d-3 !! in unit of Amp

	CALL Tbl_Chk_Pln_3DA(ny, nx, nz, 101, OTbl(1,1:101), OTbl(2,1:101), CnP(1:ny,1:nx,1:nz), OCV(1:ny,1:nx,1:nz)) 

	Kfw(1:ny,1:nx,1:nz) = ioC(1:ny,1:nx,1:nz)/(Frd*0.001              )*EXP( 0.5*Cst1*OCV(1:ny,1:nx,1:nz))
	Kbw(1:ny,1:nx,1:nz) = ioC(1:ny,1:nx,1:nz)/(Frd*CnP(1:ny,1:nx,1:nz))*EXP(-0.5*Cst1*OCV(1:ny,1:nx,1:nz))

	CvgP = 1.0
	CvgL = 1.0
	nlp = 1
	DO WHILE ( CvgP > 1.0d-8 .AND. CvgL > 1.0d-8 )
! 	DO WHILE ( CvgL > 1.0d-7 )
! 	DO nlp = 1, 1

		WHERE ( AvP(1:ny,1:nx,1:nz) > 1.0d-2 )
			!! in the unit of Volt
			dPh(1:ny,1:nx,1:nz) = phP(1:ny,1:nx,1:nz) - phE(1:ny,1:nx,1:nz)

			!! in the unit of flux
			Rxn(1:ny,1:nx,1:nz) = Kfw(1:ny,1:nx,1:nz)*CnE(1:ny,1:nx,1:nz)*EXP(-0.5*Cst1*dPh(1:ny,1:nx,1:nz)) - &
						          Kbw(1:ny,1:nx,1:nz)*CnP(1:ny,1:nx,1:nz)*EXP( 0.5*Cst1*dPh(1:ny,1:nx,1:nz))

			!! in the unit of Amp
			RxC(1:ny,1:nx,1:nz) = Rxn(1:ny,1:nx,1:nz)*Frd
		END WHERE

		!! particle potential
		ppO(1:ny,1:nx,1:nz) = phP(1:ny,1:nx,1:nz)
		RHS(1:ny,1:nx,1:nz) = -RxC(1:ny,1:nx,1:nz)*AvPx(1:ny,1:nx,1:nz)	!! in unit of Amp

		phP(0:ny+1,nx+1,0:nz+1) = BvP
		!! conjugate gradient method
		RES(1:ny,1:nx,1:nz) = RHS(1:ny,1:nx,1:nz) - 1.0/dx**2* &
			(pkN(1:ny,1:nx,1:nz,2)*(phP(2:ny+1,1:nx,1:nz)-phP(1:ny,1:nx,1:nz)) - &
			 pkN(1:ny,1:nx,1:nz,1)*(phP(1:ny,1:nx,1:nz)-phP(0:ny-1,1:nx,1:nz)) + &
			 pkN(1:ny,1:nx,1:nz,4)*(phP(1:ny,2:nx+1,1:nz)-phP(1:ny,1:nx,1:nz)) - &
			 pkN(1:ny,1:nx,1:nz,3)*(phP(1:ny,1:nx,1:nz)-phP(1:ny,0:nx-1,1:nz)) + &
			 pkN(1:ny,1:nx,1:nz,6)*(phP(1:ny,1:nx,2:nz+1)-phP(1:ny,1:nx,1:nz)) - &
			 pkN(1:ny,1:nx,1:nz,5)*(phP(1:ny,1:nx,1:nz)-phP(1:ny,1:nx,0:nz-1)))
			 
		PRJ(1:ny,1:nx,1:nz) = RES(1:ny,1:nx,1:nz)
		
		RMS = 1.0
		itlp = 1
		DO WHILE ( RMS > CrtP )
! ! 		DO itlp = 1, 10001
		
			tmp(1:ny,1:nx,1:nz) = phP(1:ny,1:nx,1:nz)
			
			CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, PRJ(0:ny+1,0:nx+1,0:nz+1) )
			PRJ(0:ny+1,nx+1,0:nz+1) = 0.0d0
			PRJ(0:ny+1,0,0:nz+1) = PRJ(0:ny+1,1,0:nz+1)
			
			WMX(1:ny,1:nx,1:nz) = 1.0/dx**2* &
				(pkN(1:ny,1:nx,1:nz,2)*(PRJ(2:ny+1,1:nx,1:nz)-PRJ(1:ny,1:nx,1:nz)) - &
				 pkN(1:ny,1:nx,1:nz,1)*(PRJ(1:ny,1:nx,1:nz)-PRJ(0:ny-1,1:nx,1:nz)) + &
				 pkN(1:ny,1:nx,1:nz,4)*(PRJ(1:ny,2:nx+1,1:nz)-PRJ(1:ny,1:nx,1:nz)) - &
				 pkN(1:ny,1:nx,1:nz,3)*(PRJ(1:ny,1:nx,1:nz)-PRJ(1:ny,0:nx-1,1:nz)) + &
				 pkN(1:ny,1:nx,1:nz,6)*(PRJ(1:ny,1:nx,2:nz+1)-PRJ(1:ny,1:nx,1:nz)) - &
				 pkN(1:ny,1:nx,1:nz,5)*(PRJ(1:ny,1:nx,1:nz)-PRJ(1:ny,1:nx,0:nz-1))) 
			
			RMS = SUM(RES(1:ny,1:nx,1:nz)**2)
			CALL MPI_ALLGATHER(RMS,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)			
			DMY = SUM(VaER(0:np-1))
		
			RMS = SUM(PRJ(1:ny,1:nx,1:nz)*WMX(1:ny,1:nx,1:nz))
			CALL MPI_ALLGATHER(RMS,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)			
			BTM = SUM(VaER(0:np-1))
									
			alp = DMY/BTM
			
			phP(1:ny,1:nx,1:nz) = phP(1:ny,1:nx,1:nz) + alp*PRJ(1:ny,1:nx,1:nz)
			
			RES(1:ny,1:nx,1:nz) = RES(1:ny,1:nx,1:nz) - alp*WMX(1:ny,1:nx,1:nz)

			RMS = SUM(RES(1:ny,1:nx,1:nz)**2)
			CALL MPI_ALLGATHER(RMS,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)			
			BTM = SUM(VaER(0:np-1))
			
			bet = BTM/DMY
			PRJ(1:ny,1:nx,1:nz) = RES(1:ny,1:nx,1:nz) + bet*PRJ(1:ny,1:nx,1:nz)
			
			TOR(1:ny,1:nx,1:nz) = phP(1:ny,1:nx,1:nz) - tmp(1:ny,1:nx,1:nz)
! 
! ! 			Lsm = SUM(RES(1:ny,1:nx,1:nz)**2*psP(1:ny,1:nx,1:nz))
			Lsm = SUM(TOR(1:ny,1:nx,1:nz)**2*psP(1:ny,1:nx,1:nz))
			
			CALL MPI_ALLGATHER(Lsm,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)
			RMS = SQRT(SUM(VaER(0:np-1))/tPP)
! 
			if ( mod(itlp,200) == 1 .AND. mod(iter,500) == 1 .AND. rank == 0 ) Print*, 'phP', iter, nlp, itlp, RMS
! 			if ( mod(itlp,20) == 1 .AND. rank == 0 ) Print*, 'phP', iter, nlp, itlp, RMS
! 			
			IF ( itlp > 30000 ) EXIT

			itlp = itlp + 1

		ENDDO
		
		!! BC
		CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, phP(0:ny+1,0:nx+1,0:nz+1) )	
		phP(0:ny+1,nx+1,0:nz+1) = BvP
		phP(0:ny+1,0,0:nz+1) = phP(0:ny+1,1,0:nz+1)			

		Lsm = SUM((phP(1:ny,1:nx,1:nz)-ppO(1:ny,1:nx,1:nz))**2*psP(1:ny,1:nx,1:nz))
		CALL MPI_ALLGATHER(Lsm,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)
		CvgP = SQRT(SUM(VaER(0:np-1))/tPP)

		!! electrolyte potential
		peO(1:ny,1:nx,1:nz) = phE(1:ny,1:nx,1:nz)
		RHS(1:ny,1:nx,1:nz) = (Rxn(1:ny,1:nx,1:nz))*AvPx(1:ny,1:nx,1:nz) + RHE(1:ny,1:nx,1:nz) !! flux

		RES(1:ny,1:nx,1:nz) = RHS(1:ny,1:nx,1:nz) - 1.0/dx**2* &
			(pcN(1:ny,1:nx,1:nz,2)*(phE(2:ny+1,1:nx,1:nz)-phE(1:ny,1:nx,1:nz)) - &
			 pcN(1:ny,1:nx,1:nz,1)*(phE(1:ny,1:nx,1:nz)-phE(0:ny-1,1:nx,1:nz)) + &
			 pcN(1:ny,1:nx,1:nz,4)*(phE(1:ny,2:nx+1,1:nz)-phE(1:ny,1:nx,1:nz)) - &
			 pcN(1:ny,1:nx,1:nz,3)*(phE(1:ny,1:nx,1:nz)-phE(1:ny,0:nx-1,1:nz)) + &
			 pcN(1:ny,1:nx,1:nz,6)*(phE(1:ny,1:nx,2:nz+1)-phE(1:ny,1:nx,1:nz)) - &
			 pcN(1:ny,1:nx,1:nz,5)*(phE(1:ny,1:nx,1:nz)-phE(1:ny,1:nx,0:nz-1)))

		PRJ(1:ny,1:nx,1:nz) = RES(1:ny,1:nx,1:nz)
		PRJ(0:ny+1,0,0:nz+1) = 0.0d0

		RMS = 1.0
		itle = 1
		DO WHILE ( RMS > CrtL )
! 		DO itlp = 1, 10001
! 
			tmp(1:ny,1:nx,1:nz) = phE(1:ny,1:nx,1:nz)

			CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, PRJ(0:ny+1,0:nx+1,0:nz+1) )
			PRJ(0:ny+1,nx+1,0:nz+1) = PRJ(0:ny+1,nx,0:nz+1)
			PRJ(0:ny+1,0,0:nz+1) = 0.0d0

			WMX(1:ny,1:nx,1:nz) = 1.0/dx**2* &
				(pcN(1:ny,1:nx,1:nz,2)*(PRJ(2:ny+1,1:nx,1:nz)-PRJ(1:ny,1:nx,1:nz)) - &
				 pcN(1:ny,1:nx,1:nz,1)*(PRJ(1:ny,1:nx,1:nz)-PRJ(0:ny-1,1:nx,1:nz)) + &
				 pcN(1:ny,1:nx,1:nz,4)*(PRJ(1:ny,2:nx+1,1:nz)-PRJ(1:ny,1:nx,1:nz)) - &
				 pcN(1:ny,1:nx,1:nz,3)*(PRJ(1:ny,1:nx,1:nz)-PRJ(1:ny,0:nx-1,1:nz)) + &
				 pcN(1:ny,1:nx,1:nz,6)*(PRJ(1:ny,1:nx,2:nz+1)-PRJ(1:ny,1:nx,1:nz)) - &
				 pcN(1:ny,1:nx,1:nz,5)*(PRJ(1:ny,1:nx,1:nz)-PRJ(1:ny,1:nx,0:nz-1))) 
			
			RMS = SUM(RES(1:ny,1:nx,1:nz)**2)
			CALL MPI_ALLGATHER(RMS,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)			
			DMY = SUM(VaER(0:np-1))
			
			RMS = SUM(PRJ(1:ny,1:nx,1:nz)*WMX(1:ny,1:nx,1:nz))
			CALL MPI_ALLGATHER(RMS,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)			
			BTM = SUM(VaER(0:np-1))
									
			alp = DMY/BTM

			phE(1:ny,1:nx,1:nz) = phE(1:ny,1:nx,1:nz) + alp*PRJ(1:ny,1:nx,1:nz)
			
			RES(1:ny,1:nx,1:nz) = RES(1:ny,1:nx,1:nz) - alp*WMX(1:ny,1:nx,1:nz)

			RMS = SUM(RES(1:ny,1:nx,1:nz)**2)
			CALL MPI_ALLGATHER(RMS,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)			
			BTM = SUM(VaER(0:np-1))
			
			bet = BTM/DMY
			PRJ(1:ny,1:nx,1:nz) = RES(1:ny,1:nx,1:nz) + bet*PRJ(1:ny,1:nx,1:nz)

			TOR(1:ny,1:nx,1:nz) = phE(1:ny,1:nx,1:nz) - tmp(1:ny,1:nx,1:nz)

! 			Lsm = SUM(RES(1:ny,1:nx,1:nz)**2*psE(1:ny,1:nx,1:nz))
			Lsm = SUM(TOR(1:ny,1:nx,1:nz)**2*psE(1:ny,1:nx,1:nz))
			
			CALL MPI_ALLGATHER(Lsm,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)
			RMS = SQRT(SUM(VaER(0:np-1))/tPE)

! 			if ( mod(itle,200) == 1 .AND. mod(iter,500) == 1 .AND. rank == 1 ) Print*, 'phE', iter, nlp, itle, RMS
			if ( rank == 1 .and. mod(itle,50) == 1 ) Print*, 'phE', iter, nlp, itle, RMS

			IF ( itle > 100000 ) EXIT

			itle = itle + 1

		ENDDO

		!! BC
		CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, phE(0:ny+1,0:nx+1,0:nz+1) )
		phE(0:ny+1,0,0:nz+1) = BvE
		phE(0:ny+1,nx+1,0:nz+1) = phE(0:ny+1,nx,0:nz+1)

		Lsm = SUM((phE(1:ny,1:nx,1:nz)-peO(1:ny,1:nx,1:nz))**2*psE(1:ny,1:nx,1:nz))
		CALL MPI_ALLGATHER(Lsm,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)
		CvgL = SQRT(SUM(VaER(0:np-1))/tPE)


		if ( rank == 3 .AND. mod(iter,50) == 1 ) print*, 'Converg', iter, nlp, itlp, itle, CvgP, CvgL
! 		if ( rank == 3 ) print*, 'Converg', iter, nlp, itlp, itle, CvgP, CvgL

		nlp = nlp + 1

	ENDDO


! 	! =============================================== 
! 	!               _ _           _                                    _   
! 	!      /\      | (_)         | |                                  | |  
! 	!     /  \   __| |_ _   _ ___| |_    ___ _   _ _ __ _ __ ___ _ __ | |_ 
! 	!    / /\ \ / _` | | | | / __| __|  / __| | | | '__| '__/ _ \ '_ \| __|
! 	!   / ____ \ (_| | | |_| \__ \ |_  | (__| |_| | |  | | |  __/ | | | |_ 
! 	!  /_/    \_\__,_| |\__,_|___/\__|  \___|\__,_|_|  |_|  \___|_| |_|\__|
! 	!               _/ |                                                   
! 	!              |__/                                                    
! 	! ===============================================             
             
	Lsm = SUM(Rxn(1:ny,1:nx,1:nz)*AvPx(1:ny,1:nx,1:nz))
	CALL MPI_ALLGATHER(Lsm,1,MPI_DOUBLE_PRECISION,VaER(0:np-1),1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,errcode)
	sCrnt = SUM(VaER(0:np-1))

	BvE = BvE + DSIGN(Vsr,trgI-sCrnt)*dt
	phE(1:ny,1:nx,1:nz) = phE(1:ny,1:nx,1:nz) + DSIGN(Vsr,trgI-sCrnt)*dt

	CALL BcMpiYZ(rank, rnb, cnb, ddR, ddC, ny, nx, nz, phE(0:ny+1,0:nx+1,0:nz+1) )
	phE(0:ny+1,0,1:nz) = BvE
	phE(0:ny+1,nx+1,1:nz) = phE(0:ny+1,nx,1:nz)

	if ( rank == 0 .AND. mod(iter,100) == 1 ) print*, 'Current', iter, trgI, sCrnt, BvE, Xrf
! 	if ( rank == 0 ) print*, 'Current', iter, trgI, sCrnt, BvE, Xrf

	tm = tm + dt

	IF ( Xrf > 0.02 + (cnt-1)*0.005 ) THEN
! 	IF ( MOD(iter,500) == 1 ) THEN

		IF ( MOD(cnt,5) == 1 ) THEN
! 		IF ( MOD(cnt,10) == 1 ) THEN
			CALL OutPut3D_A(cnt, ny, nx, nz, 'CP'//RkNb, 2, CnP(1:ny,1:nx,1:nz))
			CALL OutPut3D_A(cnt, ny, nx, nz, 'CE'//RkNb, 2, CnE(1:ny,1:nx,1:nz))
			CALL OutPut3D_A(cnt, ny, nx, nz, 'PP'//RkNb, 2, phP(1:ny,1:nx,1:nz))
			CALL OutPut3D_A(cnt, ny, nx, nz, 'PE'//RkNb, 2, phE(1:ny,1:nx,1:nz))
! 			CALL OutPut3D_A(cnt, ny, nx, nz, 'PD'//RkNb, 2, dPh(1:ny,1:nx,1:nz))
		ENDIF
	
		IF ( rank == 0 ) THEN
			CALL ReC_TmFX(cnt, OUTPF, iter, tm, Xrf, sCrnt, trgI, BvE)
		ENDIF

		cnt = cnt + 1
	ENDIF
	
	IF ( BvE >= BvP +0.005 .OR. Xrf > 0.97 ) THEN 
		IF ( rank == 0 ) THEN
			CALL ReC_TmFX(cnt, OUTPF, iter, tm, Xrf, sCrnt, trgI, BvE)
		ENDIF
		
		CALL OutPut3D_A(cnt, ny, nx, nz, 'CP'//RkNb, 2, CnP(1:ny,1:nx,1:nz))
		CALL OutPut3D_A(cnt, ny, nx, nz, 'CE'//RkNb, 2, CnE(1:ny,1:nx,1:nz))
		CALL OutPut3D_A(cnt, ny, nx, nz, 'PP'//RkNb, 2, phP(1:ny,1:nx,1:nz))
		CALL OutPut3D_A(cnt, ny, nx, nz, 'PE'//RkNb, 2, phE(1:ny,1:nx,1:nz))		
		
		CALL MPI_BARRIER(MPI_COMM_WORLD,errcode)
	 	EXIT
	ENDIF
	
	iter = iter + 1

ENDDO


CALL MPI_BARRIER(MPI_COMM_WORLD,errcode)
CALL MPI_FINALIZE(errcode)



END PROGRAM SBMES_3D_MPI_CarGr_ConjG_v01
!! ============================================================================================== !!
!! SUBROUTINES
!! ============================================================================================== !!
!! Data IOs
!! ============================================================================================== !!

SUBROUTINE OutPut3D_A(fr,py,px,pz,prfn,AB,Cn1)
IMPLICIT NONE

INTEGER,INTENT(IN) :: fr, px, py, pz, AB
REAL(KIND=8),INTENT(IN),DIMENSION(1:py,1:px,1:pz) :: Cn1
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
		OPEN(UNIT=20,FILE='/mnt/scratch/hcy/SBMES_2025/2025_0818_GraMS/output1/'//flnm, &
			FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE')
			WRITE(20) SNGL(Cn1)
! 			WRITE(20) Cn1
		CLOSE(20)
END SELECT

END SUBROUTINE OutPut3D_A

!! ==================== ==================== ==================== ==================== !!
SUBROUTINE ReC_TmFX(fr, prfn, iter, tm, frx, crnt, Tg, BvE)

IMPLICIT NONE

INTEGER, INTENT(IN) :: fr, iter
REAL(KIND=8), INTENT(In) :: tm, frx, crnt, Tg, BvE
CHARACTER(LEN=6),INTENT(In) :: prfn
CHARACTER(LEN=14) :: flname

flname = prfn//'TmFx.txt'
IF ( fr == 1 ) THEN
	OPEN(UNIT=5,FILE=flname,STATUS='REPLACE',ACTION='WRITE')
		WRITE(5,*) fr, iter, tm, frx, crnt, Tg, BvE
	CLOSE(5)
ELSE
	OPEN(UNIT=6,FILE=flname,STATUS='OLD',ACTION='WRITE',POSITION='APPEND')
		WRITE(6,*) fr, iter, tm, frx, crnt, Tg, BvE
	CLOSE(6)
ENDIF


END SUBROUTINE ReC_TmFX
!! ==================== ==================== ==================== ==================== !!
!! Domain decomposition
!! ==================== ==================== ==================== ==================== !!
SUBROUTINE MyDmnCmpn_YZ(rank,tlyg,tlzg,rnb,cnb,ddR,ddC,py,lwy,upy,pz,lwz,upz)
IMPLICIT NONE
INTEGER,INTENT(IN) :: rank, tlyg, tlzg, rnb, cnb
INTEGER,INTENT(OUT) :: ddR, ddC, py,lwy, upy, pz, lwz, upz
INTEGER :: bivy, rmvy, bivz, rmvz

ddC = INT(rank/rnb)+1
ddR = MOD(rank,rnb)+1

!! determine subdomain size along Z-axis
bivz = INT(tlzg/rnb)
rmvz = MOD(tlzg,rnb)
IF ( (ddR-1) <= rmvz-1 ) THEN
	pz = bivz+1
ELSE
	pz = bivz
ENDIF
!! global index of each subdomain
lwz = (ddR-1)*pz+1
IF ( (ddR-1) >= rmvz ) lwz = lwz+rmvz
upz = lwz+(pz-1)

!! determine subdomain size along X-axis
bivy = INT(tlyg/cnb)
rmvy = MOD(tlyg,cnb)
IF ( (ddC-1) <= rmvy-1 ) THEN
	py = bivy+1
ELSE
	py = bivy
ENDIF
!! global index of each subdomain
lwy = (ddC-1)*py+1
IF ( (ddC-1) >= rmvy ) lwy = lwy+rmvy
upy = lwy+(py-1)

END SUBROUTINE MyDmnCmpn_YZ
!! ==================== ==================== ==================== ==================== !!
!! ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
!! Look-up tables
!! ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
SUBROUTINE Tbl_Chk_Pln_3DA(py, px, pz, sz, Cord, Vals, Cn, fn)
IMPLICIT NONE
INTEGER,INTENT(IN) :: px, py, pz,sz
REAL(KIND=8),INTENT(IN),DIMENSION(sz) :: Cord, Vals
REAL(KIND=8),INTENT(IN),DIMENSION(py,px,pz) :: Cn
REAL(KIND=8),INTENT(OUT),DIMENSION(py,px,pz) :: fn
REAL(KIND=8) :: ditv
INTEGER :: i, j, k, lb
REAL(KIND=8),DIMENSION(py,px,pz) :: Ctmp


Ctmp(1:py,1:px,1:pz) = Cn(1:py,1:px,1:pz)
ditv = Cord(2)-Cord(1)

DO k = 1, pz
	DO j = 1, px
		DO i = 1, py
			
			IF ( Ctmp(i,j,k) < 0.0 ) Ctmp(i,j,k) = 1.0d-6
			IF ( Ctmp(i,j,k) > 1.0 ) Ctmp(i,j,k) = 1.0
			
			! IF ( Ctmp(i,j,k) < 0.0 ) THEN
			! 	fn(i,j,k) = Vals(1) + (Ctmp(i,j,k) - Cord(1))/ditv*(Vals(2) - Vals(1))
			! ELSEIF ( Ctmp(i,j,k) > 1.0 ) THEN
			! 	fn(i,j,k) = Vals(sz-1) + (Ctmp(i,j,k) - Cord(sz-1))/(Cord(sz) - Cord(sz-1))*(Vals(sz) - Vals(sz-1))
			! ELSE
				lb = FLOOR((Ctmp(i,j,k))/ditv) + 0
						
				IF (lb == 0) THEN
					lb = 1
				ELSEIF (lb == sz) THEN
					lb = sz-1
				ENDIF
			
				fn(i,j,k) = Vals(lb) + (Ctmp(i,j,k) - Cord(lb))/ditv*(Vals(lb+1) - Vals(lb))
			! ENDIF

		ENDDO
	ENDDO
ENDDO

END SUBROUTINE Tbl_Chk_Pln_3DA
!! ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
!! ==================== ==================== ==================== ==================== !!
SUBROUTINE BcMpiYZ(rank, rnb, cnb, ddR, ddC, py, px, pz, Udis)
IMPLICIT NONE
INCLUDE 'mpif.h'

INTEGER,INTENT(IN) :: rank, rnb, cnb, ddR, ddC, py, px, pz
REAL(KIND=8),INTENT(INOUT) :: Udis(0:py+1,0:px+1,0:pz+1)
REAL(KIND=8),DIMENSION(py, px) :: BcMpiA, BcMpiB
REAL(KIND=8),DIMENSION(px,0:pz+1) :: BcMpiC, BcMpiD
INTEGER :: errcode
INTEGER,DIMENSION(MPI_STATUS_SIZE) :: nstatus

!! apply BCs along Z-axis using MPI communication (sending up)
BcMpiA(1:py,1:px) = Udis(1:py,1:px,pz)
IF ( ddR /= rnb ) THEN
	CALL MPI_SEND(BcMpiA(1:py,1:px),py*px,MPI_DOUBLE_PRECISION,rank+1,99,MPI_COMM_WORLD,errcode)
ENDIF
IF ( ddR /= 1 ) THEN
	CALL MPI_RECV(BcMpiB(1:py,1:px),py*px,MPI_DOUBLE_PRECISION,rank-1,99,MPI_COMM_WORLD,nstatus,errcode)
ELSE
	BcMpiB(1:py,1:px) = Udis(1:py,1:px,1)
ENDIF
Udis(1:py,1:px,0) = BcMpiB(1:py,1:px)

!! apply BCs along Z-axis using MPI communication (sending down)
BcMpiA(1:py,1:px) = Udis(1:py,1:px,1)
IF ( ddR /= 1 ) THEN
	CALL MPI_SEND(BcMpiA(1:py,1:px),py*px,MPI_DOUBLE_PRECISION,rank-1,98,MPI_COMM_WORLD,errcode)
ENDIF
IF ( ddR /= rnb ) THEN
	CALL MPI_RECV(BcMpiB(1:py,1:px),py*px,MPI_DOUBLE_PRECISION,rank+1,98,MPI_COMM_WORLD,nstatus,errcode)
ELSE
	BcMpiB(1:py,1:px) = Udis(1:py,1:px,pz)
ENDIF
Udis(1:py,1:px,pz+1) = BcMpiB(1:py,1:px)

!! apply BCs along Y-axis using MPI communication (sending right)
BcMpiC(1:px,0:pz+1) = Udis(py,1:px,0:pz+1)
IF ( ddC /= cnb ) THEN
	CALL MPI_SEND(BcMpiC(1:px,0:pz+1),px*(pz+2),MPI_DOUBLE_PRECISION,rank+rnb,97,MPI_COMM_WORLD,errcode)
ENDIF
IF ( ddC /= 1 ) THEN
	CALL MPI_RECV(BcMpiD(1:px,0:pz+1),px*(pz+2),MPI_DOUBLE_PRECISION,rank-rnb,97,MPI_COMM_WORLD,nstatus,errcode)
ELSE
	BcMpiD(1:px,0:pz+1) = Udis(1,1:px,0:pz+1)
ENDIF
Udis(0,1:px,0:pz+1) = BcMpiD(1:px,0:pz+1)

!! apply BCs along X-axis using MPI communication (sending left)
BcMpiC(1:px,0:pz+1) = Udis(1,1:px,0:pz+1)
IF ( ddC /= 1 ) THEN
	CALL MPI_SEND(BcMpiC(1:px,0:pz+1),px*(pz+2),MPI_DOUBLE_PRECISION,rank-rnb,96,MPI_COMM_WORLD,errcode)
ENDIF
IF ( ddC /= cnb ) THEN
	CALL MPI_RECV(BcMpiD(1:px,0:pz+1),px*(pz+2),MPI_DOUBLE_PRECISION,rank+rnb,96,MPI_COMM_WORLD,nstatus,errcode)
ELSE
	BcMpiD(1:px,0:pz+1) = Udis(py,1:px,0:pz+1)
ENDIF
Udis(py+1,1:px,0:pz+1) = BcMpiD(1:px,0:pz+1)

END SUBROUTINE BcMpiYZ
!! ==================== ==================== ==================== ==================== !!
