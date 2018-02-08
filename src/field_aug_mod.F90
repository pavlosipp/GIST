MODULE field_aug_mod
  USE type_mod
  USE specifications_mod
  USE lagrange_mod, ONLY: indef, plag2d
  IMPLICIT NONE

  PRIVATE

  INTEGER :: nEDIT=0
  PUBLIC :: field_aug, data_aug
  
  
!================================================================================!
CONTAINS    

SUBROUTINE field_aug(rrr,ppp,zzz,Brad,Bphi,Bzet,dBraddr,dBraddp,dBraddz,&
     &dBphidr,dBphidp,dBphidz,dBzetdr,dBzetdp,dBzetdz)

  USE type_mod

  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: rrr,ppp,zzz
  REAL(DP), INTENT(OUT) :: Brad,Bphi,Bzet,dBraddr,dBraddp,dBraddz
  REAL(DP), INTENT(OUT) :: dBphidr,dBphidp,dBphidz,dBzetdr,dBzetdp,dBzetdz
  INTEGER :: ierr, LPFx=4, first_time_field=0,i,j,icp
  INTEGER, PARAMETER:: nr=129,nz=257
  REAL :: tSHOTswap
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: splpsi
  REAL(DP) :: d2psidz2,d2psidrdz,dpsidz,dpsidr,psif,d2psidr2
  INTEGER, DIMENSION(:), ALLOCATABLE :: imi,ima,jmi,jma
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: ipoint

  REAL(DP),DIMENSION(nr) :: r
  REAL(DP),DIMENSION(nz) :: z
  REAL(DP),DIMENSION(nr,nz):: psi
  REAL(DP):: dr,dz
  REAL(SP) :: Bradial, Bzett, Bt
  REAL(SP) :: Brr, Brz, Bzr,Bzz,Btr,Btz,rhoPF,rhoTF
  REAL(SP) :: fPF, fJp,dum,fTF
  REAL(SP), DIMENSION(0:4) :: PFxx
  REAL(SP), DIMENSION(0:4) :: RPFx, ZPFx
  LOGICAL:: first=.TRUE.,kk=.f.,firstrho=.TRUE.
SAVE
  tSHOTswap=tSHOT
  rmin = 0.75
  rmax = 2.67
  zmin = -1.504
  zmax = 1.504

!-------Get B field and derivatives from subroutine kdrzBrzt-------

kkswitch: IF (kk) THEN

     CALL kdrzBrzt (ierr, EXPAUG, DIAG, shot, nEDIT, tSHOT, &
          & REAL(rrr), REAL(zzz), 1, 4, Bradial, Bzett, Bt, fPF, fJp,&
          & Brr, Brz, Bzr,Bzz,Btr,Btz)

     IF (ierr.NE.0) STOP '(field_aug_mod): !Line hit boundary!'

!B components

     Brad = Bradial
     Bzet = Bzett
     Bphi = Bt
     dBraddp = 0.
     dBphidp = 0.
     dBzetdp = 0.
     dBphidr = Btr
     dBphidz = Btz
     dBraddr = Brr
     dBzetdz = Bzz
     dBraddz = Brz
     dBzetdr = Bzr

ELSE
!-------Get Psi(R,z) grid from subroutine kkrzBrzt------


firsttime: IF (first) THEN
first=.FALSE.
!Create r and z arrays
dr=dble(rmax-rmin)/dble(nr-1)
dz=dble(zmax-zmin)/dble(nz-1)
DO i=1,nr
   r(i)=dble(rmin+(i-1)*dr)
ENDDO
DO i=1,nz
   z(i)=dble(zmin+(i-1)*dz)
ENDDO

DO i=1,nr
   DO j=1,nz
      CALL kkrzBrzt(ierr,expaug,DIAG,shot,nEDIT,tSHOT, REAL(r(i)),REAL(z(j)),1,&
           &dum,dum,dum,fPF,dum)
      psi(i,j)=-dble(fPF)
   ENDDO
ENDDO


psi=psi/2.d0/Pi

DO j=1,nz
   DO i=1,nr
      psi(i,j) =psi(i,j)/r(i)  !D3D: psi+psib
   ENDDO
ENDDO



!------------------------------ Splines (Tracer)  -----------------------------------!

!Rectangular domain
       
   ALLOCATE(imi(nz),ima(nz),jmi(nr),jma(nr))    

   imi=1
   ima=nr
   jmi=1
   jma=nz

!Number of data in splpsi
  
   DO i=1,nz
      IF ( imi(i) > 0 .AND. ima(i) > 0 ) THEN
         icp = icp + ima(i) - imi(i) + 1
      ENDIF
   ENDDO
   ALLOCATE(splpsi(6,6,icp),ipoint(nr,nz))

   CALL s2dcut(nr,nz,dr,dz,psi,imi,ima,jmi,jma,icp,splpsi,ipoint)
!OPEN(55,FILE='splpsiout',status='replace')
!DO i=1,nr
!   DO j=1,nz
!      WRITE(55,*) r(i),z(j),psi(i,j)
!   ENDDO
!ENDDO


ENDIF firsttime



CALL spline(nr,nz,r,z,dr,dz,icp,splpsi,ipoint,rrr,zzz,&
     &psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2,ierr)


!B components

Brad = -dpsidz
Bzet = psif/rrr + dpsidr
Bphi =Btor*rmag/rrr
dBraddp = 0.
dBphidp = 0.
dBzetdp = 0.
dBphidr = -Btor*rmag/rrr**2
dBphidz = 0.
dBraddr = -d2psidrdz
dBzetdz = dpsidz/rrr + d2psidrdz
dBraddz = -d2psidz2
dBzetdr = -psif/rrr**2 + dpsidr/rrr + d2psidr2


ENDIF kkswitch


IF ((firstrho).AND.(swsymm)) THEN
   firstrho=.FALSE.
   !Determine actual rho of traced flux surface for first field line
   CALL kkrzPTFn(iERR,expaug,DIAG,shot,nEDIT,tSHOT,REAL(r_in(1)),REAL(zsym),1,fPF,rhoPF,fTF,rhoTF)
   print *,'Actually traced rho_tor:',rhoTF
   print *,'Actually traced rho_pol:',rhoPF
ENDIF


END SUBROUTINE field_aug




SUBROUTINE data_aug(rrr,zzz,rstart,dens1,dens2,temp1,temp2,omti,omni,omte,omne &
     &,safety,c11,alpha_aug,shat)

USE lagrange_mod, ONLY : plag1d, indef,plag1d_uneven,indef_uneven
USE type_mod
  
  INTEGER, PARAMETER :: Mp=4 !Interpolation (do !not! change)
  REAL(DP), INTENT(IN) :: rrr,zzz
  REAL(DP), INTENT(OUT) :: temp1, temp2, omti, omni, omte, omne, shat,safety
  REAL(DP), INTENT(OUT) :: c11,alpha_aug,rstart,dens1,dens2
  INTEGER, DIMENSION(Mp) :: indrho
  REAL(DP), DIMENSION(Mp) :: q_sten, rhosten,te_sten,ti_sten,ne_sten!,ni_sten
  REAL(DP), DIMENSION(Mp) :: rtot_sten
  REAL(DP),DIMENSION(4,Mp):: rhotorsten
  INTEGER, DIMENSION(4,Mp)::indrhotor
  REAL(DP) :: rmin,rmax,drho,drhotor,overdrho,overdrhotor,polylag,poly1r,rho_min
  REAL(DP) :: te_int,ti_int,ne_int,ni_int,dtidrho,dtedrho,&
       & dnedrho,dnidrho,dqdrho
  REAL(DP) :: rtot_int
  INTEGER :: nmax=50, LPFx=4,nfmax=1024,ine,ite,iti,imax
  INTEGER :: ierr=0, i,j, cnt,first_time=0
  REAL(SP) :: fPF,fJP, fTF,rhoTF,rhoPF, tSHOTswap,dum
  REAL(SP), DIMENSION(50):: r, u, ne, te, ti,qp,rtot
  REAL(SP), DIMENSION(0:4):: PFxx,RPFx,ZPFx
  REAL(SP) :: Rn, zn,Bphi
  CHARACTER(68):: line
  INTEGER,DIMENSION(4)::nrho
!  REAL(SP),DIMENSION(:),ALLOCATABLE:: rho_tor
  REAL(SP),DIMENSION(:,:),ALLOCATABLE:: fit,rho_pol,rho_tor
  CHARACTER(10)::quant
  LOGICAL:: swTi=.TRUE.
  SAVE
nEDIT=1
  tSHOTswap=tSHOT
  rmin = 0.75
  rmax = 2.67
  zmin = -1.504
  zmax = 1.504

  time=tSHOT

first_reading: IF (first_time==0) THEN

   first_time=1

!Refine EQI grid
   
   CALL kkEQintS(iERR,20)
        IF (iERR.NE.0) STOP 'error refining EQI'


!Prepare a set of rho_tor values

   DO cnt=1,nmax
      r(cnt)=1.*REAL(cnt-1)/REAL(nmax-1)
   ENDDO


!Get positions of magn. axis and x-point

   CALL kkEQpfx(iERR,expaug,DIAG,shot,nEDIT,tSHOT,LPFx,PFxx,RPFx,ZPFx)
   IF (ierr/=0) STOP 'Error opening shotfile.'
   rmag=dble(RPFx(0))  
   zmag=dble(ZPFx(0))
   Psi_Sep=dble(PFxx(1))

print *,'Rmag=',rmag

!Get Btor at pos. of magnetic axis
   CALL kkrzBrzt(iERR,expaug,DIAG,shot,nEDIT,tSHOT,RPFx(0),ZPFx(0)& 
        &,1,dum,dum,Bphi,fPF,dum)
   Psi_Ax=dble(fPF)
   Btor=-dble(Bphi)
print *,'Btor=',Btor

!Get Phi_Sep

   CALL kkrhoToP(iERR,expaug,diag,shot,nEDIT,tSHOT,1.,1,rhoPF,fPF,fTF)
Phi_Sep=abs(dble(fTF))
alpha_aug=sqrt(Phi_Sep/Pi/abs(Btor))

!Get rho_pol
   
   CALL kkrhoToP(iERR,expaug,diag,shot,nEDIT,tSHOT,REAL(rrr),1,rhoPF,fPF,fTF)
   print *,'rho_pol=',rhoPF


!Get q values

   DO cnt=1,nmax-1
      CALL kkrhoTPq(iERR,expaug,diag,shot,nEDIT,tSHOT,REAL(r(cnt)),1,qp(cnt),&
           & fPF,fTF)
   ENDDO
   print *,'tEQI=',tSHOT
   qp=abs(qp) !fix negative sign
   tSHOT=tSHOTswap


!Read n, T profiles from augped_output.dat




   ALLOCATE(fit(4,1024),rho_pol(4,1024),rho_tor(4,1024))



   OPEN (66, file='../input/AUG/augped_output.dat', ACTION='READ',iostat=iERR,status='OLD')

   IF (iERR/=0) THEN 
      dbswitch=.FALSE.
      print *,'Profile database not found; Te, Ti, ne will not be given!'
   ENDIF

read_profiles: IF (dbswitch) THEN
imax=4
   loop_data_set: DO i=1,4

      DO WHILE (INDEX(line,"Fit to:").EQ.0) 
         READ(66,'(A)',iostat=iERR) line
!         print *,line
         IF (iERR<0) THEN
            imax=imax-1
            EXIT loop_data_set
         ENDIF
         IF (iERR<0) THEN
            !print *,'err',iERR
            EXIT
         ENDIF
         IF (iERR>0) STOP 'Error reading file.'   
      ENDDO
      BACKSPACE(66)

      READ(66,'(A10)') quant
      quant=quant(9:10)
!      print *,quant
      SELECT CASE (quant)
      CASE ('Te')
         ite=i
!         print *,ite
      CASE ('Ti') 
         iti=i
!         print *,iti
      CASE ('ne') 
         ine=i
!         print *,ine
      END SELECT
         
      
      DO WHILE (INDEX(line,'Rmaj').EQ.0)
         READ(66,*) line
      ENDDO   
      !Read fit data set
      
      DO cnt=1,1024
         
         READ (66,*,iostat=iERR) dum, rho_pol(i,cnt),dum,fit(i,cnt)!,dum
         
         IF (iERR>0) THEN
            print *,dum
            STOP 'error reading profiles'
         ENDIF
         
      ENDDO
         
 

!Convert rho_pol to rho_tor

   DO cnt=1,nfmax
      IF (rho_pol(i,cnt).le.1.0d0) THEN
         CALL kkrhoPTo(iERR,expaug,diag,shot,nEDIT,tSHOT,rho_pol(i,cnt),1,rho_tor(i,cnt)&
              &,fPF,fTF)
      ELSE
         EXIT
      ENDIF
   ENDDO
   IF (iERR/=0) STOP 'error converting data'
   nrho(i)=cnt





   ENDDO loop_data_set
   CLOSE(66)
   OPEN(85,file='check.rho',status='replace')
   DO cnt=1,nrho(1)
      write (85,*) rho_tor(1,cnt),rho_pol(1,cnt)
   ENDDO
   CLOSE(85)
ENDIF read_profiles



!Create rtot array
DO cnt=1,nmax
   CALL kkrhoToP(iERR,expaug,diag,shot,nEDIT,tSHOT,r(cnt),1,rhoPF,fPF,fTF)
   CALL kkrhoRz(iERR,expaug,diag,shot,nEDIT,tSHOT,rhoPF,1,0,edgecorr,rtot(cnt),zn)
ENDDO
print *,'converting rho to R using z=',zn,'(should be close to zsym)'

OPEN(64,file='rtot-rho',status='replace')
DO cnt=1,nmax
   write (64,*) rtot(cnt),r(cnt)
ENDDO
CLOSE(64)

!Check
OPEN(64,file='neprof',status='replace')
DO cnt=1,nrho(ine)
   write (64,*) rho_pol(ine,cnt),fit(ine,cnt)
ENDDO
CLOSE(64)

OPEN(64,file='tiprof',status='replace')
DO cnt=1,nrho(iti)
   write (64,*) rho_pol(iti,cnt),fit(iti,cnt)
ENDDO
CLOSE(64)

OPEN(64,file='teprof',status='replace')
DO cnt=1,nrho(ite)
   write (64,*) rho_pol(ite,cnt),fit(ite,cnt)
ENDDO
CLOSE(64)


ENDIF first_reading




!rho_tor grid, 50 points for data from experiment

rho_min=dble(MINVAL(r))
drho=dble((1.-rho_min)/REAL(nmax-1))
overdrho=1.d0/drho



!Stencil

CALL indef(rrr,rho_min,overdrho,nmax,indrho)


DO j=1,Mp
   rhosten(j)=rho_min+(indrho(j)-1)*drho
   rtot_sten(j)=rtot(indrho(j))
   q_sten(j)=qp(indrho(j))
ENDDO


CALL plag1d(rrr,rtot_sten,overdrho,rhosten,polylag,poly1r)
rtot_int=polylag
drtotdrho=poly1r
print *,'dr/drho=',drtotdrho
CALL plag1d(rrr,q_sten,overdrho,rhosten,polylag,poly1r)
safety=polylag
dqdrho=poly1r


IF (dbswitch) THEN
!rho_tor grid points for fitted data (Te,Ti,ne)
!print *,imax
   DO i=1,imax
      CALL indef_uneven(rrr,dble(rho_tor(i,:)),nrho(i),indrhotor(i,:))
!      print *,nrho(i)
   ENDDO

      DO j=1,Mp
         DO i=1,imax
            rhotorsten(i,j)=dble(rho_tor(i,indrhotor(i,j)))
         ENDDO
         Te_sten(j)=fit(ite,indrhotor(ite,j))
         IF (swTi) Ti_sten(j)=fit(iti,indrhotor(iti,j))
         IF (.NOT.(swTi)) Ti_sten(j)=Te_sten(j)
         ne_sten(j)=fit(ine,indrhotor(ine,j))
      ENDDO



   CALL plag1d_uneven(rrr,Te_sten,rhotorsten(ite,:),polylag,poly1r)
   Te_int=polylag
   dTedrho=poly1r
!print *,'Te'
!print *,Te_int,'|',dTedrho,'|',rhotorsten(ite,:),'|',Te_sten


   CALL plag1d_uneven(rrr,Ti_sten,rhotorsten(iti,:),polylag,poly1r)
   Ti_int=polylag
   dTidrho=poly1r
!print *,'Ti'
!print *,Ti_int,'|',dTidrho,'|',rhotorsten(iti,:),'|',Ti_sten


   CALL plag1d_uneven(rrr,ne_sten,rhotorsten(ine,:),polylag,poly1r)
   ne_int=polylag
   dnedrho=poly1r
!print *,'ne'
!print *,ne_int,'|',dnedrho,'|',rhotorsten(ine,:),'|',ne_sten

   dens2=ne_int
   dens1=dens2 
   temp1=Ti_int/1000 !convert to keV
   temp2=Te_int/1000 
   omti=-1./Ti_int*dTidrho
   omte=-1./Te_int*dTedrho
   omne=-1./ne_int*dnedrho
   omni=omne


ENDIF


!Quantities

IF (.NOT.dbswitch) THEN
   dens2=1.
   dens1=1.
   omni=0.
   omne=0.
   omti=0.
   omte=0.
   temp1=1.
   temp2=1.
ENDIF

   c11=alpha_aug/drtotdrho
   shat=rrr/safety*dqdrho

   rstart=rtot_int !for tracing
!   print *,''
!   print *,'rho=',rrr,'q=',safety,'shat=',shat
!   IF (dbswitch)  print *,'Te=',Te_int,'Ti=',Ti_int,'ne=',ne_int
!   IF (dbswitcH)  print *,'omte=',omte,'omti=',omti,'omne=',omne

END SUBROUTINE data_aug


END MODULE field_aug_mod
 


