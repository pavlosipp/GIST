
MODULE field_nstx_mod !Used also for JET
  USE type_mod
  USE specifications_mod, ONLY : rmin,rmax,zmin,zmax,rmag,device,drtotdrho
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: field_nstx, data_nstx
  
  REAL(DP), PARAMETER :: per_phi=2*Pi
  INTEGER :: nwEQD,nhEQD !Grid dimensions treated as global parameters

!================================================================================!
CONTAINS           !                MODULE SUBPROGRAMS                           !
!================================================================================!

  SUBROUTINE field_nstx(rr,ppp,zz,Brad,Bphi,Bzet,dBraddr,dBraddp,dBraddz,&
                     &dBphidr,dBphidp,dBphidz,dBzetdr,dBzetdp,dBzetdz)

    REAL(DP), INTENT(IN) :: rr,ppp,zz
    REAL(DP), INTENT(OUT) :: Brad,Bphi,Bzet,dBraddr,dBraddp,dBraddz
    REAL(DP), INTENT(OUT) :: dBphidr,dBphidp,dBphidz,dBzetdr,dBzetdp,dBzetdz
 
    REAL(DP) :: ploc,rm,zm,rrr,zzz
    REAL(DP) :: Brc,Bpc,Bzc,dBrdRc,dBrdpc,dBrdZc,dBpdRc
    REAL(DP) :: dBpdpc,dBpdZc,dBzdRc,dBzdpc,dBzdZc 
    REAL(DP), PARAMETER :: Bscale=1e-4,Rscale=1e-2 

!--------------------------------- PERIODICITY ----------------------------------! 
    ploc=ppp

    IF (ploc >= per_phi) ploc = ploc - (INT(ploc/per_phi))*per_phi
    IF (ploc < 0.) ploc = ploc + (INT(ABS(ploc)/per_phi) +1)*per_phi
!--------------------------------------------------------------------------------!

!CGS units

    rrr=rr/Rscale ; zzz=zz/Rscale

    rm=rrr ; zm=zzz

!Equilibrium file

    CALL field_eq(rm,ploc,zm,Brad,Bphi,Bzet,dBraddr,dBraddp,dBraddz,&
                 &dBphidr,dBphidp,dBphidz,dBzetdr,dBzetdp,dBzetdz) 

!SI units
   
    Brad=Brad*Bscale              ; Bphi=Bphi*Bscale       ; Bzet=Bzet*Bscale 
    dBraddr=dBraddr*Bscale/Rscale ; dBraddp=dBraddp*Bscale ; dBraddz=dBraddz*Bscale/Rscale
    dBphidr=dBphidr*Bscale/Rscale ; dBphidp=dBphidp*Bscale ; dBphidz=dBphidz*Bscale/Rscale
    dBzetdr=dBzetdr*Bscale/Rscale ; dBzetdp=dBzetdp*Bscale ; dBzetdz=dBzetdz*Bscale/Rscale


  END SUBROUTINE field_nstx

!============================== READ EFIT FILE ==================================!

  SUBROUTINE field_eq(rrr,ppp,zzz,Brad,Bphi,Bzet,dBraddr,dBraddp,dBraddz,&             
                     &dBphidr,dBphidp,dBphidz,dBzetdr,dBzetdp,dBzetdz)


    REAL(DP), INTENT(IN) :: rrr,ppp,zzz
    REAL(DP), INTENT(OUT) :: Brad,Bphi,Bzet,dBraddr,dBraddp,dBraddz
    REAL(DP), INTENT(OUT) :: dBphidr,dBphidp,dBphidz,dBzetdr,dBzetdp,dBzetdz

    REAL(DP), DIMENSION(:), ALLOCATABLE :: rad,zet,xi,f
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: psi
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: splpsi
    REAL(DP) :: dummy,hrad,hzet,btf,rtf,psib
    INTEGER, DIMENSION(:), ALLOCATABLE :: imi,ima,jmi,jma
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ipoint
    INTEGER :: icp=0,i,j,ierr
    REAL(DP) :: d2psidz2,d2psidrdz,dpsidz,dpsidr,psif,d2psidr2
    LOGICAL :: first=.TRUE.
    SAVE
  
    first_reading:  IF (first) THEN
  
       first=.FALSE.

!Get grid dimensions 

       CALL dimeq(nwEQD,nhEQD)

       ALLOCATE(rad(nwEQD),zet(nhEQD))
       ALLOCATE(psi(nwEQD,nhEQD))
     
!Read equilibrium file

       CALL eqfile(psib,btf,rtf,rad,zet,psi)
   
!CGS units     

       rad = rad*100. 
       zet = zet*100. 
       rtf = rtf*100. 
       psi = psi*1.e8
       psib= psib*1.e8
       btf = btf*1.e4

! Physical $\phi$ component of vector potential:

       DO j=1,nhEQD
          DO i=1,nwEQD
             psi(i,j) = (psi(i,j) + psib)/rad(i)
          ENDDO
       ENDDO
!      OPEN(42,FILE='psiout',status='replace')
!      DO i=1,nwEQD
!         DO j=1,nhEQD
!            WRITE(42,*) rad(i),zet(j),psi(i,j)
!         ENDDO
!      ENDDO
!      CLOSE(42)


!OPEN(40,FILE='psiout',status='replace')
!DO i=1,nwEQD
!   DO j=1,nhEQD
!      WRITE(40,*) rad(i),' ',zet(j),' ',psi(i,j)
!   ENDDO
!ENDDO
!CLOSE(40)
!Resolution

       hrad = rad(2) - rad(1)
       hzet = zet(2) - zet(1)

!------------------------------ Splines -----------------------------------!

!Rectangular domain

       ALLOCATE(imi(nhEQD),ima(nhEQD),jmi(nwEQD),jma(nwEQD))    

       imi = 1
       ima = nwEQD
       jmi = 1
       jma = nhEQD

!Number of data in splpsi
  
       DO i=1,nhEQD
          IF ( imi(i) > 0 .AND. ima(i) > 0 ) THEN
             icp = icp + ima(i) - imi(i) + 1
          ENDIF
       ENDDO
   
       ALLOCATE(splpsi(6,6,icp),ipoint(nwEQD,nhEQD))

       CALL s2dcut(nwEQD,nhEQD,hrad,hzet,psi,imi,ima,jmi,jma,icp,splpsi,ipoint)
!print *,ipoint
    ENDIF first_reading

!---------------------------------------------------------------------------!
!print *,hrad,hzet
!print *,rad
!print *,zet
!print *,icp
!print *,rrr,zzz
    CALL spline(nwEQD,nhEQD,rad,zet,hrad,hzet,icp,splpsi,ipoint,rrr,zzz,&
               &psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2,ierr)

!B components

    Brad = -dpsidz
    Bzet = psif/rrr + dpsidr
    Bphi = -btf*rtf/rrr
    dBraddp = 0.
    dBphidp = 0.
    dBzetdp = 0.
    dBphidr = +btf*rtf/rrr**2
    dBphidz = 0.
    dBraddr = -d2psidrdz
    dBzetdz = dpsidz/rrr + d2psidrdz
    dBraddz = -d2psidz2
    dBzetdr = -psif/rrr**2 + dpsidr/rrr + d2psidr2

  END SUBROUTINE field_eq

!================================= READ R-Z DIMENSIONS =========================!

  SUBROUTINE dimeq(nw,nh)

    INTEGER, INTENT(OUT) :: nw,nh

    INTEGER :: idum,i
    CHARACTER(10), DIMENSION(6) ::  cdum

    !Either NSTX or JET 
    IF (device == 7) THEN
       OPEN(UNIT=11, FILE='../mesh/mesh_jet', STATUS="OLD", ACTION="READ", ERR=1000)
    ELSE
       OPEN(UNIT=11, FILE='../mesh/mesh_nstx', STATUS="OLD", ACTION="READ")
    ENDIF

    READ(11,"(6A8,3I4)")(cdum(i),i=1,6),idum,nw,nh
    CLOSE(11)

  
    RETURN
    
1000 STOP 'DIIID FILE MISSING!!'

  END SUBROUTINE dimeq

!============================ READ EQUILIBRIUM FIELD ============================!

  SUBROUTINE eqfile(psiSep,bt0,rzero,rad,zet,psiRZ)
    USE specifications_mod, ONLY: Btor
    REAL(DP), INTENT(OUT) :: psiSep,bt0,rzero
    REAL(DP), DIMENSION(nwEQD), INTENT(OUT) :: rad
    REAL(DP), DIMENSION(nhEQD), INTENT(OUT) :: zet
    REAL(DP), DIMENSION(nwEQD,nhEQD), INTENT(OUT) :: psiRZ
    REAL(DP) :: zmid,r1,zdim,xdim,psiAxis,zmaxis,rmaxis,xdum,plas_cur
    REAL(DP), DIMENSION(nwEQD) :: fpol,pres,ffprim,pprime,qpsi
    REAL(DP), DIMENSION(:), ALLOCATABLE :: LCFS, limEQD
    CHARACTER(10), DIMENSION(6) :: cdum
    INTEGER :: n_bndyxy,nlimEQD,idum
    INTEGER :: i,j

    !Either NSTX or JET 
    IF (device == 7) THEN
       OPEN(UNIT=11, FILE='../mesh/mesh_jet', STATUS="OLD", ACTION="READ")
    ELSE
       OPEN(UNIT=11, FILE='../mesh/mesh_nstx', STATUS="OLD", ACTION="READ")
    ENDIF

    READ(11,"(6A8,3I4)")(cdum(i),i=1,6),idum,idum,idum
    READ(11,"(5E16.9)") xdim,zdim,rzero,r1,zmid
    READ(11,"(5E16.9)") rmaxis,zmaxis,psiAxis,psiSep,bt0
    READ(11,"(5E16.9)") plas_cur,psiAxis,xdum,rmaxis,xdum
    READ(11,"(5E16.9)") zmaxis,xdum,psiSep,xdum,xdum
    READ(11,"(5E16.9)") (fpol(i),i=1,nwEQD)
    READ(11,"(5E16.9)") (pres(i),i=1,nwEQD)
    READ(11,"(5E16.9)") (ffprim(i),i=1,nwEQD)
    READ(11,"(5E16.9)") (pprime(i),i=1,nwEQD)
    READ(11,"(5E16.9)") ((psiRZ(i,j),i=1,nwEQD),j=1,nhEQD)
    READ(11,"(5E16.9)") (qpsi(i),i=1,nwEQD)

!Export Rmag and Btor

    Rmag=rmaxis
    Btor=-bt0

!Boundary Data

    READ(11,*) n_bndyxy,nlimEQD    

    ALLOCATE(LCFS(2*n_bndyxy),limEQD(2*nlimEQD))
               
    READ(11,"(5E16.9)") (LCFS(i),i=1,2*n_bndyxy)
    READ(11,"(5E16.9)") (limEQD(i),i=1,2*nlimEQD)

    CLOSE(11)

    CALL set_eqcoords(xdim,zdim,r1,zmid,rad,zet)           

  END SUBROUTINE eqfile

!======================== CALCULATE R-Z COORDINATES =======================!

  SUBROUTINE set_eqcoords(xdim,zdim,r1,zmid,rad,zet)
    
    REAL(DP), INTENT(IN) :: xdim,zdim,r1,zmid
    REAL(DP), DIMENSION(nwEQD), INTENT(OUT) :: rad
    REAL(DP), DIMENSION(nhEQD), INTENT(OUT) :: zet
 
    REAL(DP) :: z1
    INTEGER :: j,k
    INTRINSIC REAL

    DO j=1,nwEQD
       rad(j) = r1 + (j-1)*(xdim/REAL(nwEQD-1))
    ENDDO

    z1 = zmid - zdim/2.		

    DO k=1,nhEQD			
       zet(k) = z1 + (k-1)*(zdim/REAL(nhEQD-1))
    ENDDO

    rmin=MINVAL(rad) ; rmax=MAXVAL(rad)
    zmin=MINVAL(zet) ; zmax=MAXVAL(zet)

  END SUBROUTINE set_eqcoords

!============================================================================!

SUBROUTINE data_nstx(rrr,rstart,temp1,omti,dens1,omni,temp2,omte,dens2,&
     &omne,c11,alpha_nstx)
  USE specifications_mod, ONLY: alpha, shot, time, Phi_Sep,Btor
  USE lagrange_mod, ONLY:indef_uneven,plag1d_uneven,plag1d,indef
  REAL(DP), INTENT(IN):: rrr
  REAL(DP), INTENT(OUT):: rstart, temp1, omti, dens1, omni, temp2, omte&
       &,dens2, omne, c11, alpha_nstx
  

  INTEGER, PARAMETER:: nval=20, Mp=4
  INTEGER:: i,j
  INTEGER, DIMENSION(Mp):: indr,indrho
  CHARACTER(46):: line
  REAL(DP), DIMENSION(nval):: rho, ne, Te, Ti, Zeff, rtot
  REAL(DP), DIMENSION(Mp):: rhosten,nesten,Testen,Tisten,Zeffsten, rtotsten
  REAL(DP) :: polylag, poly1r, ne_int, dnedrho, Te_int, dTedrho, &
       &Ti_int, dTidrho, Zeff_int, dZdrho, drho, overdrho, rtot_int
  REAL(DP):: dum1,dum2,dum3
  REAL(DP), DIMENSION(:,:), ALLOCATABLE:: dum6
  REAL(DP), DIMENSION(:), ALLOCATABLE:: dum4, dum5
  LOGICAL:: first=.TRUE.
  
IF (first) THEN  
   first=.FALSE.
  OPEN(89,file='../input/NSTX/profiles.inp',status='old')
  READ(89,*) line
  READ(89,*)
  DO i=1,nval
     READ(89,*) rho(i), ne(i), Te(i), Ti(i), Zeff(i), rtot(i)
  ENDDO

CALL dimeq(nwEQD,nhEQD)
ALLOCATE(dum6(nwEQD,nhEQD),dum4(nwEQD),dum5(nwEQD))
CALL eqfile(dum1,dum2,dum3,dum4,dum5,dum6)
DEALLOCATE(dum4,dum5,dum6)

Phi_Sep=dble(1.13)
  
ENDIF 

!testing
  drho=(MAXVAL(rho)-MINVAL(rho))/dble(nval-1)
  overdrho=1./drho

!gradients wrt rho
  CALL indef(rrr,MINVAL(rho),overdrho,nval,indr)
  
  DO i=1,Mp
     rtotsten(i)=rtot(indr(i))
     rhosten(i)=rho(indr(i))
     nesten(i)=ne(indr(i))
     Testen(i)=Te(indr(i))
     Tisten(i)=Ti(indr(i))
     Zeffsten(i)=Zeff(indr(i))
  ENDDO
  
  
  CALL plag1d(rrr,nesten,overdrho,rhosten,polylag,poly1r)
  ne_int=polylag
  dnedrho=poly1r
  
  CALL plag1d(rrr,Testen,overdrho,rhosten,polylag,poly1r)
  Te_int=polylag
  dTedrho=poly1r
  
  CALL plag1d(rrr,Tisten,overdrho,rhosten,polylag,poly1r)
  Ti_int=polylag
  dTidrho=poly1r
  
  CALL plag1d(rrr,Zeffsten,overdrho,rhosten,polylag,poly1r)
  Zeff_int=polylag
  dZdrho=poly1r

  CALL plag1d(rrr,rtotsten,overdrho,rhosten,polylag,poly1r)
  rtot_int=polylag
  drtotdrho=poly1r

  dens2=ne_int*1.e19
  dens1=dens2 
  temp1=Ti_int/1000 !convert to keV
  temp2=Te_int/1000 
  omti=-1./Ti_int*dTidrho
  omte=-1./Te_int*dTedrho
  omne=-1./ne_int*dnedrho
  omni=omne
  

  time=1.0
  alpha_nstx=sqrt(Phi_Sep/Pi/abs(Btor))
  shot=12345
  c11=alpha_nstx/drtotdrho
  rstart=rtot_int

!100 STOP 'STOP:: No profiles.inp found'
END SUBROUTINE data_nstx





END MODULE field_nstx_mod


