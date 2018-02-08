MODULE field_d3d_mod
  USE type_mod
  USE specifications_mod, ONLY : rmin,rmax,zmin,zmax,rmag,Btor,ishot,time,dist,npol,code,Rref_gs2,drtotdrho
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: field_d3d,d3d_iterdb
  
  REAL(DP), PARAMETER :: per_phi=2*Pi
  INTEGER :: nwEQD,nhEQD !Grid dimensions treated as global parameters

!================================================================================!
CONTAINS           !                MODULE SUBPROGRAMS                           !
!================================================================================!

  SUBROUTINE field_d3d(rr,ppp,zz,Brad,Bphi,Bzet,dBraddr,dBraddp,dBraddz,&
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

    CALL stretch_coords(rrr,zzz,rm,zm)

!Equilibrium file

    CALL field_eq(rm,ploc,zm,Brad,Bphi,Bzet,dBraddr,dBraddp,dBraddz,&
                 &dBphidr,dBphidp,dBphidz,dBzetdr,dBzetdp,dBzetdz) 

!SI units
   
    Brad=Brad*Bscale              ; Bphi=Bphi*Bscale       ; Bzet=Bzet*Bscale 
    dBraddr=dBraddr*Bscale/Rscale ; dBraddp=dBraddp*Bscale ; dBraddz=dBraddz*Bscale/Rscale
    dBphidr=dBphidr*Bscale/Rscale ; dBphidp=dBphidp*Bscale ; dBphidz=dBphidz*Bscale/Rscale
    dBzetdr=dBzetdr*Bscale/Rscale ; dBzetdp=dBzetdp*Bscale ; dBzetdz=dBzetdz*Bscale/Rscale


  END SUBROUTINE field_d3d

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
   
!Export Btor=toroidal field at Rmaj for coord system (tracer.f90)

       Btor=-btf

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

    ENDIF first_reading

!---------------------------------------------------------------------------!

    CALL spline(nwEQD,nhEQD,rad,zet,hrad,hzet,icp,splpsi,ipoint,rrr,zzz,&
               &psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2,ierr)

!B components

    Brad = -dpsidz
    Bzet = psif/rrr + dpsidr
    Bphi =-btf*rtf/rrr
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

    INTEGER :: idum,i,ii,value=0
    CHARACTER(10), DIMENSION(6) ::  cdum
    CHARACTER (len=10), PARAMETER :: digit = '0123456789'

    OPEN(UNIT=11, FILE='../mesh/mesh_d3d_equ', STATUS="OLD", ACTION="READ", ERR=1000)
    READ(11,"(6A8,3I4)")(cdum(i),i=1,6),idum,nw,nh
    PRINT *,cdum(1),cdum(4)

    !Transform string to integer
    DO i=1,LEN(cdum(4))
       ii=INDEX(digit,cdum(4)(i:i))
       IF (ii == 0) CYCLE
       value=value*10+(ii-1)
    ENDDO

    !Check consistency with iterdb
    IF (value /= ishot) STOP '!STOP (dimeq): iterdb not consistent with EFIT file!'
    CLOSE(11)
  
    RETURN
    
1000 STOP 'STOP (dimeq): !DIII-D field file missing!'

  END SUBROUTINE dimeq

!============================ READ EQUILIBRIUM FIELD ============================!

  SUBROUTINE eqfile(psiSep,bt0,rzero,rad,zet,psiRZ)

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

    OPEN(UNIT=11, FILE='../mesh/mesh_d3d_equ', STATUS="OLD", ACTION="READ")
     
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

!Export Rmag

    Rmag=rmaxis

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

  SUBROUTINE stretch_coords(r,z,rm,zm)

    REAL(DP), INTENT(IN) :: r,z
    REAL(DP), INTENT(OUT) :: rm,zm

    INTEGER :: nrz=0,i,j
    INTEGER, PARAMETER :: nrzmx=100, nrhotht=360 
    REAL(DP), DIMENSION(100) :: rad_w=0.,zet_w=0.
    REAL(DP) :: R0,htht,a,b,Rw,Zw,R1,Z1
    REAL(DP) :: rho, tht, rho_c
    REAL(DP), DIMENSION(:), ALLOCATABLE :: rho_w,tht_w 
    REAL(DP), DIMENSION(nrhotht) :: rho_wall, tht_wall
    REAL(DP), PARAMETER :: delta=3.
    LOGICAL :: first=.TRUE.
    INTRINSIC MINVAL, MAXVAL
    SAVE

    first_reading: IF(first) THEN
       
       first=.FALSE.
 
       OPEN(UNIT=1, FILE='../mesh/convex_d3d', ACTION="READ", STATUS="OLD")
       
       DO i=1,nrzmx
          READ(1,*,END=10) rad_w(i),zet_w(i)
          nrz = nrz + 1
       ENDDO
10     CONTINUE

       CLOSE(1)

       nrz = nrz+1
       rad_w(nrz) = rad_w(1)
       zet_w(nrz) = zet_w(1)
   
       ALLOCATE(rho_w(nrz),tht_w(nrz))

       R0 = (MAXVAL(rad_w(1:nrz)) + MINVAL(rad_w(1:nrz)))*0.5

       DO i=1,nrz
          rho_w(i) = SQRT((rad_w(i)-R0)**2 + zet_w(i)**2)
          tht_w(i) = ATAN2(zet_w(i),(rad_w(i)-R0))
          IF (tht_w(i) < 0.) tht_w(i) = tht_w(i) + 2.*Pi
       ENDDO

       htht = 2.*pi/(nrhotht-1)

       DO i=2,nrhotht
          tht_wall(i) = htht*(i-1)
          DO j=1,nrz-1
             IF (tht_wall(i) >= tht_w(j) .AND. tht_wall(i) <= tht_w(j+1)) THEN
                IF (ABS((rad_w(j+1) - rad_w(j))/rad_w(j)) > 1.e-3) THEN
                   a = (zet_w(j+1) - zet_w(j))/(rad_w(j+1) - rad_w(j))
                   b = zet_w(j) - a*(rad_w(j) - R0)
                   Rw = b/(TAN(tht_wall(i)) - a) + R0
                   Zw = a*(Rw - R0) + b
                ELSE
                   a = (rad_w(j+1) - rad_w(j))/(zet_w(j+1) - zet_w(j))
                   b = rad_w(j)-R0 - a*zet_w(j)
                   Zw = b/(1./TAN(tht_wall(i)) - a)
                   Rw = a*Zw + b + R0
                ENDIF
             ENDIF
          ENDDO
          rho_wall(i) = SQRT((Rw-R0)**2 + Zw**2)
       ENDDO

       tht_wall(1) = 0.
       rho_wall(1) = rho_wall(nrhotht)

    END IF first_reading

    rm = r
    zm = z
    rho = SQRT((r-R0)**2 + z**2)
    tht = ATAN2(z,(r-R0))

    IF (tht < 0.) tht = tht + 2.*pi
 
    i = INT(tht/htht) + 1
    rho_c = (rho_wall(i+1) - rho_wall(i))/(tht_wall(i+1) - tht_wall(i))   &
         *(tht - tht_wall(i)) + rho_wall(i)
    
    IF(rho >= rho_c) THEN
       rho = rho_c + delta*ATAN2((rho-rho_c), delta)/pi*2.
       rm = rho*COS(tht) + R0
       zm = rho*SIN(tht)
    END IF

 END SUBROUTINE stretch_coords

!=========================== READ ITERDB =======================================!


 SUBROUTINE d3d_iterdb(rho_norm,rstart,temp1,omti,dens1,omni,temp2,&
      omte,dens2, omne,safety,c11,alpha_d3d,shat)

   USE lagrange_mod, ONLY : plag1d, indef

   REAL(DP), INTENT(IN) :: rho_norm !r/a
   REAL(DP), INTENT(OUT) :: rstart,temp1,temp2,omti,omni,omte,omne,safety,c11,alpha_d3d
   REAL(DP), INTENT(OUT) :: dens1, dens2   

   INTEGER, PARAMETER :: Mp=4 !Interpolation (do !not! change)
   INTEGER :: nf,resol
   INTEGER, DIMENSION(Mp) :: indrho
   INTEGER :: i,j,k, lun=11, ierr=0
   REAL(DP) :: beta,b0
   REAL(DP), DIMENSION(:), ALLOCATABLE :: rho,q,ti,te,ne,ni,rmajor,rminor,rtot
   REAL(DP), DIMENSION(Mp) :: rhosten,te_sten,ti_sten,ne_sten,ni_sten,q_sten
   REAL(DP), DIMENSION(Mp) :: rtot_sten,rmin_sten
   REAL(DP) :: rmin,rmax,drho,overdrho,polylag,poly1r,rho_min
   REAL(DP) :: te_int,ti_int,ne_int,ni_int,dtidrho,dtedrho,dnedrho,dnidrho,dqdrho,shat
   REAL(DP) :: rtot_int,rmin_int
   CHARACTER(len=80) :: line=''
   LOGICAL :: header=.TRUE.

!Read data from iterdb link to iterdb.* file
   
   OPEN(UNIT=lun, FILE='../input/D3D/iterdb', ACTION='READ', STATUS='OLD', ERR=100)  

!  Reading header information
   DO WHILE (INDEX(line,"ishot").EQ.0)
      READ(lun, "(A)",END=200,ERR=200) line
   END DO
   READ(lun,"(I12)",END=200,ERR=200) ishot   
   REWIND(lun)

   DO WHILE (INDEX(line,'time :').EQ.0)
      READ(lun, "(A)", END=200, ERR=200) line
   END DO
   READ(lun,"(E16.9)",END=200,ERR=200) time    
   REWIND(lun)

   DO WHILE (INDEX(line,'nj').EQ.0)
      READ(lun, "(A)",END=200,ERR=200) line
   END DO
   READ(lun,"(I12)",END=200,ERR=200) nf

  DO WHILE (INDEX(line,'Btor :').EQ.0)
      READ(lun, "(A)",END=200,ERR=200) line
   END DO
   READ(lun,"(E16.9)",END=200,ERR=200) b0
   b0=-b0

   DO WHILE (INDEX(line,'beta : toroidal beta').EQ.0)
      READ(lun, "(A)",END=200,ERR=200) line
   END DO
   READ(lun,"(E16.9)",END=200,ERR=200) beta

!  continue with fields
   ALLOCATE(rho(nf),te(nf),ti(nf),ne(nf),ni(nf),q(nf),rmajor(nf),rminor(nf),rtot(nf))
   
   DO WHILE (INDEX(line,'rho grid, meters').EQ.0)
      READ(lun, "(A)", END=200, ERR=200) line
   END DO
   READ(lun,*, END=200, ERR=200) rho

   DO WHILE (INDEX(line,'electron temperature, keV').EQ.0)
      READ(lun, "(A)", END=200, ERR=200) line
   END DO
   READ(lun,*, END=200, ERR=200) te

   DO WHILE ((ierr.EQ.0).AND.(INDEX(line,'ion temperatue, keV').EQ.0))
      READ(lun, "(A)", IOSTAT=ierr) line
   END DO
   IF (ierr.EQ.0) THEN
      READ(lun,*,END=200,ERR=200) ti
   ELSE
      REWIND(lun)
      DO WHILE (INDEX(line,'ion temperature, keV').EQ.0)
         READ(lun, "(A)",END=200,ERR=200) line
      END DO
      READ(lun,*,END=200,ERR=200) ti
   ENDIF

   DO WHILE (INDEX(line,'safety factor').EQ.0)
      READ(lun, "(A)",END=200,ERR=200) line
   END DO
   READ(lun,*, END=200, ERR=200) q

   DO WHILE (INDEX(line,'electron density').EQ.0)
      READ(lun, "(A)",END=200,ERR=200) line
   END DO
   READ(lun,*,END=200, ERR=200) ne

   DO WHILE (INDEX(line,'primary ion density').EQ.0)
      READ(lun, "(A)",END=200,ERR=200) line
   END DO
   READ(lun,*,END=200,ERR=200) ni

   DO WHILE (INDEX(line,'average major radius').EQ.0)
      READ(lun, "(A)",END=200,ERR=200) line
   END DO
   READ(lun,*,END=200,ERR=200) rmajor

   DO WHILE (INDEX(line,'average minor radius').EQ.0)
      READ(lun, "(A)",END=200,ERR=200) line
   END DO
   READ(lun,*,END=200,ERR=200) rminor

   CLOSE(lun)

!Print header

   IF (header) WRITE(*,"(A,I7)") 'Reading iterdb...Shot #',ishot
   header=.FALSE.

!Normalize rho

   alpha_d3d=MAXVAL(rho)
   rho=rho/alpha_d3d
  
!Total radius

   rtot=rmajor+rminor

!-------------------------------- INTERPOLATION ------------------------------------!
!We calculate the derivatives wrt to rho_norm

!rho_norm grid (equidistant)

   rho_min=MINVAL(rho)
   drho=rho(2)-rho(1)
   overdrho=1./drho

!Stencil

   CALL indef(rho_norm,rho_min,overdrho,nf,indrho)

   DO i=1,Mp
      rhosten(i)=rho_min+(indrho(i)-1)*drho
   ENDDO

   DO j=1,Mp
      rtot_sten(j)=rtot(indrho(j))
      Te_sten(j)=Te(indrho(j))
      Ti_sten(j)=Ti(indrho(j))
      ne_sten(j)=ne(indrho(j))
      ni_sten(j)=ni(indrho(j))
      q_sten(j)=q(indrho(j))
   ENDDO

   CALL plag1d(rho_norm,rtot_sten,overdrho,rhosten,polylag,poly1r)
   rtot_int=polylag
   drtotdrho=poly1r

   CALL plag1d(rho_norm,rmin_sten,overdrho,rhosten,polylag,poly1r)
   rmin_int=polylag

   CALL plag1d(rho_norm,te_sten,overdrho,rhosten,polylag,poly1r)
   te_int=polylag
   dTedrho=poly1r

   CALL plag1d(rho_norm,ti_sten,overdrho,rhosten,polylag,poly1r)
   ti_int=polylag
   dTidrho=poly1r

   CALL plag1d(rho_norm,ne_sten,overdrho,rhosten,polylag,poly1r)
   ne_int=polylag
   dnedrho=poly1r

   CALL plag1d(rho_norm,ni_sten,overdrho,rhosten,polylag,poly1r)
   ni_int=polylag
   dnidrho=poly1r

   CALL plag1d(rho_norm,q_sten,overdrho,rhosten,polylag,poly1r)
   safety=polylag
   dqdrho=poly1r

   !Quantities 
   c11=alpha_d3d/drtotdrho
!   IF (code == 4)   c11=1./drtotdrho

   rstart=rtot_int 
   temp1=ti_int ; temp2=te_int 
   dens1=ni_int ; dens2=ne_int

   !Lengths
   omti=-1./Ti_int*dTidrho
   omte=-1./Te_int*dTedrho
   omne=-1./ne_int*dnedrho
   omni=-1./ni_int*dnidrho

   !shat
   shat=rho_norm/safety*dqdrho

   !Length normalization for GS2
    Rref_gs2=rminor(nf)

   RETURN
100 STOP 'STOP (d3d_iterdb): !No iterdb link found in input directory!'
200 STOP 'STOP (d3d_iterdb): !Error while reading iterdb!'

 END SUBROUTINE d3d_iterdb


!-------------------------------------------------------------------------!

END MODULE field_d3d_mod


