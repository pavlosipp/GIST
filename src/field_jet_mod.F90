MODULE field_jet_mod
  USE type_mod
  USE specifications_mod, ONLY : rmin,rmax,zmin,zmax,rmag,Btor,ishot,time,dist,npol,code,Rref_gs2,drtotdrho
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: field_jet,data_jet
  
  REAL(DP), PARAMETER :: per_phi=2*Pi
  INTEGER :: np=65 !Grid dimensions treated as global parameters
  REAL(DP), DIMENSION(:,:), ALLOCATABLE:: psi
  REAL(DP), DIMENSION(:), ALLOCATABLE:: rad,zet
  REAL(DP):: hrad,hzet,btf,rtf
  SAVE

!================================================================================!
CONTAINS           !                MODULE SUBPROGRAMS                           !
!================================================================================!

  SUBROUTINE field_jet(rr,ppp,zz,Brad,Bphi,Bzet,dBraddr,dBraddp,dBraddz,&
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

!Equilibrium file

    CALL field_eq(rrr,ploc,zzz,Brad,Bphi,Bzet,dBraddr,dBraddp,dBraddz,&
                 &dBphidr,dBphidp,dBphidz,dBzetdr,dBzetdp,dBzetdz) 

!SI units
   
    Brad=Brad*Bscale              ; Bphi=Bphi*Bscale       ; Bzet=Bzet*Bscale 
    dBraddr=dBraddr*Bscale/Rscale ; dBraddp=dBraddp*Bscale ; dBraddz=dBraddz*Bscale/Rscale
    dBphidr=dBphidr*Bscale/Rscale ; dBphidp=dBphidp*Bscale ; dBphidz=dBphidz*Bscale/Rscale
    dBzetdr=dBzetdr*Bscale/Rscale ; dBzetdp=dBzetdp*Bscale ; dBzetdz=dBzetdz*Bscale/Rscale


  END SUBROUTINE field_jet

!============================== READ EFIT FILE ==================================!

  SUBROUTINE field_eq(rrr,ppp,zzz,Brad,Bphi,Bzet,dBraddr,dBraddp,dBraddz,&             
                     &dBphidr,dBphidp,dBphidz,dBzetdr,dBzetdp,dBzetdz)


    REAL(DP), INTENT(IN) :: rrr,ppp,zzz
    REAL(DP), INTENT(OUT) :: Brad,Bphi,Bzet,dBraddr,dBraddp,dBraddz
    REAL(DP), INTENT(OUT) :: dBphidr,dBphidp,dBphidz,dBzetdr,dBzetdp,dBzetdz

    REAL(DP), DIMENSION(:), ALLOCATABLE :: xi,f
!    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: psi
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: splpsi
    REAL(DP) :: dummy,psib
    INTEGER, DIMENSION(:), ALLOCATABLE :: imi,ima,jmi,jma
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ipoint
    INTEGER :: icp=0,i,j,ierr
    REAL(DP) :: d2psidz2,d2psidrdz,dpsidz,dpsidr,psif,d2psidr2
    LOGICAL :: first=.TRUE.
    SAVE
  
    first_reading:  IF (first) THEN
  
       first=.FALSE.

     
!Export Btor=toroidal field at Rmaj for coord system (tracer.f90)

       Btor=-btf

!CGS units     

       rad = rad*100. 
       zet = zet*100. 
       rtf = rtf*100. 
       psi = psi*1.e8
!       psib= psib*1.e8
       btf = btf*1.e4

! Physical $\phi$ component of vector potential:

       DO j=1,np
          DO i=1,np
             psi(i,j) = (psi(i,j))/rad(i) !+psib
          ENDDO
       ENDDO

OPEN(44,file='psiout',status='replace')
DO i=1,np
   DO j=1,np
      write(44,*) rad(i),zet(j),psi(i,j)
   ENDDO
ENDDO
CLOSE(44)
!Resolution

      hrad = rad(2) - rad(1)
      hzet = zet(2) - zet(1)

!------------------------------ Splines -----------------------------------!

!Rectangular domain

       ALLOCATE(imi(np),ima(np),jmi(np),jma(np))    

       imi = 1
       ima = np
       jmi = 1
       jma = np

!Number of data in splpsi
  
       DO i=1,np
          IF ( imi(i) > 0 .AND. ima(i) > 0 ) THEN
             icp = icp + ima(i) - imi(i) + 1
          ENDIF
       ENDDO
   
       ALLOCATE(splpsi(6,6,icp),ipoint(np,np))
       CALL s2dcut(np,np,hrad,hzet,psi,imi,ima,jmi,jma,icp,splpsi,ipoint)

    ENDIF first_reading

!---------------------------------------------------------------------------!

    CALL spline(np,np,rad,zet,hrad,hzet,icp,splpsi,ipoint,rrr,zzz,&
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





!=========================== READ JET DATA =====================================!


 SUBROUTINE data_jet(rrr,line,rstart,c11jet,alpha_jet,r_alt,safety_db,shat_db)
   USE specifications_mod, ONLY: num_lines,shot,TSHOT,rmag,zmag,Phi_Sep,rmin,rmax,diag
   USE lagrange_mod, ONLY : plag1d_uneven, indef_uneven

   REAL(DP), INTENT(IN) :: rrr
   INTEGER, INTENT(IN):: line 
   REAL(DP), INTENT(OUT) :: rstart,safety_db,c11jet,alpha_jet,shat_db
   
   
   !MDS+
   INTEGER, PARAMETER:: IDTYPE_LONG=8
   INTEGER, PARAMETER:: IDTYPE_FLOAT=10,IDTYPE_DOUBLE=11
   
   INTEGER:: descr, MdsValue, MdsConnect, MdsOpen, MdsClose, ierr, icp=0
   INTEGER:: MdsDisconnect
   INTEGER:: descr2, MdsValue2
   INTEGER:: IDTYPE_UNSIGNED_LONG, IDTYPE_LONG
   INTEGER:: IDTYPE_FLOAT, IDTYPE_DOUBLE, IDTYPE_COMPLEX
   INTEGER:: status,length,i,nt,j,jt,nrho,nangle
   
   !FLUSH
   REAL:: flushtmp(1000)
   REAL,DIMENSION(:), ALLOCATABLE:: r, z,rbndyallt,zbndallt
   REAL,DIMENSION(:,:), ALLOCATABLE:: psinorm, brg,bzg
   CHARACTER(LEN=200*100):: flush_command
   CHARACTER(LEN=20):: timeparameter
   CHARACTER(LEN=10):: np_str,np_str
   CHARACTER(20*100):: r_str, z_str, tempstr

   REAL:: dr,dz
   REAL, DIMENSION(:), ALLOCATABLE:: timear,Bvac,zax,rax,phit,rmjo,rmji,qrho
   REAL, DIMENSION(200):: psin_efit

   INTEGER, PARAMETER :: Mp=4 !Interpolation (do !not! change)
   INTEGER :: resol
   INTEGER, DIMENSION(Mp) :: indrho
   INTEGER :: i,j,k, ierr=0
   REAL(DP), DIMENSION(Mp) :: rmajisten,rmajosten,rhotorsten,psisten,rsten,q_sten
   REAL(DP), DIMENSION(Mp) :: rtot_sten,rmin_sten
   REAL(DP) :: polylag,poly1r,alpha_jet
   REAL(DP)::fbnd,rbnd,faxs,r_inb,r_alt
   REAL(DP) :: polylag,poly1r,dqdrho, fbnd, fcentre
   REAL(DP) :: rtot_int,rmin_int
   REAL(DP), DIMENSION(:), ALLOCATABLE::time,rmajo,rmaji,rhotor,phitor,q_rhotor,rbndy,zbnd
   CHARACTER(6)::shotstr
   LOGICAL:: first=.TRUE.
   LOGICAL,DIMENSION(:),ALLOCATABLE:: linefirst
SAVE

first_read: IF (first) THEN

   ALLOCATE(linefirst(num_lines))
   linefirst=.TRUE.
   first=.FALSE.

!|--------------------Get field data from MDS+ server---------------------------|

   !Convert shotnumber to string
   WRITE(shotstr,*) shot
IF (shot.lt.63446) STOP 'Toroidal flux signal is not available for shot no. < 63446.'
IF (shot.gt.63446.AND.shot.lt.66080) print *,'Warning: Toroidal flux signal may be missing for this shot.'

   ALLOCATE(timear(4000))
   !Connect to MDS+ server
   status=MdsConnect('mdsplus.jet.efda.org'//CHAR(0))
!   IF (mod(status,2)==0) STOP 'error connecting to MDS+ server'





   !Read EFIT timebase
   status=MdsValue('dim_of(_S=jet("ppf/'//diag//'/sspr",'//shotstr//'),1)'//CHAR(0), descr(IDTYPE_FLOAT,timear,4000,0),0,length)
   nt=length

   !Select element from EFIT-timear which is closest to input tSHOT   
   jt=1   
   DO i=1,nt
      IF (ABS((timear(i)-tSHOT)).LE.ABS((timear(jt)-tSHOT))) jt=i
   ENDDO
   PRINT *,'selected time ('//TRIM(diag)//')=',timear(jt)



   ALLOCATE(r(np),z(np),Bvac(nt),rax(nt),zax(nt),rad(np),zet(np))



!Get r and z of separatrix for ~100 angular pos. to determine size of plasma
ALLOCATE(rbndyallt(200000),zbndallt(200000))
   status=MdsValue('_S=jet("ppf/'//diag//'/rbnd",'//shotstr//')'//CHAR(0),descr(IDTYPE_FLOAT,rbndyallt,200000,0),0,length)
nangle=length/nt
   status=MdsValue('_S=jet("ppf/'//diag//'/zbnd",'//shotstr//')'//CHAR(0),descr(IDTYPE_FLOAT,zbndallt,200000,0),0,length)
   

ALLOCATE(rbndy(nangle),zbnd(nangle))
   DO i=1,nangle
      rbndy(i)=DBLE(rbndyallt((jt-1)*nangle+i))
      zbnd(i)=DBLE(zbndallt((jt-1)*nangle+i))
   ENDDO

   !Boundaries
   rmin=MINVAL(rbndy)-0.01 !1.65 in EFIT file
   rmax=MAXVAL(rbndy)+0.01 !4.05
   zmin=MINVAL(zbnd)-0.01 !-1.9
   zmax=MAXVAL(zbnd)+0.01 !1.9




   !Prepare a set of R and z values
   dr=REAL(rmax-rmin)/REAL(np-1)
   dz=REAL(zmax-zmin)/REAL(np-1)
   DO i=1,np
      r(i)=REAL(rmin)+REAL(i-1)*dr
   ENDDO
   DO i=1,np
      z(i)=REAL(zmin)+REAL(i-1)*dz
   ENDDO
   rad=dble(r)
   zet=dble(z)
   hrad=dble(dr)
   hzet=dble(dz)

! Build strings containing np, r and z	
	write(np_str,*)np
	write(r_str(1:20),'(A1,1f10.2)')'[',100*r(1)
	write(z_str(1:20),'(A1,1f10.2)')'[',100*z(1)
	
	do i=2,np
		write(tempstr,'(A1,2f10.2)')',',100*r(i)
		r_str(1+(i-1)*20:20*i)=tempstr

		write(tempstr,'(A1,2f10.2)')',',100*z(i)
		z_str(1+(i-1)*20:20*i)=tempstr
	enddo

	r_str(20*np:)=']'
	z_str(20*np:)=']'



!------------------ Read poloidal flux psi via FLUSH routines ------------------|

! Init flush
	write(timeparameter,*)tSHOT
	timeparameter='_tSHOT='//timeparameter
	status=mdsvalue(timeparameter//CHAR(0),descr(IDTYPE_FLOAT,flushtmp,1000,0),0,length)


	flush_command='flushinit(15,'//shotstr//',_tSHOT,0,"JETPPF","'//trim(diag)//'",_ierr)'

	status=mdsvalue(flush_command//CHAR(0),descr(IDTYPE_FLOAT,flushtmp,1000,0),0,length)
	IF (status.ne.1) STOP 'Error in flushinit'


! Read normalized poloidal flux at these positions
	ALLOCATE(psinorm(np,np), brg(np,np),bzg(np,np))



	flush_command='flupng('//np_str//','//np_str//','//r_str//','//z_str//',_psinorm,_brg,_bzg,_ierr)'
	!flush_command='flupng(2,2,[300,382.5],[0,25],_psinorm,_brg,_bzg,_ierr)'
	
	
	status=mdsvalue(flush_command//CHAR(0),descr(IDTYPE_FLOAT,flushtmp,1000,0),0,length)
	IF (status.ne.1) STOP 'Error in flupn'


! Fill psinorm array
	status=mdsvalue('_psinorm'//CHAR(0),descr(IDTYPE_FLOAT,psinorm,np*np,0),0,length)

	status=mdsvalue('_brg'//CHAR(0),descr(IDTYPE_FLOAT,brg,np*np,0),0,length)
	status=mdsvalue('_bzg'//CHAR(0),descr(IDTYPE_FLOAT,bzg,np*np,0),0,length)

!------------------- END of FLUSH call ----------------------------|



!------------------- Get additional equilibrium quantities --------|

   !Read toroidal field at R=2.96m
   status=MdsValue('_S=jet("ppf/'//diag//'/BVAC",'//shotstr//')'//CHAR(0),descr(IDTYPE_FLOAT,Bvac,nt,0),0,length)

   !Read coordinates of magnetic axis
   
   status=MdsValue('_S=jet("ppf/'//diag//'/RMAG",'//shotstr//')'//CHAR(0),descr(IDTYPE_FLOAT,rax,nt,0),0,length)
   status=MdsValue('_S=jet("ppf/'//diag//'/ZMAG",'//shotstr//')'//CHAR(0),descr(IDTYPE_FLOAT,zax,nt,0),0,length)


   !Get length of rho_tor array
   
   status=MdsValue('dim_of(_S=jet("ppf/'//diag//'/q",'//shotstr//'),0)'//CHAR(0), descr(IDTYPE_FLOAT,psin_efit,200,0),0,length)

   !psin_efit now contains nrho equidistant points; phit,rmjo,rmji and q
   !are given at these points 
   nrho=length

   ALLOCATE(phit(nt*nrho),qrho(nt*nrho),q_rhotor(nrho))
   
   !Get phi_normalized=rhotor and phi_tor
   status=MdsValue('_S=jet("ppf/'//diag//'/ftor",'//shotstr//')'//CHAR(0),descr(IDTYPE_FLOAT,phit,nt*nrho,0),0,length)
   ALLOCATE(rhotor(nrho),rmjo(nrho*nt),rmji(nrho*nt),rmajo(nrho),rmaji(nrho),phitor(nrho))

   !For psin_efit-array: read corresponding (outer and inner) Rmajors
   
   status=MdsValue('_S=jet("ppf/'//diag//'/RMJO",'//shotstr//')'//CHAR(0),descr(IDTYPE_FLOAT,rmjo,nrho*nt,0),0,length)
   
   status=MdsValue('_S=jet("ppf/'//diag//'/RMJI",'//shotstr//')'//CHAR(0),descr(IDTYPE_FLOAT,rmji,nrho*nt,0),0,length)
   
   status=MdsValue('_S=jet("ppf/'//diag//'/Q",'//shotstr//')'//CHAR(0), descr(IDTYPE_FLOAT,qrho,nt*nrho,0),0,length)    


! Create psi grid
   ALLOCATE(psi(np,np))

! Read FBND and FAXS signals to un-normalize psinorm
	status=MdsValue('_S=jet("ppf/'//diag//'/fbnd",'//shotstr//')'//CHAR(0),descr(IDTYPE_FLOAT,flushtmp,1000,0),0,length)
   fbnd=DBLE(flushtmp(jt))

	status=MdsValue('_S=jet("ppf/'//diag//'/faxs",'//shotstr//')'//CHAR(0),descr(IDTYPE_FLOAT,flushtmp,1000,0),0,length)
        fcentre=DBLE(flushtmp(jt))
DO i=1,np
   DO j=1,np
      psi(i,j) = DBLE(psinorm(i,j)) *DBLE(fbnd-fcentre) + DBLE(fcentre)
   ENDDO
ENDDO

   btf=DBLE(Bvac(jt))
   rtf=2.96d0
   Btor=-btf
   rmag=DBLE(rax(jt))
   zmag=DBLE(zax(jt))
!print *,rmag, zmag

   !Setup rhotor and rmaj arrays for later conversion
   DO i=1,nrho
      phitor(i)=DBLE(phit(nrho*(jt-1)+i))
   ENDDO
   
   Phi_Sep=abs(phitor(nrho))
   DO i=1,nrho
      rhotor(i)=DBLE(sqrt(abs(phitor(i))/Phi_Sep))
      q_rhotor(i)=DBLE(qrho(nrho*(jt-1)+i))
      rmajo(i)=DBLE(rmjo((jt-1)*nrho+i))
      rmaji(i)=DBLE(rmji((jt-1)*nrho+i))
   ENDDO


   !Disconnect from MDS+ server
   status=MdsDisconnect()

ENDIF first_read


!----------- Compute Phi from safety factor integration wrt Psi --------




!-------------------------------- INTERPOLATION ------------------------------------!
!We calculate the derivatives wrt to rho_tor


!-------------Conversion of rho_tor to R only once for each line----------------!

line_first: IF (linefirst(line)) THEN    
   
   linefirst(line)=.FALSE.

   !-------------Setup conversion of rho_tor to R,z--------------------------!
   
   !Stencil
   
   CALL indef_uneven(rrr,rhotor,nrho,indrho)
   
   DO j=1,Mp
      rhotorsten(j)=rhotor(indrho(j))
      rmajosten(j)=rmajo(indrho(j))
      rmajisten(j)=rmaji(indrho(j))
      q_sten(j)=q_rhotor(indrho(j)) !stencil for safety factor
   ENDDO
   
   
   !Calculate r(outboard) for given rho_tor
   
   CALL plag1d_uneven(rrr,rmajosten,rhotorsten,polylag,poly1r)
   
   rtot_int=polylag   !this r is now used for tracing
   drtotdrho=poly1r
   print *,'drdrho=',drtotdrho,'Rmag=',rmag
   !Calculate r(inboard) for given rho_tor
   CALL plag1d_uneven(rrr,rmajisten,rhotorsten,polylag,poly1r)
   
   r_inb=polylag 
   
   !Calculate q and dqdrho for given rho_tor
   CALL plag1d_uneven(rrr,q_sten,rhotorsten,polylag,poly1r)
   safety_db=polylag
   dqdrho=poly1r
   shat_db=rrr/safety_db*dqdrho

   !print *,'rho=',rrr
!   print *,'q=',safety_db,'dqdrho=',dqdrho,'shat=',shat_db
   
   r_alt=(rstart-r_inb)/2.
!   print *,'r_alt=',r_alt

   alpha_jet=sqrt(Phi_Sep/Pi/abs(Btor))


   !Quantities 
   c11jet=alpha_jet/drtotdrho
   rstart=rtot_int 

   !Length normalization for GS2
!    Rref_gs2=rminor(nf)


ENDIF line_first


END SUBROUTINE data_jet


!-------------------------------------------------------------------------!

END MODULE field_jet_mod

