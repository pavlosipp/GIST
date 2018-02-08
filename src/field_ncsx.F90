
 
SUBROUTINE field_ncsx(rrr,ppp,zzz,Brad,Bphi,Bzet,dBraddr,dBraddp,dBraddz,&
     &dBphidr,dBphidp,dBphidz,dBzetdr,dBzetdp,dBzetdz)

  USE util_mod
  USE lagrange_mod
  USE specifications_mod, ONLY : rmin,rmax,zmin,zmax,read_mode,vmec,b0

  IMPLICIT NONE



!!===================================================================!! 
!! This routine reads the physical components of the NCSX field 
!! on a cylindrical mesh and performs 3-degree Lagrange interpolation
!!===================================================================!!

  REAL(DP), INTENT(IN) :: rrr,ppp,zzz
  REAL(DP), INTENT(OUT) :: Brad,Bphi,Bzet,dBraddr,dBraddp,dBraddz
  REAL(DP), INTENT(OUT) :: dBphidr,dBphidp,dBphidz,dBzetdr,dBzetdp,dBzetdz
  REAL(DP) :: r,z,phi,philoc,sinp,cosp,rr,pp,zz,dr,dphi,dz,overdr,overdp,overdz
  REAL(DP) :: pmin,pmax
  REAL(DP) :: polylag, poly1r, poly1p,poly1z
  REAL(DP), PARAMETER :: per_phi=2.*Pi/3. !Device dependent
  REAL(DP), DIMENSION(3,3):: c
  INTEGER, PARAMETER :: Mp=4 
  REAL(DP), DIMENSION(Mp) :: rsten,psten,zsten
  REAL(DP), DIMENSION(Mp,Mp,Mp) :: Brsten,Bpsten,Bzsten
  INTEGER, DIMENSION(Mp) :: indr,indp,indz
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Br,Bz,Bp
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Bx,By
  REAL(DP), DIMENSION(:), ALLOCATABLE :: xi,f
  REAL(DP), DIMENSION(:), ALLOCATABLE :: phi_real 
  REAL(DP) :: dum1,dum2,dum3
  REAL(DP),PARAMETER :: eps0=1.256e-6 
  INTEGER :: i,j,k,l,ll,m,ir,iz,ip,count=0
  INTEGER :: first_time=0
  INTRINSIC SIN,COS,SQRT

!---------------------------------- MESH PARAMETERS ---------------------------!
INTEGER :: nr=150, np=289, nz=151, values=6545850 
!------------------------------------------------------------------------------!

SAVE 

loading: IF (first_time == 0) THEN
   first_time=first_time+1

   ALLOCATE(bx(nr,np,nz),by(nr,np,nz))
   ALLOCATE(br(nr,0:np+2,nz),bz(nr,0:np+2,nz),bp(nr,0:np+2,nz))

!Read mesh
   read_mod: IF(read_mode == 0) THEN

     OPEN(UNIT=10, ACTION="READ", FILE='../mesh/mesh_ncsx', STATUS='OLD', ERR=1000)
     
     WRITE(*, FMT="(A)") 'Loading NCSX field file ...'

     DO ir=1,nr
        DO ip=1,np
           DO iz=1,nz
              count=count+1
              READ(10,*) dum1,rr,pp,zz,Br(ir,ip,iz),Bp(ir,ip,iz),Bz(ir,ip,iz),dum2,dum3
!Extrema

              IF (count == 1) THEN
                 rmin=rr
                 pmin=pp
                 zmin=zz
              ELSEIF (count == values) THEN
                 rmax=rr
                 pmax=pp
                 zmax=zz
              ENDIF
           ENDDO
        ENDDO
     ENDDO

     CLOSE(10)
   ELSE

     OPEN(UNIT=10, ACTION="READ", FILE='../mesh/mesh_ncsx_bin', &
          STATUS='OLD',FORM='unformatted', ERR=1000)
     
     WRITE(*, FMT="(A)") 'Loading NCSX field file (bin)...'

     READ(10) nr,np,nz

     DO ir=1,nr
        DO ip=1,np
           DO iz=1,nz
              count=count+1
              READ(10) dum1,rr,pp,zz,Br(ir,ip,iz),Bp(ir,ip,iz),Bz(ir,ip,iz),dum2,dum3

!Extrema

              IF (count == 1) THEN
                 rmin=rr
                 pmin=pp
                 zmin=zz
              ELSEIF (count == values) THEN
                 rmax=rr
                 pmax=pp
                 zmax=zz
              ENDIF
           ENDDO
        ENDDO
     ENDDO

     CLOSE(10)
   ENDIF read_mod
!Convert to Tesla

   Br=Br*eps0 ; Bp=Bp*eps0 ; Bz=Bz*eps0

!Consistency check 

   IF (count /= values) STOP 'STOP (field_ncsx): !NCSX field file not read properly!'

   dphi=(pmax-pmin)/(Np-1)
   dr=(rmax-rmin)/(Nr-1)
   dz=(zmax-zmin)/(Nz-1)

!Enforce periodicity on guard cells

Br(:,np+1,:)=Br(:,1,:)
Bp(:,np+1,:)=Bp(:,1,:)
Bz(:,np+1,:)=Bz(:,1,:)
Br(:,np+2,:)=Br(:,2,:)
Bp(:,np+2,:)=Bp(:,2,:)
Bz(:,np+2,:)=Bz(:,2,:)
Br(:,0,:)=Br(:,np,:)
Bp(:,0,:)=Bp(:,np,:)
Bz(:,0,:)=Bz(:,np,:)

ENDIF loading

!---------------------------- Interpolation --------------------------------!

!------------------ ***** Periodic constraint ***** ------------------------!

philoc=ppp  
IF (philoc >= per_phi) philoc=philoc-(INT(philoc/per_phi))*per_phi
IF (philoc < 0.) philoc = philoc + (INT(ABS(philoc)/per_phi) +1)*per_phi

overdr=1./dr
overdz=1./dz
overdp=1./dphi

!Produce 4x4x4-point stencil for each (r,p,z) 

CALL indef(rrr,rmin,overdr,Nr,indr)
CALL indef(zzz,zmin,overdz,Nz,indz)
CALL indef_per(philoc,pmin,overdp,Np,indp) !periodic grid

!Coordinates of stencil points 

DO i=1,Mp
   rsten(i)=rmin+(indr(i)-1)*dr
   psten(i)=pmin+(indp(i)-1)*dphi
   zsten(i)=zmin+(indz(i)-1)*dz
ENDDO

!Interpolate Br,Bp,Bz 

DO k=1,mp
   DO l=1,mp
      DO m=1,mp
         Brsten(m,l,k)=Br(indr(m),indp(l),indz(k))
      ENDDO
   ENDDO
ENDDO

CALL plag3d(rrr,philoc,zzz,Brsten,overdr,overdp,overdz,&
     &rsten,psten,zsten,polylag,poly1r,poly1p,poly1z)

Brad=polylag
dBraddr=poly1r
dBraddp=poly1p
dBraddz=poly1z

DO k=1,mp
   DO l=1,mp
      DO m=1,mp
         Bzsten(m,l,k)=Bz(indr(m),indp(l),indz(k))
      ENDDO
   ENDDO
ENDDO

CALL plag3d(rrr,philoc,zzz,Bzsten,overdr,overdp,overdz,&
        &rsten,psten,zsten,polylag,poly1r,poly1p,poly1z)

Bzet=polylag
dBzetdr=poly1r
dBzetdp=poly1p
dBzetdz=poly1z  

DO k=1,mp
   DO l=1,mp
      DO m=1,mp
         Bpsten(m,l,k)=Bp(indr(m),indp(l),indz(k))
      ENDDO
   ENDDO
ENDDO

CALL plag3d(rrr,philoc,zzz,Bpsten,overdr,overdp,overdz,&
           &rsten,psten,zsten,polylag,poly1r,poly1p,poly1z)

Bphi=polylag
dBphidr=poly1r
dBphidp=poly1p
dBphidz=poly1z      

RETURN
1000 STOP 'STOP(field_ncsx.f90): !NCSX field file missing!'


END SUBROUTINE field_ncsx







