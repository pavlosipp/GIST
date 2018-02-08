!!===================================================================!! 
!! This routine reads the cartesian components of the W7-X field 
!! on a cylindrical mesh and performs 3-degree Lagrange interpolation
!!===================================================================!!

SUBROUTINE field_w7x(rrr,ppp,zzz,Brad,Bphi,Bzet,dBraddr,dBraddp,dBraddz,&
     &dBphidr,dBphidp,dBphidz,dBzetdr,dBzetdp,dBzetdz)

  USE util_mod
  USE lagrange_mod
  USE specifications_mod, ONLY : rmin,rmax,zmin,zmax,read_mode
  USE UtilCDF
  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: rrr,ppp,zzz
  REAL(DP), INTENT(OUT) :: Brad,Bphi,Bzet,dBraddr,dBraddp,dBraddz
  REAL(DP), INTENT(OUT) :: dBphidr,dBphidp,dBphidz,dBzetdr,dBzetdp,dBzetdz

  REAL(DP) :: r,z,phi,philoc,sinp,cosp
  REAL(DP) :: rr,pp,zz,dr,dphi,dz,overdr,overdp,overdp_sp,overdz
  REAL(DP) :: pmin,pmax,brf,zloc
  REAL(DP) :: polylag, poly1r, poly1p,poly1z
  REAL(DP), PARAMETER :: per_phi=2.*Pi/5. !device dependent
  REAL(DP), DIMENSION(3,3):: c
  INTEGER, PARAMETER :: Mp=4 !for 3rd degree interpolation
  REAL(DP), DIMENSION(Mp) :: rsten,psten,zsten
  REAL(DP), DIMENSION(Mp,Mp,Mp) :: Brsten,Bpsten,Bzsten,slope
  REAL(DP), DIMENSION(Mp,Mp,Mp) :: aux1,aux2,aux3,aux4,aux5
  INTEGER, DIMENSION(Mp) :: indr,indp,indz
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: bxt,byt,bzt
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: br,bz,bp
  REAL(DP), DIMENSION(:), ALLOCATABLE :: xi,f
  REAL(DP), DIMENSION(:), ALLOCATABLE :: phi_real 
  INTEGER :: i,j,j1,j2,j2_ind,j3,k,l,ll,m,ir,iz,ip,count
  INTEGER :: first_time=0
  INTRINSIC SIN,COS,SQRT

  INTEGER :: ndim,nvar,nread,ncid,rid,pid,zid
  INTEGER :: icount(3),istart(3)
  CHARACTER(LEN=30) :: c_dummy

!-------------------------------------------------------------------------!
  INTEGER :: nr, np, nz, values
!-------------------------------------------------------------------------!

  SAVE 

  loading: IF (first_time == 0) THEN
     first_time=first_time+1

     count=0

!Read mesh

     read_mod: IF(read_mode==0) THEN
       OPEN(UNIT=10, ACTION="READ", FILE='../mesh/mesh_w7x',&
          &STATUS='OLD', ERR=1000)
       ! PRINT *, 'Loading W7X field data...'

       READ(10,*)
       READ(10,*)
!       READ(10,"(A15,3(2x,I3))") c_dummy, nr, np, nz
       READ(10,"(A11,3(1x,I3))") c_dummy, nr, np, nz 
       READ(10,*)
       READ(10,*)
       READ(10,*) ndim

       ALLOCATE(bxt(nr,np,nz),byt(nr,np,nz),bzt(nr,np,nz))

!Guard cells in the phi direction: 0,np+1,np+2
!Later we enforce periodicity. This has the benefit that in the interpolation
!the point is always centered

       ALLOCATE(br(nr,0:np+2,nz),bz(nr,0:np+2,nz),bp(nr,0:np+2,nz))
      
       READ(10,*)
       READ(10,*)
       READ(10,*) nvar
      
       READ(10,*)
       READ(10,*)
       DO i=1,ndim
         READ(10,*)
       ENDDO
       DO i=1,nvar
         READ(10,*)
       ENDDO
      
       READ(10,*)
       READ(10,*)
       READ(10,*) values
       READ(10,*)
       READ(10,*)

       DO ir=1,nr
          DO ip=1,np
             DO iz=1,nz
                count=count+1
                READ(10,*) rr,pp,zz,bxt(ir,ip,iz),byt(ir,ip,iz),bzt(ir,ip,iz)
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
     ELSEIF(read_mode==1) THEN
        OPEN(UNIT=10, ACTION="READ", FILE='../mesh/mesh_w7x_unf', STATUS='OLD',&
             FORM='UNFORMATTED', ERR=1000)
       ! PRINT *, 'Loading W7X field data (binary)...'

       READ(10) nr,np,nz

       ALLOCATE(bxt(nr,np,nz),byt(nr,np,nz),bzt(nr,np,nz))

!Guard cells in the phi direction: 0,np+1,np+2
!Later we enforce periodicity. This has the benefit that in the interpolation
!the point is always centered

       ALLOCATE(br(nr,0:np+2,nz),bz(nr,0:np+2,nz),bp(nr,0:np+2,nz))

       DO ir=1,nr
          DO ip=1,np
             DO iz=1,nz
                count=count+1
                READ(10) rr,pp,zz,bxt(ir,ip,iz),byt(ir,ip,iz),bzt(ir,ip,iz)

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

     ELSEIF(read_mode==2) THEN
#ifdef NETCDF
        ! PRINT *, 'Loading W7X field data (cdf)...'
        CALL open_cdf("../mesh/mesh_w7x_cdf",READ_ACC_CDF,ncid)
        CALL get_dim_cdf(ncid,'r',rid,nr)
        CALL get_dim_cdf(ncid,'phi',pid,np)
        CALL get_dim_cdf(ncid,'z',zid,nz)
        count=nr*np*nz
        values=count
        ALLOCATE(bxt(nr,np,nz))
        ALLOCATE(byt(nr,np,nz))
        ALLOCATE(bzt(nr,np,nz))
        ALLOCATE(br(nr,0:np+2,nz))
        ALLOCATE(bp(nr,0:np+2,nz))
        ALLOCATE(bz(nr,0:np+2,nz))
        istart(1)=1;istart(2)=1;istart(3)=1
        icount(1)=1;icount(2)=1;icount(3)=1
        CALL get_var_cdf(ncid,'r',rmin,istart,icount)
        CALL get_var_cdf(ncid,'phi',pmin,istart,icount)
        CALL get_var_cdf(ncid,'z',zmin,istart,icount)
        ! IF(g_rank==0) PRINT *,'cdf rmin ',rmin
        ! IF(g_rank==0) PRINT *,'cdf pmin ',pmin
        ! IF(g_rank==0) PRINT *,'cdf zmin ',zmin
        istart(1)=nr;istart(2)=np;istart(3)=nz
        CALL get_var_cdf(ncid,'r',rmax,istart,icount)
        CALL get_var_cdf(ncid,'phi',pmax,istart,icount)
        CALL get_var_cdf(ncid,'z',zmax,istart,icount)
        ! IF(g_rank==0) PRINT *,'cdf rmax ',rmax
        ! IF(g_rank==0) PRINT *,'cdf pmax ',pmax
        ! IF(g_rank==0) PRINT *,'cdf zmax ',zmax
        istart(1)=1;istart(2)=1;istart(3)=1
        icount(1)=nr;icount(2)=np;icount(3)=nz
        CALL get_mdvar_cdf(ncid,'Bx',bxt,istart,icount)
        CALL get_mdvar_cdf(ncid,'By',byt,istart,icount)
        CALL get_mdvar_cdf(ncid,'Bz',bzt,istart,icount)
        CALL close_cdf(ncid)
#else
        STOP 'STOP (field_w7x.F90): &
             &The code was compiled without netCDF support.'
#endif
     ELSE
        STOP 'STOP (field_w7x.F90): !Unknown read mode for W7-X.'
     ENDIF read_mod

!Consistency check 

     IF (count /= values) STOP 'STOP (field_w7x.f90): !W7-X file not read properly!'

     dphi=(pmax-pmin)/(np-1)
     dr=(rmax-rmin)/(nr-1)
     dz=(zmax-zmin)/(nz-1)

!Transform cartesian components to cylindrical

     DO ip=1,np
        phi=pmin+(ip-1)*dphi
        cosp=COS(phi)
        sinp=SIN(phi)
        br(:,ip,:)=bxt(:,ip,:)*cosp + byt(:,ip,:)*sinp
        bp(:,ip,:)=-bxt(:,ip,:)*sinp + byt(:,ip,:)*cosp
        bz(:,ip,:)=bzt(:,ip,:)
     ENDDO

!Enforce periodicity on guard cells

br(:,np+1,:)=br(:,1,:)
bp(:,np+1,:)=bp(:,1,:)
bz(:,np+1,:)=bz(:,1,:)
br(:,np+2,:)=br(:,2,:)
bp(:,np+2,:)=bp(:,2,:)
bz(:,np+2,:)=bz(:,2,:)
br(:,0,:)=br(:,np,:)
bp(:,0,:)=bp(:,np,:)
bz(:,0,:)=bz(:,np,:)

  ENDIF loading


!--------------------------------- Interpolation -------------------------------!

!------------------------ ***** Periodic constraint ***** ----------------------!

  philoc=ppp 
  IF (philoc >= per_phi) philoc=philoc-(INT(philoc/per_phi))*per_phi
  IF (philoc < 0.) philoc = philoc + (INT(ABS(philoc)/per_phi) +1)*per_phi

!----------------------- W7X antisymmetry -------------------

  brf=1.
  zloc=zzz  
  IF (philoc >= per_phi/2) THEN
    philoc = per_phi - philoc
    zloc = -zzz
    brf=-1.0
  ENDIF

!--------------------------------------------------------
 
  overdr=1./dr
  overdz=1./dz
  overdp=1./dphi

!Produce 4x4x4-point stencil for each (r,p,z)

  CALL indef(rrr,rmin,overdr,nr,indr)
  CALL indef(zloc,zmin,overdz,nz,indz)
  CALL indef_per(philoc,pmin,overdp,np,indp) !periodic stencil

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

CALL plag3d(rrr,philoc,zloc,Brsten,overdr,overdp,overdz,&
            rsten,psten,zsten,polylag,poly1r,poly1p,poly1z)
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

CALL plag3d(rrr,philoc,zloc,Bzsten,overdr,overdp,overdz,&
     rsten,psten,zsten,polylag,poly1r,poly1p,poly1z)

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

CALL plag3d(rrr,philoc,zloc,Bpsten,overdr,overdp,overdz,&
     rsten,psten,zsten,polylag,poly1r,poly1p,poly1z)

Bphi=polylag
dBphidr=poly1r
dBphidp=poly1p
dBphidz=poly1z

!---------------------------- W7X symmetry --------------------!

dBraddr = brf*dBraddr
dBphidz = brf*dBphidz
dBzetdz = brf*dBzetdz
dBphidp = brf*dBphidp
dBzetdp = brf*dBzetdp
Brad = brf*Brad

!---------------------------------------------------------------!






!------------- 2D Lagrange + 1D Hermite interpolation for divB=0 ---------------------

!Compute slopes on the stencil points
!!$slope=0.
!!$do j1=1,Mp
!!$   do j2=2,3
!!$j2_ind=j2-1
!!$      do j3=1,Mp
!!$CALL plag3d(rsten(j1),psten(j2),zsten(j3),Brsten,overdr,overdp,overdz,&
!!$     rsten,psten,zsten,polylag,poly1r,poly1p,poly1z)
!!$aux1(j1,j2,j3)=polylag !Br
!!$aux2(j1,j2,j3)=poly1r  !dBrdr
!!$aux3(j1,j2,j3)=poly1p  !dBrdp
!!$CALL plag3d(rsten(j1),psten(j2),zsten(j3),Bzsten,overdr,overdp,overdz,&
!!$     rsten,psten,zsten,polylag,poly1r,poly1p,poly1z)
!!$aux4(j1,j2,j3)=poly1z  !dBzdz
!!$aux5(j1,j2,j3)=poly1p  !dBzdp
!!$slope(j1,j2_ind,j3)=-(rsten(j1)*aux2(j1,j2,j3)+aux1(j1,j2,j3)+rsten(j1)*aux4(j1,j2,j3))
!!$enddo
!!$enddo
!!$enddo
!!$
!!$
!!$CALL plag2d_herm1d_3(rrr,philoc,zloc,Bpsten,slope,overdr,overdp,overdz,&
!!$                  &rsten,psten,zsten,polylag,poly1r,poly1p,poly1z)
!!$Bphi=polylag
!!$dBphidr=poly1r
!!$dBphidp=poly1p
!!$dBphidz=poly1z


!!$dBraddr = brf*dBraddr
!!$dBphidz = brf*dBphidz
!!$dBzetdz = brf*dBzetdz
!!$dBphidp = brf*dBphidp
!!$dBzetdp = brf*dBzetdp
!!$Brad = brf*Brad
 
RETURN
1000 STOP 'STOP (field_w7x.f90): !W7-X field file missing!'

END SUBROUTINE field_w7x


 


