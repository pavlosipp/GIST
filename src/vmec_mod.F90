   MODULE vmec_mod


  character(len=256),PRIVATE::fname
  integer(8),PRIVATE::mcObj=0     ! mcObj array holds the address of C++ object

  real(8),PRIVATE :: dR=0.02,dZ=0.02,dFi=2,epsTrunc=1d-8,epsA=1d-5,pi
  integer(4),PRIVATE :: FALSE=0,TRUE=1

  CONTAINS
  
  !************************************************************
  ! returns .TRUE. if loaded 
  SUBROUTINE bc_setMeshParm(dR_,dZ_,dFi_,epsA_,epsTrunc_)
  real(8),INTENT(in):: dR_,dZ_,dFi_,epsA_,epsTrunc_
    dR = dR_   ! in m
    dZ = dZ_   ! in m
    dFi = dFi_ ! in degree
    epsA = epsA_    ! in m
    epsTrunc = epsTrunc_
    pi = 4*datan(1.d0)
  END SUBROUTINE bc_setMeshParm

  !************************************************************
  ! returns .TRUE. if loaded 
  FUNCTION bc_load(filename) RESULT(ok) ! load Boozer file
  CHARACTER(*), INTENT(in) :: filename
    real(8) mcM3Dcyl2s
    integer(4) mcload,mcCreateMeshUsingSymmetry,MCwriteMesh,mcIsMeshOK
    external mcload,mcfree,mctruncate,mcsetAccuracy
    external mcCreateMeshUsingSymmetry,mcM3DgetdBcyl,MCwriteMesh,mcM3Dcyl2s
    integer(4) i
    logical ok

    pi = 4*datan(1.d0)

    ok = .FALSE.
    if(mcload (mcObj,filename)==FALSE ) return ! return .FALSE. if error

    !************************************************************
    !Create 3dMesh  
    call mctruncate   (mcObj,epsTrunc)      ! truncate spectrum
    call mcsetAccuracy(mcObj,epsA)          ! set accuracy of coordinate transformation

    IF (mcIsMeshOK(mcObj) /= TRUE) THEN
       call mcm3dsetsmax(mcObj,1.1d0)
      i=mcCreateMeshUsingSymmetry(mcObj,dR,dZ,dFi*pi/180)   ! create mesh for whole machine
      if(i==FALSE) return                    !return .FALSE. if error
      i=MCwriteMesh(mcObj,'savedmesh.bin4')  ! save for future use
  !    print '(a, I2// )','MCwriteMesh =', i
    ENDIF
    ok = .TRUE.
    fname = filename
  END FUNCTION bc_load

!iota
  FUNCTION bc_iota(s) RESULT(iota)
    REAL(8),INTENT(in)  :: s
    REAL(8):: iota,mciota
    external mciota
    iota  = mciota(mcObj,s)   
   END FUNCTION bc_iota

!d(iota)/ds
 FUNCTION bc_iotaPrime(s) RESULT(iotaP)
    REAL(8),INTENT(in)  :: s
    REAL(8):: iota,mcIotaprime,iotaP
    external mcIota
    iotaP = mcIotaPrime(mcObj,s)
   END FUNCTION bc_iotaPrime

!toroidal flux 
  FUNCTION bc_torflux(s) RESULT(torflux)
    REAL(8),INTENT(in)  :: s
    REAL(8):: torflux,mctorflux
    external mctorflux
    torflux  = mctorflux(mcObj,s)   
   END FUNCTION bc_torflux

!minor radius
  FUNCTION bc_rminor() RESULT(a)
    REAL(8):: a,mcrminor
    external mcrminor
    a = mcrminor(mcObj)
  END FUNCTION bc_rminor

!major radius = R_00 of LCFS
FUNCTION bc_rmajor() RESULT(a)
    REAL(8):: a,mcR0
    external mcR0
    a = mcR0(mcObj)   
  END FUNCTION bc_rmajor

 !************************************************************
  ! Function retuns the fraction of trapped particles  at flux surface  s 
  ! Input:
  !  s - normalized toroidal flux, 
  ! Output: f_trap(s)
 
  FUNCTION bc_ft(s) RESULT(ft)
    REAL(8),INTENT(in)  :: s
    REAL(8):: ft,MCFTRAPPED
    external MCFTRAPPED
    ft = MCFTRAPPED(mcObj,s)   
   END FUNCTION bc_ft  

  !************************************************************
  ! Function retuns the flux surface label s for the given point cyl
  ! Input:
  !  cyl  - array of size 3 containing cylindrical coordinate of the point
  ! Output:
  !  s - normalized toroidal flux, 
  !       s > 1 if the point cyl lies outside the last closed surface 
  FUNCTION bc_cyl2s(cyl) RESULT(s)
    REAL(8),INTENT(in)  :: cyl(3)
    REAL(8):: s,mcM3Dcyl2s
    external mcM3Dcyl2s
    s  = mcM3Dcyl2s(mcObj,cyl)   ! function returns flux label at the point cyl
   END FUNCTION bc_cyl2s

  !*********************************************************************
  ! Function retuns the cylindrical coordinate for the given point in
  ! coordinate (s,theta_boozer,phi_cylindrical)
  ! Input:
  !  boozer  - array of size 3 containing (s,theta_boozer,phi_cylindrical)
  !            coordinate of the point
  ! Output:
  !  cyl  - array of size 3 containing cylindrical coordinate of the point
  !  s - normalized toroidal flux, 
  !       s > 1 if the point cyl lies outside the last closed surface 
 ! FUNCTION bc_booz2cyl(boozer) RESULT(cyl)
 !   REAL(8),INTENT(in)  :: boozer(3)
 !   REAL(8):: cyl(3)
 !   external mcsemiboozer2cyl
 !   call mcsemiboozer2cyl(mcObj,boozer,cyl)    ! cyl is output 
 ! END FUNCTION bc_booz2cyl


  !*********************************************************************
  ! Function retuns the magnetic field vector B and it derivatives 
  ! for the given point cyl
  ! Input:
  !  cyl  - array of size 3 containing cylindrical coordinate of the point
  ! Output:
  !  B - array with the components of magnetic field in Tesla
  SUBROUTINE bc_Bfield(cyl,B,dBdr,dBdfir,dBdz)
    REAL(8), INTENT(in)  :: cyl(3)
    REAL(8), INTENT(out) :: B(3)      ! cylindrical coordinates of B = (Br, Bfi, Bz)
    REAL(8), INTENT(out) :: dBdr(3)   ! dB/dr, where dBdr(1)=dBr/dr, dBdr(2)=dBfi/dr, dBdr(3)=dBz/dr
    REAL(8), INTENT(out) :: dBdfir(3) ! dB/dfi/r
    REAL(8), INTENT(out) :: dBdz(3)   ! dB/dz
    external mcM3DgetdBcyl

    call mcM3DgetdBcyl(mcObj,cyl,B,dBdr,dBdfir,dBdz) ! B,dBdr,dBdfir,dBdz are at the point cyl
  END SUBROUTINE bc_Bfield

  !***************************************************************************************
  ! Function retuns the cylindrical coordinate of the point given in
  ! magnetic coordinate (s,theta,phi)
  ! Input:
  !  magn  - array of size 3 containing magnetic
  !            coordinate (s,theta,phi) of the point
  ! Output:
  !  cyl  - array of size 3 containing cylindrical coordinate of the point
  FUNCTION bc_magn2cyl(magn) RESULT(cyl)
    REAL(8),INTENT(in) :: magn(3)
    REAL(8):: cyl(3)
    external mcmag2cyl
    call mcmag2cyl(mcObj,magn,cyl)    ! cyl is output 
  END FUNCTION bc_magn2cyl

 !***************************************************************************************
  ! Retuns BOOZER coors (s,theta,zeta) given CYLINDRICAL coords (r,phi,zeta)
  ! Input:
  !  magn  - array of size 3 containing magnetic
  !            coordinate (s,theta,phi) of the point
  ! Output:
  !  cyl  - array of size 3 containing cylindrical coordinate of the point
  FUNCTION bc_cyl2mag(cyl) RESULT(boo)
    REAL(8),INTENT(in) :: cyl(3)
    REAL(8):: boo(3),dum(3)
    external mccyl2mag
    call mccyl2mag(mcObj,cyl,boo,dum)    
  END FUNCTION bc_cyl2mag

  !*********************************************************************
  SUBROUTINE bc_free()
    external mcfree
    call mcfree(mcObj)         ! destroy C++ object and free memory
  END SUBROUTINE bc_free

!=======================================================================================

  SUBROUTINE field_vmec(rrr,ppp,zzz,Brad,Bphi,Bzet,dBraddr,dBraddp,dBraddz,&
     &dBphidr,dBphidp,dBphidz,dBzetdr,dBzetdp,dBzetdz)

  USE type_mod  
  USE specifications_mod, ONLY:sign_of_Bp
   
  REAL(DP), INTENT(IN) :: rrr,ppp,zzz
  REAL(DP), INTENT(OUT) :: Brad,Bphi,Bzet,dBraddr,dBraddp,dBraddz
  REAL(DP), INTENT(OUT) :: dBphidr,dBphidp,dBphidz,dBzetdr,dBzetdp,dBzetdz
  REAL(DP), DIMENSION(3) :: cyl,Bout,dBdrout,dBdfiout,dBdzout 
  INTEGER :: s_ 

    cyl(1)=rrr ; cyl(2)=ppp ; cyl(3)=zzz

!print *,'test',bc_cyl2s(cyl)

    CALL bc_Bfield(cyl,Bout,dBdrout,dBdfiout,dBdzout)
    sign_of_Bp=INT(Bout(2)/ABS(Bout(2)))
    s_=1

    Brad=s_*Bout(1) ; Bphi=s_*Bout(2) ; Bzet=s_*Bout(3)
    dBraddr=s_*dBdrout(1)  ; dBphidr=s_*dBdrout(2)  ; dBzetdr=s_*dBdrout(3)
    dBraddp=s_*dBdfiout(1)*rrr ; dBphidp=s_*dBdfiout(2)*rrr ; dBzetdp=s_*dBdfiout(3)*rrr
    dBraddz=s_*dBdzout(1)  ; dBphidz=s_*dBdzout(2)  ; dBzetdz=s_*dBdzout(3) 

   END SUBROUTINE field_vmec

!==========================================================================================

SUBROUTINE arclength(sss,ttt,ppp,cyl)
USE type_mod
USE specifications_mod,ONLY:vmec_file,vmec_dir
REAL(DP), INTENT(IN) :: sss, ttt, ppp
REAL(DP), DIMENSION(3), INTENT(OUT) :: cyl
REAL(DP), DIMENSION(3) :: in
CHARACTER (len=256) fname
LOGICAL :: first=.TRUE.

fname = trim(vmec_dir)//'/'//trim(vmec_file)

if (first) then
    IF (bc_load(fname) .EQV. .FALSE. ) STOP 'VMEC file not found'
        CALL bc_setMeshParm(0.01d0,0.01d0,.4d0,1d-5,1d-15)
        first=.false.
    endif

    in(1)=sss ; in(2)=ttt ; in(3)=ppp
    cyl=bc_magn2cyl(in)
    END SUBROUTINE arclength

!==================================================================    

SUBROUTINE vmec_data(sss,ttt,ppp,c11_vmec,q_vmec,alpha_vmec,shat_vmec,ftrap,lgrads,lrefs)

    USE type_mod  
    USE lagrange_mod
    USE specifications_mod,ONLY:b0,sign_of_Bp,qpr,profw7,grads,refs,&
                          &vmec_file, table_tag, vmec_dir, R0_vmec

    REAL(DP), INTENT(INOUT) :: sss,ttt,ppp !(rho,theta,phi)
    REAL(DP), INTENT(OUT) :: alpha_vmec,shat_vmec,c11_vmec,q_vmec,ftrap   
    REAL(DP), DIMENSION(3), INTENT(OUT) :: lgrads, lrefs  

    CHARACTER (len=256) fname
    REAL(DP),DIMENSION(:,:),ALLOCATABLE :: vmec_in,cyl
    REAL(DP),DIMENSION(:),ALLOCATABLE :: s,vmec_in_aux,torflux
    REAL(DP) :: sj,qj,qnew,fj
    REAL(DP), DIMENSION(3) :: cylj,in,aux
    INTEGER :: j
    INTEGER, PARAMETER :: Ntable=101
    LOGICAL :: first=.TRUE.


    fname = trim(vmec_dir)//'/'//trim(vmec_file) 

 IF (first) THEN
!    WRITE(*,"(/,2X,A,/)") '=========================== Preparing GENE data ============================'
    IF (bc_load(fname) .EQV. .FALSE. ) STOP 'VMEC file not found'
    CALL bc_setMeshParm(0.01d0,0.01d0,.4d0,1d-5,1d-15) !grid spacing


!Produce table of s,r,q for interpolation
 
    OPEN(88,file='../int/shear_'//trim(vmec_file)//'_' &
     &    //trim(table_tag),ACTION='write')

      DO j=1,Ntable
      sj=(j-1)/100.
      in(1)=sj ; in(2)=0.d0 ; in(3)=0.d0
      cylj=bc_magn2cyl(in)
      qj=1./bc_iota(sj)
      fj=bc_ft(sj)
      WRITE(88,'(4F13.8)') sj,cylj(1),qj,fj
      ENDDO

  CLOSE(88)

first=.FALSE.
ENDIF

!Define normalizing length

     alpha_vmec=bc_rminor() 
     R0_vmec=bc_rmajor()

!Calculate fraction of trapped particles
     ftrap=bc_ft(sss)

!Define B0

     b0=bc_torflux(1.d0)/Pi/alpha_vmec**2

!Define q

    q_vmec=1./bc_iota(sss)

CALL vmec_shat(sss,shat_vmec)

CALL vmec_c11(sss,c11_vmec)

if (profw7) then
CALL gradw7(sss,lgrads,lrefs)
else
lgrads=0.0
endif

grads=lgrads ; refs=lrefs

!Redefine sss to be the r-cylindrical
!ppp to be the phi-cylindrical
!and ttt to be the z-cylindrical 

    in(1)=sss ; in(2)=ttt ; in(3)=ppp
    aux=bc_magn2cyl(in)
    sss=aux(1)
    ppp=aux(2)
    ttt=aux(3)

!------------
CONTAINS
!-----------


SUBROUTINE gradw7(sin,omt,ref)

    REAL(DP), INTENT(IN) :: sin
    REAL(DP), DIMENSION(3), INTENT(OUT) :: omt,ref
    integer,parameter::ntable=51
    REAL(DP),DIMENSION(Ntable) :: sdat,cdat,tidat,tedat,ndat
    INTEGER :: stat=0,k
    INTEGER, PARAMETER :: mp=4
    REAL(DP), DIMENSION(mp) :: ssten,csten,tisten,testen,nsten
    INTEGER, DIMENSION(mp) :: ind
    REAL(DP) :: polylag,poly1x,iota_int,ds,dsmin1,sqrtsin,slope


    OPEN(UNIT=62, FILE='../inp/prof_w7.dat', STATUS="OLD", ACTION="READ", ERR=10)
    DO j=1,Ntable
       READ(62,*) sdat(j),tidat(j),tedat(j),ndat(j)
    ENDDO
    CLOSE(62)

!NB: The radial label must be sqrt(s) 

    ds=sdat(2)-sdat(1)
    dsmin1=1./ds
    sqrtsin=sqrt(sin)

    CALL indef(sqrtsin,0.d0,dsmin1,Ntable,ind)

    DO k=1,mp
       ssten(k)=(ind(k)-1)*ds
       tisten(k)=tidat(ind(k))
       testen(k)=tedat(ind(k))
       nsten(k)=ndat(ind(k))
    ENDDO

!Determine Ti gradient
CALL plag1d(sqrtsin,tisten,dsmin1,ssten,polylag,poly1x)
slope=-poly1x/polylag
ref(1)=polylag
omt(1)=slope

!Determine Te gradient
CALL plag1d(sqrtsin,testen,dsmin1,ssten,polylag,poly1x)
slope=-poly1x/polylag
ref(2)=polylag
omt(2)=slope

!Determine n gradient
CALL plag1d(sqrtsin,nsten,dsmin1,ssten,polylag,poly1x)
slope=-poly1x/polylag
ref(3)=polylag
omt(3)=slope

RETURN
10 PRINT *, 'STOP! No profile file prof_w7.dat found in the inp directory!'
stop

END SUBROUTINE	gradw7

!---------------------------------------------------------------------------!
!===========================================================================!

SUBROUTINE vmec_shat(sin,shat)

    REAL(DP), INTENT(IN) :: sin
    REAL(DP), INTENT(OUT) :: shat
    REAL(DP),DIMENSION(Ntable) :: sdat,cdat,qdat
    INTEGER :: stat=0,k
    INTEGER, PARAMETER :: mp=4
    REAL(DP), DIMENSION(mp) :: ssten,csten,qsten
    INTEGER, DIMENSION(mp) :: ind
    REAL(DP) :: polylag,poly1x,iota_int,ds,dsmin1

    OPEN(UNIT=32, FILE='../int/shear_'//trim(vmec_file)//'_'//trim(table_tag),&
               &STATUS="OLD", ACTION="READ", ERR=10)
    DO j=1,Ntable
       READ(32,*) sdat(j),cdat(j),qdat(j)
    ENDDO

    CLOSE(32)

!Determine q

    ds=sdat(2)-sdat(1)
    dsmin1=1./ds

    CALL indef(sin,0.d0,dsmin1,Ntable,ind)

    DO k=1,mp
       ssten(k)=(ind(k)-1)*ds
       qsten(k)=qdat(ind(k))
    ENDDO

    CALL plag1d(sin,qsten,dsmin1,ssten,polylag,poly1x)
 
!Determine shat

!    shat=2.*sin/polylag*poly1x
!    qpr=poly1x
    
    shat=-2*sin*bc_iotaprime(sin)/bc_iota(sin) !yury
    
   
    RETURN
10 PRINT *, '!No iota.dat link found in input directory. q0,shat will not be calculated!'
  
 
END SUBROUTINE vmec_shat

!----------------------------------------------------------------------

SUBROUTINE vmec_c11(sin_,c11)

    REAL(DP), INTENT(IN) :: sin_
    REAL(DP), INTENT(OUT) :: c11
    REAL(DP),DIMENSION(Ntable) :: sdat,rdat,qdat
    INTEGER :: stat=0,k
    INTEGER, PARAMETER :: mp=4
    REAL(DP), DIMENSION(mp) :: ssten,rsten,qsten
    INTEGER, DIMENSION(mp) :: ind
    REAL(DP) :: polylag,poly1x,iota_int,ds,dsmin1

   OPEN(UNIT=32, FILE='../int/shear_'//trim(vmec_file)//'_'//trim(table_tag),&
              & STATUS="OLD", ACTION="READ", ERR=10)
    DO j=1,Ntable
       READ(32,*) sdat(j),rdat(j),qdat(j)
    ENDDO

    CLOSE(32)

    ds=sdat(2)-sdat(1)
    dsmin1=1./ds

    CALL indef(sin_,0.d0,dsmin1,Ntable,ind)

    DO k=1,mp
       ssten(k)=(ind(k)-1)*ds
       rsten(k)=rdat(ind(k))
    ENDDO

    CALL plag1d(sin_,rsten,dsmin1,ssten,polylag,poly1x)
 
!Determine c11

    c11=.5*alpha_vmec/SQRT(sin_)/poly1x

    RETURN
10 PRINT *, '!No shear file in directory int found!'

  
 
END SUBROUTINE vmec_c11

!------------------------

 END SUBROUTINE vmec_data

!============================================================================


 END MODULE vmec_mod

