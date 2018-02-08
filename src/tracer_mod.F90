
MODULE tracer_mod
  USE type_mod
  USE specifications_mod
  IMPLICIT NONE

!=============================================================================!
CONTAINS               !                MODULE SUBPROGRAMS                    !
!=============================================================================!

 !============================= SET DPHI ==================================!

  SUBROUTINE set_dphi(code_,line_num_,trace_p_num_)
    INTEGER, INTENT(IN) :: code_,line_num_,trace_p_num_
    INTEGER :: n
    !Initial set (only used for all codes exept FINDIF)
     IF(trace_p_num_==0) THEN
       dphi = dist/REAL(resol-1)
       IF(trace_type(line_num_) == 0 .AND. (code_==1 .OR. code_==4)&
            & .AND. MOD(resol,2) == 0) dphi=dist/REAL(resol)
     ELSE
       IF( code_ == 2) THEN !Only for FINDIF
         ! trace_type > 1 means that we have to go
         !   trace_type - 1 time in forward direction to
         !   reach plates (the last point in a line)
         ! trace_type < -1 means that we have to go
         !   -trace_type - 1 time in backward direction to
         !   reach plates (the first point in a line)
         IF(trace_type(line_num_)>1.OR.    &
              & trace_type(line_num_)<-1) THEN
           ! Trace type values begins from |2|
           ! that's why we have to subtract |1|
           IF( trace_p_num_==(trace_type(line_num_)- &
                & SIGN(1,trace_type(line_num_))) ) THEN
                !Did we reach the first/last point???
             n=line_num_+trace_p_num_
             dphi = p1_in(n)*Pi !Here we have already dphi in rad
           ENDIF
           ! We have to reset correct dphi for the case
           !   of trace_type(line_num_) < -1
           IF(trace_type(line_num_) < -1.AND.trace_p_num_>0) THEN
             dphi = dist/REAL(resol-1)
           ENDIF
         ENDIF
       ENDIF
     ENDIF
  END SUBROUTINE set_dphi
!============================= SET RESOLUTION ============================!

  SUBROUTINE set_resolution(line_num_)
    INTEGER, INTENT(IN) :: line_num_
    INTEGER :: res_mod

    res_mod=MOD(resol,2)
    resol_b=(resol-res_mod)/2
    resol_f=resol_b - 1 + res_mod
    IF(trace_type(line_num_) == 1) THEN 
       resol_b=0 ; resol_f=resol-1
    ELSEIF(trace_type(line_num_) == -1) THEN
       resol_f=0 ; resol_b=resol-1
    ELSEIF(trace_type(line_num_)>1) THEN
      resol_f=trace_type(line_num_)-1
    ELSEIF(trace_type(line_num_)<-1)THEN
      resol_b=-trace_type(line_num_)-1
    ENDIF

  END SUBROUTINE set_resolution

!============================= PREPARE CUTS ==============================!

  SUBROUTINE prep_cuts(resol)

    INTEGER, INTENT(IN) :: resol 
    REAL :: dist
    INTEGER :: k
    
    cut(1)=1 ; cut(cuts+1)=resol+3 !for ghost cells
   
!Multiple equidistant cuts
 
    IF (cuts > 1) THEN
       dist=INT((resol-1)/REAL(cuts))
       DO k=2,cuts
          cut(k)=(k-1)*dist
       ENDDO
    ENDIF
   
  END SUBROUTINE prep_cuts


!=========================== MATRIX D=INV(C) =============================!

  SUBROUTINE matrix_d
    INTRINSIC TINY

    detc=c(1,1)*c(2,2)-c(1,2)*c(2,1)
    IF (detc < TINY(1.0)) STOP &
         '!Fatal error: (r0,z0) out of bounds or line hit boundary!'
            
    d(1,1)=c(2,2)/detc
    d(1,2)=-c(1,2)/detc
    d(1,3)=(c(1,2)*c(2,3)-c(1,3)*c(2,2))/detc
    d(2,1)=-c(2,1)/detc
    d(2,2)=c(1,1)/detc
    d(2,3)=(c(1,3)*c(2,1)-c(1,1)*c(2,3))/detc
    d(3,1)=0.
    d(3,2)=0.
    d(3,3)=1.

  END SUBROUTINE matrix_d

!============================= CLEBSCH CONTRAVARIANT METRICS ====================!

  SUBROUTINE gij_contravar

    gij(1,1) = c(1,1)*c(1,1) + c(1,2)*c(1,2) + c(1,3)*c(1,3)/r**2
    gij(1,2) = c(1,1)*c(2,1) + c(1,2)*c(2,2) + c(1,3)*c(2,3)/r**2
    gij(1,3) = c(1,1)*c(3,1) + c(1,2)*c(3,2) + c(1,3)*c(3,3)/r**2
    gij(2,1) = gij(1,2) 
    gij(2,2) = c(2,1)*c(2,1) + c(2,2)*c(2,2) + c(2,3)*c(2,3)/r**2
    gij(2,3) = c(2,1)*c(3,1) + c(2,2)*c(3,2) + c(2,3)*c(3,3)/r**2
    gij(3,1) = gij(1,3)
    gij(3,2) = gij(2,3)
    gij(3,3) = c(3,1)*c(3,1) + c(3,2)*c(3,2) + c(3,3)*c(3,3)/r**2

  END SUBROUTINE gij_contravar

!=========================== CLEBSCH COVARIANT METRICS ===========================!

  SUBROUTINE gij_covar 

    g_ij(1,1) = d(1,1)*d(1,1) + d(2,1)*d(2,1) + d(3,1)*d(3,1)*r**2
    g_ij(1,2) = d(1,1)*d(1,2) + d(2,1)*d(2,2) + d(3,1)*d(3,2)*r**2
    g_ij(1,3) = d(1,1)*d(1,3) + d(2,1)*d(2,3) + d(3,1)*d(3,3)*r**2
    g_ij(2,1) = g_ij(1,2) 
    g_ij(2,2) = d(1,2)*d(1,2) + d(2,2)*d(2,2) + d(3,2)*d(3,2)*r**2
    g_ij(2,3) = d(1,2)*d(1,3) + d(2,2)*d(2,3) + d(3,2)*d(3,3)*r**2
    g_ij(3,1) = g_ij(1,3)
    g_ij(3,2) = g_ij(2,3)
    g_ij(3,3) = d(1,3)*d(1,3) + d(2,3)*d(2,3) + d(3,3)*d(3,3)*r**2

  END SUBROUTINE gij_covar

!=============================== CALCULATE QUANTITIES  ===========================!

  SUBROUTINE quantities

    INTRINSIC SQRT,ABS,MATMUL

!Field and derivatives 

    Bfield=SQRT(Br**2+Bz**2+Bp**2)
    divB=Br/r+dBrdr+dBpdp/r+dBzdz 
    b3=Bp/(Bfield*r)  !b^/tau
    Bt=Bz*(r-rmag)/((r-rmag)**2+(z-z0)**2)-Br*(z-z0)/((r-rmag)**2+(z-z0)**2) !B^theta

!Derivs wrt cylindrical system

    dBdr=Br*dBrdr/Bfield+Bz*dBzdr/Bfield+Bp*dBpdr/Bfield
    dBdp=Br*dBrdp/Bfield+Bz*dBzdp/Bfield+Bp*dBpdp/Bfield
    dBdz=Br*dBrdz/Bfield+Bz*dBzdz/Bfield+Bp*dBpdz/Bfield

!Derivs wrt Clebsch system
               
    dBdv1=d(1,1)*dBdr+d(2,1)*dBdz
    dBdv2=d(1,2)*dBdr+d(2,2)*dBdz 
    dBdtau=d(1,3)*dBdr+d(2,3)*dBdz+dBdp 




!Jacobian (using 3 methods)
 
    jac_chain=-r/detc !Chain rule *gives also sign*
    jac_contravar = SQRT(1./detgij(gij)) !Via contravariant elements
    jac_covar = SQRT(detgij(g_ij)) !Via covariant elements
    jac=jac_contravar


!Stream function JB^phi

    inv=jac*Bfield*b3

!Cartesian coordinates

    cart_x=r*COS(phi)
    cart_y=r*SIN(phi)

!Poloidal angle

    IF (z > z0 .AND. r > rmag)   theta=ATAN((z-z0)/(r-rmag))
    IF (z > z0 .AND. r < rmag)   theta=Pi-ATAN((z-z0)/(rmag-r))
    IF (z < z0 .AND. r > rmag)   theta=ATAN((z-z0)/(r-rmag))
    IF (z <= z0 .AND. r < rmag)   theta=-Pi+ATAN((z-z0)/(r-rmag))
    
   

!----------------------------- Inhomogeneity ----------------------------------!

!Auxilliaries

    fac1=(gij(1,3)*gij(2,2)-gij(1,2)*gij(2,3))/&
         &(gij(1,1)*gij(2,2)-gij(1,2)*gij(1,2))
    fac2=(gij(1,1)*gij(2,3)-gij(1,2)*gij(1,3))/&
         &(gij(1,1)*gij(2,2)-gij(1,2)*gij(1,2))
    fac3=SQRT(gij(1,1))*Bfield*jac*b3

!Covariant components

    k1=dBdv1/Bfield+fac1*dBdtau/Bfield
    k2=dBdv2/Bfield+fac2*dBdtau/Bfield

!Normal curvature
 
    k_norm=(gij(1,1)*dBdv1+gij(1,2)*dBdv2&
         &+gij(1,3)*dBdtau)/SQRT(gij(1,1))/Bfield

!Geodesic curvature

    k_geo=-(dBdv2+fac2*dBdtau)/fac3              

!square of Total curvature (explicit expression)

    kappa_sq=(dBdv1**2*gij(1,1)+dBdv2**2*gij(2,2)&
         &+dBdtau**2*(gij(3,3)-b3**2)&
         &+2*dBdv1*dBdv2*gij(1,2)+2*dBdv1*dBdtau*gij(1,3)&
         &+2*dBdv2*dBdtau*gij(2,3))/Bfield**2  

!square of Total curvature (implicit expression) 

    kappa_sq_imp=k_norm**2+k_geo**2          

!-------------------------- Consistency checks ----------------------------!

!Curvature
    zero_kappa=ABS(kappa_sq-kappa_sq_imp)

!Clebsch system

    zero_cl1=Br*c(1,1)+Bz*c(1,2)+Bp/r*c(1,3)
    zero_cl2=Br*c(2,1)+Bz*c(2,2)+Bp/r*c(2,3)
    zero_cl3=ABS(Bfield**2-inv**2*(gij(1,1)*gij(2,2)-gij(1,2)**2))

!---------------------- Determinant of metrics -----------------------------!
  CONTAINS
!---------------------------------------------------------------------------!
 
    FUNCTION detgij(x)
 
      REAL(DP) :: detgij
      REAL(DP), DIMENSION(3,3) :: x
      
      detgij=x(1,1)*x(2,2)*x(3,3)+2*x(1,2)*x(2,3)*x(1,3)&
           &-x(2,2)*x(1,3)**2-x(1,1)*x(2,3)**2-x(3,3)*x(1,2)**2 
            
    END FUNCTION detgij

!----------------------------------------------------------------------------!


  END SUBROUTINE quantities


!=========================== LOCATE SYMMETRY PLANE ==========================!

  SUBROUTINE symm_plane(r0,p0,zmin,zmax,z0,ierr)
    USE specifications_mod, ONLY : device
    USE field_mod

    REAL(DP), INTENT(IN) :: r0,p0,zmin,zmax 
    REAL(DP), INTENT(OUT) :: z0 
    INTEGER, INTENT(OUT) :: ierr
    REAL(DP) :: Br_min,Br_max 
    REAL(DP) :: left,right,br,bp,bz,dbrdr,dbrdz,dbpdr,dbpdz,dbzdr,dbzdz
    REAL(DP) :: dbrdp,dbpdp,dbzdp,ztry_min,ztry_max,zdiff

!!$!---------------------- scan Br ----------------------------------- 
!!$
!!$    integer :: j,nz=500
!!$    real(dp) :: dz,zin
!!$
!!$open(1,file='br.dat')
!!$
!!$dz=(zmax-zmin)/real(nz-1)
!!$
!!$do j=1,nz-1
!!$zin=zmin+dz*(j-1)
!!$CALL field(r0,p0,zin,br,bp,bz,dbrdr,dbrdp,dbrdz,dbpdr,dbpdp,&
!!$            dbpdz,dbzdr,dbzdp,dbzdz,device) 
!!$    write(1,*) zin,br,sqrt(br**2+bp**2+bz**2)
!!$ enddo
!!$
!!$
!!$close(1)
!---------------------------------------------------------------------

!bracket zero

   ztry_min=0.3*zmin ; ztry_max=0.3*zmax
    DO
       CALL field(r0,p0,ztry_min,br,bp,bz,dbrdr,dbrdp,dbrdz,dbpdr,dbpdp,&
            dbpdz,dbzdr,dbzdp,dbzdz,device) 
       Br_min=br

       CALL field(r0,p0,ztry_max,br,bp,bz,dbrdr,dbrdp,dbrdz,dbpdr,dbpdp,&
            dbpdz,dbzdr,dbzdp,dbzdz,device) 
       Br_max=br
       IF (Br_min < 0 .AND. Br_max > 0) THEN
          left=zmin ; right=zmax
          EXIT
       ELSEIF (Br_min > 0 .AND. Br_max < 0) THEN
          left=zmax ; right=zmin
          EXIT
       ELSE
          zdiff=ABS(ztry_min-ztry_max)
          IF (zdiff > 2.*SPACING(zdiff)) THEN
             ztry_min=ztry_min/2. ; ztry_max=ztry_max/2.
          ELSE
             PRINT *, 'Symmetry plane could not be located!'
             ierr=1
             RETURN
          ENDIF
       ENDIF
    ENDDO

!bisect

    DO
       z0=(left+right)/2. 
       
       CALL field(r0,p0,z0,br,bp,bz,dbrdr,dbrdp,dbrdz,dbpdr,dbpdp,&
            dbpdz,dbzdr,dbzdp,dbzdz,device) 
       
       IF (br < 0.) THEN
          left=z0
       ELSEIF (br > 0.) THEN
          right=z0
       ELSE
          EXIT !root found
       ENDIF

       IF (ABS(right-left) < 2.*SPACING(z0)) EXIT !root found with best possible accuracy

    ENDDO

!Success

    ierr=0
    swsymm=.TRUE.
  END SUBROUTINE symm_plane
  
!=============================== CALCULATE SAFETY FACTOR ==============================!

  SUBROUTINE safety_factor
    USE lagrange_mod
    REAL(DP) :: ftor_pr,fpol_pr
    CHARACTER(80) :: iota_file    

    SELECT CASE (device)
           
    CASE(1) !W7X
       iota_file='../input/W7X/iota.dat'
       IF (.NOT. vmec) CALL extender_iota(r_norm)

    CASE(2) !NCSX
       iota_file='../input/NCSX/iota.dat'
       IF (.NOT. vmec) CALL extender_iota(r_norm)

!    CASE(3) !D3D

    CASE DEFAULT
 
       !q0=dF_tor/dF_pol
       ftor_pr=SUM(jac_out*Bp_out/r_out)
       fpol_pr=SUM(jac_out*Bt_out)
       safety=ABS(ftor_pr/fpol_pr)

    END SELECT
!
CONTAINS
!
  SUBROUTINE extender_iota(r)
    REAL(DP), INTENT(IN) :: r

    INTEGER :: count=-1,stat=0,j,k
    INTEGER, PARAMETER :: mp=4
    REAL(DP), DIMENSION(mp) :: rsten,iotasten
    INTEGER, DIMENSION(mp) :: ind
    REAL(DP), DIMENSION(:), ALLOCATABLE :: iota,rout
    REAL(DP) :: polylag,poly1x,iota_int

    OPEN(UNIT=32, FILE=iota_file, STATUS="OLD", ACTION="READ", ERR=10)
    !count entries
    DO WHILE (stat == 0) 
    READ(32,*,IOSTAT=stat)
    count=count+1
    ENDDO

    ALLOCATE(iota(count),rout(count))
    REWIND(32)

    DO j=1,count
       READ(32,*) iota(j),rout(j)
    ENDDO

    CLOSE(32)

!Interpolation

    CALL indef_uneven(r,rout,count,ind)

    DO k=1,mp
       rsten(k)=rout(ind(k))
       iotasten(k)=iota(ind(k))
    ENDDO

    CALL plag1d_uneven(r,iotasten,rsten,polylag,poly1x)
    iota_int=polylag

    safety=1./iota_int

    RETURN
10 IF(.NOT.silent_mode) PRINT *, &
  &'!No iota.dat link found in input directory. q0,shat will not be calculated!'

  END SUBROUTINE extender_iota

END SUBROUTINE safety_factor


!========================= CALCULATE Shat =====================================!

SUBROUTINE shear
!REAL(DP), DIMENSION(SIZE(gij_out(1,2,:))) :: y does not work with sun compiler
REAL(DP), DIMENSION(:), ALLOCATABLE :: y
REAL(DP) :: intercept

  ALLOCATE(y(SIZE(gij_out(1,2,:))))
  y(:)=gij_out(1,2,:)/gij_out(1,1,:)
  CALL least_sq(p_out,y,slope,intercept)
  DEALLOCATE(y)
!
CONTAINS
!
  SUBROUTINE least_sq(x,y,m,c)
    REAL(DP), DIMENSION(:), INTENT(IN) :: x, y
    REAL(DP), INTENT(OUT) :: m, c

    INTEGER :: n
    REAL :: sum_x, sum_xsq, sum_xy, sum_y

    n=SIZE(x)

    !Evaluate sums

    sum_x=SUM(x)
    sum_y=SUM(y)
    sum_xy=DOT_PRODUCT(x,y)
    sum_xsq=DOT_PRODUCT(x,x)

    !Evaluate coefficients m & c

    m=(n*sum_xy-sum_x*sum_y)/  &
         (n*sum_xsq-sum_x*sum_x)
    c=(sum_y-m*sum_x)/n
print *,'m=',m,'c=',c
  END SUBROUTINE least_sq

END SUBROUTINE shear

  SUBROUTINE local_shear(x,m)
    USE lagrange_mod,ONLY: plag1d_uneven,indef_uneven
    REAL(DP), DIMENSION(:), INTENT(IN) :: x
    REAL(DP), DIMENSION(:), INTENT(OUT) :: m
    INTEGER, PARAMETER:: Mp=4
    INTEGER:: ind,j,nphi
    INTEGER,DIMENSION(Mp):: indphi
    REAL(DP), DIMENSION(Mp):: phisten, ysten
    REAL(DP), DIMENSION(:), ALLOCATABLE:: y
    REAL(DP):: polylag, poly1r
    nphi=SIZE(gij_out(1,2,:))
    ALLOCATE(y(nphi))

    y(:)=gij_out(1,2,:)/gij_out(1,1,:)
    DO ind=1,nphi
       
       CALL indef_uneven(x(ind),x,nphi,indphi)
       DO j=1,Mp
          phisten(j)=x(indphi(j))
          ysten(j)=y(indphi(j))
       ENDDO
       CALL plag1d_uneven(x(ind),ysten,phisten,polylag,poly1r)
       m(ind)=poly1r*safety !this is local shear 
    ENDDO

  END SUBROUTINE local_shear
  

  SUBROUTINE Check_periodicity(curr,phi_period)
    INTEGER, INTENT(IN):: curr
    REAL(DP), INTENT(OUT):: phi_period

    INTEGER:: i,ref,min,mincnt,ptsperiod=0
    REAL(DP), DIMENSION(:),ALLOCATABLE:: sumerr,errors
    REAL(DP):: err,relerr
    LOGICAL:: firstrun1=.TRUE.,firstrun2=.TRUE.,first=.TRUE.
    SAVE
    IF (first) THEN 
       ALLOCATE(errors(5))
       first=.FALSE.
    ENDIF

    IF (genecnt==1) THEN
       IF (firstrun1) THEN
          ALLOCATE(sumerr(resol))
          phi_period=0
          firstrun1=.FALSE.
          ref=1
          min=1
          mincnt=0
       ENDIF
    ENDIF

    IF (genecnt==2) THEN
       IF (firstrun2) THEN
          ALLOCATE(sumerr(resol))
          firstrun2=.FALSE.
          ref=1
       ENDIF
    ENDIF

IF (curr.gt.ref+1) THEN
    
       !Calculate relative deviations for g11,g23,g33,jac and Bfield, 
       !which are all assumed periodic.
       err=gij_out(1,1,ref)-gij_out(1,1,curr)
       if (gij_out(1,1,ref)/=0) then 
          relerr=err/gij_out(1,1,ref)
       else
          relerr=0
       endif
       errors(1)=relerr

       err=gij_out(2,3,ref)-gij_out(2,3,curr)
       if (gij_out(2,3,ref)/=0) then 
          relerr=err/gij_out(2,3,ref)
       else
          relerr=0
       endif
       errors(2)=relerr

       
       err=gij_out(3,3,ref)-gij_out(3,3,curr)
       if (gij_out(3,3,ref)/=0) then 
          relerr=err/gij_out(3,3,ref)
       else
          relerr=0
       endif
       errors(3)=relerr
       
       err=jac_out(ref)-jac_out(curr)
       if (jac_out(ref)/=0) then 
          relerr=err/jac_out(ref)
       else
          relerr=0
       endif
       errors(4)=relerr
       
       err=Bfield_out(ref)-Bfield_out(curr)
       if (Bfield_out(ref)/=0) then 
          relerr=err/Bfield_out(ref)
       else
          relerr=0
       endif
       errors(5)=relerr
       
       !Add relative errors:
       sumerr(curr)=abs(errors(1))+abs(errors(2))+abs(errors(3))+abs(errors(4))+abs(errors(5))
       IF ((genecnt==1).AND.(curr.gt.60)) THEN
          IF ((sumerr(curr).gt.sumerr(curr-1)).AND.(sumerr(curr-1).lt.sumerr(curr-2))) THEN
!             print *,'curr', curr

             IF (mincnt==0.AND.sumerr(curr-1).lt.0.4) THEN 
                mincnt=mincnt+1
                IF (sumerr(curr-1).lt.sumerr(min).AND.(sumerr(curr-1).lt.0.4)) THEN
                   !                print *,'mincnt',mincnt
                   IF (curr-1-min.lt.int(0.1*REAL(min)/REAL(mincnt))) THEN
                      min=(curr-1+min)/2
                   ELSE
                      min=curr-1
                   ENDIF
                   phi_period=p_out(curr-1)/REAL(mincnt)           
                   ptsperiod=min
!                   print *,phi_period,mincnt
                ENDIF
             ELSE IF (mincnt.ge.1.AND.sumerr(curr-1).lt.0.4.AND. &
                  &curr-1-min.gt.int(0.1*ptsperiod)) THEN
                mincnt=mincnt+1
 !               print *,curr,0.1*ptsperiod,mincnt               
                IF (sumerr(curr-1).lt.sumerr(min).AND.(sumerr(curr-1).lt.0.4)) THEN
                   !                print *,'mincnt',mincnt
                   IF (curr-1-min.lt.int(0.1*ptsperiod)) THEN
                      min=(curr-1+min)/2
                   ELSE
                      min=curr-1
                   ENDIF
                   phi_period=p_out(curr-1)/REAL(mincnt)           
                ENDIF
             ENDIF

          ENDIF
       ENDIF

ELSE 
   sumerr(min)=1
   sumerr(curr)=1
ENDIF

!      IF (genecnt==1.AND.curr==resol) THEN
!          open(55,file='err',status='replace')
!          write(55,*) '# phi |  periodic quantities: sum of relative deviations'
!          DO i=1,resol
!             write(55,*) p_out(i),sumerr(i)
!          ENDDO
!          close(55)  
!       ENDIF
       IF (genecnt==1.AND.curr==resol) DEALLOCATE(sumerr)

       IF ((genecnt==2).AND.(curr==resol)) THEN
          IF (mod(resol,2)/=0) THEN
             print *,'Periodicity check (should be <~ 1E-4)=',sumerr(resol)
             IF (sumerr(resol).gt.1e-2) STOP 'Forcing periodicity failed. Please set force_period=.FALSE. for this flux tube.'
          ENDIF
       ENDIF

       IF ((genecnt==2).AND.(curr==resol)) THEN 
          shatgij=(gij_out(1,2,resol)/gij_out(1,1,resol)-&
            &gij_out(1,2,ref)/gij_out(1,1,ref))/(2.*Pi)
!          print *,shatgij
!necessary for runs with multiple field lines
          firstrun1=.TRUE.
          firstrun2=.TRUE.
          DEALLOCATE(sumerr)
       ENDIF

  END SUBROUTINE Check_Periodicity
  
  

  !------------------------------------------------------------------------------!
  
END MODULE tracer_mod
