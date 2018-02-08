!-------------------------------------------------------------------------------------------!
!                    3rd degree LAGRANGE INTERPOLATION on UNIFORM mesh                      !
!Features:                                                                                  ! 
!          Calculation of 1D 4-point stencil : indef                                        !
!          1D interpolation on 4-point stencil : plag1d          
!          2D interpolation on 4x4-point stencil : plag2d
!          3D interpolation on 4x4x4-point stencil using product method : plag3d            !
!-------------------------------------------------------------------------------------------!

MODULE lagrange_mod
  USE type_mod
  IMPLICIT NONE
!----------------------------------------------------------
  PRIVATE
  PUBLIC :: indef, plag1d, plag2d, plag3d, herm1d_7, plag2d_herm1d_7, plag2d_herm1d_3
  PUBLIC :: indef_uneven, plag1d_uneven, indef_per
!----------------------------------------------------------

  INTEGER, PARAMETER :: Mp=4
  REAL(DP), PARAMETER :: One6=0.1666666666666667D0
  INTEGER :: i,j,k

!================================================================================!
CONTAINS           !                MODULE SUBPROGRAMS                           !
!================================================================================!

!========================== CALCULATE 4-POINT UNIFORM STENCIL ===================!

   SUBROUTINE indef(u,umin,dm1,nup,sten)

     ! u        : coordinate of the point 
     ! umin     : lower bound of mesh 
     ! dm1      : reciprocal of mesh interval
     ! nup      : number of mesh points
     ! sten     : array of stencil indices (centered around u)

     REAL(DP), INTENT(IN) :: u,umin,dm1
     INTEGER, INTENT(IN) :: nup
     INTEGER, DIMENSION(Mp), INTENT(OUT) :: sten
     INTRINSIC INT

!First index of stencil

     sten(1) = INT((u-umin)*dm1)

!Check for lower boundary

     IF (sten(1) <= 0) sten(1) = 1

!Last index of stencil
     
     sten(mp) = sten(1) + Mp - 1

!Check for upper boundary

     IF (sten(Mp) > nup ) THEN
        sten(Mp) = nup
        !Correct first index
        sten(1) = sten(Mp) - Mp + 1
     ENDIF

!Indermediate stencil indices
    
     DO i=2,Mp-1
        sten(i) = sten(i-1) + 1
     ENDDO

   END SUBROUTINE indef

!---------------------------------------------------------------------------------!
!                        Stencil on a PERIODIC grid                               !         
!---------------------------------------------------------------------------------!

   SUBROUTINE indef_per(u,umin,dm1,nup,sten)

     ! u        : coordinate of the point 
     ! umin     : lower bound of mesh 
     ! dm1      : reciprocal of mesh interval
     ! nup      : number of mesh points
     ! sten     : array of stencil indices (centered around u)


     REAL(DP), INTENT(IN) :: u,umin,dm1
     INTEGER, INTENT(IN) :: nup
     INTEGER, DIMENSION(Mp), INTENT(OUT) :: sten
     INTRINSIC INT

!First index of stencil

     sten(1) = INT((u-umin)*dm1)

!Other stencil indices
    
     DO i=2,Mp
        sten(i) = sten(i-1) + 1
     ENDDO

   END SUBROUTINE indef_per

!---------------------------------------------------------------------------!
!                     Stencil on a NON-uniform grid                         !
!---------------------------------------------------------------------------!

  SUBROUTINE indef_uneven(u,ord,nup,sten)

  ! u        : coordinate of the point 
  ! ord      : ordinates of uneven mesh
  ! nup      : number of mesh points
  ! sten     : array of stencil indices (centered around u)

      REAL(DP), INTENT(IN) :: u
      REAL(DP), DIMENSION(nup), INTENT(IN) :: ord
      REAL(DP), DIMENSION(nup) :: du
      INTEGER, INTENT(IN) :: nup
      INTEGER, DIMENSION(Mp), INTENT(OUT) :: sten
      INTEGER :: i
      INTEGER, DIMENSION(nup) :: ind_temp
      INTRINSIC INT

!First index of stencil

      du=ABS(u-ord)
      ind_temp=MINLOC(du)

      IF(u >= ord(ind_temp(1))) THEN
         sten(1)=ind_temp(1)-1
      ELSE
         sten(1)=ind_temp(1)-2
      ENDIF

!Check for lower boundary
      IF (sten(1) <= 0) sten(1) = 1
!Last index of stencil
      sten(Mp) = sten(1) + Mp - 1
!Check for upper boundary
      IF (sten(Mp) > nup ) THEN
         sten(Mp) = nup
      !Correct first index
         sten(1) = sten(Mp) - Mp + 1
      ENDIF

!Indermediate stencil indices
      DO i=2,Mp-1
         sten(i) = sten(i-1) + 1
      ENDDO

    END SUBROUTINE indef_uneven


!=========================== 1D LAGRANGE INTERPOLATION =============================!

   SUBROUTINE plag1d(x,fp,dxm1,xp,poly1d,p1x1d)

     ! x      : coordinate of the point 
     ! dxm1   : reciprocal of mesh interval
     ! xp     : stencil indices
     ! poly1d : polynomial evaluated at x
     ! poly1x  : first derivative evaluated at x
 
     REAL(DP), INTENT(IN) :: x,dxm1
     REAL(DP), DIMENSION(Mp), INTENT(IN) :: xp,fp
     REAL(DP), INTENT(OUT) :: poly1d,p1x1d
     REAL(DP), DIMENSION(Mp) :: cx,cx1

!Determine Lagrange polynomial
     
     CALL cardinals(x,xp,dxm1,cx)
     poly1d = 0.
     DO i=1,Mp
        poly1d = poly1d + fp(i)*cx(i)
     ENDDO

!Determine first derivative
     
     CALL cardinals1(x,xp,dxm1,cx1)
     p1x1d = 0.
     DO i=1,Mp
        p1x1d = p1x1d + fp(i)*cx1(i)
     ENDDO

   END SUBROUTINE plag1d

!--------------------------------------------------------------------------!
!                   1D interpolation on NON-uniform grid                   !
!--------------------------------------------------------------------------!

  SUBROUTINE plag1d_uneven(x,fp,xp,poly1d,p1x1d)

  ! x      : coordinate of the point 
  ! xp     : stencil indices

  ! poly1d : polynomial evaluated at x
  ! poly1x  : first derivative evaluated at x
 
      REAL(DP), INTENT(IN) :: x
      REAL(DP), DIMENSION(Mp), INTENT(IN) :: xp,fp
      REAL(DP), INTENT(OUT) :: poly1d,p1x1d
      REAL(DP), DIMENSION(Mp) :: cx,cx1

!Determine Lagrange polynomial

      CALL cardinals_uneven(x,xp,cx)
      poly1d = 0.
      DO i=1,Mp
         poly1d = poly1d + fp(i)*cx(i)
      ENDDO

!Determine first derivative

      CALL cardinals1_uneven(x,xp,cx1)
      p1x1d = 0.
      DO i=1,Mp
         p1x1d = p1x1d + fp(i)*cx1(i)
      ENDDO

    END SUBROUTINE plag1d_uneven

!=========================== 2D LAGRANGE INTERPOLATION =============================!
  
   SUBROUTINE plag2d(x,y,fp,dxm1,dym1,xp,yp,poly2d,poly1x,poly1y)

     ! (x,y)            : coordinate of the point 
     ! dxm1,dym1        : reciprocals of mesh intervals for each direction
     ! xp,yp            : stencil indices for each direction
     ! poly2d : polynomial evaluated at (x,y)
     ! poly1x,poly1y : first partial derivatives evaluated at (x,y)

     REAL(DP), INTENT(IN) :: x,y,dxm1,dym1
     REAL(DP), DIMENSION(Mp,Mp), INTENT(IN) :: fp
     REAL(DP), DIMENSION(Mp), INTENT(IN) :: xp,yp
     REAL(DP), INTENT(OUT) :: poly2d,poly1x,poly1y
     REAL(DP), DIMENSION(Mp) :: cx,cy,cx1,cy1

!Independent interpolation in each direction      
     
     CALL cardinals(x,xp,dxm1,cx)
     CALL cardinals(y,yp,dym1,cy)

!Calculate inteprolating polynomial using product rule
     
     poly2d = 0.
     DO j=1,Mp
        DO i=1,Mp
           poly2d = poly2d + fp(i,j)*cx(i)*cy(j)
        ENDDO
     ENDDO

     CALL cardinals1(x,xp,dxm1,cx1)
     CALL cardinals1(y,yp,dym1,cy1)

!Calculate first partial derivatives
     
     poly1x=0.
     DO j=1,Mp
        DO i=1,Mp
           poly1x = poly1x + fp(i,j)*cx1(i)*cy(j)
        ENDDO
     ENDDO

     poly1y=0.
     DO j=1,Mp
        DO i=1,Mp
           poly1y = poly1y + fp(i,j)*cx(i)*cy1(j)
        ENDDO
     ENDDO

   END SUBROUTINE plag2d

!=========================== 3D LAGRANGE INTERPOLATION =============================!

   SUBROUTINE plag3d(x,y,z,fp,dxm1,dym1,dzm1,xp,yp,zp,&
                     &poly3d,poly1x,poly1y,poly1z)

     ! (x,y,z)          : coordinate of the point 
     ! dxm1,dym1,dzm1   : reciprocals of mesh intervals for each direction
     ! xp,yp,zp         : stencil indices for each direction
     ! poly3d : polynomial evaluated at (x,y,z)
     ! poly1x,poly1y,poly1z  : first partial derivatives evaluated at (x,y,z)

     REAL(DP), INTENT(IN) :: x,y,z,dxm1,dym1,dzm1
     REAL(DP), DIMENSION(Mp,Mp,Mp), INTENT(IN) :: fp
     REAL(DP), DIMENSION(Mp), INTENT(IN) :: xp,yp,zp
     REAL(DP), INTENT(OUT) :: poly3d,poly1x,poly1y,poly1z
     REAL(DP), DIMENSION(Mp) :: cx,cy,cz,cx1,cy1,cz1

!Independent interpolation in each direction      
     
     CALL cardinals(x,xp,dxm1,cx)
     CALL cardinals(y,yp,dym1,cy)
     CALL cardinals(z,zp,dzm1,cz)

!Calculate inteprolating polynomial using product rule
     
     poly3d = 0.
     DO k=1,Mp
        DO j=1,Mp
           DO i=1,Mp
              poly3d = poly3d + fp(i,j,k)*cx(i)*cy(j)*cz(k)
           ENDDO
        ENDDO
     ENDDO

     CALL cardinals1(x,xp,dxm1,cx1)
     CALL cardinals1(y,yp,dym1,cy1)
     CALL cardinals1(z,zp,dzm1,cz1)

!Calculate first partial derivatives
     
     poly1x=0.
     DO k=1,Mp
        DO j=1,Mp
           DO i=1,Mp
              poly1x = poly1x + fp(i,j,k)*cx1(i)*cy(j)*cz(k)
           ENDDO
        ENDDO
     ENDDO

     poly1y=0.
     DO k=1,Mp
        DO j=1,Mp
           DO i=1,Mp
              poly1y = poly1y + fp(i,j,k)*cx(i)*cy1(j)*cz(k)
           ENDDO
        ENDDO
     ENDDO

     poly1z=0.
     DO k=1,Mp
        DO j=1,Mp
           DO i=1,Mp
              poly1z = poly1z + fp(i,j,k)*cx(i)*cy(j)*cz1(k)
           ENDDO
        ENDDO
     ENDDO
      

   END SUBROUTINE plag3d

!============================ COMPUTE CARDINALS =================================!

   SUBROUTINE cardinals(u,up,dm1,cu)

     !Compute cardinal functions for Lagrange interpolation at u
       
     REAL(DP), DIMENSION(Mp), INTENT(IN) :: up
     REAL(DP), INTENT(IN) :: u,dm1
     REAL(DP), DIMENSION(Mp), INTENT(OUT) :: cu
     REAL(DP) :: du3
      
     du3 = dm1**3
 
     cu(1) = (u - up(2)) * (u - up(3)) * (u - up(4)) * (-One6*du3)
     cu(2) = (u - up(1)) * (u - up(3)) * (u - up(4)) * (0.5*du3)
     cu(3) = (u - up(1)) * (u - up(2)) * (u - up(4)) * (-0.5*du3)
     cu(4) = (u - up(1)) * (u - up(2)) * (u - up(3)) * (One6*du3)
 
   END SUBROUTINE cardinals

!---------------------------------------------------------------------

    SUBROUTINE cardinals_uneven(u,up,cu)

      REAL(DP), DIMENSION(Mp), INTENT(IN) :: up
      REAL(DP), INTENT(IN) :: u
      REAL(DP), DIMENSION(Mp), INTENT(OUT) :: cu
      REAL(DP) :: fac1,fac2,fac3,fac4
 
      fac1=(up(1)-up(2))*(up(1)-up(3))*(up(1)-up(4))
      fac2=(up(2)-up(1))*(up(2)-up(3))*(up(2)-up(4))
      fac3=(up(3)-up(1))*(up(3)-up(2))*(up(3)-up(4))
      fac4=(up(4)-up(1))*(up(4)-up(2))*(up(4)-up(3))
      
      cu(1) = (u - up(2)) * (u - up(3)) * (u - up(4)) / fac1
      cu(2) = (u - up(1)) * (u - up(3)) * (u - up(4)) / fac2
      cu(3) = (u - up(1)) * (u - up(2)) * (u - up(4)) / fac3
      cu(4) = (u - up(1)) * (u - up(2)) * (u - up(3)) / fac4
  
    END SUBROUTINE cardinals_uneven

!======================== COMPUTE FIRST DERIVS OF CARDINALS =======================!
 
   SUBROUTINE cardinals1(u,up,dm1,cu1)

     !Compute first derivatives of cardinal functions at u

     REAL(DP), DIMENSION(Mp), INTENT(IN) :: up
     REAL(DP), INTENT(IN) :: u,dm1
     REAL(DP), DIMENSION(Mp), INTENT(OUT) :: cu1
     REAL(DP) :: du3

     du3 = dm1**3

     cu1(1) = (  (u - up(3))*(u - up(4)) &
          &          + (u - up(2))*(u - up(4)) &
          &          + (u - up(2))*(u - up(3)) &
          &         )* (-One6*du3)
     cu1(2) = (  (u - up(3))*(u - up(4)) &
          &          + (u - up(1))*(u - up(4)) &
          &          + (u - up(1))*(u - up(3)) &
          &         ) * (0.5*du3)
     cu1(3) = (  (u - up(2))*(u - up(4)) &
          &          + (u - up(1))*(u - up(4)) &
          &          + (u - up(1))*(u - up(2)) &
          &         )* (-0.5*du3)
     cu1(4) = (  (u - up(2))*(u - up(3)) &
          &          + (u - up(1))*(u - up(3)) &
          &          + (u - up(1))*(u - up(2)) &
          &         )* (One6*du3)

   END SUBROUTINE cardinals1

!---------------------------------------------------------------------

  SUBROUTINE cardinals1_uneven(u,up,cu1)

      REAL(DP), DIMENSION(Mp), INTENT(IN) :: up
      REAL(DP), INTENT(IN) :: u
      REAL(DP), DIMENSION(Mp), INTENT(OUT) :: cu1
      REAL(DP) :: fac1,fac2,fac3,fac4

      fac1=(up(1)-up(2))*(up(1)-up(3))*(up(1)-up(4))
      fac2=(up(2)-up(1))*(up(2)-up(3))*(up(2)-up(4))
      fac3=(up(3)-up(1))*(up(3)-up(2))*(up(3)-up(4))
      fac4=(up(4)-up(1))*(up(4)-up(2))*(up(4)-up(3))


      cu1(1) = (  (u - up(3))*(u - up(4)) &
     &          + (u - up(2))*(u - up(4)) &
     &          + (u - up(2))*(u - up(3)) &
     &         )/ fac1
      cu1(2) = (  (u - up(3))*(u - up(4)) &
     &          + (u - up(1))*(u - up(4)) &
     &          + (u - up(1))*(u - up(3)) &
     &         )/ fac2
      cu1(3) = (  (u - up(2))*(u - up(4)) &
     &          + (u - up(1))*(u - up(4)) &
     &          + (u - up(1))*(u - up(2)) &
     &         )/ fac3
      cu1(4) = (  (u - up(2))*(u - up(3)) &
     &          + (u - up(1))*(u - up(3)) &
     &          + (u - up(1))*(u - up(2)) &
     &         )/ fac4

    END SUBROUTINE cardinals1_uneven

!==========================================================================================!
!                          1-d HERMITE INTERPOLATION 7th degree                            !
!==========================================================================================!

      SUBROUTINE herm1d_7(x,fp,dfp,dxm1,xp,poly1d,p1x1d)

  ! x      : coordinate of the point 
  ! dxm1   : reciprocal of mesh interval
  ! xp     : ordinate at stencil points
  ! fp     : value at stencil points
  ! dfp    : slope at stencil points

  ! poly1d : polynomial evaluated at x
  ! poly1x  : first derivative evaluated at x
 
      REAL(DP), INTENT(IN) :: x,dxm1
      REAL(DP), DIMENSION(Mp), INTENT(IN) :: xp,fp,dfp
      REAL(DP), INTENT(OUT) :: poly1d,p1x1d
      REAL(DP), DIMENSION(Mp) :: cx,cx1,cx2,cx3

CALL card_herm_7(x,xp,dxm1,cx) 
CALL card1_herm_7(x,xp,dxm1,cx1)
CALL der_card_herm_7(x,xp,dxm1,cx2)
CALL der_card1_herm_7(x,xp,dxm1,cx3)


!Determine Hermite polynomial
      
      poly1d = 0.
      DO i=1,Mp
         poly1d = poly1d + fp(i)*cx(i) + dfp(i)*cx1(i)
      ENDDO

!Determine first derivative
      
      p1x1d = 0.

      DO i=1,Mp
         p1x1d = p1x1d + fp(i)*cx2(i) + dfp(i)*cx3(i)
      ENDDO

    END SUBROUTINE herm1d_7

!---------------------------------------------------------------------
    
    SUBROUTINE card_herm_7(u,up,dm1,cu)

      !Compute cardinal functions (type H) for Hermite interpolation at u
       
      REAL(DP), DIMENSION(Mp), INTENT(IN) :: up
      REAL(DP), INTENT(IN) :: u,dm1
      REAL(DP), DIMENSION(Mp), INTENT(OUT) :: cu
      REAL(DP) :: du7
      
      du7 = dm1**7
 
      cu(1) = 11/108.*du7*(u-up(1)+3/11./dm1)*(u-up(2))**2*(u-up(3))**2*(u-up(4))**2
      cu(2) = .25*du7*(u-up(2)+1/dm1)*(u-up(1))**2*(u-up(3))**2*(u-up(4))**2
      cu(3) = -.25*du7*(u-up(3)-1/dm1)*(u-up(1))**2*(u-up(2))**2*(u-up(4))**2
      cu(4) = -11/108.*du7*(u-up(4)-3/11./dm1)*(u-up(1))**2*(u-up(2))**2*(u-up(3))**2
        
 
    END SUBROUTINE card_herm_7


!---------------------------------------------------------------------

       SUBROUTINE card1_herm_7(u,up,dm1,cu)

      !Compute cardinal functions (type M) for Hermite interpolation at u
       
      REAL(DP), DIMENSION(Mp), INTENT(IN) :: up
      REAL(DP), INTENT(IN) :: u,dm1
      REAL(DP), DIMENSION(Mp), INTENT(OUT) :: cu
      REAL(DP) :: du6
      
      du6 = dm1**6
 
      cu(1) = 1./36.*du6*(u-up(1))*(u-up(2))**2*(u-up(3))**2*(u-up(4))**2
      cu(2) = .25*du6*(u-up(2))*(u-up(1))**2*(u-up(3))**2*(u-up(4))**2
      cu(3) = .25*du6*(u-up(3))*(u-up(1))**2*(u-up(2))**2*(u-up(4))**2
      cu(4) = 1./36.*du6*(u-up(4))*(u-up(1))**2*(u-up(2))**2*(u-up(3))**2
 
    END SUBROUTINE card1_herm_7


!---------------------------------------------------------------------

       SUBROUTINE der_card_herm_7(u,up,dm1,cu)

 !Compute first derivative of cardinal functions (type H) 
       
      REAL(DP), DIMENSION(Mp), INTENT(IN) :: up
      REAL(DP), INTENT(IN) :: u,dm1
      REAL(DP), DIMENSION(Mp), INTENT(OUT) :: cu
      REAL(DP) :: du7
      
      du7 = dm1**7
 
      cu(1) = 11./108.*du7*(u-up(2))**2*(u-up(3))**2*(u-up(4))**2&
      &+11./54.*du7*(u-up(1)+3./11./dm1)*(u-up(2))*(u-up(3))**2*(u-up(4))**2&
      &+11./54.*du7*(u-up(1)+3./11./dm1)*(u-up(3))*(u-up(2))**2*(u-up(4))**2&
      &+11./54.*du7*(u-up(1)+3./11./dm1)*(u-up(4))*(u-up(2))**2*(u-up(3))**2
      
      cu(2) = .25*du7*(u-up(1))**2*(u-up(3))**2*(u-up(4))**2&
      &+.5*du7*(u-up(2)+1./dm1)*(u-up(1))*(u-up(3))**2*(u-up(4))**2&
      &+.5*du7*(u-up(2)+1./dm1)*(u-up(3))*(u-up(1))**2*(u-up(4))**2&
      &+.5*du7*(u-up(2)+1./dm1)*(u-up(4))*(u-up(1))**2*(u-up(3))**2

      cu(3) = -.25*du7*(u-up(1))**2*(u-up(2))**2*(u-up(4))**2&
      &-.5*du7*(u-up(3)-1./dm1)*(u-up(1))*(u-up(2))**2*(u-up(4))**2&
      &-.5*du7*(u-up(3)-1./dm1)*(u-up(2))*(u-up(1))**2*(u-up(4))**2&
      &-.5*du7*(u-up(3)-1./dm1)*(u-up(4))*(u-up(1))**2*(u-up(2))**2

      cu(4) = -11./108.*du7*(u-up(1))**2*(u-up(2))**2*(u-up(3))**2&
      &-11./54.*du7*(u-up(4)-3./11./dm1)*(u-up(1))*(u-up(2))**2*(u-up(3))**2&
      &-11./54.*du7*(u-up(4)-3./11./dm1)*(u-up(2))*(u-up(1))**2*(u-up(3))**2&
      &-11./54.*du7*(u-up(4)-3./11./dm1)*(u-up(3))*(u-up(1))**2*(u-up(2))**2

    END SUBROUTINE der_card_herm_7



!---------------------------------------------------------------------

       SUBROUTINE der_card1_herm_7(u,up,dm1,cu)

 !Compute first derivative of cardinal functions (type M) 
       
      REAL(DP), DIMENSION(Mp), INTENT(IN) :: up
      REAL(DP), INTENT(IN) :: u,dm1
      REAL(DP), DIMENSION(Mp), INTENT(OUT) :: cu
      REAL(DP) :: du6
      
      du6 = dm1**6
 
      cu(1) = 1./36.*du6*(u-up(2))**2*(u-up(3))**2*(u-up(4))**2&
      &+1./18.*du6*(u-up(1))*(u-up(2))*(u-up(3))**2*(u-up(4))**2&
      &+1./18.*du6*(u-up(1))*(u-up(3))*(u-up(2))**2*(u-up(4))**2&
      &+1./18.*du6*(u-up(1))*(u-up(4))*(u-up(2))**2*(u-up(3))**2

      cu(2) = .25*du6*(u-up(1))**2*(u-up(3))**2*(u-up(4))**2&
      &+.5*du6*(u-up(2))*(u-up(1))*(u-up(3))**2*(u-up(4))**2&
      &+.5*du6*(u-up(2))*(u-up(3))*(u-up(1))**2*(u-up(4))**2&
      &+.5*du6*(u-up(2))*(u-up(4))*(u-up(1))**2*(u-up(3))**2

      cu(3) = .25*du6*(u-up(1))**2*(u-up(2))**2*(u-up(4))**2&
      &+.5*du6*(u-up(3))*(u-up(1))*(u-up(2))**2*(u-up(4))**2&
      &+.5*du6*(u-up(3))*(u-up(2))*(u-up(1))**2*(u-up(4))**2&
      &+.5*du6*(u-up(3))*(u-up(4))*(u-up(1))**2*(u-up(2))**2

      cu(4) = 1./36.*du6*(u-up(1))**2*(u-up(2))**2*(u-up(3))**2&
      &+1./18.*du6*(u-up(4))*(u-up(1))*(u-up(2))**2*(u-up(3))**2&
      &+1./18.*du6*(u-up(4))*(u-up(2))*(u-up(1))**2*(u-up(3))**2&
      &+1./18.*du6*(u-up(4))*(u-up(3))*(u-up(1))**2*(u-up(2))**2

 
    END SUBROUTINE der_card1_herm_7


!=========================== 2D LAGRANGE + 1D HERMITE 7th degree =============================!

   SUBROUTINE plag2d_herm1d_7(x,y,z,fp,dfp,dxm1,dym1,dzm1,xp,yp,zp,&
                     &poly3d,poly1x,poly1y,poly1z)

     ! (x,y,z)          : coordinate of the point 
     ! dxm1,dym1,dzm1   : reciprocals of mesh intervals for each direction
     ! xp,yp,zp         : stencil indices for each direction
     ! poly3d : polynomial evaluated at (x,y,z)
     ! poly1x,poly1y,poly1z  : first partial derivatives evaluated at (x,y,z)

     REAL(DP), INTENT(IN) :: x,y,z,dxm1,dym1,dzm1
     REAL(DP), DIMENSION(Mp,Mp,Mp), INTENT(IN) :: fp,dfp
     REAL(DP), DIMENSION(Mp), INTENT(IN) :: xp,yp,zp
     REAL(DP), INTENT(OUT) :: poly3d,poly1x,poly1y,poly1z
     REAL(DP), DIMENSION(Mp) :: cx,cy,cz,cx1,cy1,cy11,cy2,cy21,cz1

!Independent interpolation in r,z,phi - directions      
     
     CALL cardinals(x,xp,dxm1,cx)
     CALL card_herm_7(y,yp,dym1,cy1)
     CALL card1_herm_7(y,yp,dym1,cy2)
     CALL cardinals(z,zp,dzm1,cz)

!Calculate inteprolating polynomial using product rule
     
     poly3d = 0.
     DO k=1,Mp
        DO j=1,Mp
           DO i=1,Mp
              poly3d = poly3d + fp(i,j,k)*cx(i)*cy1(j)*cz(k)&
                             &+ dfp(i,j,k)*cx(i)*cy2(j)*cz(k)
           ENDDO
        ENDDO
     ENDDO

     CALL cardinals1(x,xp,dxm1,cx1)
     CALL der_card_herm_7(y,yp,dym1,cy11)
     CALL der_card1_herm_7(y,yp,dym1,cy21) 
     CALL cardinals1(z,zp,dzm1,cz1)

!Calculate first partial derivatives
     
     poly1x=0.
     DO k=1,Mp
        DO j=1,Mp
           DO i=1,Mp
              poly1x = poly1x + fp(i,j,k)*cx1(i)*cy1(j)*cz(k)&
                             &+ dfp(i,j,k)*cx1(i)*cy2(j)*cz(k)
           ENDDO
        ENDDO
     ENDDO

     poly1y=0.
     DO k=1,Mp
        DO j=1,Mp
           DO i=1,Mp
              poly1y = poly1y + fp(i,j,k)*cx(i)*cy11(j)*cz(k)&
                             &+ dfp(i,j,k)*cx(i)*cy21(j)*cz(k)
           ENDDO
        ENDDO
     ENDDO

     poly1z=0.
     DO k=1,Mp
        DO j=1,Mp
           DO i=1,Mp
              poly1z = poly1z + fp(i,j,k)*cx(i)*cy1(j)*cz1(k)&
                             &+ dfp(i,j,k)*cx(i)*cy2(j)*cz1(k)
           ENDDO
        ENDDO
     ENDDO
      

   END SUBROUTINE plag2d_herm1d_7

!==========================================================================================!
!                            HERMITE interpolation 3rd degree
!==========================================================================================!    

    SUBROUTINE card_herm_3(u,up,dm1,cu)

      !Compute cardinal functions (type H) for Hermite interpolation at u
       
      REAL(DP), DIMENSION(Mp), INTENT(IN) :: up
      REAL(DP), INTENT(IN) :: u,dm1
      REAL(DP), DIMENSION(2), INTENT(OUT) :: cu
      REAL(DP) :: du
      
      du = dm1**3
 
      cu(1) = 2*du*(u-up(1)+0.5/dm1)*(u-up(2))**2
      cu(2) =-2*du*(u-up(2)-0.5/dm1)*(u-up(1))**2
 
    END SUBROUTINE card_herm_3


!---------------------------------------------------------------------

       SUBROUTINE card1_herm_3(u,up,dm1,cu)

      !Compute cardinal functions (type M) for Hermite interpolation at u
       
      REAL(DP), DIMENSION(Mp), INTENT(IN) :: up
      REAL(DP), INTENT(IN) :: u,dm1
      REAL(DP), DIMENSION(2), INTENT(OUT) :: cu
      REAL(DP) :: du
      
      du = dm1**2
 
      cu(1) = du*(u-up(1))*(u-up(2))**2
      cu(2) = du*(u-up(2))*(u-up(1))**2
 
    END SUBROUTINE card1_herm_3


!---------------------------------------------------------------------


    SUBROUTINE der_card_herm_3(u,up,dm1,cu)

      !Derivative of  cardinal functions (type H) for Hermite interpolation at u
       
      REAL(DP), DIMENSION(Mp), INTENT(IN) :: up
      REAL(DP), INTENT(IN) :: u,dm1
      REAL(DP), DIMENSION(2), INTENT(OUT) :: cu
      REAL(DP) :: du
      
      du = dm1**3
 
      cu(1) = 2*du*(u-up(2))**2+4*du*(u-up(1)+0.5/dm1)*(u-up(2)) 
      cu(2) =-2*du*(u-up(1))**2-4*du*(u-up(2)-0.5/dm1)*(u-up(1)) 
 
    END SUBROUTINE der_card_herm_3


!---------------------------------------------------------------------

       SUBROUTINE der_card1_herm_3(u,up,dm1,cu)

      !Compute cardinal functions (type M) for Hermite interpolation at u
       
      REAL(DP), DIMENSION(Mp), INTENT(IN) :: up
      REAL(DP), INTENT(IN) :: u,dm1
      REAL(DP), DIMENSION(2), INTENT(OUT) :: cu
      REAL(DP) :: du
      
      du = dm1**2
 
      cu(1) = du*(u-up(2))**2+2*du*(u-up(1))*(u-up(2))
      cu(2) = du*(u-up(1))**2+2*du*(u-up(1))*(u-up(2)) 
 
    END SUBROUTINE der_card1_herm_3

!=========================== 2D LAGRANGE + 1D HERMITE 3th degree =============================!

   SUBROUTINE plag2d_herm1d_3(x,y,z,fp,dfp,dxm1,dym1,dzm1,xp,yp,zp,&
                     &poly3d,poly1x,poly1y,poly1z)

     ! (x,y,z)          : coordinate of the point 
     ! dxm1,dym1,dzm1   : reciprocals of mesh intervals for each direction
     ! xp,yp,zp         : stencil indices for each direction
     ! poly3d : polynomial evaluated at (x,y,z)
     ! poly1x,poly1y,poly1z  : first partial derivatives evaluated at (x,y,z)

     REAL(DP), INTENT(IN) :: x,y,z,dxm1,dym1,dzm1
     REAL(DP), DIMENSION(Mp,2,Mp), INTENT(IN) :: fp,dfp
     REAL(DP), DIMENSION(Mp), INTENT(IN) :: xp,yp,zp
     REAL(DP), INTENT(OUT) :: poly3d,poly1x,poly1y,poly1z
     REAL(DP), DIMENSION(Mp) :: cx,cz,cx1,cz1
     REAL(DP), DIMENSION(2) :: cy1,cy11,cy2,cy21

!Independent interpolation in r,z,phi - directions      
     
     CALL cardinals(x,xp,dxm1,cx)
     CALL card_herm_3(y,yp,dym1,cy1)
     CALL card1_herm_3(y,yp,dym1,cy2)
     CALL cardinals(z,zp,dzm1,cz)

!Calculate inteprolating polynomial using product rule
     
     poly3d = 0.
     DO k=1,Mp
        DO j=1,2
           DO i=1,Mp
              poly3d = poly3d + fp(i,j,k)*cx(i)*cy1(j)*cz(k)&
                             &+ dfp(i,j,k)*cx(i)*cy2(j)*cz(k)
           ENDDO
        ENDDO
     ENDDO

     CALL cardinals1(x,xp,dxm1,cx1)
     CALL der_card_herm_3(y,yp,dym1,cy11)
     CALL der_card1_herm_3(y,yp,dym1,cy21) 
     CALL cardinals1(z,zp,dzm1,cz1)

!Calculate first partial derivatives
     
     poly1x=0.
     DO k=1,Mp
        DO j=1,2
           DO i=1,Mp
              poly1x = poly1x + fp(i,j,k)*cx1(i)*cy1(j)*cz(k)&
                             &+ dfp(i,j,k)*cx1(i)*cy2(j)*cz(k)
           ENDDO
        ENDDO
     ENDDO

     poly1y=0.
     DO k=1,Mp
        DO j=1,2
           DO i=1,Mp
              poly1y = poly1y + fp(i,j,k)*cx(i)*cy11(j)*cz(k)&
                             &+ dfp(i,j,k)*cx(i)*cy21(j)*cz(k)
           ENDDO
        ENDDO
     ENDDO

     poly1z=0.
     DO k=1,Mp
        DO j=1,2
           DO i=1,Mp
              poly1z = poly1z + fp(i,j,k)*cx(i)*cy1(j)*cz1(k)&
                             &+ dfp(i,j,k)*cx(i)*cy2(j)*cz1(k)
           ENDDO
        ENDDO
     ENDDO
      

   END SUBROUTINE plag2d_herm1d_3
















 END MODULE lagrange_mod
