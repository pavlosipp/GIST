
MODULE rk5_mod
  USE type_mod
  USE util_mod
  IMPLICIT NONE
!================================================================================!
CONTAINS           !                MODULE SUBPROGRAMS                           !
!================================================================================!
  SUBROUTINE odeint(ystart,x1,x2,eps,h1,hmin,derivs)
    IMPLICIT NONE

    INTERFACE derivs_interface
       SUBROUTINE derivs(x,y,dydx)
         USE type_mod
         REAL(DP), INTENT(IN) :: x
         REAL(DP), DIMENSION(:), INTENT(IN) :: y
         REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx
       END SUBROUTINE derivs
    END INTERFACE derivs_interface

    REAL(DP), DIMENSION(:), INTENT(INOUT) :: ystart
    REAL(DP), INTENT(IN) :: x1,x2,eps,h1,hmin
    REAL(DP), PARAMETER :: TINY=1.0e-30_sp
    INTEGER(I4B), PARAMETER :: MAXSTP=1000000
    INTEGER(I4B) :: nstp
    LOGICAL :: first=.TRUE.
    REAL(DP) :: h,hdid,hnext,x
    REAL(DP), DIMENSION(size(ystart)) :: dydx,y,yscal

    x=x1
    h=sign(h1,x2-x1)
    y(:)=ystart(:)
    
    do nstp=1,MAXSTP
       CALL derivs(x,y,dydx)
       yscal(:)=abs(y(:))+abs(h*dydx(:))+TINY
       
       if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x
       CALL rkqs(y,dydx,x,h,eps,yscal,hdid,hnext,derivs)
       
       if ((x-x2)*(x2-x1) >= 0.0) then
          ystart(:)=y(:)
          RETURN
       end if
       if (abs(hnext) < hmin)&
            call nrerror('stepsize smaller than minimum in odeint')
       h=hnext
    end do
    call nrerror('too many steps in odeint')
    
  END SUBROUTINE odeint

!================================================================================!

  SUBROUTINE rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
    
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: y
    REAL(DP), DIMENSION(:), INTENT(IN) :: dydx,yscal
    REAL(DP), INTENT(INOUT) :: x
    REAL(DP), INTENT(IN) :: htry,eps
    REAL(DP), INTENT(OUT) :: hdid,hnext
    
    INTERFACE myintfc
       SUBROUTINE derivs(x,y,dydx)
         USE type_mod
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: x
         REAL(DP), DIMENSION(:), INTENT(IN) :: y
         REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx
       END SUBROUTINE derivs
    END INTERFACE
    
    INTEGER(I4B) :: ndum
    REAL(DP) :: errmax,h,htemp,xnew
    REAL(DP), DIMENSION(size(y)) :: yerr,ytemp,dydx_sign
    REAL(DP), PARAMETER :: SAFETY=0.9_DP,PGROW=-0.2_DP,PSHRNK=-0.25_DP,&
         ERRCON=1.89e-4

    ndum=assert_eq(size(y),size(dydx),size(yscal),'rkqs')
    h=htry
    do
       call rkck(y,dydx,x,h,ytemp,yerr,derivs)
       errmax=0.
       errmax=maxval(abs(yerr(:)/yscal(:)))/eps
       if (errmax <= 1.0) exit
       htemp=SAFETY*h*(errmax**PSHRNK)
       h=sign(max(abs(htemp),0.1_DP*abs(h)),h)
       xnew=x+h
       if (xnew == x) call nrerror('stepsize underflow in rkqs')
    end do
    if (errmax > ERRCON) then
       hnext=SAFETY*h*(errmax**PGROW)
    else
       hnext=5.0_DP*h
    end if
    hdid=h
    x=x+h
    y(:)=ytemp(:)

  END SUBROUTINE rkqs
  
!================================================================================!

  SUBROUTINE rkck(y,dydx,x,h,yout,yerr,derivs)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: y,dydx
    REAL(DP), INTENT(IN) :: x,h
    REAL(DP), DIMENSION(:), INTENT(OUT) :: yout,yerr
    INTERFACE myintfc2
       SUBROUTINE derivs(x,y,dydx)
         USE type_mod
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: x
         REAL(DP), DIMENSION(:), INTENT(IN) :: y
         REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx
       END SUBROUTINE derivs
    END INTERFACE
    INTEGER(I4B) :: ndum
    REAL(DP), DIMENSION(size(y)) :: ak2,ak3,ak4,ak5,ak6,ytemp
    REAL(DP), PARAMETER :: A2=0.2_DP,A3=0.3_DP,A4=0.6_DP,A5=1.0_DP,&
         A6=0.875_DP,B21=0.2_DP,B31=3.0_DP/40.0_DP,B32=9.0_DP/40.0_DP,&
         B41=0.3_DP,B42=-0.9_DP,B43=1.2_DP,B51=-11.0_DP/54.0_DP,&
         B52=2.5_DP,B53=-70.0_DP/27.0_DP,B54=35.0_DP/27.0_DP,&
         B61=1631.0_DP/55296.0_DP,B62=175.0_DP/512.0_DP,&
         B63=575.0_DP/13824.0_DP,B64=44275.0_DP/110592.0_DP,&
         B65=253.0_DP/4096.0_DP,C1=37.0_DP/378.0_DP,&
         C3=250.0_DP/621.0_DP,C4=125.0_DP/594.0_DP,&
         C6=512.0_DP/1771.0_DP,DC1=C1-2825.0_DP/27648.0_DP,&
         DC3=C3-18575.0_DP/48384.0_DP,DC4=C4-13525.0_DP/55296.0_DP,&
         DC5=-277.0_DP/14336.0_DP,DC6=C6-0.25_DP
    ndum=assert_eq(size(y),size(dydx),size(yout),size(yerr),'rkck')
    ytemp=y+B21*h*dydx
    call derivs(x+A2*h,ytemp,ak2)
    ytemp=y+h*(B31*dydx+B32*ak2)
    call derivs(x+A3*h,ytemp,ak3)
    ytemp=y+h*(B41*dydx+B42*ak2+B43*ak3)
    call derivs(x+A4*h,ytemp,ak4)
    ytemp=y+h*(B51*dydx+B52*ak2+B53*ak3+B54*ak4)
    call derivs(x+A5*h,ytemp,ak5)
    ytemp=y+h*(B61*dydx+B62*ak2+B63*ak3+B64*ak4+B65*ak5)
    call derivs(x+A6*h,ytemp,ak6)
    yout=y+h*(C1*dydx+C3*ak3+C4*ak4+C6*ak6)
    yerr=h*(DC1*dydx+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6)
  END SUBROUTINE rkck

END MODULE rk5_mod
