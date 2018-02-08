!-----------------------------------------------------------------------VBA1490
        subroutine equin                                                VBA1500
!
      use vbal_param
      use cidata
      use crdata
      use cfoura
!
!.. Implicits ..
      implicit none
!       include 'combal.inc'                                            VBA1510
!.. Local Scalars ..
      integer :: mn
!                                                                       VBA1520
!..formatted equilibrium input                                          VBA1530
!cc                                                                     VBA1540
!--READ-IN FOURIER AMPLITUDES OF THE BALLOONING COEFFICIENTS AND THE    VBA1550
!--VALUES OF s, q, p' and V' ON THE FLUX SURFACE UNDER CONSIDERATION.   VBA1560
        read (25,98) rs, qs, pp, vp, qppsip, sqgbsq, ftp ,pth, qprime
        if (lmodel.eq.2)qppsip=qppsip*qppsip
        if (lmodel.eq.3)qppsip=0.
        read (25,99) (a0mn(mn),mn=1,mpnt)                               VBA1580
        read (25,99) (a1mn(mn),mn=1,mpnt)                               VBA1590
        read (25,99) (a2mn(mn),mn=1,mpnt)                               VBA1600
        read (25,99) (rgmn(mn),mn=1,mpnt)                               VBA1610
        read (25,99) (d0mn(mn),mn=1,mpnt)                               VBA1620
        read (25,99) (d1mn(mn),mn=1,mpnt)                               VBA1630
        read (25,99) (gssmn(mn),mn=1,mpnt)
        read (25,99) (gstmn(mn),mn=1,mpnt)
        read (25,99) (bsmn(mn),mn=1,mpnt)
!        read (25,99) (d2mn(mn),mn=1,mpnt)                              VBA1640
!        read (25,99) (d2mn(mn),mn=1,mpnt)                              VBA1650
!        read (25,99) (d2mn(mn),mn=1,mpnt)                              VBA1660
!  98   format (1p2e22.14,1p2e16.8,1p1e22.14,1p2e16.8)                  VBA1670
   98   format (1p9e16.8)                                               VBA1670
   99   format (1p6e22.14)                                              VBA1680
!                                                                       VBA1690
        return                                                          VBA1700
        end subroutine equin                                            VBA1710
