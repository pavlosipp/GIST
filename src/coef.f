!-----------------------------------------------------------------------VBA5800
        subroutine coef (nt,sl,ca,cb)                                   VBA5810
!
      use vbal_param
      use cidata
      use ceigen
      use ccoef
!
!.. Implicits ..
      implicit none
!       include 'combal.inc'                                            VBA5820
       integer :: nt
       double precision :: sl, ca, cb
!                                                                       VBA5830
!..model dependent coefficients to be used in the ode                   VBA5840
!                                                                       VBA5850
!..ideal incompressible mhd                                             VBA5860
!                                                                       VBA5870
        sl   = s(nt)                                                    VBA5880
!                                                                       VBA5890
        if (lmodel .ne. 0) then                                         VBA5900
!                                                                       VBA5910
        ca = -c2(nt)                                                    VBA5920
        cb = -c1(nt) - c0(nt) * g                                       VBA5930
!                                                                       VBA5940
        return                                                          VBA5950
        end if                                                          VBA5960
!                                                                       VBA5970
!..ideal mhd marginal point analysis                                    VBA5980
!                                                                       VBA5990
        ca = -c2(nt)                                                    VBA6000
        cb = -c1(nt) * (1. - g)                                         VBA6010
!                                                                       VBA6020
        return                                                          VBA6030
        end subroutine coef                                             VBA6040
