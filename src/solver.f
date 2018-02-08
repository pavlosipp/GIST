!-----------------------------------------------------------------------VBA3040
        subroutine solver                                               VBA3050
!                                                                       VBA3070
!  the ordinary differential equation ballooning mode equation is integrVBA3080
!  zero crossings are searched for                                      VBA3090
!  the growth rate is iterated until an eigenvalue is approximated      VBA3100
!  error checking is done                                               VBA3110
!x                                                                      VBA3120
      use vbal_param
      use cidata
      use crdata
      use ceigen
      use cresis
      use ccoef
!
C.. Implicits ..
      implicit none
!       include 'combal.inc'                                            VBA3130
!
! ... LOCAL SCALARS
      integer :: lode, ngu, ngv, ns, nconv, j
      double precision ::
     &               su, gu, gv, x1, x2, x3, x4, f2, tol, f1, g1, slop
!.. External Calls ..
      external odeint
!
! .. Intrinsic Functions ..
      intrinsic  sqrt, abs
!                                                                       VBA3130
!..initialize                                                           VBA3140
                                                                        VBA3150
        nskip  = 1                                                      VBA3160
        stend  = s(jstop)                                               VBA3170
        sbegin = s(jstart)                                              VBA3180
        su     = 0.                                                     VBA3190
        lode   = 0                                                      VBA3200
        g      = ginit                                                  VBA3210
        ncount = 0                                                      VBA3220
        ngu = 0                                                         VBA3230
        gu = gmax                                                       VBA3240
        ngv = 0                                                         VBA3250
        gv = gmin                                                       VBA3260
        nconv = 0                                                       VBA3270
        ns    = 0                                                       VBA3280
                                                                        VBA3290
!        write (nout,101) sbegin,stend,jtot,nskip                        VBA3300
! 101  format (/5x,"integrating ode over arclength from "                VBA3310
!     c ,e13.5," to ",e13.5                                              VBA3320
!     c ," divided into ",i5," intervals, skipping every "               VBA3330
!     c ,i5," intervals")                                                VBA3340
!        write (nout,110)                                                VBA3350
! 110  format(//5x,"ncount  growth rate    arclength",6x,"x1",11x,"x2"   VBA3360
!     c ,11x,"x3",11x,"x4",11x,"f2",7x,"lode",/)                         VBA3370
! 112  format(5x,i6,7e13.5,2x,i4)                                        VBA3380
                                                                        VBA3390
c..coarse tuning of iteration                                           VBA3400
                                                                        VBA3410
  10  continue                                                          VBA3420
        if (g .lt. gmin) go to 90                                       VBA3430
                                                                        VBA3440
        rgs     = 0.                                                    VBA3450
        if (lresis .eq. 1) rgs = sresis * eta / sqrt(g)                 VBA3460
                                                                        VBA3470
        ncount = ncount + 1                                             VBA3480
        call odeint (su,x1,x2,x3,x4,lode)                               VBA3490
        x2 = x2 / x1                                                    VBA3500
        x4 = x4 / x3                                                    VBA3510
        f2 = x2 - x4                                                    VBA3520
!        write (nout,112) ncount, g, su, x1, x2, x3, x4, f2, lode        VBA3530
        x1end = x2                                                      VBA3540
        x2end = f2                                                      VBA3550
                                                                        VBA3560
c..check convergence                                                    VBA3570
                                                                        VBA3580
        tol = abs (relerr * g) + abs (abserr)                           VBA3590
        if (abs (gu-g) .lt. tol) nconv = nconv + 1                      VBA3600
        if (abs (g-gv) .lt. tol) nconv = nconv + 1                      VBA3610
        if (nconv .ge. 2 .and. lode .ne. -1) go to 40                   VBA3620
        if (ncount .gt. itrmax) go to 92                                VBA3630
        if (lode) 20,14,12                                              VBA3640
                                                                        VBA3650
c  lode = +1,  ode solution exceedingly large, lower g and try again    VBA3660
                                                                        VBA3670
  12  continue                                                          VBA3680
        gu = g                                                          VBA3690
        ngu = ngu + 1                                                   VBA3700
        if (ngv .ge. 1) go to 30                                        VBA3710
        g = (g / gscal1) - dg1                                          VBA3720
        go to 10                                                        VBA3730
                                                                        VBA3740
c..medium tuning                                                        VBA3750
                                                                        VBA3760
c  lode = 0, ode solution never passed through zero, try again          VBA3770
                                                                        VBA3780
  14  continue                                                          VBA3790
        ns = ns + 1                                                     VBA3800
        if (f2 .lt. 0.) go to 16                                        VBA3810
c..f2 positive at the origin; lower g.                                  VBA3820
        ngu = ngu + 1                                                   VBA3830
        gu = g                                                          VBA3840
        if (ns .gt. 2) go to 35                                         VBA3850
        if (ngv .ge. 1) go to 30                                        VBA3860
        f1 = f2                                                         VBA3870
        g1 = g                                                          VBA3880
        g = (g / gscal2) - dg2                                          VBA3890
        go to 10                                                        VBA3900
                                                                        VBA3910
  16  continue                                                          VBA3920
c..f2 negative at the origin; increase g.                               VBA3930
        ngv = ngv + 1                                                   VBA3940
        gv = g                                                          VBA3950
        if (ns .gt. 2) go to 35                                         VBA3960
        if (ngu .ge. 1) go to 30                                        VBA3970
        f1 = f2                                                         VBA3980
        g1 = g                                                          VBA3990
        g = (g * gscal2) + dg2                                          VBA4000
        go to 10                                                        VBA4010
                                                                        VBA4020
c  lode = -1,  ode solution crossed zero,  raise g, try again           VBA4030
                                                                        VBA4040
  20  continue                                                          VBA4050
        gv = g                                                          VBA4060
        ngv = ngv + 1                                                   VBA4070
        if (ngu .ge. 1) go to 30                                        VBA4080
        g = (g * gscal2) + dg2                                          VBA4090
        go to 10                                                        VBA4100
                                                                        VBA4110
c..fine tuning, both gv and gu filled in, interpolate g                 VBA4120
c..bisection of growth rate method used                                 VBA4130
                                                                        VBA4140
  30  continue                                                          VBA4150
        g = gv + 0.5 * (gu - gv)                                        VBA4160
        go to 10                                                        VBA4170
c                                                                       VBA4180
c..converge growth rate using newton's method.                          VBA4190
c                                                                       VBA4200
  35  continue                                                          VBA4210
        slop = (g - g1) / (f2 - f1)                                     VBA4220
        g1   = g                                                        VBA4230
        g    = g - f2 * slop                                            VBA4240
        f1   = f2                                                       VBA4250
        if (g .lt. gv .or. g .gt. gu) go to 30                          VBA4260
        go to 10                                                        VBA4270
                                                                        VBA4280
                                                                        VBA4290
c..growth rate converged                                                VBA4300
                                                                        VBA4310
  40  continue                                                          VBA4320
!        write (nout,103) g                                              VBA4330
 103  format (/5x,"growth rate converged to ",e13.5)                    VBA4340
        ih = 1                                                          VBA4350
        if (lode .eq. -1) ih = -1                                       VBA4360
c                                                                       VBA4370
c..normalise the eigenfunction.                                         VBA4380
   43   x1 = 1. / x1                                                    VBA4390
        x3 = 1. / x3                                                    VBA4400
CDIR$ IVDEP                                                             VBA4410
        do 44 j=1,jh                                                    VBA4420
           y(jtot+1-j) = y(jtot+1-j) * x3                               VBA4430
           y(j)        = y(j)        * x1                               VBA4440
   44   end do     
        y(jh)       = 1.                                                VBA4450
                                                                        VBA4460
        return                                                          VBA4470
                                                                        VBA4480
  90  continue !write (nout,190) g,gmin                                           VBA4490
 190  format (/5x,"error: growth rate ",e13.5                           VBA4500
     c ," less than gmin ",e13.5,"; leaving subrtn solver")             VBA4510
        ih = 3                                                          VBA4520
        if (lode .eq. -1) ih = -3                                       VBA4530
        return                                                          VBA4540
                                                                        VBA4550
  92  continue !write (nout,192) ncount,itrmax                                    VBA4560
 192  format (/5x,"error: ncount ",i5," exceeding itrmax ",i5           VBA4570
     c ,"; too many iterations; leaving subrtn solver")                 VBA4580
        ih = 2                                                          VBA4590
        if (lode .eq. -1) ih = -2                                       VBA4600
        go to 43                                                        VBA4610
                                                                        VBA4620
!  94  write (nout,194)                                                  VBA4630
 194  format (/5x,"error: exceedingly large value of ode solution"      VBA4640
     c ," (lode=+1) during fine tuning part of iteration in solver")    VBA4650
        ih = 4                                                          VBA4660
                                                                        VBA4670
                                                                        VBA4680
        return                                                          VBA4690
        end subroutine solver                                           VBA4700
