!-----------------------------------------------------------------------VBA4710
        subroutine odeint (su,xu1,xu2,xu3,xu4,lode)                     VBA4720
!
      use vbal_param
      use cidata
      use crdata
      use ceigen
      use ccoef
!
!.. Implicits ..
      implicit none
!       include 'combal.inc'                                            VBA4730
!                                                                       VBA4740
!.. External Calls ..
      external coef
!
!..Scalars
      integer :: lode
      double precision :: su, xu1, xu2, xu3, xu4
!
!.. Local Scalars ..
      integer :: nt, n1, n
      double precision ::
     &        xl1, xl2, xl3, xl4,  st,  sv,  sw, cu1, cu2, cu3, cu4,
     &        cl1, cl2, cl3, cl4, au1, au2, au3, au4, al1, al2, al3,
     &        al4, dsh, recip
!..ode integrator using trapezoidal rule                                VBA4750
!       this particular version is set up to integrate                  VBA4760
!       simultaneously forward and backward towards zero                VBA4770
!       in the independent variable and to test the solution            VBA4780
!       at every point in the integration.                              VBA4790
!                                                                       VBA4800
!..initial values                                                       VBA4810
!                                                                       VBA4820
        xl1 = x0                                                        VBA4830
        xl2 = dlta                                                      VBA4840
        xl3 = x0                                                        VBA4850
        xl4 = - dlta                                                    VBA4860
        y(1)= xl1                                                       VBA4870
        y(jtot)= xl3                                                    VBA4880
        nt  = jstart                                                    VBA4890
        st  = sbegin                                                    VBA4900
        sv  = stend                                                     VBA4910
        cu1 = 0.                                                        VBA4920
        cu2 = 0.                                                        VBA4930
        cl1 = 0.                                                        VBA4940
        cl2 = 0.                                                        VBA4950
        cu3 = 0.                                                        VBA4960
        cu4 = 0.                                                        VBA4970
        cl3 = 0.                                                        VBA4980
        cl4 = 0.                                                        VBA4990
                                                                        VBA5000
        call coef (nt,st,cl1,cl2)                                       VBA5010
        nt  = jstop                                                     VBA5020
        call coef (nt,sv,cl3,cl4)                                       VBA5030
                                                                        VBA5040
c..step over unequally spaced intervals                                 VBA5050
                                                                        VBA5060
        n1 = 1 + nskip                                                  VBA5070
        do 10 n=n1,jh,nskip                                             VBA5080
!                                                                       VBA5090
           nt = jstart + n + nskip - n1                                 VBA5100
           call coef (nt,su,cu1,cu2)                                    VBA5110
           nt = jstop  - n - nskip + n1                                 VBA5120
           call coef (nt,sw,cu3,cu4)                                    VBA5130
                                                                        VBA5140
!  trapezoidal rule over one step                                       VBA5150
                                                                        VBA5160
           dsh = 0.5 * (su - st)                                        VBA5170
           al1 = dsh * cl1                                              VBA5180
           al2 = dsh * cl2                                              VBA5190
           au1 = dsh * cu1                                              VBA5200
           au2 = dsh * cu2                                              VBA5210
           recip = 1. / (1. - au1 * au2)                                VBA5220
!                                                                       VBA5230
           xu1 = (xl1*(1.+au1*al2) + xl2*(au1+al1)) * recip             VBA5240
           xu2 = (xl2*(1.+au2*al1) + xl1*(au2+al2)) * recip             VBA5250
!                                                                       VBA5260
           dsh = 0.5 * (sw - sv)                                        VBA5270
           al3 = dsh * cl3                                              VBA5280
           al4 = dsh * cl4                                              VBA5290
           au3 = dsh * cu3                                              VBA5300
           au4 = dsh * cu4                                              VBA5310
           recip = 1. / (1. - au3 * au4)                                VBA5320
!                                                                       VBA5330
           xu3 = (xl3*(1.+au3*al4) + xl4*(au3+al3)) * recip             VBA5340
           xu4 = (xl4*(1.+au4*al3) + xl3*(au4+al4)) * recip             VBA5350
                                                                        VBA5360
             if (ldiag(1) .eq. 2) then                                  VBA5370
               y(jtot+1-n)= xu3                                         VBA5380
               y(n)       = xu1                                         VBA5390
             end if                                                     VBA5400
                                                                        VBA5410
!  test for  zero crossing or exceedingly large values                  VBA5420
                                                                        VBA5430
             if (xu1 .gt. xmax .or. xu3 .gt. xmax) then                 VBA5440
!                                                                       VBA5730
!  exceedingly large value encountered                                  VBA5740
             lode = + 1                                                 VBA5770
             return                                                     VBA5780
             end if
!
             if (xu1 .lt. 0.0 .or. xu3 .lt. 0.) then                    VBA5450
!  zero crossing                                                        VBA5640
             lode = - 1                                                 VBA5670
             if (xu1 - xl1 .ne. 0.)                                     VBA5680
     &               su = (st*xu1 - su*xl1) / (xu1 - xl1)               VBA5690
             if (xu3 - xl3 .ne. 0.)                                     VBA5700
     &               sw = (sv*xu3 - sw*xl3) / (xu3 - xl3)               VBA5710
             return                                                     VBA5720
             end if
!                                                                       VBA5460
          st  = su                                                      VBA5470
          cl1 = cu1                                                     VBA5480
          cl2 = cu2                                                     VBA5490
          xl1 = xu1                                                     VBA5500
          xl2 = xu2                                                     VBA5510
          sv  = sw                                                      VBA5520
          cl3 = cu3                                                     VBA5530
          cl4 = cu4                                                     VBA5540
          xl3 = xu3                                                     VBA5550
          xl4 = xu4                                                     VBA5560
  10  end do                                                            VBA5570
                                                                        VBA5580
!  normal ending -- no zero crossing                                    VBA5590
                                                                        VBA5600
        lode = 0                                                        VBA5610
        return                                                          VBA5620
                                                                        VBA5630
        end subroutine odeint                                           VBA5790
