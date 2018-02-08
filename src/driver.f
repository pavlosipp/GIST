!-----------------------------------------------------------------------VBA1720
        subroutine driver                                               VBA1730
!
       use vbal_param
       use cidata
       use crdata
       use ceigen
       use cresis
       use cmodes
       use cfoura
       use ccoef
       use specifications_mod, ONLY:vmec_file,table_tag
!
!       include 'combal.inc'                                            VBA1740
C.. Implicits ..
!      implicit none
!
! ... LOCAL SCALARS
      integer :: mn, nscr, kf, n1, n2, n, jmax, l, i, ibg
      real(4) :: tim
!
!.. External Calls ..
      external equin, prepar, solver
!.. External Functions ..
      integer :: isamax
      external  isamax
!
!.. Intrinsic Functions ..
      intrinsic MIN
!                                                                       VBA1750
!  this subroutine calls subrtns to                                     VBA1760
!  sweep through a sequence of flux surfaces to test stability          VBA1770
!  iterate on the growth rate or eigenvalue until convergence is obtaineVBA1780
!  check for errors and react accordingly                               VBA1790
!x                                                                      VBA1800
      tim=0.0
!--READ-IN EQUILIBRIUM MODES.    
      OPEN(25,file='../int/tpr_'//trim(vmec_file)//'_'//trim(table_tag),
     &POSITION="REWIND",STATUS="OLD",ACTION="READWRITE", ERR=216)
        read (25,410) kf, mpnt, b0r0, vp 
        read (25,413) (mb(mn),mn=1,mpnt),(nb(mn),mn=1,mpnt)             VBA1830
  410   format (2i5,1p2e15.6)                                           VBA1840
  413   format (12i6)    
       
        if (lmodel .eq. 0) write (nout,103)                             VBA1870
  103   format(//,8x,"ideal mhd marginal point analysis$")              VBA1880
!                                                                       VBA1890
!..plots of eigenfunctions along arclength                              VBA1900
!                                                                       VBA1910
!..resistive coefficients                                               VBA1920
!                                                                       VBA1930
        if (lresis .eq. 1) then                                         VBA1940
        if (gmin .le. 0.) gmin = .000001                                VBA1950
        if (ginit .le. 0.)  ginit = 11. * gmin                          VBA1960
!                                                                       VBA1970
!        write (nout,104)                                                VBA1980
!  104   format (//,8x,"resistive form of baloon code (lmodel = 1)$")    VBA1990
!        write (nout,105) sresis                                         VBA2000
!  105   format (8x,e13.5,"s = magnetic reynolds number$")               VBA2010
!        write (nout,106) gmin                                           VBA2020
!  106   format (8x,e13.5,"gmin = min value allowed for eigenfn$")       VBA2030
        end if                                                          VBA2040
!                                                                       VBA2050
!..choose contours                                                      VBA2060
!                                                                       VBA2070
        nscr=6                                                          VBA2080
        modes=MIN(modes,mpnt)                                           VBA2090
!pax        write (ntype,124) modes,mpnt,jtot,sk,alf,spk                    VBA2100
      
!        write (nscr,124) modes,mpnt,jtot,sk,alf,spk                     VBA2110
  124   format(3x,"modes=",i4,"  mpnt=",i4," jtot=",i6,"  sk=",f7.4,    VBA2120
     >         "  alfa=",f7.4," spk=",f7.3)
!pax        write (ntype,125) vp                                            VBA2140
 125    format (1p1e13.5)                                               VBA2150
!pax        write (ntype,120)                                               VBA2160
! 120  format (7x,"rho",9x,"p-prime",6x,"v-prime",9x,"q",6x,"nps",5x,    VBA2170
!     c"sbegin",7x,"s--max",4x,"d(ln y)/ds",4x,"diffr-log"               VBA2180
!     c,4x,"eigenvalue",1x,"ih")                                         VBA2190
!        write (nscr,123)                                                VBA2200
 123  format (7x,"rho",6x,"p-prime",7x,"q",7x,"s--max",4x,              VBA2210
     c"diffr-log",1x,"eigenvalue",1x,"ih")                              VBA2220
!
!        write (24,128)
! 128    format(2x,"n",2x,"mb",2x,"nb",5x,"a0mn",8x,"a1mn",8x,"a2mn",8x
!     >        ,"d0mn",8x,"d1mn")
!                                                                       VBA2230
        n1 = 1                                                          VBA2240
        n2 = kf                                                         VBA2250
        if (lxcont .ne. 0) n1 = lxcont                                  VBA2260
        if (lxcont .ne. 0) n2 = lxcont                                  VBA2270
        if (lxcont .eq. -1) n1=nstart                                   VBA2280
        if (lxcont .eq. -1) n2=nstop                                    VBA2290
!                                                                       VBA2300
!       write (14,133) n2                                               VBA2310
        do 5 n=1,n1-1                                                   VBA2320
           call equin                                                   VBA2330
    5   end do   
                                                                        VBA2340
        do 10 n=n1,n2                                                   VBA2350
                                                                        VBA2360
          call equin                                                    VBA2450
!          write (nout,100)                                              VBA2370
!       rs = (n-0.5)/real(kf)                                           VBA2380
!          write (nout,101) rs                                           VBA2390
!  100     format (//                                                    VBA2400
!     &    "-----------------------------------------------------------" VBA2410
!     &    ,"----------------------------------------------------------" VBA2420
!     &    //)                                                           VBA2430
 101      format (5x,"new contour with rho = ",e13.5)                   VBA2440
           call prepar    
           call solver                   
          CLOSE(25)
!                                                                       VBA2480
!--output for post-processed graphics                                   VBA2490
!                                                                       VBA2500
           jmax = jstart - 1 + isamax(jtot,y,nskip)


!            do 33 l=1,modes
!             write(24,129) n,mb(l),nb(l),a0mn(l),a1mn(l),a2mn(l),d0mn(l)
!     &                                  ,d1mn(l)
!   33       end do  
!          if (lxcont .gt. 0) then
!       write (14,134) jtot,rs,pp,qs,g,jmax,s(jmax),x2end,alf,sk        VBA2510
!       write (14,135) (s(i),i=jstart,jstop,nskip)                      VBA2520
!       write (14,135) (y(i),i=1,jtot,nskip)                            VBA2530
!           do 34 i=1,jtot,nskip
!              ibg=i+jstart-1
!pax               write(14,135) s(ibg),y(i),c1(ibg),1./c2(ibg),c2(ibg)
!   34      end do
!          endif
!                                                                       VBA2540
!..intermediate output                                                  VBA2550
!                                                                       VBA2560
!           write(ntype,121)n,rs,pp,vp,qs,jtot,sbegin,s(jmax),x1end,x2endVBA2570
!     &                    ,g,ih                                         VBA2580
!      !     write(nscr,122)n,rs,pp,qs,s(jmax),x2end,g,ih   
                                                                        VBA2600
  10  end do                                                            VBA2610
!                                                                       VBA2620
!        write (nout,100)                                                VBA2630
      call second(tim)
!      write(nout,136) tim
!      write(ntype,136) tim
!      write(nscr,136) tim
 121  format (i2,1p4e13.5,i5,1p5e13.5,1x,i2)                            VBA2640
 122  format (i2,1p6e11.3,1x,i2)                                        VBA2650
 129  format (3i4,1p5e12.4)
 133     format(i2)
 134  format (i6,1p4e12.4,/i6,1p4e12.4)
 135  format (1p5e14.6)
 136  format("total running time =",1p1e15.6," seconds")
!                        
        return
 216    WRITE(*,"(2X,A,/)") '-------- The equilibrium has not been yet
     . initialized. Try initialize=.t. -------' 
       STOP 
        end subroutine driver  
