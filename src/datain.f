!-----------------------------------------------------------------------VBA0670
        subroutine datain                                               VBA0680
!
      use vbal_param
      use cidata
      use crdata
      use ceigen
      use cresis
      use specifications_mod, only:surf,nz0,pol_turns,theta0
      use specifications_mod, only:alpha0,gs2
!
!.. Implicits ..
      implicit none
!
!       include 'combal.inc'                                            VBA0690
!
C.. Local Scalars ..
      integer :: n



C.. Intrinsic Functions ..
      intrinsic atan
!
!        namelist        /baloin/                                        VBA0710
!     &   pol_turns     ,theta_k    ,itrmax ,nz0  ,ci     ,mercie                 VBA0720
!     &  ,xmax   ,abserr ,relerr ,gmin   ,gmax   ,ginit                  VBA0730
!     &  ,dg1    ,dg2    ,gscal1 ,gscal2 ,jstart ,jstop                  VBA0740
!     &  ,lresis ,sresis ,dlta   ,sk     ,x0     ,alfa                    VBA0750
!     &  ,modes  ,surf ,lmodel ,ldiag  ,nstart ,nstop                    VBA0760
                                                                        VBA0770

!pax: Parameters which normally belong to the namelist baloin 
!     and usually do not change 
!                                                                       VBA0700

    

      LDIAG(1)=2
      LMODEL=2 !0
      MODES=527 !137
      NSTART=1 
      NSTOP=94 !49
      JSTART=1
      LRESIS=0 
      GINIT=0.0000
      DG2=0.08
      GMAX=0.995
      SK=0.0
      DLTA=1.E-250
      XMAX=1.E250

 
!------------------------------------------------------------------
                                                                        VBA0780
        nout = 16                                                       VBA0790
        ntype = 59                                                      VBA0800
        twopi = 8. * atan(1.)                                           VBA0810
        twopin = 1. / twopi                                             VBA0820
!        write (nout,101)                                                VBA0830
 101  format(5x,"CRPP Ballooning Code for 3-D Equilibria by  W. Anthony VBA0840
     & Cooper (October 1986)")                                          VBA0850
!        write (nout,102)                                                VBA0860
 102  format(5x,"shooting from plus and minus infinity to the origin"/  VBA0870
     & ,/)                                                              VBA0880
                                                                        VBA0890
                                                                        VBA0900
!..initialize data                                                      VBA0910
                                                                        VBA0920
        do 2 n=1,16                                                     VBA0930
            ldiag(n) = 0                                                VBA0940
 2      end do                                                          VBA0950
        lmodel = 1                                                      VBA0960
        lxcont = 0                                                      VBA0970
                                                                        VBA0980
        itrmax = 30                                                     VBA0990
        ntmax  = 131                                                    VBA1000
        mercie = 1                                                      VBA1010
                                                                        VBA1020
        s0in     = 13.0                                                   VBA1030
        spk    = -0.                                                    VBA1040
        alf    = 0.                                                     VBA1050
        jstart = 1                                                      VBA1060
        jstop  = 131                                                    VBA1070
        nstart = 2                                                      VBA1080
        nstop  = 10                                                     VBA1090
        xmax   = 1.e90                                                  VBA1100
        abserr = .0                                                     VBA1110
        relerr = .001                                                   VBA1120
        gmin   = -100.                                                  VBA1130
        gmax   =  100.                                                  VBA1140
        ginit  = 0.                                                     VBA1150
        dg1    = .008                                                   VBA1160
        dg2    = .003                                                   VBA1170
        gscal1 = 1.                                                     VBA1180
        gscal2 = 1.                                                     VBA1190
        x0     = 0.                                                     VBA1200
        dlta   = 1.e-90                                                 VBA1210
                                                                        VBA1220
        lresis = 0                                                      VBA1230
        sresis = .0001                                                  VBA1240
        eta   = 1.                                                      VBA1250
        modes  = ntmp                                                   VBA1260
        sk     = 0.                                                     VBA1270
        ci     = 1.                                                     VBA1280
!                                                                       VBA1290
!  read in data                                                         VBA1300
!                                                                       VBA1310
!       input data baloin,5,6     
!        OPEN(15, file='../inp/trabal.inp', action="read")
!        read  (15,baloin)   
!        CLOSE(15)

!pax:set parameters

        lxcont=surf
        if (gs2) then
        ntmax=2*nz0+1
        else
        ntmax=nz0+1 !add extra point for dBdz
        endif
        jstop=ntmax
        spk=theta0
        s0in=pol_turns
        alf=alpha0 !alf=zeta - q (theta-theta_k)
                                                                        VBA1340
!        write (nout,baloin)                                             VBA1350
!        write (ntype,103) ntmax, lxcont, s0, ginit, dg1, dg2,           VBA1360
!     &  ci, alf, jstart, jstop, sk, sresis, lresis, lmodel, spk, dlta,  VBA1370
!     &  itrmax,x0                                                       VBA1380
  103 format(2x,"ntmax=",i5,2x,"lxcont=",i3,2x                          VBA1390
     &,"s0=",f6.1,2x,"ginit=",e15.5,2x,"dg1=",e15.5,2x,"dg2=",e15.5,    VBA1400
     &2x,"ci=",f9.4,2x,"alf=",f9.4,2x,"jstart=",i5,/2x,"jstop=",i5,     VBA1410
     &2x,"sk=",f8.4,2x,"sresis=",e15.5,2x,"lresis=",i1,2x,"lmodel=",i1, VBA1420
     &2x,"spk=",f9.4,2x,"dlta=",e10.3,2x,"itrmax=",i2,2x,"x0=",e10.3)   VBA1430
        jtot = 1 + jstop - jstart                                       VBA1440
        jh   = (jtot + 1) / 2                                           VBA1450
                                                                        VBA1460
        return                                                          VBA1470
        end subroutine datain                                           VBA1480
