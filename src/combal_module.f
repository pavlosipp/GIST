!CIPP    CLICHE  COMBAL                                                 COM0010
        module vbal_param                                               COM0020
        integer, parameter :: NTMP=1151
        integer, parameter :: NTURN=15
        integer, parameter :: MXX=2000*NTURN+2       
        end module vbal_param                                           COM0050
!
        module cidata
        integer :: ntmax,nout   ,ntype  ,mpnt   ,modes                  COM0060
     &  ,mercie ,lxcont ,itrmax ,nmod   ,ih     ,nstart ,nstop          COM0070
     &  ,nskip  ,ncount ,lmodel ,jstart ,jstop  ,jtot   ,jh             COM0080
        integer, dimension(16) :: ldiag
!       common  /cidata/                                                COM0090
!    &   ntmax  ,nout   ,ntype  ,mpnt   ,modes  ,ldiag                  COM0100
!    &  ,mercie ,lxcont ,itrmax ,nmod   ,ih     ,nstart ,nstop          COM0110
!    &  ,nskip  ,ncount ,lmodel ,jstart ,jstop  ,jtot   ,jh             COM0120
        end module cidata
                                                                        COM0130
        module crdata
        real ::  dlta   ,x0    ,sk     ,alf    ,spk    ,s0in     ,ci     COM0140
     &  ,rs     ,twopi  ,twopin ,qs     ,eta    ,b0r0   ,vp     ,pp     COM0150
     &  ,qppsip ,sqgbsq ,ftp    ,pth,   qprime
!        common  /crdata/                                               COM0170
!     &           dlta   ,x0     ,sk     ,alf    ,spk    ,s0     ,ci    COM0180
!     &  ,rs     ,twopi  ,twopin ,qs     ,eta    ,b0r0   ,vp     ,pp    COM0190
!     &  ,qppsip ,sqgbsq ,ftp    ,pth                                   COM0200
        end module crdata
!                                                                       COM0210
        module cmodes
        integer, allocatable ::       mb(:)    ,nb(:)                   COM0220
!        common /cmodes/  mb          ,nb                               COM0230
        end module cmodes
                                                                        COM0240
        module cfoura
        real, allocatable ::  a0mn(:)  ,a1mn(:)  ,a2mn(:)  ,rgmn(:)     COM0250
     &                       ,d0mn(:)  ,d1mn(:)  ,d2mn(:), gssmn(:)     COM0260
     &                       ,gstmn(:) ,bsmn(:)
!        common  /cfoura/ a0mn  ,a1mn  ,a2mn  ,rgmn  ,d0mn  ,d1mn  ,d2mnCOM0270
        end module cfoura
                                                                        COM0280
        module ceigen
        real ::  g      ,gmax   ,gmin   ,ginit  ,dg1    ,dg2    ,gscal1 COM0290
     &  ,gscal2 ,xmax   ,relerr ,abserr ,stend  ,sbegin ,x1end  ,x2end  COM0300
!        common  /ceigen/                                               COM0310
!     &           g      ,gmax   ,gmin   ,ginit  ,dg1    ,dg2    ,gscal1COM0320
!     &  ,gscal2 ,xmax   ,relerr ,abserr ,stend  ,sbegin ,x1end  ,x2end COM0330
        end module ceigen
                                                                        COM0340
        module ccoef
        real, allocatable :: s(:),c0(:),c1(:),c2(:),fs(:),fq(:),y(:)
     &      ,gss(:),gst(:),pax1(:),bs(:),c3(:),c111(:),cbel(:)
!        common  /ccoef/ s     ,c0    ,c1    ,c2    ,fs    ,fq    ,y    COM0360
        end module ccoef
                                                                        COM0370
        module cresis
        integer :: lresis                                               COM0380
        real    :: sresis ,rgs                                          COM0390
!        common  /cresis/                                               COM0400
!     &     lresis ,sresis ,rgs                                         COM0410
        end module cresis
!                                                                       COM0420
!CIPP    ENDCLICHE                                                      COM0430
