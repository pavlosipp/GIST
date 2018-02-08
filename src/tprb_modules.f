c--------0---------0---------0---------0---------0---------0---------0-c
C  18.11.05       LAST MODIFICATION 18.11.05        WAC: tprb_modules  D
c--------0---------0---------0---------0---------0---------0---------0-c
!
!
      module tprb_param
!
!      integer, parameter :: NI=196, MLMNV=238
      integer, parameter :: NJ=445
      integer, parameter :: NK=65
      integer, parameter :: NJK=NJ*NK
      integer, parameter :: MLMNB=527
      integer, parameter :: LSSL=3184
      integer, parameter :: MLMNL=526
      integer, parameter :: IVAC=0
      integer, parameter :: NPITCH=100
!
      end module tprb_param
!

!
C..   ALLOCATE BLOCKS
!
      module PHYSEQ
      use tprb_param
C...FORMERLY Variables in Common Block /PHYSEQ/ ...
      integer :: NPER,nprocs,mms,nsmin,nsmax,mm,nmin,nmax ,nbalmn, lmnb0
      real :: GAM,PI,TPI,EPS,EL,TPONP,parfac,pvac,qvac,dsvac,qn
     &       ,rplmin,xplo,dsta,dth,dph,dtdp,      alpha,djp
!      COMMON /PHYSEQ/ GAM, NPER, PI, TPI, EPS, EL, TPONP, parfac, nprocs
!     &               ,pvac, qvac, dsvac,  mms, nsmin, nsmax, qn, nowall
!     &               ,awall, ewall, dwall, gwall, drwal, dzwal, npwall
!     &               ,mm, nmin, nmax, rplmin, xplo, dsta, dth, dph, dtdp
!     &               ,nbalmn  ,lmnb0,alpha,    djp
      end module PHYSEQ
C
      module LAMBDA
      use tprb_param
!C...FORMERLY Variables in Common Block /LAMBDA/
      real, allocatable ::
     &     VL(:,:) , vlt(:  ), vlp(:  ), PSIVT(:)
     &    ,PSIVP(:), PGPPV(:), PGTTV(:), PGTPV(:)
!      COMMON /LAMBDA/ VL,vlt,vlp,PSIVT,PSIVP,PGPPV,PGTTV,PGTPV
      end module LAMBDA
C
      module COLAMN
      use tprb_param
C...FORMERLY Variables in Common Block /COLAMN/
      real, allocatable ::
     &     AL(:,:   ), gla(:,:  )
     &    ,dna(:,:,:), dnb(:,:,:)
!      COMMON /COLAMN/ AL,gla,dna,dnb
      end module COLAMN
C
      module VMCFOU
      use tprb_param
C...FORMERLY Variables in Common Block /VMCFOU/
      real, allocatable ::
     &    FRV(:,: ), FZV(:,:  ), FVL(:,:)
     &   ,FVLI(:,:), FVJAC(:,:), FVGAM(: )
     &   ,FVALF(: ), FVB(:    ), fvq(:   )
!      COMMON /VMCFOU/ FRV,FZV,FVL,FVLI,FVJAC,FVGAM,FVALF,FVB,fvq
      end module VMCFOU
C
      module TRIGGS
      use tprb_param
C...FORMERLY Variables in Common Block /TRIGGS/
      integer LMNV, LMNL, LMNB, lss
      integer, allocatable ::
     &       MX(:)   ,NX(:)       , ML(:)    , NL(:) , MB(:)
     &      ,NB(:)   , lfrz(:,:  ), lfx(:,: )
     &      ,lsx(:,:), msx(:)     , nsx(:)
      real, allocatable ::
     &     TCOS(:,:) ,TSIN(:,:)   , tsss(:,:)
     &    ,POL(:   ) ,TOR(:      ), TSC(:,:)
c     &              ,TSS(:S,MLMNS),TSC(:S,MLMNS)
chsx     &              ,TSS(200,njk),TSC(2094,njk)
cw7as     &              ,TSS(200,njk),TSC(2691,njk)
cmhh     &              ,TSS(200,njk),TSC(2976,njk)
ch1     &              ,TSS(200,njk),TSC(2463,njk)
clhd     &              ,TSS(200,njk),TSC(1804,njk)
!      COMMON /TRIGGS/TCOS,TSIN,tsss,POL,TOR,TSC,MX,NX,ML,NL,MB,NB
!     &              ,LMNV,LMNL,LMNB,lfrz,lfx,lsx,msx,nsx,lss
      end module TRIGGS
C
      module VMCRSP
      use tprb_param
C...FORMERLY Variables in Common Block /VMCRSP/
      real, allocatable ::
     &     RV(:,:  ), ZV(:,:  ), RVT(:,:), ZVT(:,: )
     &    ,RVP(:,: ), ZVP(:,: ), RVQ(:,:), RVQT(:,:)
     &    ,GTTV(:,:), GTPV(:,:),GPPV(:,:), VJAC(:,:)
     &    ,VBT(:,: ), VBP(:,: ),VBSQ(:,:)
!      COMMON /VMCRSP/ RV,ZV,RVT,ZVT,RVP,ZVP,RVQ,RVQT,GTTV,GTPV,GPPV,VJAC
!     &               ,VBT,VBP,VBSQ
      end module VMCRSP
C
      module MAPPIN
      use tprb_param
C...FORMERLY Variables in Common Block /MAPPIN/
      real, allocatable ::
     &     VQ(:   ), sqbbsq(:), VALF(: ), VGAM(:), THV(:)
     &    ,VALFT(:), VGAMT(: ),VJTOBJ(:)
!      COMMON /MAPPIN/ VQ,sqbbsq,VALF,VGAM,THV,VALFT,VGAMT,VJTOBJ
      end module MAPPIN
C
      module BOOFOU
      use tprb_param
C...FORMERLY Variables in Common Block /BOOFOU/
      real, allocatable ::
     &     FR(:,:    ) , FZ(:,:  ), FPHV(:,:   )
     &   , FRI(:,: )   ,FZI(:,:)  , FPHVI(:,:  )
     &   ,FPHVT(:, : ) ,FPHVP(:,:), fbsq(:,:   )
     &   ,FBJAC(:,:)   ,FBS(:,:)
!      COMMON /BOOFOU/ FR,FZ,FPHV,FRI,FZI,FPHVI,FPHVT,FPHVP,fbsq,FBJAC
!     &               ,FBS
      end module BOOFOU
C
      module BOORSP
      use tprb_param
C...FORMERLY Variables in Common Block /BOORSP/
      real, allocatable ::
     &      R (:,:)  , Z (:,:   )  , phv (:,:) ,  RI (:,:)
     &    ,ZI (:,:)  , phvi (:,:)  , RT (:,: ) ,  ZT (:,:  )
     &    ,rs (:,:  ), RP (:,:    ), ZP (:,: ) ,  zs (:,:  )
     &    ,RSq(:,:)  , phvp(:     ), phvt(:   )
!      COMMON /BOORSP/ R,Z,phv,RI,ZI,phvi,RT,ZT,rs,RP,ZP,zs,RSq,phvp,phvt
      end module BOORSP
C
      module COMERC
      use tprb_param
C...FORMERLY Variables in Common Block /COMERC/
      real, allocatable :: 
     &     BST(:)    , baldp(:)    , curvds(:) , shearl(:)
!      COMMON /COMERC/ BST, baldp, curvds, shearl
      end module COMERC
C
      module MERFOU
      use tprb_param
C...FORMERLY Variables in Common Block /MERFOU/
      real, allocatable ::
     &     fcurvs(:) , flshea(:)   , fjpar(:)  , fgssu(:)
!      COMMON /MERFOU/ fcurvs,flshea,fjpar,fgssu
      end module MERFOU
C
      module BOOMET
      use tprb_param
C...FORMERLY Variables in Common Block /BOOMET/
      real, allocatable ::
     &     BJAC(:,:) ,BJACS(:,:)   , GTTL(:,:) ,GTPL(:,:)
     &    ,GPPL(:,:) ,GSSL(:,: )   , GSTL(:,:) ,GSsu(:,:)
!      COMMON /BOOMET/ BJAC,BJACS,GTTL,GTPL,GPPL,GSSL,GSTL,GSSU
      end module BOOMET
C
      module BDERIV
      use tprb_param
C...FORMERLY Variables in Common Block /BDERIV/
      real, allocatable ::
     &     BT(:,:)   ,   BP(:,:)   ,  BSQ(:,:) ,    BS(:,:)
!      COMMON /BDERIV/ BT,   BP,  BSQ,  BS
      end module BDERIV
C
      module RADFCT
      use tprb_param
C...FORMERLY Variables in Common Block /RADFCT/
      real, allocatable ::
     &                   AM(  :),  PTH(  :),  PVP(  :),PVPI(  :)
     &              ,   FPP(:  ),  FTP(:  ),  PPI(  :),  PP(  :)
     &              ,  FPPP(:  ), FTPP(:  ),AIOTA(  :), VVP(  :)
     &              ,    f0(  :),   f1(  :),   f2(  :),  f3(  :)
     &              , fpppp(  :),ftppp(  :),WMAGV(  :),WMAG(  :)
     &              ,   CIV(  :),  CJV(  :), CIVP(  :),CJVP(  :)
     &              ,   CI (  :),  CJ (  :), CIP (  :),CJP (  :)
     &              ,  CIPI(  :), CJPI(  :), ppp (  :),glm (:,:)
     &              ,    VP ( :),  VPP(  :),EQUIV(  :),EQUI(  :)
     &              ,    s (:)
!      COMMON /RADFCT/ AM,PTH,PVP,PVPI,FPP,FTP,PPI,PP,FPPP,FTPP,AIOTA
!     &               ,VVP,f0,f1,f2,f3,fpppp,ftppp,WMAGV,WMAG,CIV,CJV
!     &               ,CIVP,CJVP,CI,CJ,CIP,CJP,CIPI,CJPI,ppp,glm,VP,VPP
!     &               ,EQUIV,EQUI,s
      end module RADFCT
C
      module COBALL
      use tprb_param
C...FORMERLY Variables in Common Block /COBALL/
       real, allocatable ::
     &                fbaldp(:), fbalds(:), fbalcp(:)
     &              , fbalcs(:), fbalcq(:)
!     COMMON /COBALL/ fbaldp, fbalds, fbalcp, fbalcs, fbalcq
      end module COBALL
C
      module COBOOT
      use tprb_param
C...FORMERLY Variables in Common Block /COBOOT/
      integer, allocatable ::  jkmax(:)
      real, allocatable ::
     &        pitch(:)  ,bnorm(:, :)  ,g2(:)      ,g4(:)     ,gb(:)
     &       ,g1mnin(:) ,g4mn(:)      ,rl31(:)    ,rl32e(:)  ,gc(:)
     &       ,boot(:)   ,bsqav(:)     ,bsqmax(:)  ,rl32i(:)  ,ft(:)
     &       ,g2av(:)   ,tempp(:)     ,dboot(:)   ,bjav(:)   ,fc(:)
     &       ,g1av(:,:)               ,g4av(:,:)
     &       ,rlampl(:) ,gbplat(:)    ,rnue(:)    ,bgradgmn(:)
     &       ,bdgradb(:),bdgradg(:)   ,bdglnb(:)  ,bgradbmn(:)

!      COMMON /COBOOT/ pitch,bnorm,g2,g4,gb,g1mnin,g4mn,rl31,rl32e,rl32i
!     &               ,gc,boot,bsqav,bsqmax,ft,fc,g2av,tempp,dboot,bjav
!     &               ,g1av,g4av,jkmax,gbplat,rlampl,rnue,bdgradg,bdgradb
!     &               ,bdglnb,bgradbmn,bgradgmn
      end module COBOOT
C
      module PAPRPL
C...FORMERLY Variables in Common Block /PAPRPL/
C     PAPRPL: PRINT & PLOT SWITCHES
C
      integer ::
     &            LLAMPR,     LVMTPR,     LMETPR,     LRIPPL
     &          , LRADMN,     LRADMX,     LNBALN,     LNBALX
     &          , LXYZPR,     LIOTPL,     LSKEWN,     LMXBAL
     &          , LCURRF,     LMESHP,     LMESHV,     LITERS
     &          , LXYZPL,     LEFPLS,     LBOOTJ,     LPRESS
!      COMMON /PAPRPL/
!     &            LLAMPR,     LVMTPR,     LMETPR,     LRIPPL
!     &          , LRADMN,     LRADMX,     LNBALN,     LNBALX
!     &          , LXYZPR,     LIOTPL,     LSKEWN,     LMXBAL
!     &          , LCURRF,     LMESHP,     LMESHV,     LITERS
!     &          , LXYZPL,     LEFPLS,     LBOOTJ,     LPRESS
      end module PAPRPL
