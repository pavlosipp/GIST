C-----------------------------------------------------------------------
C  18.11.05        LAST MODIFICATION 25.11.05   UHS: tprall_bal DATA   D
C-----------------------------------------------------------------------
C
      subroutine TPRALL_BAL
C
C     FORMERLY  TPRCOM DATA D
C     COMMON BLOCKS FOR TERPSICHORE BALLOONING --LAST MODIFIED 20.11.05
c--------0---------0---------0---------0---------0---------0---------0-c
!
!
      use tprb_param
!
      use specifications_mod, ONLY:ns_vmec,mpnt_vmec
      use physeq
      use radfct
      use triggs
      use paprpl
      use vmcfou
      use vmcrsp
      use colamn
      use boofou
      use boorsp
      use lambda
      use mappin
      use boomet
      use bderiv
      use comerc
      use merfou
      use coball
      use coboot
!
C..   ALLOCATE BLOCKS  

       NI=ns_vmec-1 
       MLMNV=mpnt_vmec

      
!
!C...FORMERLY Variables in Common Block /LAMBDA/
      allocate
     &    (VL(NJK,NI), vlt(NJK  ), vlp(NJK  ), PSIVT(NJK)
     &    ,PSIVP(NJK), PGPPV(NJK), PGTTV(NJK), PGTPV(NJK),stat=ierr1)
!      COMMON /LAMBDA/ VL,vlt,vlp,PSIVT,PSIVP,PGPPV,PGTTV,PGTPV
C
C...FORMERLY Variables in Common Block /COLAMN/
      allocate
     &    (AL(mlmnl,mlmnl+3 ), gla(3,lssl       )
     &    ,dna(mlmnl,2,mlmnl), dnb(mlmnl,2,mlmnl),stat=ierr2)
!      COMMON /COLAMN/ AL,gla,dna,dnb
C
C...FORMERLY Variables in Common Block /VMCFOU/
      allocate
     &   (FRV(MLMNV,0:NI ), FZV(MLMNV,  0:NI), FVL(MLMNL,NI)
     &   ,FVLI(MLMNL,0:NI), FVJAC(MLMNB,0:NI), FVGAM(MLMNB )
     &   ,FVALF(MLMNB    ), FVB(MLMNB       ), fvq(mlmnb   ),stat=ierr3)
!      COMMON /VMCFOU/ FRV,FZV,FVL,FVLI,FVJAC,FVGAM,FVALF,FVB,fvq
C
C...FORMERLY Variables in Common Block /TRIGGS/
!      integer LMNV, LMNL, LMNB, lss
      allocate
     &      (MX(MLMNV),NX(MLMNV), ML(MLMNB), NL(MLMNB) , MB(MLMNB)
     &      ,NB(MLMNB), lfrz(0:36,-36:36  ), lfx(-36:72,-72:72 )
     &      ,lsx(-36:72,-72:72)            , msx(lssl) , nsx(lssl)
     &      ,stat=ierr4)
      allocate
     &    (TCOS(NJK,MLMNB),TSIN(NJK,MLMNB), tsss(njk,mlmnb-1)
     &    ,POL(NJK       ),TOR(NJK       ), TSC(lssl,njk),stat=ierr5)
c     &              ,TSS(NJKS,MLMNS),TSC(NJKS,MLMNS)
chsx     &              ,TSS(200,njk),TSC(2094,njk)
cw7as     &              ,TSS(200,njk),TSC(2691,njk)
cmhh     &              ,TSS(200,njk),TSC(2976,njk)
ch1     &              ,TSS(200,njk),TSC(2463,njk)
clhd     &              ,TSS(200,njk),TSC(1804,njk)
!      COMMON /TRIGGS/TCOS,TSIN,tsss,POL,TOR,TSC,MX,NX,ML,NL,MB,NB
!     &              ,LMNV,LMNL,LMNB,lfrz,lfx,lsx,msx,nsx,lss
C
C...FORMERLY Variables in Common Block /VMCRSP/
      allocate
     &    (RV(NJK,0:NI  ), ZV(NJK,0:NI  ), RVT(NJK,0:NI), ZVT(NJK,0:NI )
     &    ,RVP(NJK,0:NI ), ZVP(NJK,0:NI ), RVQ(NJK,0:NI), RVQT(NJK,0:NI)
     &    ,GTTV(NJK,0:NI), GTPV(NJK,0:NI),GPPV(NJK,0:NI), VJAC(NJK,0:NI)
     &    ,VBT(NJK,0:NI ), VBP(NJK,0:NI ),VBSQ(NJK,0:NI), stat=ierr6   )
!      COMMON /VMCRSP/ RV,ZV,RVT,ZVT,RVP,ZVP,RVQ,RVQT,GTTV,GTPV,GPPV,VJAC
!     &               ,VBT,VBP,VBSQ
C
C...FORMERLY Variables in Common Block /MAPPIN/
      allocate
     &    (VQ(NJK   ), sqbbsq(ni), VALF(NJK ), VGAM(NJK), THV(NJK)
     &    ,VALFT(NJK), VGAMT(NJK),VJTOBJ(NJK), stat=ierr7)
!      COMMON /MAPPIN/ VQ,sqbbsq,VALF,VGAM,THV,VALFT,VGAMT,VJTOBJ
C
C...FORMERLY Variables in Common Block /BOOFOU/
      allocate
     &   ( FR(MLMNB,NI    ), FZ(MLMNB,NI  ), FPHV(MLMNB,NI   )
     &   , FRI(MLMNB,0:NI ),FZI(MLMNB,0:NI), FPHVI(MLMNB,0:NI)
     &   ,FPHVT(MLMNB, NI ),FPHVP(MLMNB,NI), fbsq(mlmnb,ni   )
     &   ,FBJAC(MLMNB,0:NI),FBS(MLMNB,0:NI), stat=ierr8      )
!      COMMON /BOOFOU/ FR,FZ,FPHV,FRI,FZI,FPHVI,FPHVT,FPHVP,fbsq,FBJAC
!     &               ,FBS
C
C...FORMERLY Variables in Common Block /BOORSP/
      allocate
     &    ( R (NJK,0:NI), Z (NJK,0:NI   ), phv (njk,ni),  RI (NJK,0:NI)
     &    ,ZI (NJK,0:NI), phvi (njk,0:ni), RT (NJK,NI ),  ZT (NJK,NI  )
     &    ,rs (njk,ni  ), RP (NJK,NI    ), ZP (NJK,NI ),  zs (njk,ni  )
     &    ,RSq(NJK,0:NI), phvp(njk      ), phvt(njk   ),  stat=ierr9)
!      COMMON /BOORSP/ R,Z,phv,RI,ZI,phvi,RT,ZT,rs,RP,ZP,zs,RSq,phvp,phvt
C
C...FORMERLY Variables in Common Block /COMERC/
      allocate 
     &    (BST(njk),  baldp(njk), curvds(njk), shearl(njk),stat=ierr10)
!      COMMON /COMERC/ BST, baldp, curvds, shearl
C
C...FORMERLY Variables in Common Block /MERFOU/
      allocate
     &    (fcurvs(mlmnb), flshea(mlmnb), fjpar(mlmnb), fgssu(mlmnb)
     &    ,stat=ierr11)
!      COMMON /MERFOU/ fcurvs,flshea,fjpar,fgssu
C
C...FORMERLY Variables in Common Block /BOOMET/
      allocate
     &    (BJAC(NJK,0:Ni),BJACS(NJK,0:NI), GTTL(NJK,0:Ni),GTPL(NJK,0:Ni)
     &    ,GPPL(NJK,0:Ni),GSSL(NJK,0:Ni ), GSTL(NJK,0:Ni),GSsu(NJK,0:Ni)
     &    ,stat=ierr12)
!      COMMON /BOOMET/ BJAC,BJACS,GTTL,GTPL,GPPL,GSSL,GSTL,GSSU
C
C...FORMERLY Variables in Common Block /BDERIV/
      allocate
     &    (BT(NJK,0:NI),   BP(NJK,0:NI),  BSQ(NJK,0:NI),    BS(NJK,0:NI)
     &    ,stat=ierr13)
!      COMMON /BDERIV/ BT,   BP,  BSQ,  BS
C
C...FORMERLY Variables in Common Block /RADFCT/
      allocate
     &              (    AM(  NI),  PTH(  NI),  PVP(  NI),PVPI(  NI)
     &              ,   FPP(NI+1),  FTP(NI+1),  PPI(  NI),  PP(  NI)
     &              ,  FPPP(NI+1), FTPP(NI+1),AIOTA(  NI), VVP(  NI)
     &              ,    f0(  ni),   f1(  ni),   f2(  ni),  f3(  ni)
     &              , fpppp(  ni),ftppp(  ni),WMAGV(  NI),WMAG(  NI)
     &              ,   CIV(  NI),  CJV(  NI), CIVP(  NI),CJVP(  NI)
     &              ,   CI (  NI),  CJ (  NI), CIP (  NI),CJP (  NI)
     &              ,  CIPI(  NI), CJPI(  NI), ppp (  ni),glm (ni,2)
     &              ,    VP ( NI),  VPP(  NI),EQUIV(  NI),EQUI(  NI)
     &              ,    s (0:ni),stat=ierr14)
!      COMMON /RADFCT/ AM,PTH,PVP,PVPI,FPP,FTP,PPI,PP,FPPP,FTPP,AIOTA
!     &               ,VVP,f0,f1,f2,f3,fpppp,ftppp,WMAGV,WMAG,CIV,CJV
!     &               ,CIVP,CJVP,CI,CJ,CIP,CJP,CIPI,CJPI,ppp,glm,VP,VPP
!     &               ,EQUIV,EQUI,s
C
C...FORMERLY Variables in Common Block /COBALL/
       allocate
     &    ( fbaldp(mlmnb), fbalds(mlmnb), fbalcp(mlmnb)
     &    , fbalcs(mlmnb), fbalcq(mlmnb), stat=ierr15)
!     COMMON /COBALL/ fbaldp, fbalds, fbalcp, fbalcs, fbalcq
C
C...FORMERLY Variables in Common Block /COBOOT/
      allocate ( jkmax(NI), stat=ierr16)
      allocate
     &       (pitch(npitch),bnorm(NJK, NI),g2(NJK)    ,g4(NJK)   ,gb(NI)
     &       ,g1mnin(MLMNB),g4mn(MLMNB)   ,rl31(NI)   ,rl32e(NI) ,gc(NI)
     &       ,boot(0:NI)   ,bsqav(NI)     ,bsqmax(NI) ,rl32i(NI) ,ft(NI)
     &       ,g2av(NI)     ,tempp(NI)     ,dboot(NI)  ,bjav(NI)  ,fc(NI)
     &       ,g1av(NPITCH,NI)             ,g4av(NPITCH,NI)
     &       ,rlampl(NI)   ,gbplat(NI)    ,rnue(NI)   ,bgradgmn(MLMNB)
     &       ,bdgradb(NJK) ,bdgradg(NJK)  ,bdglnb(NJK),bgradbmn(MLMNB)
     &       ,stat=ierr17)
!      COMMON /COBOOT/ pitch,bnorm,g2,g4,gb,g1mnin,g4mn,rl31,rl32e,rl32i
!     &               ,gc,boot,bsqav,bsqmax,ft,fc,g2av,tempp,dboot,bjav
!     &               ,g1av,g4av,jkmax,gbplat,rlampl,rnue,bdgradg,bdgradb
!     &               ,bdglnb,bgradbmn,bgradgmn
C
C...FORMERLY Variables in Common Block /PAPRPL/
C     PAPRPL: PRINT & PLOT SWITCHES
C
      end subroutine TPRALL_BAL
