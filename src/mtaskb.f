c----------------------------------------------------------------------c
C  18.12.89        LAST MODIFICATION 28.11.05     UHS: mtaskb DATA   D
c-----------------------------------------------------------------------
C
      subroutine mtaskb(tim)
C
      use specifications_mod, ONLY:ns_vmec,mpnt_vmec
      use tprb_param
      use physeq
      use radfct
      use triggs
      use paprpl
      use vmcfou
      use vmcrsp
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
C
C.. Implicits ..
      implicit none
C
!      include 'parbal.inc'
C
C.. External Calls ..
      external BOPHYS, extint, metric, vmtobo, mercie, second
!      include 'tprcom.bal'
C
C.. Internal scalars
      integer :: i, l,NI1, MLMNV1
      real :: tim(*), shalf
C ... Executable Statements ...
C
C
      NI1=ns_vmec-1 
      MLMNV1=mpnt_vmec

C     PERFORM MAPPING TO BOOZER COORDINATES.
c
c      do 10 i = 1,ni1
      call vmtobo (ni1,njk,nj,lmnb,lmnb0,mlmnb,lmnl,nper,lvmtpr,
     &             lrippl,tpi,dtdp,dph,dth,fpp,ftp,cjv,civ,mb,nb,ml,nl,
     &             tsin,tcos,vbp,vbt,vbsq,vjac,r,z,vl,phv,fvl,fr,fz,
     &             fphv,fphvt,fphvp,fbsq,fbjac,fvalf,fvgam,fvb,fvq,
     &             valft,vgamt,valf,vgam,vjtobj,thv,vq,sqbbsq)
C 10   continue
      call second (tim(5))
c
C     COMPUTE R/Z/PHI ON THE HALF INTEGER GRID.
c
c      do 20 l = 1,lmnb
      call extint (ni1,njk,lmnb,mm,nmin,nmax,nper,lmetpr,rplmin,mb,nb,
     &             lfrz,tsin,tcos,ri,zi,phvi,fr,fz,fphv,fri,fzi,fphvi,
     &             nbalmn)
C 20   continue
c
c     calculate metric elements, jacobian and geometry (r/z/phi) on
c     the half integer grid in real boozer coordinate space.
c
c      do 30 i = 1,ni1
      call metric(ni1,njk,lmnb,s,fpp,ftp,cjv,civ,mb,nb,tsin,
     &  tcos,r,z,phv,rt,zt,rp,zp,rs,zs,rsq,ri,zi,phvi,bsq,
     &bjac,vjac,gttl,gtpl,gppl,gssl,gstl,fr,fz,fbsq,fphv,phvp,phvt)
C 30   continue
      call second (tim(6))
c
      call bophys(ni1,njk,lmnb,lmnb0,mm,nmin,nmax,lcurrf,lvmtpr,lpress,
     &            lmetpr,rplmin,djp,dtdp,s,ftp,fpp,ftpp,fppp,ci,cj,cip,
     &            cjp,pp,wmag,vp,cipi,cjpi,pvpi,pth,ppi,equi,civ,cjv,
     &            vvp,mb,nb,tsin,tcos,lfrz,gttl,gtpl,gppl,gssu,bjac,
     &            bjacs,vjac,bs,bp,bt,bsq,fr,fbjac,fbs,fphv,fbsq)
      call second(tim(7))
!           print *,' lmnb0=',lmnb0,'   fr(lmnb0,ni1)=',fr(lmnb0,ni1)
C
      if (lpress.ne.-9)
     &call mercie(ni1,njk,lmnb,lmnb0,dtdp,s,ftp,fpp,ftpp,fppp,ci,cj,pp,
     &            pth,vp,cjp,cip,mb,nb,tcos,tsin,gssu,gssl,gstl,bjac,
     &            bjacs,bsq,glm,bs,bst,fbs,fbjac,baldp,fbaldp,fbalds,
     &            fbalcp,fbalcs,fbalcq,djp,nper,fr,fz,fcurvs,flshea,
     &            fjpar,fgssu,curvds,shearl,fbsq,sqbbsq)
c
C..OUTPUT IN UNIT 37 FOR GUIDING CENTRE CALCULATIONS
!       write(37,820) ni1, lmnb
      do i=1,ni1
      do l=1,lmnb
!       write(37,821) mb(l),nb(l),fr(l,i),fz(l,i),fphv(l,i),fbsq(l,i)
      enddo
      enddo
      do i=1,ni1
       shalf=0.5*(s(i)+s(i-1))
!       write(37,822)shalf,-ftp(i),-fpp(i),-ci(i),-cj(i),-cip(i),-cjp(i)
      enddo
 820  format(1x,2i6)
 821  format(1x,2i8,1p4e14.6)
 822  format(1x,1p7e14.6)
!.. BOOTSTRAP CURRENT CALCULATION
      call second(tim(8))
      if (lbootj.eq.1) then
      call bootsj(ni1,njk,lmnb,npitch,dtdp,pitch,bnorm,fbjac,bjac,
     &bsq,mb,nb,tcos,tsin,fc,ft,gb,g1mnin,g4mn,g1av,g4av,g2,g4,g2av,
     &s,vp,ftp,fpp,ci,cj,pth,pp,ppi,tempp,tpi,boot,bjav,sqbbsq,
     &dboot,rl31,rl32e,rl32i,jkmax,bsqav,bsqmax,djp,bdgradg,bdgradb,
     &bgradbmn,bgradgmn,bdglnb,gbplat,rlampl,rnue)
      end if
!
      return
      end subroutine mtaskb
