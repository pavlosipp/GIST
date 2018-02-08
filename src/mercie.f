!--------0---------0---------0---------0---------0---------0---------0-c
!  03.09.90          LAST MODIFICATION 28.11.05       UHS: MERCIE DATA D   
!--------0---------0---------0---------0---------0---------0---------0-c
      subroutine mercie(ni,njk,lmnb,lmnb0,dtdp,s,ftp,fpp,ftpp,fppp,ci,cj
     &,pp,pth,vp,cjp,cip,mb,nb,tcos,tsin,gssu,gssl,gstl,bjac,bjacs,bsq
     &,glm,bs,bst,fbs,fbjac,baldp,fbaldp,fbalds,fbalcp,fbalcs,fbalcq,djp
     &,nper,fr,fz,fcurvs,flshea,fjpar,fgssu,curvds,shearl,fbsq,sqgbsq)

          use specifications_mod, ONLY:ns_vmec, mpnt_vmec, vmec_file,
     &                                     table_tag

! 
! .. IMPLICITS
      implicit none
!.. Parameters ..
!      real gbgmin
!      parameter (gbgmin = 1e-300)
!
!
! .. FORMAL ARGUMENTS
      real :: dtdp,s(0:*),ftp(*),fpp(*),ftpp(*),fppp(*),ci(*),cj(*)
     &,pp(*),cjp(*),cip(*),vp(*),pth(*),tsin(njk,*),tcos(njk,*)
     &,bjac(njk,0:*),bjacs(njk,0:*),bsq(njk,0:*)
     &,fbsq(lmnb,ni),fbs(lmnb,0:ni),fbjac(lmnb,0:ni),glm(ni,*)
     &,gssu(njk,0:*),gssl(njk,0:*),gstl(njk,0:*),bs(njk,0:*)
     &,djp
      integer mb(*),nb(*),ni,njk,lmnb,lmnb0,nper
!           
      real ::  bst(*)   , baldp(*),shearl(*),curvds(*) ,sqgbsq(*)
     &,fbaldp(*),fbalds(*),fbalcp(*),fbalcs(*),fbalcs1(lmnb),fbalcq(*)
     &,fcurvs(*),flshea(*),fjpar(*) ,fgssu(*) ,fr(lmnb,*),fz(lmnb,*)
     &,fbalbs(lmnb)
!
! ... LOCAL SCALARS
      integer :: i, jk, l
      real :: t0, t1, t2, t3, t4, t5,t6, t7, t8, t9,hq, dper, dsec, 
     &        hdsc, hdsq, dsca, qsurf, gmer, gmerc, shamp, cjogb2,
     &        ppogb2, cmnjac, qmnjac, sqrbsq, c1bamp, c2bamp, c3bamp,
     &        gbgrad, tbalin, fintls, djpboo

      real, dimension(ni) :: t_0

! .. Intrinsic Functions ..
      intrinsic  sqrt, max
!
      data djpboo /1.0e-16/
      djpboo = max(djp,djpboo)
!
!      write(66,858)nper,lmnb,ni-2
!
      do 100 i=2,ni-1
!..surface quantities
      sqgbsq(i)=fpp(i)*cj(i) - ftp(i)*ci(i)
      t0 = ftpp(i) - ftp(i) * fppp(i) / fpp(i)
      t1 = sqgbsq(i) / (t0 * t0)
      t2 = pp(i) / (fpp(i) * fpp(i))
      t3 = cj(i) * fppp(i) - ci(i) * ftpp(i)
      t4 = t0 / (fpp(i) * ftp(i))
      t5 = cj(i) * pp(i)
      t6 = t1 * t4
      t_0= t0/fpp(i) !qprime

c..initialisation
      hq   = 0.
      dper = 0.
      dsec = 0.
      hdsc = 0.
      hdsq = 0.
      do 5 jk=1,njk
         baldp(jk)  = 0.
         curvds(jk) = 0.
 5       bst(jk)    = 0.
c
c..the derivative of b-sub-s with respect to theta required for Mercier.
      do 10 l=1,lmnb
      do 10 jk=1,njk
 10      bst(jk) = bst(jk) + mb(l) * fbs(l,i) * tcos(jk,l)
c
c..the mercier coefficients.
      do 20 jk=1,njk
         hq   = hq   + 1. / gssu(jk,i)
         dper = dper - bjacs(jk,i) + (pp(i)*bjac(jk,i) + t3) / bsq(jk,i)
         dsca = cjp(i) - bst(jk) + t5 / bsq(jk,i)
         dsec = dsec + dsca
         hdsc = hdsc + dsca / gssu(jk,i)
 20      hdsq = hdsq + dsca * dsca / gssu(jk,i)
      hq   = dtdp * t1 * hq
      dper = dtdp * t2 * dper
      dsec = dtdp * t4 * dsec
      hdsc = dtdp * t6 * hdsc
      hdsq = dtdp * t6 * t4 * hdsq
      dsca = hq * (dper + dsec + hdsq)
      dper = hq * (dper + dsec) - hdsc
      hdsq = 1. / (hq * hdsq - hdsc * hdsc)
      glm(i,2) = sqrt((dper*dper*hdsq + 1.)*hdsq)
      glm(i,2) = 1. + 0.5 * (dper*hdsq - glm(i,2))
      glm(i,1) = (0.5 + hdsc) * (0.5 + hdsc) - dsca      
      qsurf = ftp(i)/fpp(i)
      gmerc = glm(i,1)
      gmer  = - glm(i,2)
!      write (6,890) i,qsurf,gmerc,gmer
!  810 write(16,890) i,qsurf,gmerc,gmer
  890 format (i4,1p1e16.8,1p2e12.4)
C
C --- DETERMINATION OF THE COEFFICIENTS FOR BALLOONING CALCULATION
C
c..surface quantities
      t8 = fpp(i) * cjp(i) - ftp(i) * cip(i) + t3
      t9 = - t2 / (t1 * t0)
      t3 = 2.*dtdp
      t7 = t3 / t1
      t6 = 1. / sqgbsq(i)
      t5 = 2. * t3 * t0 / ftp(i)
      t4 = 2. * t0 * t6 * cj(i) / ftp(i)
      t2 = - t3 * t2
      t1 = 2. * pp(i)
      shamp  = 1. / (t6*ftp(i))
      cjogb2 = cj(i) * t6
      ppogb2 = pp(i) * t6
c
c..the derivatives of the jacobian with respect to theta and phi required for
c..normal curvature.
      do 50 l=1,lmnb
      qmnjac = - (mb(l) * fpp(i) - nb(l) * ftp(i)) * fbjac(l,i)
      cmnjac = - (mb(l) * ci(i)  - nb(l) * cj(i) ) * fbjac(l,i)
      do 50 jk=1,njk
      curvds(jk) = curvds(jk) + cmnjac * tsin(jk,l)
 50   baldp(jk)  = baldp(jk)  + qmnjac * tsin(jk,l)
C
c..evaluation of the integrated local shear contribution in array shearl.
c..determination of the normal curvature(2p'(s)*jacobian*(kappa.grad(s))
c..array gssl used to store the coefficients of field line bending.
      do 55 jk=1,njk
      sqrbsq = sqrt(bsq(jk,i))
      shearl(jk)= shamp * (gstl(jk,i)-cjogb2*bs(jk,i))/gssu(jk,i)
c      curvds(jk)= ppogb2*gssu(jk,i)*(bjac(jk,i)*(t1*bjac(jk,i) + t8)
c     >    - bjacs(jk,i)/t6 - shearl(jk)*curvds(jk) + bs(jk,i)*baldp(jk))
      baldp(jk) = (t1*bjac(jk,i) + t8)/bsq(jk,i) + t6*bs(jk,i)*baldp(jk)
     >           - bjacs(jk,i)     
      curvds(jk)= pp(i)*gssu(jk,i) *(baldp(jk)-t6*shearl(jk)*curvds(jk))
!00      gssl(jk,1) = 1. / (bjac(jk,i)*gssu(jk,i))
!00      gssl(jk,2) = (cj(i)*bs(jk,i) - sqgbsq(i)*gstl(jk,i))/sqrbsq
!00 55   gssl(jk,3) = gssu(jk,i)/sqrbsq
 55   bst(jk) = gssl(jk,i) - t6 * bs(jk,i) * bs(jk,i)
C
      do 60 l=1,lmnb
      fcurvs(l) = 0.
      flshea(l) = 0.
      fbaldp(l) = 0.
      fbalcp(l) = 0.
      fbalcs(l) = 0.
      fbalbs(l) = 0.
      fbalcs1(l) = 0.
 60   fbalcq(l) = 0.
      do 65 l=1,lmnb
      do 65 jk=1,njk
      fcurvs(l) = fcurvs(l) + curvds(jk)  * tcos(jk,l)
      fbaldp(l) = fbaldp(l) + baldp(jk)   * tcos(jk,l)
      fbalcp(l) = fbalcp(l) +   bst(jk)   * tcos(jk,l)
      fbalcq(l) = fbalcq(l) +  gssu(jk,i) * tcos(jk,l)
      fbalbs(l) = fbalbs(l) +  bs(jk,i) * tsin(jk,l) !pax
!00      fbalcq(l) = fbalcq(l) + gssl(jk,3)  * tcos(jk,l)
!00      fbalcp(l) = fbalcp(l) + gssl(jk,1)  * tcos(jk,l)
      flshea(l) = flshea(l) + shearl(jk)  * tsin(jk,l)
!00 65   fbalcs(l) = fbalcs(l) + gssl(jk,2)  * tsin(jk,l)
      fbalcs1(l) = fbalcs1(l) +  gstl(jk,i) * tsin(jk,l)
 65   fbalcs(l) = fbalcs(l) +  gstl(jk,i) * tsin(jk,l)
     
      
C
      c1bamp = 2.0*dtdp
      c2bamp = c1bamp/ftp(i)
      c3bamp = c1bamp*t0
      do 70 l=1,lmnb
      fcurvs(l) = t3 * fcurvs(l)
      flshea(l) = -t3 * (mb(l)*fpp(i) - nb(l)*ftp(i))*flshea(l)
      gbgrad    = mb(l)*fpp(i) - nb(l)*ftp(i)
      gbgrad    = gbgrad / (gbgrad*gbgrad + djpboo * (gbgrad+fpp(i)) *
     &                      (gbgrad+fpp(i)))
      fjpar(l)  = ppogb2 * fbjac(l,i) * (mb(l)*ci(i) - nb(l)*cj(i))
     &                  * gbgrad
!      gbgrad    = cvmgn(gbgrad,gbgmin,gbgrad)
!      fjpar(l)  = ppogb2 * fbjac(l,i) * (mb(l)*ci(i) - nb(l)*cj(i))
!     &                 /gbgrad
      fgssu(l)  = t3 * fbalcq(l)
      fbaldp(l) = t2 * fbaldp(l)
      fbalds(l) = t9 * (mb(l)*ci(i) - nb(l)*cj(i)) * fbjac(l,i)
!00      fbalcp(l) = c1bamp * fbalcp(l)
!00      fbalcs(l) = c2bamp * fbalcs(l)
!00 70   fbalcq(l) = c3bamp * fbalcq(l)
      fbalcp(l) = t3 * fbalcp(l)
      fbalcs(l) = t4 * fbs(l,i) - t5 * fbalcs(l)
 70   fbalcq(l) = t7 * fbalcq(l)
C
      fcurvs(lmnb0) = 0.5 * fcurvs(lmnb0)
      flshea(lmnb0) = fpp(i)*ftpp(i) - ftp(i)*fppp(i)
      fjpar(lmnb0)  = t6 * (cj(i)*cip(i) - ci(i)*cjp(i))
      fgssu(lmnb0)  = 0.5 * fgssu(lmnb0)
      fbaldp(lmnb0) = 0.5 * fbaldp(lmnb0)
      fbalcs(lmnb0) = 0.
      fbalcs1(lmnb0) = 0.
      fbalbs(lmnb0) = 0.
      fbalcp(lmnb0) = 0.5 * fbalcp(lmnb0)
      fbalcq(lmnb0) = 0.5 * fbalcq(lmnb0)
C
      qsurf = ftp(i) / fpp(i)
      t1 = 0.5 * (s(i) + s(i-1))
      tbalin = t0 / (2.*ftp(i))
      
      OPEN(25,file='../int/tpr_'//trim(vmec_file)//'_'//trim(table_tag),
     &        action="READWRITE",STATUS="OLD",POSITION="APPEND")
      write (25,98)t1,qsurf, pp(i),vp(i),tbalin,sqgbsq(i), ftp(i),pth(i)
     &             ,t_0(i)
      write (25,99) (fbalcp(l),l=1,lmnb)
      write (25,99) (fbalcs(l),l=1,lmnb)
      write (25,99) (fbalcq(l),l=1,lmnb)
      t1 = 1. / (fpp(i) * fpp(i))
      write (25,99) (t1*fbjac(l,i),l=1,lmnb)
      write (25,99) (fbaldp(l),l=1,lmnb)
      write (25,99) (fbalds(l),l=1,lmnb)
      write (25,99) (fgssu(l),l=1,lmnb)
      write (25,99) (fbalcs1(l),l=1,lmnb)
      write (25,99) (fbalbs(l),l=1,lmnb)
      CLOSE(25)
C
c..print the normal curvature, the local shear, the parallel currrent
c..density, |grad s|^2, B^2, the integrated local shear and the Jacobian
c..fourier components for plotting purposes.
      do 95 l=1,lmnb
      gbgrad    = mb(l)*fpp(i) - nb(l)*ftp(i)
      fintls = flshea(l)
      if(l.ne.lmnb0)fintls = fintls/gbgrad
!      write(66,859) mb(l),nb(l),fr(l,i),fz(l,i),fcurvs(l),flshea(l)
!     >             ,fjpar(l),fgssu(l),fbsq(l,i),fintls,fbjac(l,i)
 95   continue
C
 100  continue
c 98   format (1p2e22.14,1p2e16.8,1p1e22.14,1p2e16.8)
 98   format (1p9e16.8)
 99   format (1p6e22.14)
C
!     write(66,860) (sqgbsq(i),i=2,ni-1)
 858  format(3i6)
 859  format(1x,2i6,1p9e13.6)
 860  format(1x,1p6e13.6)
C
      return
      end
