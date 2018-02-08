c-----------------------------------------------------------------------
C  29.03.88        LAST MODIFICATION 04.04.96     UHS: VMTOBO DATA D
c-----------------------------------------------------------------------
C
      subroutine vmtobo (ni,njk,nj,lmnb,lmnb0,mlmnb,lmnl,nper,lvmtpr
     &,lrippl,tpi,dtdp,dph,dth,fpp,ftp,cjv,civ,mb,nb,ml,nl,tsin,tcos,vbp
     &,vbt,vbsq,vjac,r,z,vl,phv,fvl,fr,fz,fphv,fphvt,fphvp,fbsq,fbjac
     &,fvalf,fvgam,fvb,fvq,valft,vgamt,valf,vgam,vjtobj,thv,vq,sqbbsq)
C
C     COMPUTES VALF, FVALF, VGAM, FVGAM, VQ, FVQ AND VJTOBJ FOR MAPPING
C        GIVES FR, FZ, FPHV, R, Z, PHV IN BC
C              BY DIRECT INTEGRATION
C
C.. Implicits ..
      implicit none
C
C.. Formal Arguments ..
      integer njk,lmnb,lmnb0,ni,nj,mlmnb,lmnl,nper,lvmtpr,mb(*),nb(*),
     &        ml(*),nl(*),lrippl
      real tpi,dtdp,dph,dth,fpp(*),ftp(*),cjv(*),civ(*),tsin(njk,*),
     &     tcos(njk,*),vbp(njk,0:*),vbt(njk,0:*),vbsq(njk,0:*),
     &     vjac(njk,0:*),r(njk,0:*),z(njk,0:*),vl(njk,*),phv(njk,*),
     &     fvl(lmnl,*),fr(lmnb,*),fz(lmnb,*),fphv(lmnb,*),fphvt(lmnb,*),
     &     fphvp(lmnb,*),fbsq(lmnb,*),fbjac(lmnb,0:*),fvalf(*),fvgam(*),
     &     fvb(*),fvq(*),valft(*),vgamt(*),valf(*),vgam(*),vjtobj(*),
     &     thv(*),vq(*),sqbbsq(*)
C
C.. Local Scalars ..
      integer I,J,JK,K,L,LB,LL,NBPFP
      real AT,cosat,FPPI,FTPI,TBZT,TWO
C
C.. Intrinsic Functions ..
!      intrinsic MOD, SIN, cos, cvmgm
      intrinsic MOD, SIN, cos, merge
C
C ... Executable Statements ...
C
C     if (i.eq.0) go to 300
C     if (i.eq.1)WRITE( 6,1001)
C     if (i.eq.1)WRITE(16,1001)
 1001 format(/1x,'vmtobo: i,j,vl,vq,valf,vgam'/)
C
c..DETERMINE THE INDEX FOR WHICH mb=0 AND nb=0
      do l=1,lmnb
        if (mb(l).eq.0 .and. nb(l).eq.0) then
        lmnb0 = l
        endif
      end do
C
C     COMPUTE FVQ AND VQ FROM VBT AND VBP
C
      do 300 i = 1,ni
C
        do  l = 1,lmnb
         fvb(l) = 0.
         fvq(l) = 0.
C
         if (nb(l).ne.0) then
             do  jk = 1,njk
             fvb(l) = fvb(l) + 2.* vbp(jk,i) * tcos(jk,l)
             end do
cperf   fvq(l) = - dtdp / nb(l) * fvb(l) / tpi
             fvq(l) = - dtdp / nb(l) * fvb(l)
         elseif (nb(l).eq.0.and.mb(l).ne.0) then
             do  jk = 1,njk
             fvb(l) = fvb(l) + 2.* vbt(jk,i) * tcos(jk,l)
             end do
cperf   fvq(l) =   dtdp / mb(l) * fvb(l) / tpi
        fvq(l) =   dtdp / mb(l) * fvb(l)
C
        endif
C
         end do  
cc141   continue
cc  140   continue
C
cpara      do 150 i=1,ni
           do  jk=1,njk
C
            vq(jk) = 0.
C
            do  l=1,lmnb
C
               vq(jk) = vq(jk) + fvq(l) * tsin(jk,l)
            end do
           end do
c
C     ALL MAPPING FUNCTIONS COMPUTED ON HALF GRID
C
cpara      do 200 i=1,ni
C
         fppi = fpp(i)
         ftpi = ftp(i)
         sqbbsq(i)  = fppi * cjv(i) - ftpi * civ(i)
C
C..NOTE THAT THE SUPPLEMENTARY MODES FOR BALLOONING DONT CONTRIBUTE
         do  ll = 1,lmnl
c         do  lb = 1,lmnb
c            if (ml(ll).eq.mb(lb).and.nl(ll).eq.nb(lb)) then
!            lb=cvmgm(ll,ll+1,ll-lmnb0)
            lb = merge(ll,ll+1,ll<lmnb0)
            fvalf(lb ) = ( fppi*fvq(lb)- civ(i)*fvl(ll,i))/sqbbsq(i)
            fvgam(lb ) = ( ftpi*fvq(lb)- cjv(i)*fvl(ll,i))/sqbbsq(i)
C      IF(I.LE.NI-3) WRITE(16,1003) LL,LB,ML(LL),MB(LB),NL(LL),NB(LB)
C    $              ,FVQ(LB),FVL(LL,I),FVLI(LL,I)
 1003 format(' vmtobo:',6i3,3x,1p3e16.6)
C
c            end if
         end do
            fvalf(lmnb0) = 0.
            fvgam(lmnb0) = 0.
C
         do  jk=1,njk
C
            valf(jk  ) = ( fppi*vq(jk) - civ(i)*vl(jk,i) )
     $                   / sqbbsq(i)
            vgam(jk  ) = ( ftpi*vq(jk) - cjv(i)*vl(jk,i) )
     $                   / sqbbsq(i)
          vjtobj(jk  ) = vjac(jk,i) * vbsq(jk,i) / sqbbsq(i)
C        WRITE (16,1002)I,JK,VL(JK,I),VQ(JK),VALF(JK),VGAM(JK)
        end do
 1002 format (1x,2i3,1p4e15.6)
C
C        IF (I.LE.8) THEN
C        DO 208 L = 1,LMNB
C 208    WRITE (16, *  ) I,ML(L),MB(L),FVL(L),FVQ(L)
C
C
C        END IF
      do  jk = 1,njk
C
C     PRESET PHV ETC
C
            k    =      ((jk-1)/nj)
           j     =       mod((jk-1),nj)
cperf  phv(jk,i) = dph * k
cperf  thv(jk  ) = dth * j
       phv(jk,i) = dph * k * tpi / nper
       thv(jk  ) = dth * j * tpi
      vgamt(jk ) = 0.
      valft(jk ) = 0.
C
      end do
C
C     COMPUTE VALFT, VGAMT
C
      do  l  = 1,lmnb
        do  jk = 1,njk
cperf    tbzt = tpi * mb(l) * tcos(jk,l)
         tbzt =       mb(l) * tcos(jk,l)
         vgamt(jk  ) = vgamt(jk  ) + fvgam(l  ) * tbzt
         valft(jk  ) = valft(jk  ) + fvalf(l  ) * tbzt
        end do
      end do
C
C
C     DIRECT FOURIER INTEGRALS FOR BC COEFFICIENTS
C
      do l=1,lmnb
C
           fr(l,i) = 0.
           fz(l,i) = 0.
         fbsq(l,i) = 0.
        fbjac(l,i) = 0.
         fphv(l,i) = 0.
        fphvt(l,i) = 0.
        fphvp(l,i) = 0.
C
         two = 2.
         if (mb(l).eq.0.and.nb(l).eq.0) then
           two = 1.
         end if
            do  jk = 1,njk
C
cperf       at = tpi * (mb(l)*( thv(jk  ) + valf(jk  ) )
            at =       (mb(l)*( thv(jk  ) + valf(jk  ) )
     $                - nb(l)*( phv(jk,i) + vgam(jk  ) ) )
            cosat = cos(at)
            fr(l,i) = fr(l,i) + r (jk,i) * vjtobj(jk  ) * cosat
            fz(l,i) = fz(l,i) + z (jk,i) * vjtobj(jk  ) * sin(at)
            fbsq(l,i) = fbsq(l,i)+ vbsq(jk,i) * vjtobj(jk) * cosat
            fbjac(l,i) = fbjac(l,i) + vjac(jk,i) * cosat
            fphvt(l,i) = fphvt(l,i) - vgamt(jk  ) * cosat
            fphvp(l,i) = fphvp(l,i) + (1.+ valft(jk) ) * cosat
C
            end do
C
            fr(l,i) = dtdp * fr(l,i) * two
            fz(l,i) = dtdp * fz(l,i) * two
          fbsq(l,i) = dtdp * fbsq(l,i) * two
         fbjac(l,i) = dtdp *fbjac(l,i) * two
         fphvt(l,i) = dtdp *fphvt(l,i) * two
         fphvp(l,i) = dtdp *fphvp(l,i) * two
C
         end do
C
         do  l = 1,lmnb
C
          if (mb(l).ne.0 .or. nb(l).ne.0) then
           if (mb(l).ne.0) then
cperf       fphv(l,i) = + 1./(tpi*mb(l)) * fphvt(l,i)
            fphv(l,i) = + 1./(    mb(l)) * fphvt(l,i)
           else
cperf       fphv(l,i) = - 1./(tpi*nb(l)) * fphvp(l,i)
            fphv(l,i) = - 1./(    nb(l)) * fphvp(l,i)
           end if
          end if
         end do
C
         if((i.le.16.or.i.gt.ni-16).and.lvmtpr.eq.1) then
C         IF(I.GT.NI-5          ) THEN
          do  l = 1,lmnb
          nbpfp = nb(l)/nper
cperf     write (16,1000) i,l,mb(l),nb(l),fr(l,i),fz(l,i),fphv(l,i)
!          write (16,1000) i,l,mb(l),nbpfp,fr(l,i),fz(l,i),fphv(l,i)
!     $                                   ,fbsq(l,i)
cperf     write ( 6,1000) i,l,mb(l),nb(l),fr(l,i),fz(l,i),fphv(l,i)
          write ( 6,1000) i,l,mb(l),nbpfp,fr(l,i),fz(l,i),fphv(l,i)
     $                                   ,fbsq(l,i)
          end do
          end if
 1000     format(1x'bofou',4i3,1p4e15.6)
  200     continue
C
C          if (lrippl.ne.1) return
C
C          if (i.eq.1) write (6,1010)
C 1010     format(//,'    i         s          rms-ripple')
C          bmn0 = 0.
C          bmnp = 0.
C          do 300 l=1,lmnb
C          nbpfp = nb(l)/nper
C          bmn0 = bmn0 + fbsq(l,i)
C          bmnp = bmnp + fbsq(l,i)*(-1)**nbpfp
  300      end do
C          rmsrip1 = 2.*(bmn0-bmnp)/(bmn0+bmnp)
C          bmn0 = sqrt(bmn0)
C          bmnp = sqrt(bmnp)
C          rmsrip2 = 2.*(bmn0-bmnp)/(bmn0+bmnp)
C          if (rmsrip1 .ge. 0.) then
C          rmsrip1 = sqrt(rmsrip1)
C          t1 = 0.5 * (s(i)+s(i-1))
C          write (6,1011) i, t1, rmsrip1, rmsrip2
C          else
C          write (6,1012) i
C 1012   format(i6,3x,'attention: the ripple at the outside is negative')
C         endif
C 1011     format (i6,1p3e15.6)
C
      return
      end
