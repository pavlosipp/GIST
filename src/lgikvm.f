C-----------------------------------------------------------------------
C  27.03.88        LAST MODIFICATION 01.02.94     UHS: LGIKVM DATA D
C-----------------------------------------------------------------------
C
      subroutine LGIKVM(ni,nj,njk,lmnv,nper,tpi,dph,mx,nx,tcos,tsin,frv,
     &                  fzv,fvjac,r,z,ri,zi,phvi,vjac,rv,zv,rvt,zvt,rvp,
     &                  zvp,rvq,rvqt,gttv,gtpv,gppv,lmnb)
C
C 
C.. Implicits .. 
      implicit none
C 
C.. Formal Arguments .. 
      integer njk,lmnv,lmnb,ni,nj,nper,mx(*),nx(*)
      real tpi,dph,tcos(njk,*),tsin(njk,*),frv(lmnv,0:*),fzv(lmnv,0:*),
     &     fvjac(lmnb,0:*),r(njk,0:*),z(njk,0:*),ri(njk,0:*),
     &     zi(njk,0:*),phvi(njk,0:*),vjac(njk,0:*),rv(njk,0:*),
     &     zv(njk,0:*),rvt(njk,0:*),zvt(njk,0:*),rvp(njk,0:*),
     &     zvp(njk,0:*),rvq(njk,0:*),rvqt(njk,0:*),gttv(njk,0:*),
     &     gtpv(njk,0:*),gppv(njk,0:*)
C 
C.. Local Scalars .. 
      integer I,JK,L
C 
C ... Executable Statements ...
C 
C
C     COMPUTES R,Z,     LOWER METRIC ELEMENTS AND EQUILIBRIUM
C                         FUNCTIONS IN VMEC/FIT COORDINATES
C     COMPUTE RV,ZV, ... ON INTEGER GRID
C
      do I = 0,NI
C
        do JK = 1,NJK
C
          RI(JK,I) = 0.
          ZI(JK,I) = 0.
          PHVI(JK,I) = DPH * ((JK-1)/NJ) * TPI / NPER
          VJAC(JK,I) = 0.
          RV(JK,I) = 0.
          ZV(JK,I) = 0.
          RVT(JK,I) = 0.
          ZVT(JK,I) = 0.
          RVP(JK,I) = 0.
          ZVP(JK,I) = 0.
        end do
C
        do L = 1,LMNV
C
          do JK = 1,NJK
C
            VJAC(JK,I) = VJAC(JK,I) + FVJAC(L,I)*TCOS(JK,L)
            RV(JK,I) = RV(JK,I) + FRV(L,I)*TCOS(JK,L)
            ZV(JK,I) = ZV(JK,I) + FZV(L,I)*TSIN(JK,L)
CPERF          RVT(JK,I) = RVT(JK,I) - FRV (L,I) * TPI*MX(L)*TSIN(JK,L)
CPERF          ZVT(JK,I) = ZVT(JK,I) + FZV (L,I) * TPI*MX(L)*TCOS(JK,L)
CPERF          RVP(JK,I) = RVP(JK,I) + FRV (L,I) * TPI*NX(L)*TSIN(JK,L)
CPERF          ZVP(JK,I) = ZVP(JK,I) - FZV (L,I) * TPI*NX(L)*TCOS(JK,L)
            RVT(JK,I) = RVT(JK,I) - FRV(L,I)*MX(L)*TSIN(JK,L)
            ZVT(JK,I) = ZVT(JK,I) + FZV(L,I)*MX(L)*TCOS(JK,L)
            RVP(JK,I) = RVP(JK,I) + FRV(L,I)*NX(L)*TSIN(JK,L)
            ZVP(JK,I) = ZVP(JK,I) - FZV(L,I)*NX(L)*TCOS(JK,L)
          end do
C
        end do
C     IF (I.EQ.0) WRITE(16,*)((JK,RV(JK,I),ZV(JK,I)),JK=1,NJK)
C
C
Cpara      DO 140 I = 0,NI
        do JK = 1,NJK
C
          RVQ(JK,I) = RV(JK,I)**2
          RVQT(JK,I) = 2. * RV(JK,I) * RVT(JK,I)
          GTTV(JK,I) = (RVT(JK,I)**2+ZVT(JK,I)**2)
          GPPV(JK,I) = (RVP(JK,I)**2+ZVP(JK,I)**2+RVQ(JK,I))
CPERF$                  + TPONPQ *    RVQ(JK,I)          )
          GTPV(JK,I) = (RVT(JK,I)*RVP(JK,I)+ZVT(JK,I)*ZVP(JK,I))
CPERF$                  + TPONPQ *    RVQ(JK,I)          )
        end do
C
C     INITIALIZE R/Z ON HALF-INTEGER GRID
        if (i .gt. 0) then
          do jk = 1,njk
            r(jk,i) = 0.
            z(jk,i) = 0.
          end do
        end if
C
      end do
      end
