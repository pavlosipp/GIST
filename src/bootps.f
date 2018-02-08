!--------0---------0---------0---------0---------0---------0---------0-c
!  27.04.01          LAST MODIFICATION 05.23.03       UHS: BOOTSJ DATA D   
!--------0---------0---------0---------0---------0---------0---------0-c
      subroutine bootsj(ni,njk,lmnb,npitch,dtdp,pitch,bnorm,fbjac,bjac,
     &bsq,mb,nb,tcos,tsin,fc,ft,gb,g1mnin,g4mn,g1av,g4av,g2,g4,g2av,
     &s,vp,ftp,fpp,ci,cj,pth,pp,ppi,tempp,tpi,boot,bjav,sqgbsq,
     &dboot,rl31,rl32e,rl32i,jkmax,bsqav,bsqmax,djp,bdgradg,bdgradb,
     &bgradbmn,bgradgmn,bdglnb,gbplat,rlampl,rnue)
! 
! .. IMPLICITS
      implicit none
!
! .. FORMAL ARGUMENTS
      real :: dtdp,tpi,djp,vp(*),s(0:*),bnorm(njk,*),tcos(njk,*)
     &,tsin(njk,*),fbjac(lmnb,0:*),bjac(njk,0:*),bsq(njk,0:*),pitch(*)
     &,ftp(*),fpp(*),ci(*),cj(*),fc(*),ft(*),gb(*),gbplat(*)
     &,g1av(npitch,*),g4av(npitch,*),g2av(*),g2(*),g4(*),g1mnin(*)
     &,g4mn(*),boot(0:*),bjav(*),dboot(*),rl31(*),rl32e(*),rl32i(*)
     &,bsqav(*),bsqmax(*),tempp(*),pp(*),pth(*),ppi(*),sqgbsq(*)
     &,rlampl(*),rnue(*),bdgradb(*),bdgradg(*),bdglnb(*),bgradbmn(*)
     &,bgradgmn(*)
      integer :: mb(*),nb(*),ni,njk,lmnb,npitch,jkmax(*)
!           
! ... LOCAL SCALARS
      integer :: i, jk, l, lpitch, lden, jbiter, nine
      real :: zero,qden,zeff,delpch,rmue1a,rmui1a,rmue2a,rmui2a,rmue3a
     &       ,rmui3a,rlee11,rlee12,rlee22,rlii22,shalf,qlmn,g2mn,djpboo
     &       ,cn, djpb, bdgrgmn, bxgrgmn, ftpogb, dpitdx(npitch)
     &       ,boos(0:ni),dboos(ni),bjos(ni),osl31,osl32,alpha0,f32ee
     &       ,f32ei,ftsq,ftcb,ftqt,ftpsgn,twopi,sqrtpi,sqrtbsq
     &       ,densty,rlame,tonax,dn20,pax,bsqavax,betax,xlmax,pdpit
! .. Intrinsic Functions ..
      intrinsic abs, max, sqrt, min, atan, tanh, cosh
!
      data zero,qden,lden,cn,zeff,djpboo / 0.0,1.0,2,0.999,1.0,1.e-16 /
      data nine / 9 /
      data xlmax / 6.0 /
!
! ... READ INPUT DATA
      read (5,1001) tonax, zeff, cn, qden, lden, nine
 1001 format (//,1p4e12.4,6x,2i3)
! ... SET UP SUCH THAT IF THE MAGNETIC FIELD IS IN TESLA, THEN tonax
! ... IS THE TEMPERATURE ON AXIS IN kev, WITH INTERNAL PRESSURE CORRESPONDING
! ... TO mu_0 p WITH p in (N/m^2) TO YIELD dn20 AS DENSITY ON AXIS IN UNITS
! ... OF 10^(20)/m^3.
! ... MESH ACCUMULATION OF PITCH ANGLE INTEGRALS NEAR BOUNCE POINT ACCOMPLISHED
! ... WITH USE OF LAMBDA = TANH(x)
!
! ... INITIALISE
      twopi = 8.0 * atan(1.0)
      sqrtpi = sqrt(twopi/2.0)
      ftpsgn = tpi*ftp(ni)/abs(ftp(ni))
      pax    = 1.5 * pth(1) - 0.5 * pth(2)
      dn20   = 49.67 * pax / tonax
      pp(1)  = 1.5 * ppi(1)    - 0.5 * ppi(2)
      pp(ni) = 1.5 * ppi(ni-1) - 0.5 * ppi(ni-2)
      djpboo = max(djp,djpboo)
!      delpch = 1. / npitch
      delpch = xlmax / npitch
      do 5 lpitch=1,npitch
        pitch(lpitch)  = tanh((lpitch - 0.5) * delpch)
        dpitdx(lpitch) = 1.0 / cosh((lpitch - 0.5) * delpch)**2
 5    end do   
! ... INITIALISE AND EVALUATE MODEL p(0) N^hat(s) T^hat-prime(s)
! ... ASSUMES DENSITY N^hat NORMALISED TO VALUE ON AXIS AS (1-cn*s**lden)**qden
      do 10 i=1,ni
        bsqmax(i) = zero
        fc(i)     = zero
        g2av(i)   = zero
        gb(i)     = zero
        rlampl(i) = zero
        gbplat(i) = zero
        shalf     = 0.5 * (s(i) + s(i-1))
        tempp(i)  = pp(i) + pth(i) * cn*real(lden)*qden*shalf**(lden-1)
     &                / (1.0 - cn*shalf**lden)
 10   end do
! ... DETERMINE BSQ-max, B/B-max and  B-average ON EACH FLUX SURFACE
! ... DETERMINE jkmax, THE INDEX FOR WHICH BSQ=BSQ-max
      do 17 i=1,ni
         do 12 jk=1,njk
            if (bsq(jk,i).gt.bsqmax(i)) then
              bsqmax(i) = bsq(jk,i)
              jkmax(i)  = jk
            end if
 12      end do
!
         do 15 jk=1,njk
!            bsqav(i)    = bsqav(i) + bjac(jk,i) * bsq(jk,i)
            bnorm(jk,i) = sqrt(bsq(jk,i) / bsqmax(i))
 15      end do
 17   end do
      do 19 i=1,ni
!         bmax(i)  = sqrt(bsqmax(i))
!         bsqav(i) = - dtdp * bsqav(i) / vp(i)
          bsqav(i) = - sqgbsq(i) / vp(i)
 19   end do
!
! ... DETERMINE THE AVERAGE OF g1 REQUIRED FOR TRAPPED PARTICLE FRACTION
!
      do 29 i=1,ni
         do 27 lpitch=1,npitch
           g1av(lpitch,i)  = zero
           do 25 jk=1,njk
             g1av(lpitch,i) = g1av(lpitch,i)
     &            + bjac(jk,i) * sqrt(1.0 - pitch(lpitch) * bnorm(jk,i))
 25        end do
 27      end do
         do 28 lpitch=1,npitch
           g1av(lpitch,i) = - dtdp * g1av(lpitch,i) / vp(i)
 28      end do
 29   end do
!
! ... MAIN LOOP OVER THE FLUX SURFACES
      do 80 i=1,ni
!
! ... DETERMINE FLUX SURFACE AVERAGE OF g_2
!
! ... INITIALISE g_2 FOR EACH FLUX SURFACE
! ... INITIALISE sqrt(g)B.grad(sqrt(g)) AND sqrt(g)Bxgrad(s).grad(sqrt(g))
! ... WHICH ARE STORED IN ARRAYS bdgradb AND bdgradg, RESPECTIVELY.
!
        do 32 jk=1,njk
           g2(jk) = zero
           bdgradg(jk)=zero
           bdgradb(jk)=zero
 32     enddo
        do 36 l=1,lmnb 
          qlmn = mb(l)*fpp(i) - nb(l)*ftp(i)
          djpb = djpboo*(qlmn + fpp(i))*(qlmn + fpp(i))
          bdgrgmn = - qlmn * fbjac(l,i)
          bxgrgmn = (mb(l) * ci(i) - nb(l) * cj(i)) * fbjac(l,i)
!          g2mn = fbjac(l,i) * (mb(l)*ci(i) - nb(l)*cj(i)) * qlmn
!     &                      / (qlmn * qlmn + djpb)           
          g2mn = bxgrgmn * qlmn     / (qlmn * qlmn + djpb)           
            do 34 jk=1,njk
              g2(jk) = g2(jk) + g2mn * (tcos(jkmax(i),l) - tcos(jk,l))
              bdgradb(jk) = bdgradb(jk) + bdgrgmn * tsin(jk,l)
              bdgradg(jk) = bdgradg(jk) + bxgrgmn * tsin(jk,l)
 34        end do
 36     end do
! ... NOTE THAT AT THIS POINT THE ARRAY g2(jk) STORES sqrt(g)g_2/Phi'(s)
! ... DETERMINE b.grad(B)/B AND b.grad(g_2)/(2B^2)
! ... STORED IN ARRAYS bdglnb AND bdgradg, RESPECTIVELY. THE ARRAY bdgradb
! ... IS REDEFINED TO CONTAIN sqrt(g)b.grad(B) IN LOOP 38.
        ftpogb = 0.5*ftp(i)/(sqgbsq(i)*sqgbsq(i))
        do 38 jk=1,njk
          sqrtbsq = sqrt(bsq(jk,i))
          bdgradb(jk)=-0.5*bsq(jk,i)*bdgradb(jk)/sqgbsq(i)
          bdglnb(jk) = bdgradb(jk)*sqrtbsq/sqgbsq(i)
          bdgradg(jk)=ftpogb*(sqrtbsq*bdgradg(jk)
     &               + 2.0*bdglnb(jk)*g2(jk)*sqgbsq(i))
 38     end do
!
! ... FOURIER ANALYSIS TO GET [b.grad(B)/B]_mn AND [b.grad(g_2)/(2B^2)]_mn
!
          do 42 l=1,lmnb
            bgradbmn(l) = zero
            bgradgmn(l) = zero
                  do 41 jk=1,njk
                    bgradbmn(l) = bgradbmn(l) + bdglnb(jk) * tsin(jk,l)
                    bgradgmn(l) = bgradgmn(l) + bdgradg(jk)* tsin(jk,l)
 41               end do
 42       end do
!
! ... CORRECTLY NORMALISE AND DIVIDE THESE FOURIER TERMS BY  |m_bPsi'-n_bPhi'|
        do 43 l=1,lmnb
          qlmn = mb(l)*fpp(i) - nb(l)*ftp(i)
          djpb = djpboo*(qlmn + fpp(i))*(qlmn + fpp(i))
          qlmn = qlmn / (qlmn * qlmn + djpb)
          bgradbmn(l) = 2.0 * dtdp * bgradbmn(l) * qlmn
          bgradgmn(l) = 2.0 * dtdp * bgradgmn(l) * qlmn
 43     end do
!
! ... CONVERSION OF RESULTING FOURIER COMPONENTS TO REAL SPACE.
!
        do 45 jk=1,njk
          bdglnb(jk) =zero
          bdgradg(jk)=zero
 45    end do
        do 47 l=1,lmnb
          do 46 jk=1,njk
            bdglnb(jk)  = bdglnb(jk)  + bgradbmn(l) * tsin(jk,l)
            bdgradg(jk) = bdgradg(jk) + bgradgmn(l) * tsin(jk,l)
 46      end do
 47    end do
!
! ... DETERMINATION OF THE AVERAGE OF g_2
! ... EVALUATION OF SCALE LENGTH AND GEOMETRIC FACTOR IN PLATEAU REGIME
        do 49 jk=1,njk
          g2av(i)   = g2av(i)   + g2(jk)
          rlampl(i) = rlampl(i) + bdgradb(jk) * bdglnb(jk)
          gbplat(i) = gbplat(i) + bdgradb(jk) * bdgradg(jk)
 49    end do
        g2av(i)   = - dtdp * ftp(i) * g2av(i) / vp(i)
        rlampl(i) = - dtdp * rlampl(i) / vp(i)
        gbplat(i) = - dtdp * gbplat(i) / vp(i)
        gbplat(i) = g2av(i) - bsqav(i) * gbplat(i) / rlampl(i)
        rlampl(i) = - 2.0 * bsqav(i) / (sqrtpi * sqgbsq(i) * rlampl(i))
!
! ... DETERMINE FLUX SURFACE AVERAGE OF g4
!
        do 68 lpitch=1,npitch
          g4av(lpitch,i) = zero
!
! ... CALCULATE FOURIER SPECTRUM OF 1/g1
          do 52 l=1,lmnb
            g1mnin(l) = zero
              do 51 jk=1,njk
                g1mnin(l) = g1mnin(l) + tcos(jk,l)
     &                     / sqrt(1.0 - pitch(lpitch) * bnorm(jk,i))
 51           end do
 52       end do
! ... INITIALISE g4 FOR EACH FLUX SURFACE AND PITCH VALUE
          do 53 jk=1,njk
            g4(jk) = zero
 53       end do
! ... NORMALISE FOURIER SPECTRUM OF 1/g1 AND EVALUATE FOURIER SPECTRUM OF g4
          do 54 l=1,lmnb
            g1mnin(l) = 2. * dtdp * g1mnin(l)
            qlmn = mb(l)*fpp(i) - nb(l)*ftp(i)
            djpb = djpboo*(qlmn + fpp(i))*(qlmn + fpp(i))
            g4mn(l) = g1mnin(l) * ((mb(l)*ci(i) - nb(l)*cj(i)) * qlmn)
     &                            / (qlmn * qlmn + djpb)
 54       end do
! ... DETERMINE g4(s,theta,phi,pitch) (REPEATED FOR EACH s AND PITCH VALUE)
          do 56 l=1,lmnb
              do 55 jk=1,njk
             g4(jk) = g4(jk) + g4mn(l) * (tcos(jkmax(i),l) - tcos(jk,l))
 55           end do
 56        end do
! ... MULTIPLY g4 BY MISSING AMPLITUDE (Phi'(s)*v-parallel/v = Phi'(s) * g1)
! ... EVALUATE FLUX SURFACE AVERAGE OF g4
           do 59 jk=1,njk
             g4av(lpitch,i) = g4av(lpitch,i) + g4(jk) *
     &               bjac(jk,i) * sqrt(1.0- pitch(lpitch) * bnorm(jk,i))
 59        end do
          g4av(lpitch,i) = - ftp(i) * dtdp * g4av(lpitch,i) / vp(i)
!
 68     end do
!
! ... DETERMINE GEOMETRICAL FACTOR FOR BOOSTRAP CURRENT gb AND TRAPPED/
! ... CIRCULATING PARTICLE FRACTION
!
        do 70 lpitch=1,npitch
          pdpit = delpch * dpitdx(lpitch) * pitch(lpitch)
          gb(i) = gb(i) + pdpit * g4av(lpitch,i) / g1av(lpitch,i)
          fc(i) = fc(i) + pdpit / g1av(lpitch,i)
! ... OLD FORMS WITH UNIFORM PITCH ANGLE MESH
!          gb(i) = gb(i) + delpch * pitch(lpitch) * g4av(lpitch,i)
!     &                  / g1av(lpitch,i)
!          fc(i) = fc(i) + delpch * pitch(lpitch) / g1av(lpitch,i)
 70     end do
!
      write(93,113) i, jkmax(i), bsqmax(i), gb(i), g2av(i)/(1.0-fc(i))
      write(93,108)(bsq(jk,i),jk=1,njk)
      write(93,108)(g2(jk),jk=1,njk)
 80   end do
!      write(44,103)(g4av(lpitch,ni-1),lpitch=1,npitch)
!      write(44,103)(g4av(lpitch,ni),lpitch=1,npitch)
!      write(45,103)(g1av(lpitch,ni-1),lpitch=1,npitch)
!      write(45,103)(g1av(lpitch,ni),lpitch=1,npitch)
!      write(43,103)(fc(i),i=1,ni)
!      write(43,103)(gb(i),i=1,ni)
 103  format( 1p6e14.6)
! ... TRAPPED/PASSING PARTICLE FRACTIONS, 1/nu GEOMETRIC FACTOR AND
! ... NORMALISED COLLISION FREQUENCY
      do 85 i=1,ni
        fc(i) = 0.75 * bsqav(i) * fc(i) / bsqmax(i)
        ft(i) = 1.0 - fc(i)
        gb(i) = (g2av(i) - 0.75 * bsqav(i) * gb(i) / bsqmax(i)) / ft(i)
        densty= dn20*(1.0-0.5*cn*(s(i)+s(i-1))**lden)**qden
        rlame = 0.25e6*(1+zeff)*pth(i)/(zeff*densty)**3
        rnue(i)= 4.0 * ft(i) * rlampl(i) / (3.0*sqrtpi * fc(i) * rlame)
 85   end do
!      write(43,103)(fc(i),i=1,ni)
!      write(43,103)(gb(i),i=1,ni)
!      write(47,103)(g2av(i),i=1,ni)
!
! ... THE mu AMPLITUDES, Eqn.(3) of WATANABE et al., Nucl. Fusion 32(1992)1500
! ... THE TRAPPED TO CIRCULATING PARTICLE RATIO INCLUDED IN LOOP 90
      rmui1a = sqrt(2.0) - log(1.0+sqrt(2.0))
      rmue1a = rmui1a + zeff
      rmui2a = 2.0*sqrt(2.0) - 5.0*log(1.0+sqrt(2.0))/2.0
      rmue2a = rmui2a + 3.*zeff/2.0
      rmui3a = 39.0*sqrt(2.0)/8.0 - 25.0*log(1.0+sqrt(2.0))/4.0
      rmue3a = rmui3a + 13.0*zeff/4.0
!
! ... THE L^ij_IJ COEFICIENTS
      rlee11 = zeff
      rlee12 = 1.5 * zeff
      rlee22 = sqrt(2.0) + 13.0 * zeff / 4.
      rlii22 = sqrt(2.0)
!
! ... THE BOOTSTRAP COEFFICIENTS, Eqn.(2) of WATANABE et al.
! ... THE BOOTSTRAP CURRENT Eqn.(1) and Eqn.(4)
      boot(0) = zero
      boos(0) = zero
      do 90 i=1,ni
        dboot(i)=(rlee11+ft(i)*rmue1a/fc(i))*(rlee22+ft(i)*rmue3a/fc(i))
     &          -(rlee12+ft(i)*rmue2a/fc(i))*(rlee12+ft(i)*rmue2a/fc(i))
        rl31(i) =(rmue1a * (rlee22 + ft(i) * rmue3a / fc(i))
     &          - rmue2a * (rlee12 + ft(i) * rmue2a / fc(i))) * ft(i)
     &          / (dboot(i) * fc(i))
        rl32e(i)=(rmue3a*rlee12 - rmue2a*rlee22) *ft(i)/(dboot(i)*fc(i))
        rl32i(i)=-rl31(i) * rmui2a * rlii22
     &   /(rmui1a*(rlii22+ft(i)*rmui3a/fc(i))-rmui2a*rmui2a*ft(i)/fc(i))
        bjav(i)= - gb(i)*(rl31(i)*pp(i) + (zeff*rl32e(i)+rl32i(i))
     &                                  * tempp(i)/(1.0+zeff)) / ftp(i)
        boot(i) = boot(i-1)+ftpsgn*ftp(i)*(s(i)-s(i-1))*bjav(i)/bsqav(i)
        ftsq    = ft(i) * ft(i)
        ftcb    = ftsq  * ft(i)
        ftqt    = ftcb  * ft(i)
        osl31   =  (1.0+1.4/(zeff+1.0))*ft(i) - 1.9*ftsq/(zeff+1.0)
     &           + 0.3*ftcb/(zeff+1.0)        + 0.2*ftqt/(zeff+1.0)
        f32ee   =  (0.05+0.62*zeff)*(ft(i)-ftqt)/(zeff*(1.0+0.44*zeff))
     &           + (ftsq-ftqt-1.2*(ftcb-ftqt))/(1.0+0.22*zeff)
     &           + 1.2*ftqt/(1.0+0.5*zeff)
        f32ei   =- (0.56+1.93*zeff)*(ft(i)-ftqt)/(zeff*(1.0+0.44*zeff))
     &           + 4.95 * (ftsq-ftqt-0.55*(ftcb-ftqt))/(1.0+2.48*zeff)
     &           - 1.2*ftqt/(1.0+0.5*zeff)
        osl32   =  f32ee + f32ei
        alpha0  =- 1.17* (1.0 - ft(i)) / (1.0 - 0.22*ft(i) - 0.19*ftsq)
        bjos(i) = ci(i)*(osl31*pp(i) + zeff*(osl32+osl31*alpha0)
     &                               * tempp(i)/(1.0+zeff)) / fpp(i)
        boos(i) = boos(i-1)+ftpsgn*ftp(i)*(s(i)-s(i-1))*bjos(i)/bsqav(i)
 90   end do
!      
! ...  STORE BOOTSTRAP CURRENT ON THE HALF GRID IN ARRAY dboot
       write(6,102)
 102   format(/,'  I',7x,'S',7x,'<J_boot.B>',3x,'<J_b.B>_tok',1x,
     & 'GEOM FACTR GB','-I(s)/IOTA(s)',' PLATEAU GB','  TRAPPED FRAC.',
     & ' COLL. FREQ.')
!     & 'GEOMETRIC FACTR',' Tok. Geo. Fac.',' TRAPPED FRAC.')
       write(41,104)
 104   format(' I ',10x,'S',9x,'2pi J_boot-toroidal -- 2pi J_boot-prime'
     &,5x,'<j_boot.B>',4x,'GEOMETRIC FACTOR GB',3x,'TRAPPED FRACTION')
       write(46,111)
      do 95 i=1,ni
       dboot(i) = 0.5 * (boot(i) + boot(i-1))
       dboos(i) = 0.5 * (boos(i) + boos(i-1))
       shalf = 0.5 * (s(i) + s(i-1))
       write(41,105) i, shalf, dboot(i), ftpsgn*ftp(i)*bjav(i)/bsqav(i)
     &                , bjav(i) ,gb(i), ft(i)
!       write(42,106) i, shalf, pp(i), tempp(i),-ci(i)*ftp(i)/fpp(i)
!     &                , gb(i), fc(i)
       write(42,110)i, shalf, pp(i) ,tempp(i), gbplat(i), gb(i)
     &               , g2av(i)/ft(i),rnue(i)            , rlampl(i)
!       write(44,106)i, shalf, pp(i), tempp(i), ft(i), dboot(i), dboos(i)
!       write(44,106)i, shalf,bjav(i)/gb(i),bjos(i)*fpp(i)/(ftp(i)*ci(i))
       write(44,106)i, shalf,bjav(i),bjos(i), ft(i), dboot(i), dboos(i)
       write(6,110) i, shalf,bjav(i),bjos(i),gb(i),-ci(i)*ftp(i)/fpp(i)
     &              ,gbplat(i),ft(i),rnue(i)
       write(45,113)i, jkmax(i), bsqmax(i), gb(i), g2av(i)/ft(i)
       write(46,112)i, shalf, rl31(i),rl32e(i),rl32i(i)
     &         ,ft(i)/fc(i)*rmue1a,ft(i)/fc(i)*rmui1a,ft(i)/fc(i)*rmue2a
     &         ,ft(i)/fc(i)*rmui2a,ft(i)/fc(i)*rmue3a,ft(i)/fc(i)*rmui3a
 95   end do
!       print *, '  BOOTSTRAP CURRENT = ', boot(ni)
       bsqavax = 1.5 * bsqav(1) - 0.5 * bsqav(2)
       betax   = 2. * pax / bsqavax
       write(6,107) boot(ni), dn20, tonax, pax, bsqavax, betax
! ...  READ PREVIOUSLY COMPUTED TOROIDAL BOOTSTRAP CURRENT
       read(43,108) (dboos(i),i=1,ni)
       read(43,109) jbiter
       jbiter = min(nine,jbiter)
! ... COMBINE A FRACTION OF PREVIOUS BOOTSTRAP CURRENT WITH NEW
! ... BOOTSTRAP CURRENT FOR INPUT TO VMEC
       do i = 1,ni
        dboot(i) = (jbiter * dboos(i) + dboot(i))/(jbiter+1)
       end do
       jbiter = jbiter + 1
       rewind(43)
       write(43,108)(dboot(i),i=1,ni)
       write(43,109) jbiter
 105   format(i3,1p6e20.12)
 106   format(i3,1p6e14.6)
 107   format(' BOOTSTRAP CURRENT: = ',1p1e20.12,' N_20 = ',1p1e12.4,
     & ' T(0) = ',1p1e12.4,' mu_0 p(0) = ',1p1e12.4,/' B^2(0) = ',
     & 1p1e12.4,' beta ON AXIS = ',1p1e12.4) 

 108   format(1p5e16.8)
 109   format(i3)
 110   format(i3,1p8e13.5)
 111   format(' I ',5x,'S',9x,'RL31',6x,'RL32e',6x,'RL32i',7x,'MUe1',7x
     &        ,'MUi1',7x,'MUe2',7x,'MUi2',7x,'MUe3',7x,'MUi3')
 112   format(i3,1p10e11.3)
 113   format(2i5,1p3e14.6)
!--------0---------0---------0---------0---------0---------0---------0-c
      return
      end subroutine bootsj
