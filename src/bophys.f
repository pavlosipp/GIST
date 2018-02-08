!--------0---------0---------0---------0---------0---------0---------0-c
!  23.09.88          LAST MODIFICATION 28.11.05  UHS: BOPHYST DATA D
!-----------------------------------------------------------------------
!
      subroutine bophys (ni, njk, lmnb, lmnb0, mm, nmin, nmax, lcurrf,
     &lvmtpr, lpress, lmetpr, rplmin, djp, dtdp, s, ftp, fpp, ftpp,
     &fppp, ci, cj, cip, cjp, pp, wmag, vp, cipi, cjpi, pvpi, pth, ppi,
     &equi, civ, cjv, vvp, mb, nb, tsin, tcos, lfrz, gttl, gtpl, gppl,
     &gssu, bjac, bjacs, vjac, bs, bp, bt, bsq, fr, fbjac, fbs, fphv,
     &fbsq)


      use specifications_mod, ONLY:vmec_file, table_tag
!
!-----------------------------------------------------------------------
!     RECONSTRUCTION OF EQUILIBRIUM IN BC
!-----------------------------------------------------------------------
!
!.. Implicits ..
      implicit none  
!
!.. Formal Arguments ..
      integer :: njk, nmin, nmax, lmnb, ni, mm, lcurrf, lvmtpr, lpress,
     &lmetpr, mb ( * ), nb ( * ), lfrz (0:36, nmin:nmax)
      real :: rplmin, djp, dtdp, s (0: * ), ftp ( * ), fpp ( * ),
     &ftpp ( * ), fppp ( * ), ci ( * ), cj ( * ), cip ( * ), cjp ( * ),
     &pp ( * ), wmag ( * ), vp ( * ), cipi ( * ), cjpi ( * ), pvpi ( * )
     &, pth ( * ), ppi ( * ), equi ( * ), civ ( * ), cjv ( * ), vvp ( * 
     &), tsin (njk, * ), tcos (njk, * ), gttl (njk, 0: * ), gtpl (njk,
     &0: * ), gppl (njk, 0: * ), gssu (njk, 0: * ), bjac (njk, 0: * ),
     &bjacs (njk, 0: * ), vjac (njk, 0: * ), bs (njk, 0: * ), bp (njk,
     &0: * ), bt (njk, 0: * ), bsq (njk, 0: * ), fr (lmnb, * ),
     &fbjac (lmnb, 0: * ), fbs (lmnb, 0: * )
      real :: fphv (lmnb, * ), fbsq (lmnb, * )  
!
!.. Local Scalars ..
      integer :: I, JK, L, lmnb0, m, mboo, n, lx1, lx2, lx3, lx4, lx5,
     &lmx1, lmx2, lmx3, lmx4, lmx5, lnx1, lnx2, lnx3, lnx4, lnx5
      real :: DIFJAC, dvol, EQUIN, FPPI, FTPI, PMN, t7, t8, tds, tsi,
     &VPI, b2max, b2min, ft1, ft2, vp0, b0r0, t1, blarge, djpboo
!      parameter (blarge=1.e300)
!
!.. Intrinsic Functions ..
      intrinsic ABS , MIN, MAX 
!
      data  djpboo, blarge  / 1.0e-16, 1.0e300 /
!
! ... Executable Statements ...
!
      djpboo = max(djp,djpboo)
!
 !     if (lcurrf.eq.1) write (16, 2000)  
 2000 format(/'************** ci/j from bp/t *******************')  
 !     if (lcurrf.eq.2) write (16, 2001)  
 2001 format(/'************** ci/j from vmec *******************')  
 !     write (6, 2025)  
 2025 format   (/,14x,'s',14x,'ci',14x,'cj')  
      do 115 i = 1, ni  
!
       ci (i) = 0.  
       cj (i) = 0.  
       vp (i) = 0.  
       pvpi (i) = 0.  
       wmag (i) = 0.  
       dvol = (s (i) - s (i - 1) ) * dtdp  
!
      do 110 jk = 1, njk  
!
         bp (jk, i) = (gppl (jk, i) * ftp (i) + gtpl (jk, i) * fpp (i) )
!
         bt (jk, i) = (gtpl (jk, i) * ftp (i) + gttl (jk, i) * fpp (i) )
!
!     FBSQ(L,I) AND BSQ(JK,I) IS COMPUTED DIRECTLY IN VMTOBO (LIKE R/Z)
!
!          BSQ(JK,I) = (  GPPL(JK,I) * FTP(I)**2
!    $                +2.*GTPL(JK,I) * FPP(I) * FTP(I)
!    $                  + GTTL(JK,I) * FPP(I)**2  ) / BJAC(JK,I)
!
      if (lcurrf.eq.1) then  
       ci (i) = ci (i) - dtdp * bp (jk, i)  
       cj (i) = cj (i) + dtdp * bt (jk, i)  
      endif  
      if (lcurrf.eq.2) then  
       ci (i) = civ (i)  
       cj (i) = cjv (i)  
      endif  
       vp (i) = vp (i) - dtdp * bjac (jk, i)  
       pvpi (i) = pvpi (i) - dtdp * vjac (jk, i)  
!
  110   end do  
!
      t1 = 0.5 * (s (i) + s (i - 1) )  
!      write (6, 2026) i, t1, ci (i), cj (i)  
 2026 format     (i6,1p3e15.6)  
!
        do 112 jk = 1, njk  
          wmag (i) = wmag (i) + 0.5 * bsq (jk, i) * bjac (jk, i) * dvol
!
  112   end do
  115 end do
!
!      write (6,1033)
!      write (16, 1033)  
 1033 format ("differential volumes from geometric and magnetic boozer j
     &acobians")
      do 122 i = 1, ni  
!      write (6,1034) pvpi(i),vp(i)
!       write (16, 1034) pvpi (i), vp (i)  
  122 end do  
 1034 format (23x,1p1e13.6,5x,1p1e13.6)  
!
      mboo = 0  
      do 124 m = 0, mm  
       do 123 n = nmin, nmax  
         mboo = mboo + lfrz (m, n)  
  123  end do
  124 end do
!      write (16, 1036) rplmin, mboo  
 1036 format(/"modes in the boozer table in which r,z,phi fourier amplit
     &ude exceeds:",/25x,1p1e13.5,3x,i4,/)
!      do 125 n = nmin, nmax  
!  125 write (16, 1037) (lfrz (m, n), m = 0, 36), n  
! 1037 format(5x,31i2,i3)  
!
      do 129 i = 1, ni - 1  
!
       tsi = 2.0 / (s (i + 1) - s (i - 1) )  
       vpi = 0.5 * (vp (i + 1) + vp (i) )  
       ftpi = 0.5 * (ftp (i + 1) + ftp (i) )  
       fppi = 0.5 * (fpp (i + 1) + fpp (i) )  
       cipi (i) = tsi * (ci (i + 1) - ci (i) )  
       cjpi (i) = tsi * (cj (i + 1) - cj (i) )  
       pvpi (i) = tsi * (pth (i + 1) - pth (i) )  
       equin = abs (cjpi (i) * fppi) + abs (cipi (i) * ftpi) + abs (
     &   vpi * pvpi (i) )
       equi (i) = ( (cjpi (i) * fppi - cipi (i) * ftpi) - vpi * pvpi (
     &   i) ) / equin
  129 ppi (i) = (cjpi (i) * fppi - cipi (i) * ftpi) / vpi  
      pvpi (ni) = 0.  
!
!      write (6, 1002)  
!     WRITE(16,1002)
      do 130 i = 1, ni  
!      write (6, 1003) i, vp (i), ppi (i), pvpi (i), cjpi (i), cipi (i)
!     &   , equi (i)
!        WRITE(16,1003)I,VP(I),PPI(I),PVPI(I),CJPI(I),CIPI(I),EQUI(I)
  130 end do  
 1002 format(//1x,'***** bophys, reconstructed equilibrium *****',//3x
     &          ,'i     vp(i)    ppi(i)     pvpi(i)    cjpi(i)    cipi(i
     &)     equi',/)
!-----------------------------------------------------------------------
 1003 format(1x,i3,1p6e11.4)  
!
!     INTER/EXTRAPOLATION TO HALF GRID
      if (lpress.eq.2) then  
        do 135 i = 1, ni - 1  
  135   ppi (i) = pvpi (i)  
      endif  
!
      do 145 i = 2, ni - 1  
!
        tds = 2.0 / (s (i + 1) + s (i) - s (i - 1) - s (i - 2) )  
        fppp (i) = tds * (fpp (i + 1) - fpp (i - 1) )  
        ftpp (i) = tds * (ftp (i + 1) - ftp (i - 1) )  
        pp (i) = 0.5 * (ppi (i) + ppi (i - 1) )  
        cip (i) = 0.5 * (cipi (i) + cipi (i - 1) )  
        cjp (i) = 0.5 * (cjpi (i) + cjpi (i - 1) )  
!
         do 140 jk = 1, njk  
          bjacs (jk, i) = tds * (bjac (jk, i + 1) - bjac (jk, i - 1) )
  140    end do
  145 end do
!
      b0r0 = fbsq (lmnb0, ni) / (fr (lmnb0, ni) * fr (lmnb0, ni) )  
      vp0 = 1.5 * vp (1) - 0.5 * vp (2)  
!
! --- OUTPUT FOR BALLOONING CALCULATION

      OPEN(25,file='../int/tpr_'//trim(vmec_file)//'_'//trim(table_tag), 
     .          action="READWRITE",STATUS="REPLACE")
      write (25, 1655) ni, lmnb, b0r0, vp0  
      write (25, 1656) (mb (l), l = 1, lmnb), (nb (l), l = 1, lmnb)  
 1655 format (2i5,1p2e15.6)  
 1656 format (12i6)  
      CLOSE(25)
      write (6, 889)  
  889 format (2x,"i",8x,"q(si)",6x,"mercier criterion"/)  
!
      do 200 i = 2, ni - 1  
!
        do 205 jk = 1, njk  
  205   vjac (jk, i) = (fpp (i) * cj (i) - ftp (i) * ci (i) ) / bsq (jk,
     &   i)
!
      do 210 l = 1, lmnb  
!
!bjac         fbjac(l,i) = 0.
        fphv (l, i) = 0.  
        fbs (l, i) = 0.  
!
!
      do 220 jk = 1, njk  
!bjac         fbjac(l,i) = fbjac(l,i)
!bjac     $              + bjac(jk,i) *tcos(jk,l)
          fphv (l, i) = fphv (l, i) + vjac (jk, i) * tcos (jk, l)  
  220     end do  
!
!bjac      fbjac(l,i) = 2. * dtdp * fbjac(l,i)
        fphv (l, i) = 2. * dtdp * fphv (l, i)  
        if (mb (l) .eq.0.and.nb (l) .eq.0) fphv (l, i) = 0.5 * fphv (
     &     l, i)
!
!     COMPARE (DISCRETIZED) BJAC COEFFICIENTS WITH PHV(=BJAC) COMPUTED
!     BY DIRECT INTEGRATION OVER BSQ
!
        difjac = 2. * (fbjac (l, i) - fphv (l, i) ) / (vp (i) + vvp (
     &     i) )
!
        if (abs (difjac) .ge.1.e-06.and.lvmtpr.eq.1) then  
!      write (16, 1006) i, mb (l), nb (l), vp (i), vvp (i),
!     &       fbjac (l, i), fphv (l, i), difjac
 1006 format(1x,3i4,1p2e20.8,1p3e10.2)  
      endif  
!
      if (mb (l) .ne.0.or.nb (l) .ne.0) then  
       pmn = mb (l) * fpp (i) - nb (l) * ftp (i)  
       pmn = pmn / (pmn*pmn + djpboo * (pmn+fpp(i)) * (pmn+fpp(i)))
       fbs (l, i) = pp (i) * fbjac (l, i) * pmn  
!       fbs (l, i) = pp (i) * fbjac (l, i) * pmn / (pmn * pmn + djp)  
!perf  fbs(l,i) =  fbs(l,i) / tpi
!      WRITE(16,1234) I,MB(L),NB(L),FBS(L,I), -NB(L)*FTP(I)+MB(L)*FPP(I)
!    $ ,AIOTA(I)
! 1234  format(1x,'fbset',3i4,1p3e15.7)
      endif  
  210   end do  
!
!
      t7 = (cj (i) * fpp (i) - ci (i) * ftp (i) ) / (ftp (i) * ftp (i) )
      t8 = cj (i) * cj (i) / (ftp (i) * ftp (i) )  
      do 225 jk = 1, njk  
        gssu (jk, i) = t7 * gttl (jk, i) - t8  
  225   bs (jk, i) = 0.  
!
      do 230 l = 1, lmnb  
        do 230 jk = 1, njk  
  230   bs (jk, i) = bs (jk, i) + fbs (l, i) * tsin (jk, l)  
!
!      WRITE (16,1004)
!1004 format(/1x,'fbs(1...4),cjp,cip'/)
!      WRITE (16,1005) I,(FBS(L,I),L=1,4),CJP(I),CIP(I)
!1005 format(i4,1p4e11.4,1x,1p2e11.4)
!
  200 end do  
!
!     CALL RADPLO(NI,   CI,'   CI$',1,0)
!     CALL RADPLO(NI,  CIP,'  CIP$',1,0)
!     CALL RADPLO(NI,   CJ,'   CJ$',1,0)
!     CALL RADPLO(NI,  CJP,'  CJP$',1,0)
!     CALL RADPLO(NI,  FPP,'  FPP$',1,0)
!     CALL RADPLO(NI, FPPP,' FPPP$',1,0)
!     CALL RADPLO(NI,  FBS,'   VP$',1,0)
!      do 290 jk = 1,njk
!  290    bs(jk,1) = 2.* bs(jk,2) - bs(jk,3)
!
!     WRITE (16,1007)
!1007 format(/'  i, fppp, ftpp, cip, cjp, pp',/)
!     DO 300 I = 1,6
! 300    WRITE (16,1008) I,FPPP(I),FTPP(I),CIP(I),CJP(I),PP(I)
!1008 format(1x,i3,1p5e14.5)
!     WRITE (16,1009)
!1009 format(/'  i, bjac(1...5)              '/)
!     DO 310 I = NI/3-2,NI/3+2
! 310    WRITE (16,1008) I,(BJAC(JK,I),JK=1,5)
!     WRITE (16,1010)
!1010 format(/'  i,bjacs(1...5)              '/)
!     DO 320 I = NI/3-2,NI/3+2
! 320    WRITE (16,1008)I,(BJACS(JK,I),JK=1,5)
!     WRITE (16,1011)
!1011 format(/'  i, gttl(1...5)              '/)
!     DO 330 I = 1,6
! 330    WRITE (16,1008)I,( GTTL(JK,I),JK=1,5)
!     WRITE (16,1012)
!1012 format(/'  i, gssl(1...5)              '/)
!     DO 340 I = 1,6
! 340    WRITE (16,1008)I,( GSSL(JK,I),JK=1,5)
!     WRITE (16,1013)
!1013 format(/'  i,   bs(1...5)              '/)
!     DO 350 I = 1,6
! 350    WRITE (16,1008)I,(   BS(JK,I),JK=1,5)
!     WRITE (16,1014)
!1014 format(/'  i, gstl(1...5)              '/)
!     DO 360 I = 1,6
! 360    WRITE (16,1008)I,( GSTL(JK,I),JK=1,5)
!     WRITE (16,1015)
!1015 format(/'  i,fbjac(1...5)              '/)
!     DO 370 I = 1,5
! 370    WRITE (16,1008)I,(FBJAC( L,I), L=1,5)
!     WRITE (16,1016)
!1016 format(/'  i,fbs  (1...5)              '/)
!     DO 380 I = 1,5
! 380    WRITE (16,1008)I,(FBS  ( L,I), L=1,5)
!
!...calculate and print main 4 fourier amplitudes of b**2
!...and trapped particle fraction.
      if (lmetpr.ne.4) return  
      b2max = - blarge  
      do 421 i = 1, ni  
       do 421 l = 1, lmnb  
        if (abs (fbsq (l, i) ) .le.b2max) goto 421  
        lx1 = l  
        lmx1 = mb (l)  
        lnx1 = nb (l)  
        b2max = abs (fbsq (l, i) )  
  421 continue  
      b2max = - blarge  
      do 422 i = 1, ni  
       do 422 l = 1, lmnb  
        if (l.eq.lx1) goto 422  
        if (abs (fbsq (l, i) ) .le.b2max) goto 422  
        lx2 = l  
        lmx2 = mb (l)  
        lnx2 = nb (l)  
        b2max = abs (fbsq (l, i) )  
  422 continue  
      b2max = - blarge  
      do 423 i = 1, ni  
       do 423 l = 1, lmnb  
        if (l.eq.lx1.or.l.eq.lx2) goto 423  
        if (abs (fbsq (l, i) ) .le.b2max) goto 423  
        lx3 = l  
        lmx3 = mb (l)  
        lnx3 = nb (l)  
        b2max = abs (fbsq (l, i) )  
  423 continue  
      b2max = - blarge  
      do 424 i = 1, ni  
       do 424 l = 1, lmnb  
        if (l.eq.lx1.or.l.eq.lx2.or.l.eq.lx3) goto 424  
        if (abs (fbsq (l, i) ) .le.b2max) goto 424  
        lx4 = l  
        lmx4 = mb (l)  
        lnx4 = nb (l)  
        b2max = abs (fbsq (l, i) )  
  424 continue  
      b2max = - blarge  
      do 425 i = 1, ni  
       do 425 l = 1, lmnb  
        if (l.eq.lx1.or.l.eq.lx2.or.l.eq.lx3.or.l.eq.lx4) goto 425  
        if (abs (fbsq (l, i) ) .le.b2max) goto 425  
        lx5 = l  
        lmx5 = mb (l)  
        lnx5 = nb (l)  
        b2max = abs (fbsq (l, i) )  
  425 continue  
!
!...trapped particle fraction
      write (26, 1456) ni, lmx1, lnx1, lmx2, lnx2, lmx3, lnx3, lmx4,
     & lnx4, lmx5, lnx5
 1456 format(1x,11i6)  
      do 429 i = 1, ni  
       b2max = - blarge
       b2min = blarge
       do 426 jk = 1, njk  
        b2max = MAX (bsq (jk, i), b2max)  
  426   b2min = MIN (bsq (jk, i), b2min)  
      ft1 = sqrt (sqrt (b2max) / sqrt (b2min) - 1.)  
      ft2 = sqrt (b2min) / sqrt (b2max)  
      write (26, 1457) fbsq (lx1, i), fbsq (lx2, i), fbsq (lx3, i),
     &   fbsq (lx4, i), fbsq (lx5, i), ft2
  429 end do  
 1457 format(1x,1p6e11.4)  
!
      return  
      end subroutine bophys
