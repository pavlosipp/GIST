!-----------------------------------------------------------------------
!  15.12.89        LAST MODIFICATION 10.02.93     UHS: vforce DATA D
!-----------------------------------------------------------------------
!
      subroutine vforce (ni, njk, lmnb, nper, lnbaln, lnbalx, mm,
     &lmxbal, nmin, nmax, lskewn, s, vvp, ftp, fpp, civ, cjv, civp,
     &cjvp, pth, pvpi, equiv, pvp, mb, nb, pol, tor, tcos, tsin, lfrz)
!
!     CALCULATE THE EQUILIBRIUM FORCE BALANCE IN VMEC COORDINATES.
!     PRINT THESE AND OTHER FLUX SURFACE QUANTITIES.
!

      real :: s (0:*), vvp (*), ftp (*), fpp (*), civ (*), cjv (*),
     &civp (*), cjvp (*), pth (*), pvpi (*), equiv (*), pvp (*),
     &pol (*), tor (*), tcos (njk, *), tsin (njk, *)
      integer :: ni, njkm, lmnb, nper, lnbaln, lnbalx, mm, lmxbal, nmin,
     &nmax, lskewn, lfrz (0:36, nmin:nmax), mb (*), nb (*)
!
      do 320 i = 1, ni - 1  
!
       tsi = 2. / (s (i + 1) - s (i - 1) )  
       vvpi = 0.5 * (vvp (i + 1) + vvp (i) )  
       ftpi = 0.5 * (ftp (i + 1) + ftp (i) )  
       fppi = 0.5 * (fpp (i + 1) + fpp (i) )  
       civp (i) = tsi * (civ (i + 1) - civ (i) )  
       cjvp (i) = tsi * (cjv (i + 1) - cjv (i) )  
       pvpi (i) = tsi * (pth (i + 1) - pth (i) )  
       equin = abs (cjvp (i) * fppi) + abs (civp (i) * ftpi) + abs (
     &    vvpi * pvpi (i) )
       equiv (i) = ( (cjvp (i) * fppi - civp (i) * ftpi) - vvpi *
     &    pvpi (i) ) / equin
  320 pvp (i) = (cjvp (i) * fppi - civp (i) * ftpi) / vvpi  
!
C      write (6, 1002)  
!      write (16, 1002)  
      do 330 i = 1, ni  
C        write (6, 1003) i, vvp (i), pvp (i), pvpi (i), cjvp (i),
C     &    civp (i), equiv (i)
!        write (16, 1003) i, vvp (i), pvp (i), pvpi (i), cjvp (i),
!     &    civp (i), equiv (i)
  330 end do  
 1002 format(//1x,'***** lgikvm, reconstructed equilibrium *****',//3x
     &          ,'i',11x,'vvp(i)',7x,'pvp(i)',4x,'pvpi(i)',4x,'cjvp(i)'
     &,4x,'civp(i)',5x,'equi',/)
 1003 format(1x,i3,1pe20.9,1p5e11.4)  
!
!     COMPUTE COS/SIN FOR BC SET
!
!     extra modes for ballooning
      nminb = lnbaln  
      nmaxb = lnbalx  
      do 112 m = mm + 1, lmxbal  
        if (lskewn.gt.0) then  
          lskewn = min0 (lskewn, ni)  
          rmcent = m * fpp (lskewn) / (ftp (lskewn) * nper)  
          mcent = nint (rmcent)  
          nminb = mcent + lnbaln  
          nmaxb = mcent + lnbalx  
        endif  
      do 112 n = nminb, nmaxb  
         lmnb = lmnb + 1  
         mb (lmnb) = m  
  112    nb (lmnb) = n * nper  
!
      do 130 l = 1, lmnb  
       do 130 jk = 1, njk  
         tcos (jk, l) = cos (mb (l) * pol (jk) - nb (l) * tor (jk) )  
  130    tsin (jk, l) = sin (mb (l) * pol (jk) - nb (l) * tor (jk) )  
!
!     reinitialise lfrz array
!
      do 140 m = 0, 36  
       do 140 n = nmin, nmax  
  140 lfrz (m, n) = 0  
!
      return  
      end subroutine vforce
