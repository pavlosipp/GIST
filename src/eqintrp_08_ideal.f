        subroutine EQINTRP_08_i

      USE specifications_mod, ONLY:ns_vmec,mpnt_vmec,vmec_6_90,
     &vmec_file,table_tag, vmec_dir

C      USE strings
C
!..  INTERPOLATES OUTPUT FROM VMEC UNIT=8 TO DIFFERENT RADIAL MESH SIZE
!..  USING A SPLINE FIT AND PREPARES TERPSICHORE INPUT IN UNIT=18
!..  VERSION ADAPTED FOR 'ANIMEC'
!
!      integer, parameter :: rprec = SELECTED_REAL_KIND(12,100)
!      integer, parameter :: dp = rprec
!      integer, parameter :: NNS=197, NRHO=197, MNMAX=399
!      integer, parameter :: NDS=max(NNS,NRHO)
      INTEGER, PARAMETER :: IBASPL=73
      REAL, PARAMETER:: mu_0=4*3.1415926e-7
C
C.. Local Scalars ..
      character(30) :: dumchar,vmver
      character(60) :: mgrid_file
      integer :: nfp, mpol, ntor, kpres
     &          ,it0, it1, it2, it3, it4, it5
     &          ,J  , MN , i, iasym,ns,mpnt

      real ::    wb, wp, gamma, ofp, helsym, xm, xn, rmu0, fourpi         
     &          ,pt1, pt2, pt3, pwr,dsc,dsf,shd,yp1, ypn
 

      REAL, DIMENSION(:,:), ALLOCATABLE :: xrc,xzc,xjc,
     &                                     xrf,xzf,xjf,bcf,bcc 

      REAL, DIMENSION(:), ALLOCATABLE :: aiota, amass, pres, PHIP, VP,
     &               sc,sf,rg,f0,shalf,f1,f2,f3,shalm


      INTEGER, DIMENSION(:), ALLOCATABLE :: ixm,ixn

C.. Intrinsic Functions ..
      intrinsic nint, min, max, ATAN
!
!.. EXTERNAL FUNCTIONS
      external spline1, seval

!.. DATA
      data yp1, ypn / 1.0e30, 1.0e30 /
!
!        VOLI=0.
        fourpi = 16 * ATAN(1.0)
        rmu0 = fourpi * 1.0e-7
        GAMMA =0
        OFP=0
        MPNT=0
        MPOL=0
        NTOR=0
        helsym = 0
        kpres = 0
!        ns = nns
        

CPKM------------------------------------------------------------------
        OPEN(8,file=trim(vmec_dir)//'/'//trim(vmec_file),action="read", 
     &       err=931)
        OPEN(18,file='../int/eqin_'//trim(vmec_file)//'_' 
     &        //trim(table_tag),action="write")

!=================  VMEC VERSION ==================================

        read (8,*) dumchar,dumchar,dumchar,vmver !read version
        
        read (8,*) wb, wp, gamma, pt1, pt2, pt3, pt4 !pax added pt4

        IF (vmver == '8.00') then 

        read (8,*) nfp, ns_vmec, mpol, ntor, mpnt_vmec, it1, it2, iasym,
     1              it4,it5 !pax added it5
        
        ELSE

!  Remove trailing dummy variables; a file from Nunami is too short
        read (8,*) nfp, ns_vmec, mpol, ntor, mpnt_vmec, it0, it1, 
     &                it2, iasym  ! iasym is the last needed input
! original     &      it2, iasym,it4!,it5 !pax added it5, extra it0

        ENDIF


! shorten next read since all elements are dummies
        read (8,*) it1
!Old:        read (8,*) it1, it2, it3, it4, it5, it6 !pax added it6
      read (8,'(a)') mgrid_file
        print *, ' VMEC_file = ', trim(vmec_dir)//'/'//trim(vmec_file)          !mgrid_file
        print *, ' mgrid file = ',mgrid_file
        
        ns=ns_vmec ; mpnt=mpnt_vmec

        ALLOCATE(xrc(MPNT,NS) , xzc(MPNT,NS) , xjc(MPNT,NS)
     &          ,xrf(NS,MPNT)  , xzf(NS,MPNT)  , xjf(NS,MPNT)
     &          ,bcf(NS,MPNT)  , bcc(MPNT,NS))


        ALLOCATE(aiota(NS), amass(NS), pres(NS), PHIP(NS), VP(NS)
     &          ,sc(NS),rg(NS+1),f0(NS+1),shalf(NS+1),sf(NS)
     &          ,f1(NS+1),f2(NS+1),f3(NS+1),shalm(NS-1))


       ALLOCATE(ixm(MPNT),ixn(MPNT))

        NSF = NS

!        if (ns .gt. NDS) then
!        print *, 'STOP (eqintrp_08_i): Increase NNS,NRHO to', ns
!        stop
!        endif
!        if (ns .gt. NRHO) then
!        print *, '  DIMENSIONS ARE INCORRECT'
!        print *, '  INCREASE NRHO TO: ',ns
!        stop
!        endif
!        if (mpnt .gt. MNMAX) then
!        print *, '  DIMENSIONS ARE INCORRECT'
!        print *, '  INCREASE MNMAX TO: ',mpnt
!        stop
!        endif
!        READ (8)     GAM,NFP,NRHO,MPOL,NTOR,MPNT,IA,IB
        ofp  = 1./REAL(NFP)
!        NPERIODS = NFP
CPKME-----------------------------------------------------------------
        WRITE(18,840)WB,GAMMA,ofp,helsym,mpnt,NSF,MPOL,
     &                   NTOR,1,1,IT1,IT2,kpres
 840        FORMAT (1X,1PE22.12,1X,3(1PE12.5),7I4,I9,I2)
! Original format brok for NFIS case with large NITER (IT2)
! Orig  840        FORMAT(1x,1pe22.12,1x,1p3e12.5,8i4,i1)
! Pathscale generic fields need separating commas
! 840        FORMAT(13(G,', ')) ! Intel adds blanks
! Output the results to compare with the read values;
!   use the variable names from eqinvm.f:
        write(*,'(" >>> EQINTRP_08_IDEAL is writing:")')
        write(*,'("VOLI, GAM, ENM1, alpha=",4E11.4)')
     x   WB,GAMMA,ofp,helsym
        write(*,'("MPNT, NRHO, MPOL, NTOR, NTOR0, MN0,",
     &  " NIT, IFSQ, KPRES=",/,7I4,I9,I2)')
     x   mpnt,NSF,MPOL,NTOR,1,1,IT1,IT2,kpres
c
  720 format(2i12)
!  730 format(3d24.15)
  730 format(5e20.13)
       dsc = 1.0 / real(ns-1)
!======================================================
!  READ IN FOURIER COMPONENTS (!! Version Dependend)
!======================================================

      DO 30 J=1,ns
          LJ=(J-1)*MPNT

           IF (vmver == '8.00' .or. vmec_6_90) then 

           DO 15 MN=1,MPNT
             if (j .eq. 1) then 
             READ(8,*)  ixm(mn), ixn(mn)
             end if
           READ(8,*) XRC(mn,j), XZC(mn,j), PT1, BCC(mn,j), XJC(mn,j)
     &                ,PT2, PT3, PT1, PT2, PT3, PT1 
 15     ENDDO
           
          ELSE  !Version .ne. 8.00

!First block

           DO 16 MN=1,MPNT
             if (j .eq. 1) then 
             READ(8,*)  ixm(mn), ixn(mn)
             end if
           READ(8,*) XRC(mn,j), XZC(mn,j), PT1
 16     ENDDO

!Second block

           DO 17 MN=1,MPNT
             if (j .eq. 1) then 
             READ(8,*)  ixm(mn), ixn(mn)
             end if
           READ(8,*) BCC(mn,j), XJC(mn,j)
     &                ,PT2, PT3, PT1, PT2, PT3 
 17     ENDDO
          
          ENDIF

!======================================================
!            COMPONENTS READ
!======================================================

CPKME6---01--------01--------01--------01--------01--------01--------01>
        sc(j) = real(J-1)/(ns-1)
        shalf(j)=sc(j)-0.5*dsc
!.. 
        if (j.gt.1)then
          DO 25 MN=1,MPNT
              PWR = 0.5 * ixm(MN)
              XRC(MN,J)   = XRC(MN,J)/sc(J)**PWR
              XZC(MN,J)   = XZC(MN,J)/sc(J)**PWR
!              XJC(MN,J)   = XJC(MN,J)/(sc(J)-0.5*dsc)**PWR
 25       end do
        end if
 30    end do
       shalm(1:ns-1)=shalf(2:ns)
!.. FINE MESH GRIDS
       dsf = 1.0 / real(nsf-1)
       do 40 j=1,nsf
         sf(j)=real(j-1)/(nsf-1)
 40    end do
       shalf(1)    = 0.0
       shalf(NS+1) = 1.0
       DO 81 mn=1,mpnt
              PWR = 0.5 * ixm(MN)
!.. INTERPOLATE R_mn ONTO FINE GRID INTEGER MESH
         do 42 j=1,ns
          rg(j)=xrc(mn,j)
 42      end do
       call spline1(ns,sc,rg,f0,f1,f2)
!.. THE AXIS AND THE PVI OF COARSE AND FINE MESHES COINCIDE
         XRF(NSF,mn)=XRC(mn,NS) * sf(NSF)**pwr
         do 43 i=2,nsf-1
           xrf(i,mn)=seval(ns,sf(i),sc,rg,f0,f1,f2)
           xrf(i,mn)=xrf(i,mn)*sf(i)**pwr
 43      end do
         XRF(1,mn) = XRC(mn,1)
!.. INTERPOLATE Z_mn ONTO FINE GRID INTEGER MESH
         do 44 j=1,ns
          rg(j)=xzc(mn,j)
 44    end do
       call spline1(ns,sc,rg,f0,f1,f2)
!.. THE AXIS AND THE PVI OF COARSE AND FINE MESHES COINCIDE
         XZF(NSF,mn)=XZC(mn,NS) * sf(NSF)**pwr
         do 45 i=2,nsf-1
           xzf(i,mn)=seval(ns,sf(i),sc,rg,f0,f1,f2)
           xzf(i,mn)=xzf(i,mn)*sf(i)**pwr
 45      end do
         XZF(1,mn) = XZC(mn,1)
!.. INTERPOLATE JACOBIAN_mn ONTO FINE GRID HALF-INTEGER MESH
       do 46 j=1,ns
        rg(j)=xjc(mn,j)
 46    end do
!.. EXTRAPOLATION TO EDGE
        rg(ns+1)=1.5*rg(ns)-0.5*rg(ns-1)
!.. EXTRAPOLATION TO ORIGIN
       if (ixm(mn) /= 0) then
        rg(1)=0.0
!       call spline1(ns,shalf,rg,f0,f1,f2)
       call spline1(ns+1,shalf,rg,f0,f1,f2)
         do 47 i=2,nsf
           xjf(i,mn)=seval(ns+1,sf(i)-0.5*dsf,shalf,rg,f0,f1,f2)
!           xjf(i,mn)=xjf(i,mn)*(sf(i)-0.5*dsf)**pwr
 47      end do
       else
        if (ixm(mn) == 0)rg(1)=1.5*rg(2)-0.5*rg(3)
!        rg(1:ns-1)=rg(2:ns)
!       call spline1(ns-1,shalm,rg,f0,f1,f2)
       call spline1(ns+1,shalf,rg,f0,f1,f2)
         do 48 i=2,nsf
!           xjf(i,mn)=seval(ns-1,sf(i)-0.5*dsf,shalm,rg,f0,f1,f2)
!           xjf(i,mn)=seval(ns,sf(i)-0.5*dsf,shalf,rg,f0,f1,f2)
           xjf(i,mn)=seval(ns+1,sf(i)-0.5*dsf,shalf,rg,f0,f1,f2)
!           xjf(i,mn)=xjf(i,mn)*(sf(i)-0.5*dsf)**pwr
 48     end do
       end if
!
!.. INTERPOLATE BC_mn (MOD-B) ONTO FINE GRID HALF-INTEGER MESH
       do 69 j=1,ns
        rg(j)=bcc(mn,j)
 69   end do
!.. EXTRAPOLATION TO EDGE
        rg(ns+1)=1.5*rg(ns)-0.5*rg(ns-1)
!.. EXTRAPOLATION TO ORIGIN
       IF (ixm(mn) /= 0) then
        rg(1)=0.
!       call spline1(ns,shalf,rg,f0,f1,f2)
       call spline1(ns+1,shalf,rg,f0,f1,f2)
         do 70 i=2,nsf
           bcf(i,mn)=seval(ns+1,sf(i)-0.5*dsf,shalf,rg,f0,f1,f2)
 70     end do
       ELSE
        IF (ixm(mn) == 0)rg(1)=1.5*rg(2)-0.5*rg(3)
!        rg(1:ns-1)=rg(2:ns)
!       call spline1(ns-1,shalm,rg,f0,f1,f2)
       call spline1(ns+1,shalf,rg,f0,f1,f2)
         do 71 i=2,nsf
           bcf(i,mn)=seval(ns+1,sf(i)-0.5*dsf,shalf,rg,f0,f1,f2)
 71     end do
       END IF
!
!
 81   END DO
!
!.. PRINT GEOMETRY TO UNIT 18
!       WRITE(IBASPL,711) NFP, MPNT, NSF

!!            print *,'NSF=',nsf

        do 83 i=1,nsf
         do 82 mn=1,mpnt
            xm=real(ixm(mn))
            xn=real(ixn(mn))
            WRITE(18,850)XM, XN, XRF(i,mn), XZF(i,mn), XJF(i,mn)
!           WRITE(IBASPL,850)XM, XN, XRF(i,mn), XZF(i,mn), BCF(i,mn)
!           WRITE(73,710)XM, XN, XRF(i,mn), XZF(i,mn), XJF(i,mn),shd,shd
 82     end do
 83   end do
  850            FORMAT(1x,5es20.10)
  710            FORMAT(1x,1p7e13.6)
  711            format(3i6)
!
CPKM--------------------------------------------------------------------
!        READ (8,730)(AIOTA(J),AMASS(J),PRES(J),PHIP(J),PT1,PT2,PT3,
!     >  VP(J),PT1,PT2,PT3,PT1, J=2,NS)


         if (vmec_6_90) then


      	READ(8,*)(aiota(j),amass(j),pres(j),PT1,phip(j),PT2, PT3, !pax PTI<->phip
!NP     1            vp(j), PT2,PT3, j=2,ns) !pax added pt4, * statt 730
     1            PT4,vp(j),PT1,PT2,PT3,PT4, j=2,ns) !  NP

        else 

        READ(8,*)(PT1, PT2, PT3, PT1, PT2, PT3, j=1,ns)

	READ(8,*)(aiota(j),amass(j),pres(j),PT1,phip(j),PT2, PT3, !pax PTI<->phip
     1            vp(j), PT2,PT3, j=2,ns) !pax added pt4, * statt 730

       end if


        aiota(1) = 0.
        pres(1)  = 0.
!       AIOTA(1)  = 2.*AIOTA(2)-AIOTA(3)
!       pres(1)   = pres(2)
       amass(1)  = 0.
       phip(1)   = 0.
       vp(1)     = 0.
      do 84 j=2,ns
!      write(78,710)sc(j),aiota(j),amass(j),pres(j),-phip(j),vp(j)
 84   end do
!..  INTERPOLATION OF FLUX SURFACE FUNCTIONS ONTO FINER MESH
!..  ALL SURFACE FUNCTIONS ARE IN HALF-INTEGER MESH (THE JACOBIAN ABOVE ALSO)
!..  AND STORED IN THE RANGE 2 ---> NUMBER OF RADIAL GRID POINTS
!..  INTERPOLATION OF IOTA
      do 85 J=2,NS
!        rg(J-1) = aiota(J)
        rg(J) = aiota(J)
 85   end do
!.. EXTRAPOLATION TO ORIGIN AND EDGE
        rg(1)    = 1.5 * rg(2)  - 0.5 * rg(3)
        rg(ns+1) = 1.5 * rg(ns) - 0.5 * rg(ns-1)
       call spline1(ns+1,shalf,rg,f0,f1,f2)
         do 87 i=2,nsf
           AIOTA(i)=seval(ns+1,sf(i)-0.5*dsf,shalf,rg,f0,f1,f2)
 87      end do
!
!..  INTERPOLATION OF MASS
      do 88 J=2,NS
!        rg(J-1) = amass(J)
        rg(J) = amass(J)
 88   end do
!.. EXTRAPOLATION TO ORIGIN AND EDGE
        rg(1)    = 1.5 * rg(2)  - 0.5 * rg(3)
        rg(ns+1) = 1.5 * rg(ns) - 0.5 * rg(ns-1)
       call spline1(ns+1,shalf,rg,f0,f1,f2)
         do 89 i=2,nsf
           AMASS(i)=seval(ns+1,sf(i)-0.5*dsf,shalf,rg,f0,f1,f2)
 89     end do
!
!..  INTERPOLATION OF PRESSURE
      do 91 J=2,NS
!        rg(J) = PRES(J) * rmu0
        rg(J) = PRES(J)
 91   end do
!.. EXTRAPOLATION TO ORIGIN AND EDGE
        rg(1)    = 1.5 * rg(2)  - 0.5* rg(3)
        rg(ns+1) = 1.5 * rg(ns) - 0.5 * rg(ns-1)
       call spline1(ns+1,shalf,rg,f0,f1,f2)
         do 92 i=2,nsf
           PRES(i)=seval(ns+1,sf(i)-0.5*dsf,shalf,rg,f0,f1,f2)
 92     end do
!
!..  INTERPOLATION OF DERIVATIVE OF THE TOROIDAL FLUX FUNCTION
      do 95 J=2,NS
!        rg(J-1) = PHIP(J)
        rg(J) = PHIP(J)
 95   end do
!.. EXTRAPOLATION TO ORIGIN AND EDGE
        rg(1)    = 1.5 * rg(2)  - 0.5 * rg(3)
        rg(ns+1) = 1.5 * rg(ns) - 0.5 * rg(ns-1)
       call spline1(ns+1,shalf,rg,f0,f1,f2)
         do 96 i=2,nsf
           PHIP(i)=seval(ns+1,sf(i)-0.5*dsf,shalf,rg,f0,f1,f2)
 96     end do
!
!..  INTERPOLATION OF DERIVATIVE OF THE DIFFERENTIAL VOLUME
      do 97 J=2,NS
!        rg(J-1) = VP(J)
        rg(J) = VP(J)
 97   end do
!.. EXTRAPOLATION TO ORIGIN AND EDGE
        rg(1)    = 1.5 * rg(2)  - 0.5 * rg(3)
        rg(ns+1) = 1.5 * rg(ns) - 0.5 * rg(ns-1)
       call spline1(ns+1,shalf,rg,f0,f1,f2)
         do 98 i=2,nsf
           VP(i)=seval(ns+1,sf(i)-0.5*dsf,shalf,rg,f0,f1,f2)
 98     end do
!.. PRINT FLUX SURFACE QUANTITIES
       WRITE(18,850)(AIOTA(js),AMASS(js),PRES(js)*mu_0,-PHIP(js),VP(js),
     >               js=1,NSF)
!
      DO 99 i=2,nsf
!      write(77,710)sf(i)-0.5*dsf,aiota(i),amass(i),pres(i),-phip(i)
!     &                          ,vp(i)
 99   end do
CPKME-------------------------------------------------------------------

        CLOSE(18)
        RETURN
 931  WRITE(*,"(/,2X,A,/)") "--------------------- VMEC file not found. 
     .GIST will stop. -----------------------"
        STOP 
        END subroutine EQINTRP_08_i
!********0*********0*********0*********0*********0*********0*********0*>
ccc---------------------------------------------------------------------
      subroutine spline1 (n, x, y, b, c, d)                              00002230
      integer n                                                         00002240
      real x(n), y(n), b(n), c(n), d(n)                                 00002250
c                                                                       00002260
c  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed      00002270
c  for a cubic interpolating spline                                     00002280
c                                                                       00002290
c    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3  00002300
c                                                                       00002310
c    for  x(i) .le. x .le. x(i+1)                                       00002320
c                                                                       00002330
c  input..                                                              00002340
c                                                                       00002350
c    n = the number of data points or knots (n.ge.2)                    00002360
c    x = the abscissas of the knots in strictly increasing order        00002370
c    y = the ordinates of the knots                                     00002380
c                                                                       00002390
c  output..                                                             00002400
c                                                                       00002410
c    b, c, d  = arrays of spline coefficients as defined above.         00002420
c                                                                       00002430
c  using  p  to denote differentiation,                                 00002440
c                                                                       00002450
c    y(i) = s(x(i))                                                     00002460
c    b(i) = sp(x(i))                                                    00002470
c    c(i) = spp(x(i))/2                                                 00002480
c    d(i) = sppp(x(i))/6  (derivative from the right)                   00002490
c                                                                       00002500
c  the accompanying function subprogram  seval  can be used             00002510
c  to evaluate the spline.                                              00002520
c                                                                       00002530
c                                                                       00002540
      integer nm1, ib, i                                                00002550
      real t                                                            00002560
c                                                                       00002570
      nm1 = n-1                                                         00002580
      if ( n .lt. 2 ) return                                            00002590
      if ( n .lt. 3 ) go to 50                                          00002600
c                                                                       00002610
c  set up tridiagonal system                                            00002620
c                                                                       00002630
c  b = diagonal, d = offdiagonal, c = right hand side.                  00002640
c                                                                       00002650
      d(1) = x(2) - x(1)                                                00002660
      c(2) = (y(2) - y(1))/d(1)                                         00002670
      do 10 i = 2, nm1                                                  00002680
         d(i) = x(i+1) - x(i)                                           00002690
         b(i) = 2.*(d(i-1) + d(i))                                      00002700
         c(i+1) = (y(i+1) - y(i))/d(i)                                  00002710
         c(i) = c(i+1) - c(i)                                           00002720
   10 continue                                                          00002730
c                                                                       00002740
c  end conditions.  third derivatives at  x(1)  and  x(n)               00002750
c  obtained from divided differences                                    00002760
c                                                                       00002770
      b(1) = -d(1)                                                      00002780
      b(n) = -d(n-1)                                                    00002790
      c(1) = 0.                                                         00002800
      c(n) = 0.                                                         00002810
      if ( n .eq. 3 ) go to 15                                          00002820
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))                        00002830
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))              00002840
      c(1) = c(1)*d(1)**2/(x(4)-x(1))                                   00002850
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))                              00002860
c                                                                       00002870
c  forward elimination                                                  00002880
c                                                                       00002890
   15 do 20 i = 2, n                                                    00002900
         t = d(i-1)/b(i-1)                                              00002910
         b(i) = b(i) - t*d(i-1)                                         00002920
         c(i) = c(i) - t*c(i-1)                                         00002930
   20 continue                                                          00002940
c                                                                       00002950
c  back substitution                                                    00002960
c                                                                       00002970
      c(n) = c(n)/b(n)                                                  00002980
      do 30 ib = 1, nm1                                                 00002990
         i = n-ib                                                       00003000
         c(i) = (c(i) - d(i)*c(i+1))/b(i)                               00003010
   30 continue                                                          00003020
c                                                                       00003030
c  c(i) is now the sigma(i) of the text                                 00003040
c                                                                       00003050
c  compute polynomial coefficients                                      00003060
c                                                                       00003070
      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n))         00003080
      do 40 i = 1, nm1                                                  00003090
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.*c(i))          00003100
         d(i) = (c(i+1) - c(i))/d(i)                                    00003110
         c(i) = 3.*c(i)                                                 00003120
   40 continue                                                          00003130
      c(n) = 3.*c(n)                                                    00003140
      d(n) = d(n-1)                                                     00003150
      return                                                            00003160
c                                                                       00003170
   50 b(1) = (y(2)-y(1))/(x(2)-x(1))                                    00003180
      c(1) = 0.                                                         00003190
      d(1) = 0.                                                         00003200
      b(2) = b(1)                                                       00003210
      c(2) = 0.                                                         00003220
      d(2) = 0.                                                         00003230
      return                                                            00003240
      end                                                               00003250
      real function seval(n, u, x, y, b, c, d)                          00003260
      integer n                                                         00003270
      real  u, x(n), y(n), b(n), c(n), d(n)                             00003280
c                                                                       00003290
c  this subroutine evaluates the cubic spline function                  00003300
c                                                                       00003310
c    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3 00003320
c                                                                       00003330
c    where  x(i) .lt. u .lt. x(i+1), using horner's rule                00003340
c                                                                       00003350
c  if  u .lt. x(1) then  i = 1  is used.                                00003360
c  if  u .ge. x(n) then  i = n  is used.                                00003370
c                                                                       00003380
c  input..                                                              00003390
c                                                                       00003400
c    n = the number of data points                                      00003410
c    u = the abscissa at which the spline is to be evaluated            00003420
c    x,y = the arrays of data abscissas and ordinates                   00003430
c    b,c,d = arrays of spline coefficients computed by spline           00003440
c                                                                       00003450
c  if  u  is not in the same interval as the previous call, then a      00003460
c  binary search is performed to determine the proper interval.         00003470
c                                                                       00003480
      integer i, j, k                                                   00003490
      real dx                                                           00003500
      data i/1/                                                         00003510
      if ( i .ge. n ) i = 1                                             00003520
      if ( u .lt. x(i) ) go to 10                                       00003530
      if ( u .le. x(i+1) ) go to 30                                     00003540
c                                                                       00003550
c  binary search                                                        00003560
c                                                                       00003570
   10 i = 1                                                             00003580
      j = n+1                                                           00003590
   20 k = (i+j)/2                                                       00003600
      if ( u .lt. x(k) ) j = k                                          00003610
      if ( u .ge. x(k) ) i = k                                          00003620
      if ( j .gt. i+1 ) go to 20                                        00003630
c                                                                       00003640
c  evaluate spline                                                      00003650
c                                                                       00003660
   30 dx = u - x(i)                                                     00003670
      seval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))                    00003680
      return                                                            00003690
      end                                                               00003700
ccc---------------------------------------------------------------------

