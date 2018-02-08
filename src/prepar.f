        subroutine prepar                                               VBA2700
!
      use vbal_param
      use cidata
      use crdata
      use ceigen
      use cmodes
      use cfoura
      use ccoef
      use specifications_mod, only:s_tracer,q0_vvbal,out_tag,s_vvbal, 
     &  aref_bel,bref_bel,parity,vmec_file, table_tag, sloc_unif
      use radfct, only:ftpp,fppp
!
C.. Implicits ..
      implicit none
!       include 'combal.inc'                                            VBA2710
!
! ... LOCAL SCALARS
      integer :: j, mn, iarg, mn0
!       double precision ::
      real :: 
     &   ds, dh, st, zet, arg, cs, sn, signjac, qsq, ftpsq, rkpamp,
     &   dpdpsi, crvamp, bmod, bdgrad, rkperp, curvdr, curvgeo, gradbd
     &  ,g12,quasipax,bavg,wgt
!
! .. Intrinsic Functions ..
      intrinsic  sqrt, abs
!x                                                                      VBA2720
        sbegin = -0.5 * twopi * s0in  + spk                               VBA2730
        st = s0in * twopi                                                 VBA2740
        ds = st / real(ntmax-1)                                         VBA2750
      !!  ds=st*qs/(ntmax-2) !pax   
        dh = .5 * ds                                                    VBA2760

!INITIALIZATION
        do 2 j=1,ntmax                                                  VBA2780
          s(j)=float(j-1)*ds+sbegin                                     VBA2790
          fs(j)=0.
          fq(j)=0.
          c0(j)=0.                                                      VBA2800
          c1(j)=0.                                                      VBA2810
          c111(j)=0.
          c2(j)=0.
          c3(j)=0.
          cbel(j)=0.
          pax1(j)=0.                          
          gss(j)=0.
          gst(j)=0.                                                     VBA2820
 2      end do
!                                                                       VBA2830
        do 10 mn=1,modes                                                VBA2840
          if(mb(mn).eq.0 .and. nb(mn).eq.0)mn0=mn
CDIR$ IVDEP                

             do 5 j=1,ntmax 
             zet = alf + qs * s(j)                                      VBA2870
             arg=twopin*(mb(mn)*s(j)-nb(mn)*zet)                        VBA2880
             iarg=arg                                                   VBA2890
             arg=twopi*(arg-iarg)                                       VBA2900
             cs = cos(arg)                                              VBA2910
             sn = sin(arg)                                              VBA2920
             c2(j)=c2(j)+(a0mn(mn)+a2mn(mn)*(s(j)-sk)*(s(j)-sk))*cs     VBA2930
     &             +a1mn(mn)*(s(j)-sk)*sn                               VBA2940
       c3(j)=c3(j)+a2mn(mn)*(s(j)-sk)*cs+.5*a1mn(mn)*sn
!       cbel(j)=cbel(j)+((a0mn(mn)+a2mn(mn)*(s(j)-sk)*(s(j)-sk))*cs
!     &        +a1mn(mn)*(s(j)-sk)*sn)
        cbel(j)=cbel(j)+ ((a0mn(mn) + a2mn(mn)*s(j)*s(j))*cs         
     .             + a1mn(mn)*s(j)*sn)

!00             fs(j)=fs(j)+a1mn(mn)*sn
!00             fq(j)=fq(j)+a2mn(mn)*cs
!00             c2(j)=c2(j)+a0mn(mn)*cs
             c1(j)=c1(j)+d0mn(mn)*cs+d1mn(mn)*(s(j)-sk)*sn 
             c111(j)=c111(j)+d0mn(mn)*cs
             c0(j)=c0(j)+rgmn(mn)*cs
             pax1(j)=pax1(j)+d1mn(mn)*sn!*(s(j)-sk)
             gss(j)=gss(j)+gssmn(mn)*cs !g^ss
             gst(j)=gst(j)+gstmn(mn)*cs !g_st/jac
             bs(j)=bs(j)+bsmn(mn)*cs !pax
 5        end do
 10     end do
        signjac=rgmn(mn0)/abs(rgmn(mn0))

!--FIELD LINE BENDING TERM
!00        do j=1,ntmax
!00          c2(j) = c2(j)*(1.0+(fs(j)+fq(j)*(s(j)-sk))**2)
!00        end do
!--COEFFICIENTS REQUIRED FOR KINETIC BALLOONING

        open(23,file='../int/vvbal_'//trim(vmec_file)//'_'
     &        //trim(table_tag)//'_'//trim(out_tag),action='write')

        if (lxcont .gt. 0) then
        qsq    = qs * qs !qs=safety factor
        ftpsq  = ftp * ftp
        rkpamp = sqgbsq * qsq / ftpsq !sqgbsq=jac*B^2, ftpsq=|Phi'|^2 
        dpdpsi = pp * qs / ftp
        crvamp = qs / (2. * pp * ftp) !pp=P'
        bavg=0. ; wgt=0.
        do j=1,ntmax 
        bmod   = sqrt(rkpamp / c0(j))
        bdgrad = 1. / sqrt(sqgbsq * c0(j))
        rkperp = rkpamp * signjac * abs(c2(j))
        curvdr = crvamp * bmod * c1(j) / c0(j)
        curvgeo=.5*ftp/qs/pp/sqgbsq * Bmod**3 * pax1(j)/qprime 
        gradbd = bmod * curvdr + dpdpsi
        quasipax=.5*ftp/qs*Bmod**3/sqgbsq/pp*c111(j)
        g12 = c3(j)*sqgbsq*qsq/ftpsq/qprime
        zet = alf + qs*s(j)
!       write(15,666) zet,s(j),bmod,bdgrad,rkperp,curvdr,gradbd,dpdpsi, !8
!     & pth,gss(j),gst(j),pp/ftp/sqgbsq,curvgeo,curvdr-quasipax,quasipax
        bavg=bavg+sqrt(rkpamp/c0(j))
        wgt=wgt+1.

!Calculate parity using flux
        
        parity=ftp/abs(ftp)
   
!!FILE FOR GENE
                                             
!pass variables
         q0_vvbal=qs ; s_tracer=rs; s_vvbal=rs
         
!uniform local shear if required
         if (sloc_unif) then
         g12=gss(j)*s(j)*qprime
         endif
        
         write(23,"(14ES20.10)") s(j),gss(j),g12,rkperp,bmod,bdgrad,
     &   curvdr,gradbd/bmod,curvgeo,ftp,qprime,alf,sqgbsq/bmod**2,pp
       enddo

 666    format(1p20e17.8)

!!BELLI's 
                aref_bel=1./sqrt(rkpamp*cbel((ntmax+1)/2))
                bref_bel=bavg/wgt

        endif

       
!!        write(*,*) 'aref_bel,bref_bel',aref_bel,bref_bel


       close(23)

!      print *,'qprime=',qprime
       
!--BALLOONING MODE COEFFICIENTS ON A FIELD LINE IN REAL SPACE.          VBA2970
        if (lmodel .le. 1) then
        do 15 j=1,ntmax                                                 VBA2980
          c0(j)=b0r0*c0(j)*c0(j)*c2(j)                                  VBA2990
          c2(j)=signjac * abs(c2(j))
          c2(j)=1./c2(j)                                                VBA3000
 15     end do
        else
        do 16 j=1,ntmax
          c0(j)=c0(j)*(1.+qppsip*(s(j)-sk)*(s(j)-sk))
          c2(j)=signjac * abs(c2(j))
          c2(j)=1./c2(j)
 16   end do
        endif
!                                                                       VBA3010
        return                                                          VBA3020
        end subroutine prepar       
