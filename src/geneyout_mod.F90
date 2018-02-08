
   MODULE geneyout_mod
   implicit none

CONTAINS

SUBROUTINE geneyout
  USE type_mod
  USE specifications_mod, only:pol_turns,q0_vvbal,s_vvbal,nz0,name,b0,s_vvbal,theta0,nalpha0, &
                             aref_bel,bref_bel,surf,parity,alpha0,vmec_file,table_tag,out_tag, &
                             R0_vmec,ftrap,profw7,refs,gviu,alpha0_out,global
  USE cidata, only : ntmax 
  USE vmec_mod, ONLY : vmec_data
  IMPLICIT NONE

INTEGER :: k,n,nt

REAL(DP), DIMENSION(:), ALLOCATABLE  :: g11_gene,g12_gene,g22_gene,metr_ratio,&
                               &bmod_gene,dBdx_vv,dBdxvac_vv,dBdy_vv,dbmod_gene,jac_gene,&
                               &sloc_vv
REAL(DP), DIMENSION(:), ALLOCATABLE  :: g11_vv,g12_vv,g22_vv,curvdr_vv,&
                               &jac_vv,bmod_vv,dbmod_vv,dBdv1_vv,dBdv2_vv,s_vv,&
                               &bdgrad_vv

REAL(DP) :: ftp_vv,qpr_vv,alf_vv,safety,s0,shat,alpha,ba,dz,d1,&
drhodpsi,shat_gs2,shat_vvbal,R_norm
REAL(DP), DIMENSION(3) :: d2
REAL(DP) :: my_dpdx,pp_vv
INTEGER :: safety_sign,sign_,sign,res,pol_turns_gene
INTEGER, PARAMETER :: ncols=6
LOGICAL :: first=.true., first1=.true.,first2=.true.
!====================================================================

res=ntmax !augmented by 1 extra point

ALLOCATE(g11_vv(res+1),g12_vv(res+1),g22_vv(res+1),bmod_vv(res+1),dbmod_vv(res),sloc_vv(res))
    ALLOCATE(dBdv1_vv(res+1),dBdv2_vv(res+1),jac_vv(res+1))
    ALLOCATE(curvdr_vv(res+1),s_vv(res+1),bdgrad_vv(res+1))      
    ALLOCATE(g11_gene(res+1),g12_gene(res+1),g22_gene(res+1),bmod_gene(res+1))
    ALLOCATE(dBdx_vv(res+1),dBdxvac_vv(res+1),dBdy_vv(res+1),dbmod_gene(res+1),jac_gene(res+1))
    ALLOCATE(metr_ratio(res+1))


!read prepar.f
 open(33,file='../int/vvbal_'//trim(vmec_file)//'_' &
     &        //trim(table_tag)//'_'//trim(out_tag),action='read') 

  DO n=1,res 
                 READ(33,"(14ES20.10)") s_vv(n),&
                  &g11_vv(n),g12_vv(n),g22_vv(n),&
                  &bmod_vv(n),bdgrad_vv(n),curvdr_vv(n),dBdv1_vv(n),&
                  &dBdv2_vv(n),ftp_vv,qpr_vv,alf_vv,jac_vv(n),pp_vv
          ENDDO
close(33)

!--------------- Auxilliary calculations
 CALL vmec_data(d1,d1,d1,d1,d1,alpha,d1,d1,d2,d2)

!!GENE quantities via VVBAL =================================
          s0=s_vvbal ; safety=abs(q0_vvbal) ; shat_vvbal=2*s_vvbal/safety*qpr_vv
          safety_sign=safety/abs(safety)
          R_norm=alpha !minor radius
          g11_gene=g11_vv*R_norm**2/4./s0
          g12_gene=g12_vv*R_norm**2/2./safety
          g22_gene=g22_vv*R_norm**2*s0/safety**2
          bmod_gene=bmod_vv/ABS(B0)
!          jac_gene=ABS(safety*ABS(B0)/Bp_out*r_out/R_norm) !2*safety*ftp_vv*r_out/R_norm**3/Bp_out
          jac_gene=2*safety*ABS(jac_vv)/R_norm**3 !! x^3=theta
          dBdx_vv=-R_norm**2*sqrt(s0)/safety*curvdr_vv 
          dBdxvac_vv=-R_norm**2*sqrt(s0)/safety*dbdv1_vv
          my_dpdx=-R_norm**2*sqrt(s0)*2*pp_vv/ftp_vv/ABS(B0)
          dBdy_vv=R_norm**2/2./sqrt(s0)*dBdv2_vv !-L_1 as there is a sign inversion in GENE

!!Parallel derivative for Bfield and g12/g11
    dz = (2.0 * pi * pol_turns) / REAL(res-1)
    metr_ratio=g12_gene/g11_gene    

dbmod_vv=0.

    DO k=1, res-1 
       dbmod_vv(k) = (bmod_vv(k+1)-bmod_vv(k))/dz
       sloc_vv(k) = (metr_ratio(k+1)-metr_ratio(k))/(dz)       
    END DO


!!Apply parity
      !    sign_=Bp_out(res/2)/ABS(Bp_out(res/2))
          sign_=parity
          g12_gene=g12_gene*safety_sign
          dBdx_vv=dBdx_vv*sign_*safety_sign
          dBdxvac_vv=dBdxvac_vv*sign_*safety_sign
          dBdy_vv=dBdy_vv*sign_
          my_dpdx=my_dpdx*sign_
      !    dbmod_vv=dbmod_vv*sign_

      !    write(*,"(4X,A,I5)"), 'parity=',nint(parity)

pol_turns_gene=pol_turns


              if (first) then
             WRITE(35,"(A)") "&parameters"
             WRITE(35,"(A,F12.7,I5)") "!rho0, nalpha0 = ",sqrt(s_vvbal),nalpha0
             WRITE(35,"(A,2F12.7)") "!R0,r0[m]= ",R0_vmec,r_norm
       !      WRITE(35,"(A,I3)") "!parity = ",sign_
             WRITE(35,"(A,F12.7)") "my_dpdx = ",my_dpdx
             WRITE(35,"(A,F12.7)") "q0 = ",safety
             WRITE(35,"(A,F12.7)") "shat = ", shat_vvbal*safety_sign
             WRITE(35,"(A,I5)") "gridpoints = ",res-mod(res,2)
             WRITE(35,"(A,I3)") "n_pol = ",pol_turns_gene
             WRITE(35,"(A)") "/"
             first=.false.
              endif

          DO n=1,res-1
            WRITE(35,"(8ES20.10)") & !s_vv(n)-theta0,&
                  &g11_gene(n),g12_gene(n),g22_gene(n),&
                  &bmod_gene(n),jac_gene(n),&
                  &dBdxvac_vv(n),dBdy_vv(n),&
                  &dbmod_vv(n)/ABS(B0) !,& nunami
                !  &-R_norm/(2*sqrt(s0))*dBdv2_vv(n),R_norm*sqrt(s0)/(safety)*dbdv1_vv(n)
            ENDDO


 if (global .and. gviu) then
       if (first2) then
              WRITE(95,"(A)") "&parameters"
              WRITE(95,"(A,3I7,3X,A,3X,A)") "(#lines, #points, #cols, coordinates, code) =", &
                                         & nalpha0, res-1, ncols, "BOOZER", "GIST"
              WRITE(95,*) "s, iota = ",s_vvbal,1./safety
              WRITE(95,"(A)") "/"
              WRITE(95,"(8(A,X))") "theta","y","modB","kbad","kgeo","sloc","g11","g22"       
      endif
      first2=.false.

      sign=B0/abs(B0) ! -1 for NCSX, +1 for W7X

       DO n=1,res-1
            WRITE(95,"(8ES20.8)") s_vv(n),alpha0_out,bmod_vv(n),dBdxvac_vv(n),dBdy_vv(n),sloc_vv(n),g11_gene(n),g22_gene(n)
       ENDDO
  endif
!==================================================================================!
!-------------------------- Output on the screen ----------------------------------!
!==================================================================================!

if (first1) then

 WRITE(*,"(/,2X,A,/)")  '----------------------------- Surface characteristics ---------------------------'
 WRITE(*,"(22X,A,2F9.4,//)") '    iota, shear = ',1./abs(q0_vvbal),shat_vvbal*safety_sign

 WRITE(*,"(2X,A,/)") '------------------------------- GENE normalization ------------------------------' 
 WRITE(*,"(18X,A,2F9.4,/)") ' Bref[T], Lref(minor radius)[m] =',ABS(B0),alpha


 if (profw7) then
 WRITE(*,"(4X,A,3F12.7,/)") 'Ti_ref, Te_ref, n_ref =', refs(1),refs(2),refs(3)
 endif

write (*,"(/,2X,A,/)") "-------------------- GIST file successfully generated for GENE -------------------"

first1=.false.

endif 

END subroutine geneyout



!=============================================================!

END MODULE geneyout_mod








