MODULE gs2_mod

 character(len=256),PRIVATE::fname
 
CONTAINS
!-------------------------------------------------------------------
!-------------------------------------------------------------------

SUBROUTINE gs2_out
  USE type_mod
  USE specifications_mod, only:pol_turns,q0_vvbal,s_vvbal,nz0,b0, &
         aref_bel,bref_bel,surf,parity,gs2_press,vmec_file,table_tag,out_tag, &
         terpsichore,pest
  USE vmec_mod, ONLY : vmec_data
  IMPLICIT NONE

INTEGER :: k,n,nt
REAL(DP), DIMENSION(:), ALLOCATABLE  :: g11_vv,g12_vv,g22_vv,curvdr_vv,bdgrad_vv
REAL(DP), DIMENSION(:), ALLOCATABLE  :: jac_vv,bmod_vv,dbmod_vv,gradbdoverbmod,curvgeo,s_vv
REAL(DP) :: ftp_vv,qpr_vv,alf_vv,q0,s0,shat,alpha,ba,dz,d1,drhodpsi,shat_gs2,sign_q0
REAL(DP), DIMENSION(3) :: d2
!====================================================================

nt=2*nz0
                               
ALLOCATE(g11_vv(nt+1),g12_vv(nt+1),g22_vv(nt+1),bmod_vv(nt+1),dbmod_vv(nt))
ALLOCATE(gradbdoverbmod(nt+1),curvgeo(nt+1),jac_vv(nt+1))
ALLOCATE(curvdr_vv(nt+1),s_vv(nt+1),bdgrad_vv(nt+1))     


open(35,file='../out/gist_gs2_'//trim(vmec_file)//'_' &
     &     //trim(table_tag)//'_'//trim(out_tag),action='write') 


!read prepar.f
 open(33,file='../int/vvbal_'//trim(vmec_file)//'_' &
     &     //trim(table_tag)//'_'//trim(out_tag),action='read')
 
  DO n=1,nt+1
             READ(33,"(13ES20.10)") s_vv(n),&
                  &g11_vv(n),g12_vv(n),g22_vv(n),&
                  &bmod_vv(n),bdgrad_vv(n),curvdr_vv(n),gradbdoverbmod(n),&
                  &curvgeo(n),ftp_vv,qpr_vv,alf_vv,jac_vv(n)
          ENDDO
close(33)

!--------------- Auxilliary calculations
 CALL vmec_data(d1,d1,d1,d1,d1,alpha,d1,d1,d2,d2)

 dbmod_vv=0.
 dz = (2.0 * pi * pol_turns) / REAL(nt+1)
    DO k=2, nt
       dbmod_vv(k) = (bmod_vv(k+1)-bmod_vv(k-1))/(2.0*dz)
!       sloc_vv(k) = (metr_ratio(k+1)-metr_ratio(k-1))/(2.0*dz)       
    END DO
!!
q0=q0_vvbal; s0=s_vvbal ; sign_q0=q0/abs(q0)
!   the following will give negative s-hat for rising abs(q(r))
!   tested with NCSX flip and LHD configurations
shat=2*s0/q0*qpr_vv
!  Using sign_q0 does not give correct sign of s-hat:
! DO NOT USE: shat=2*s0/q0*qpr_vv*sign_q0
!alpha=bc_rminor()  
!ba=bc_torflux(1.d0)/Pi/alpha**2

!drhodpsi=.5/sqrt(s0)*alpha**2*abs(B0)*q0/ftp_vv
shat_gs2=qpr_vv*bref_bel*aref_bel**2/abs(ftp_vv)*sign_q0
drhodpsi=1.

!  ouput normalizing length and magnetic field in GIST version 2
!  use abs() for q0 and B0 so they'll always be positive
        WRITE(35,"(I5,7ES15.6)") nz0,pol_turns,shat,drhodpsi, &
           abs(q0),s0,alpha,abs(B0)
!was:      nz0,pol_turns,shat_gs2,drhodpsi,aref_bel,bref_bel

   if (gs2_press) then

        DO n=1,nt+1 !extra point for GS2
           WRITE(35,"(10ES20.10)") &
           s_vv(n),bmod_vv(n)/abs(B0),alpha*bdgrad_vv(n),&
           alpha**2*s0/q0**2*g22_vv(n),&
           -qpr_vv*alpha**2*s0/q0**2*g12_vv(n),&                
           qpr_vv**2*alpha**2*s0/q0**2*g11_vv(n),&
           2*abs(B0)/bmod_vv(n)*alpha**2*sqrt(s0)/q0*curvdr_vv(n)*parity,&
           2*abs(B0)/bmod_vv(n)*alpha**2*qpr_vv*sqrt(s0)/q0*curvgeo(n)*parity,&
           2*abs(B0)/bmod_vv(n)*alpha**2*sqrt(s0)/q0*gradbdoverbmod(n)*parity,&
           2*abs(B0)/bmod_vv(n)*alpha**2*qpr_vv*sqrt(s0)/q0*curvgeo(n)*parity

        ENDDO  
   else

        print *,'Pressure gradient switched OFF'

        DO n=1,nt+1 !extra point for GS2
           WRITE(35,"(10ES20.10)") &
           s_vv(n),bmod_vv(n)/abs(B0),alpha*bdgrad_vv(n),&
           alpha**2*s0/q0**2*g22_vv(n),&
           -qpr_vv*alpha**2*s0/q0**2*g12_vv(n)*sign_q0,&                
           qpr_vv**2*alpha**2*s0/q0**2*g11_vv(n),&
           2*abs(B0)/bmod_vv(n)*alpha**2*sqrt(s0)/q0*gradbdoverbmod(n)*parity,&
           2*abs(B0)/bmod_vv(n)*alpha**2*qpr_vv*sqrt(s0)/q0*curvgeo(n)*parity,&
           2*abs(B0)/bmod_vv(n)*alpha**2*sqrt(s0)/q0*gradbdoverbmod(n)*parity,&
           2*abs(B0)/bmod_vv(n)*alpha**2*qpr_vv*sqrt(s0)/q0*curvgeo(n)*parity
        ENDDO

    endif

close(35)

write (*,"(/,2X,A,/)") "-------------------- GIST file successfully generated for GS2 --------------------"

END SUBROUTINE gs2_out

END MODULE gs2_mod












