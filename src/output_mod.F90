MODULE output_mod
  USE specifications_mod
  USE UtilCDF
  USE pest_mod
  USE benchmark
  USE benchmark_mod
  USE field_d3d_mod, ONLY : d3d_iterdb
  USE vmec_mod, ONLY : bc_cyl2mag
  IMPLICIT NONE

!=========================================================================!
CONTAINS           !              MODULE SUBPROGRAMS                      !
!=========================================================================!

!====================== OPEN OUTPUT FILES ================================!
  SUBROUTINE open_out

    CHARACTER(255) :: filename

!Open output files

!    OPEN(UNIT=41,FILE='../out/tracer_'//trim(vmec_file)//'_' &
!     &        //trim(table_tag)//'_'//trim(out_tag),STATUS="REPLACE",ACTION="WRITE")
!    OPEN(UNIT=14,FILE='../out/tracer.dat_'//trim(vmec_file)//'_' &
!     &        //trim(table_tag)//'_'//trim(out_tag),STATUS="REPLACE",ACTION="WRITE")
!    OPEN(UNIT=12,FILE='../out/gij.out',STATUS="REPLACE",ACTION="WRITE")
!    OPEN(UNIT=17,FILE='../out/curv.out',STATUS="REPLACE",ACTION="WRITE")

 END SUBROUTINE open_out
!====================== CLOSE OUTPUT FILES ===============================!
  SUBROUTINE close_out

!    CLOSE(41)
!    CLOSE(14)
!    CLOSE(12)
!    CLOSE(13)

  END SUBROUTINE close_out

!====================== ALLOCATE OUTPUT ARRAYS ===========================!

  SUBROUTINE alloc_out
    ALLOCATE(c_out(3,3,-2:resol+3),d_out(3,3,-2:resol+3),gij_out(3,3,-2:resol+3),g_ij_out(3,3,-2:resol+3))
    ALLOCATE(Br_out(-2:resol+3),Bp_out(-2:resol+3),Bz_out(-2:resol+3),divB_out(-2:resol+3),b3_out(-2:resol+3))
    ALLOCATE(jac_out(-2:resol+3),inv_out(-2:resol+3),k1_out(-2:resol+3),k2_out(-2:resol+3))
    ALLOCATE(k_norm_out(-2:resol+3),k_geo_out(-2:resol+3),kappa_sq_out(-2:resol+3))
    ALLOCATE(zero_kappa_out(-2:resol+3),dBdtau_out(-2:resol+3))
    ALLOCATE(dBdp_out(-2:resol+3),dBdz_out(-2:resol+3),dBdr_out(-2:resol+3),Bfield_out(-2:resol+3))
    ALLOCATE(dbdv1_out(-2:resol+3),dBdv2_out(-2:resol+3))
    ALLOCATE(zero_cl1_out(-2:resol+3),zero_cl2_out(-2:resol+3))
    ALLOCATE(zero_cl3_out(-2:resol+3),r_out(-2:resol+3),p_out(-2:resol+3),z_out(-2:resol+3))
    ALLOCATE(Bt_out(-2:resol+3),cart_x_out(-2:resol+3),cart_y_out(-2:resol+3),sloc(-2:resol+3),theta_out(-2:resol+3))
  END SUBROUTINE alloc_out

 !========================= PREPARE OUTPUT ================================!

  SUBROUTINE prep_out(ind,line)

    INTEGER, INTENT(IN) :: ind,line
    !Internal
    INTEGER :: nout_line

    r_out(ind)=r ; p_out(ind)=phi ; z_out(ind)=z
    Bfield_out(ind)=Bfield ; divb_out(ind)=divb
    Br_out(ind)=Br ; Bp_out(ind)=Bp ; Bz_out(ind)=Bz
    dBdtau_out(ind)=dBdtau
    dBdv1_out(ind)=dBdv1 ; dBdv2_out(ind)=dBdv2  
    jac_out(ind)=jac
    k1_out(ind)=k1 ; k2_out(ind)=k2
    k_norm_out(ind)=k_norm ; k_geo_out(ind)=k_geo ; kappa_sq_out(ind)=kappa_sq  
    inv_out(ind)=inv
    gij_out(:,:,ind)=gij(:,:) ; g_ij_out(:,:,ind)=g_ij(:,:)
    c_out(:,:,ind)=c(:,:) 
    zero_cl1_out(ind)=zero_cl1 ; zero_cl2_out(ind)=zero_cl2
    zero_cl3_out(ind)=zero_cl3
    zero_kappa_out(ind)=zero_kappa
    d_out(:,:,ind)=d(:,:) 
    dBdr_out(ind)=dBdr ; dBdz_out(ind)=dBdz ; dBdp_out(ind)=dBdp
    b3_out(ind)=b3 ; Bt_out(ind)=Bt 
    cart_x_out(ind)=cart_x ; cart_y_out(ind)=cart_y
    theta_out(ind)=theta

  END SUBROUTINE prep_out

!========================== PRODUCE OUTPUT ===============================!

     SUBROUTINE output(line)

       INTEGER, INTENT(IN) :: line
       INTEGER :: j,k,l,ll,l1,l2,l3,m,n,res,nout_line,turn,rind,zind,cut,ind,sign_,safety_sign
       INTEGER :: pol_turns_gene
       LOGICAL :: header=.true., header_boris=.true.
       REAL(DP) :: det_g_2d,dz, scale_sloc,scale_kgeo,my_beta=1.,my_dpdx,pp_vv
       REAL(DP) :: ppp,per,dum,ftp_vv,alf_vv,qpr_vv,s0
       REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: g_bor
       REAL(DP), DIMENSION(:), ALLOCATABLE  :: jac_bor,theta,elstpot,int
       REAL(DP), DIMENSION(:), ALLOCATABLE  :: dbmod_vv,bmodaux
       REAL(DP), DIMENSION(:), ALLOCATABLE  :: g11_ben,g12_ben,g22_ben,&
                               &bmod_ben,dBdv1_ben,dBdv2_ben
       
       REAL(DP), DIMENSION(:), ALLOCATABLE  :: metr_ratio,&
                               &bmod_gene,dBdxvac_vv,dbmod_gene,jac_gene,&
                               &sloc_vv

       REAL(DP), DIMENSION(3) :: cylin
       REAL(DP), DIMENSION(3) :: help_booz
       REAL(DP), DIMENSION(:), ALLOCATABLE  :: s_booz



!Adjust resolution

          IF(trace_type(line) == 1) THEN
             res=resol_f+1
          ELSEIF(trace_type(line) == -1) THEN
             res=resol_b+1
          ELSE
             res=resol
          ENDIF


!Compute BOOZER coordinates from cylindrical ones
allocate(s_booz(res),tht_booz(res),zet_booz(res))
do k=1,res

               cylin(1)=r_out(k);cylin(2)=p_out(k);cylin(3)=z_out(k)
               help_booz=bc_cyl2mag(cylin)
               s_booz(k)=help_booz(1);tht_booz(k)=help_booz(2);zet_booz(k)=help_booz(3)
enddo

pol_turns_gene=pol_turns

!================================================================================================!
!---------------------------------  TRACER.OUT  -------------------------------------------------!
!================================================================================================!

#if 0
       IF (header) THEN
          WRITE(41,*) 'Configuration: ',config
          WRITE(41,"(A)") '00:index 01:line 02:r 03:phi 04:z 05:Bfield 06:Br 07:Bp 08:Bz 09:jac'
          WRITE(41,"(A)") '10:k1 11:k2 12:k_norm 13:k_geo 14:kappa 15:inv 16:g11 17:g12 18:g13 19:g22'
          WRITE(41,"(A)") '20:g23 21:g33 22:d11 23:d12 24:d13 25:d21 26:d22 27:d23 28:sloc'
       ENDIF

       WRITE(41,"(I6,I6)") res,trace_type(line)

       DO k=1,res
          WRITE(41,"(I4,1X,I2,3F12.5,38ES20.10)")  &
               &k,line,r_out(k),p_out(k),z_out(k),Bfield_out(k),Br_out(k),Bp_out(k),& !8
               &Bz_out(k),jac_out(k),k1_out(k),k2_out(k),k_norm_out(k),k_geo_out(k),& !14
               &sqrt(kappa_sq_out(k)),inv_out(k) ,gij_out(1,1,k),gij_out(1,2,k),& !18
               &gij_out(1,3,k),gij_out(2,2,k),gij_out(2,3,k),gij_out(3,3,k),d_out(1,1,k),& !23
               &d_out(1,2,k),d_out(1,3,k),d_out(2,1,k),d_out(2,2,k),d_out(2,3,k),sloc(k),& !29
              &Bt_out(k),cart_x_out(k),cart_y_out(k),c_out(1,1,k),c_out(1,2,k),& !34
               &c_out(1,3,k),c_out(2,1,k),c_out(2,2,k),c_out(2,3,k),& !38
               &theta_out(k),dBdv1_out(k),dBdv2_out(k),g_ij_out(1,3,k),& !42
               &dBdtau_out(k) 
       ENDDO
#endif
!==========================================================================================
!---------------------------------  TRACER.DAT  -------------------------------------------
!==========================================================================================


          first_line: IF (line == 1) THEN


          IF (tracer) THEN
             R_norm=alpha ; safety=safetyf(line) !q0_vvbal
#if 0
             WRITE(14,"(A)") "&parameters"
             WRITE(14,"(A,F12.8)") "q0 = ",safetyf(line)
             WRITE(14,"(A,F12.8)") "shat = ",shatf(line)
             WRITE(14,"(A,F12.8)") "!rho = ",sqrt(s_in)
             WRITE(14,"(A,F12.8)") "!theta_0 = ",theta0 !t0_in(line)
             WRITE(14,"(A,F12.8)") "!R_ref= ",r_norm !t0_in(line)
             WRITE(14,"(A,I5)") "gridpoints = ",res-mod(res,2)
             WRITE(14,"(A,I3)") "n_pol = ",pol_turns_gene
             WRITE(14,"(A,F12.8)") "!fraction of trapped particles = ",ftrap
             WRITE(14,"(A)") "/"
#endif

          ELSE
#if 0
             WRITE(14,"(A)") "&parameters"
             WRITE(14,"(A,F8.4)") "q0 = ",safety
             IF (force_period) THEN
                WRITE(14,"(A,F8.4)") "shat = ",shatgij
             ELSE
                WRITE(14,"(A,F8.4)") "shat = ",slope*safety
             ENDIF
             WRITE(14,"(A,F8.4)") "!rho = ",rho_in(line)
IF (device==6) WRITE(14,"(A,F8.4)") "!edgecorr = ",edgecorr
             WRITE(14,"(A,I5)") "gridpoints = ",res-mod(res,2)
             WRITE(14,"(A,I3)") "n_pol = ",npol
             WRITE(14,"(A)") "!temp_ions,temp_electrons,beta normalized to electron temperature"

             WRITE(14,"(A,F8.4)") "temp_ions = ",temp1/temp2
             WRITE(14,"(A,F8.4)") "temp_electrons = ",1.0
             WRITE(14,"(A,F8.4)") "omt_ions = ",omti
             WRITE(14,"(A,F8.4)") "omt_electrons = ",omte
             !to match quasineutrality condition omn1 is assumed to be equal to omn2
             WRITE(14,"(A,F8.4)") "omn_ions = ",omne
             WRITE(14,"(A,F8.4)") "omn_electrons = ",omne
             !usually: beta = 2.0 * \mu_0 * n k_b T / B_ref^2
             !GENE expects: beta = 2.0 \mu *n_e0 k_b T_Norm / B_ref^2
             WRITE(14,"(A,E12.4)") "beta = ", &
                  2.0*1.2566370614E-6/btor**2*&
                  (dens2*1.60217733E-19)*1E3*temp2
             WRITE(14,"(A,E12.4)") "coll = ", &
                  2.3031E-5*(dens2*1.0E-19)/temp2**2*alpha*&
                  (24.-LOG(SQRT(dens2*1.0E-6)/temp2*0.001))
             !last line is coulomb logarithm for ee, ei collisions
             !as defined in NRL, p.35
             WRITE(14,"(A,F8.4)") "Bref = ", Btor
             WRITE(14,"(A,F8.4)") "Lref = ",alpha
             WRITE(14,"(A,F8.4)") "Tref = ",temp2
             WRITE(14,"(A,E12.4)") "nref = ",dens1/1.e19
             WRITE(14,"(A)") "!Units: Bref in T;Lref in m;Tref in keV;&
                  &nref in 1e19/m**3"
             WRITE(14,"(A)") "!additional information:"
             WRITE(14,"(A,I8)") "!shot = ", shot   
             WRITE(14,"(A,F8.4)") "!time = ", time
             WRITE(14,"(A)") "!Gradients:"
             WRITE(14,"(A,F8.4)") "!R_LTe = ",Rmag/drtotdrho*omte
             WRITE(14,"(A,F8.4)") "!R_LTi = ",Rmag/drtotdrho*omti
             WRITE(14,"(A,F8.4)") "!R_Ln = ",Rmag/drtotdrho*omne
             WRITE(14,"(A)") "!temperatures in keV"
             WRITE(14,"(A,F8.4)") "!T_ions = ",temp1
             WRITE(14,"(A,F8.4)") "!T_electrons = ",temp2
             WRITE(14,"(A)") "!densities in 1e19/meter**3"  
             WRITE(14,"(A,E12.4)") "!n_ions = ",dens1/1.e19
             WRITE(14,"(A,E12.4)") "!n_electrons = ",dens2/1.e19
             WRITE(14,"(A,E12.4)") "debye2 = ",&
                  8.8542E-12*Btor**2/dens2/3.348E-27

             WRITE(14,"(A)") "/" 
#endif
          ENDIF

!=======================================================================================!
!---------------------------------------------------------------------------------------!
!=======================================================================================!

 for_vvbal: IF (vvbal) THEN


          open(33,file='../int/vvbal_'//trim(vmec_file)//'_' &
     &     //trim(table_tag)//'_'//trim(out_tag),action='read') !from prepar.f
          open(34,file='../out/gist_genet_'//trim(vmec_file)//'_' &
     &        //trim(table_tag)//'_'//trim(out_tag),action='write')
     !     open(37,file='../out/curv_stellopt_'//trim(vmec_file)//'_' &
     !&        //trim(table_tag)//'_'//trim(out_tag),action='write')
         

    ALLOCATE(g11_vv(res+1),g12_vv(res+1),g22_vv(res+1),bmod_vv(res+1),dbmod_vv(res),sloc_vv(res))
    ALLOCATE(dBdv1_vv(res+1),dBdv2_vv(res+1),jac_vv(res+1))
    ALLOCATE(curvdr_vv(res+1),s_vv(res+1),bdgrad_vv(res+1))      
    ALLOCATE(g11_gene(res+1),g12_gene(res+1),g22_gene(res+1),bmod_gene(res+1))
    ALLOCATE(dBdx_vv(res+1),dBdxvac_vv(res+1),dBdy_vv(res+1),dbmod_gene(res+1),jac_gene(res+1))
    ALLOCATE(metr_ratio(res+1),bmodaux(res+1))

          DO n=1,res-mod(res,2)+1 !add exta point for dBdz
             READ(33,"(14ES20.10)") s_vv(n),&
                  &g11_vv(n),g12_vv(n),g22_vv(n),&
                  &bmod_vv(n),bdgrad_vv(n),curvdr_vv(n),dBdv1_vv(n),&
                  &dBdv2_vv(n),ftp_vv,qpr_vv,alf_vv,jac_vv(n),pp_vv
          ENDDO


!!===================================================================================================!!
!!                                      GIST for GENE                                                !!
!!===================================================================================================!!
!Auxilliary          
		     safety_sign=safety/abs(safety) !some equilibria might have q < 0
                     s0=s_in ; safety=ABS(q0_vvbal)                      
                     R_norm=alpha !minor radius

!Metrics
	             g11_gene=g11_vv*R_norm**2/4./s0
                     g12_gene=g12_vv*R_norm**2/2./safety*safety_sign
                     g22_gene=g22_vv*R_norm**2*s0/safety**2
                     jac_gene=2*safety*ABS(jac_vv)/R_norm**3 
                     bmod_gene=bmod_vv/ABS(B0)
             !Uniform local shear if required
                     if (sloc_unif) then
                     g22_gene=(bmod_gene**2+g12_gene**2)/g11_gene
                     endif


!Shear
		     shat_vvbal=2*s_in/safety*qpr_vv*safety_sign 

!Inhomogeneity
                     dBdxvac_vv=-R_norm**2*sqrt(s0)/safety*dbdv1_vv*safety_sign		     
                     dBdy_vv=R_norm**2/2./sqrt(s0)*dBdv2_vv !-L_1 : sign inversion in GENE/src/geometry.F90

!Bad curvature
		     dBdx_vv=-R_norm**2*sqrt(s0)/safety*curvdr_vv*safety_sign 

!Pressure gradient
                     my_dpdx=-R_norm**2*sqrt(s0)*2*pp_vv/ftp_vv/ABS(B0)

!Parallel derivative 
	  	     dz = (2.0 * pi * pol_turns) / REAL(res)
                     metr_ratio=g12_gene/g11_gene    
                     bmodaux=bmod_vv
       DO k=1, res
!       dbmod_vv(k) = (bmodaux(k+1)-bmodaux(k))/dz
       sloc_vv(k) = (metr_ratio(k+1)-metr_ratio(k))/(dz)       
       END DO

       dbmod_vv=0.
       DO k=2, res
       dbmod_vv(k) = (bmodaux(k+1)-bmodaux(k-1))/(2.0*dz)
       ENDDO

!Parity (W7X/HSX=-1, NCSX=+1)

                    sign_=parity
                    dBdx_vv=dBdx_vv*sign_
                    dBdxvac_vv=dBdxvac_vv*sign_
                    dBdy_vv=dBdy_vv*sign_
                    my_dpdx=my_dpdx*sign_

!Output file

                   WRITE(34,"(A)") "&parameters"
                   WRITE(34,"(A,3F12.7)") "!(rho0/a,alpha0,theta0) = ",sqrt(s_in),alpha0,theta0
                   WRITE(34,"(A,2F12.7)") "!R,a[m]= ",R0_vmec,r_norm
                   WRITE(34,"(A,F12.7)") "my_dpdx = ",my_dpdx
                   WRITE(34,"(A,F12.7)") "q0 = ",safety
                   WRITE(34,"(A,F12.7)") "shat = ", shat_vvbal
                   WRITE(34,"(A,I5)") "gridpoints = ",res-mod(res,2)
                   WRITE(34,"(A,I3)") "n_pol = ",pol_turns_gene
!                   WRITE(34,"(A,F12.7)") "!fraction of trapped particles = ",ftrap
if (profw7) then
                   WRITE(34,"(A)") "!GRADIENTS for W7" 
                   WRITE(34,"(A,2F12.7)") "!Ti_0, omega_Ti = ",refs(1),grads(1)
                   WRITE(34,"(A,2F12.7)") "!Te_0, omega_Te = ",refs(2),grads(2)
                   WRITE(34,"(A,2F12.7)") "!n_0, omega_n = ",refs(3),grads(3)
endif
                   WRITE(34,"(A)") "/"

                   DO n=1,res-mod(res,2)
                   WRITE(34,"(11ES20.10)") & !s_vv(n)-theta0,&
                  &g11_gene(n),g12_gene(n),g22_gene(n),&
                  &bmod_gene(n),jac_gene(n),&
                  &dBdxvac_vv(n),dBdy_vv(n),&
                  &dbmod_vv(n)/ABS(B0),sloc_vv(n),Bmod_vv(n)*jac_vv(n)/ftp_vv*safety,s_vv(n)
                   ENDDO
!Nunami
!  &-R_norm/(2*sqrt(s0))*dBdv2_vv(n),R_norm*sqrt(s0)/(safety)*dbdv1_vv(n)
     

!!===================================================================================================!!
!!------------------------------------------  BENCHMARKING  -----------------------------------------!!
!!===================================================================================================!!

	IF (benchmarks) THEN

            IF (vvbal .and. (alpha0 .eq. 0) .and. (global .eqv. .false.)) THEN
            OPEN(36,file='../out/bench_'//trim(vmec_file)//'_' &
           &//trim(table_tag)//'_'//trim(out_tag),action='write')
            
            ALLOCATE(dBdx_tr(res),dBdy_tr(res))
            CALL bench
            DEALLOCATE(dBdx_tr,dBdy_tr)

            ELSE IF (alpha0 .ne. 0) THEN
            WRITE(*,"(/,6X,A,//)") "!!! Benchmark only applicable for alpha0=0. !!!"
            ENDIF 

        ENDIF
!!==========================================  DEALLOCATE  ============================================!!
       DEALLOCATE(g11_vv,g12_vv,g22_vv,bmod_vv,dbmod_vv)
       DEALLOCATE(dBdv1_vv,dBdv2_vv,jac_vv,curvdr_vv,s_vv,bdgrad_vv)
       DEALLOCATE(g11_gene,g12_gene,g22_gene,bmod_gene)
       DEALLOCATE(dBdx_vv,dBdxvac_vv,dBdy_vv,jac_gene)
  
ENDIF for_vvbal

       ENDIF first_line

       close(33)
       close(34)

!!====================================================================================================!!
!!------------------------------------------ Output on the screen ------------------------------------!!
!!====================================================================================================!!
if (.not. pest) then
WRITE(*,"(2X,A,/)")  '----------------------------- Surface characteristics ---------------------------' 
 IF (vvbal) THEN
 WRITE(*,"(21X,A,2F8.4,/)") '       iota, shear = ',1./abs(q0_vvbal),shat_vvbal
 ELSE 
 WRITE(*,"(21X,A,2F8.4,/)") '       iota, shear = ',1./abs(q0_tracer),shat_tracer
 ENDIF

WRITE(*,"(2X,A,/)") '------------------------------- GENE normalization ------------------------------' 
 WRITE(*,"(16X,A,2F8.4,/)") '   Bref[T], Lref(minor radius)[m] =', ABS(B0),alpha
endif

if (profw7) then
 WRITE(*,"(4X,A,3F8.4,/)") 'Ti_ref, Te_ref, n_ref =', refs(1),refs(2),refs(3)
 endif

write (*,"(/,2X,A,/)") "-------------------- GIST file successfully generated for GENE -------------------"

!!=====================================================================================================!!

END SUBROUTINE output

!!======================================  DEALLOCATE  =================================================!!
     SUBROUTINE dealloc_out

       DEALLOCATE(c_out,d_out,gij_out,g_ij_out,Bfield_out,Br_out,Bp_out)
       DEALLOCATE(Bz_out,divB_out,jac_out,inv_out,k1_out,k2_out)
       DEALLOCATE(k_norm_out,k_geo_out,kappa_sq_out)
       DEALLOCATE(zero_kappa_out,dBdtau_out)
       DEALLOCATE(zero_cl1_out,zero_cl2_out,zero_cl3_out,r_out,p_out,z_out)
       DEALLOCATE(dBdr_out,dBdp_out,dBdz_out,dBdv1_out,dBdv2_out)
       DEALLOCATE(b3_out,Bt_out,cart_x_out,cart_y_out,sloc,theta_out)
      
     END SUBROUTINE dealloc_out

!!=====================================================================================================!!

   END MODULE output_mod
