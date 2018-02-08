
MODULE specifications_mod
  USE type_mod
  IMPLICIT NONE

!generic
LOGICAL :: nodrift=.false., notrap=.false.

!vmec
 INTEGER :: ns_vmec,mpnt_vmec
 REAL(DP) :: R0_vmec

!trabal
 LOGICAL :: initialize=.TRUE., vvbal=.TRUE., gs2=.FALSE., tracer=.TRUE., midplane=.false.,profw7=.false.,global=.false. 
 LOGICAL :: gviu=.false., stellopt=.false., benchmarks=.false., terpsichore=.false.
 LOGICAL :: gs2_press=.TRUE., pest=.false., boozer=.false., fake3d=.false.
 character(250) :: table_tag="highres", vmec_file="", out_tag="", vmec_dir='./'
 integer :: surf=0,nz0=-1,dir_tr
 character(30) :: name
 real(dp) :: theta0=0.,q0_vvbal,alpha0=0.,shat_vvbal,q0_tracer,shat_tracer,rho_tr=0.5,t_tr=0.,z_tr=0.,s_tr, x0ref=-1, x0min=-1, x0max=-1
 integer :: nalpha=1,nalpha0=0 
 integer :: nsurf=-1
 real(dp) :: alpha0_start,alpha0_end, alpha0_out 
 real(dp), DIMENSION(:), ALLOCATABLE :: alpha0_vec 
 real :: pol_turns=1.

!pest.f90
  LOGICAL :: sloc_unif=.FALSE., sloc_curv=.FALSE., jorge=.FALSE., wilcox=.FALSE.

!prepar.f
 real(dp)::s_tracer,s_vvbal,parity


!read_input
 real(dp) :: s_in

!common 
  !for parallelization
  INTEGER :: g_rank
  INTEGER :: g_num_proc
!read_parms

  INTEGER :: device=2,code=1,cuts=1,num_lines=1,read_mode=0,cs_gs2=1
  INTEGER :: out_mode = 0
  LOGICAL :: silent_mode
  REAL(DP) ::rmag=0.,z0,zmag=0.,Bref_gs2=1.,Rref_gs2=1.
  INTEGER, DIMENSION(:), ALLOCATABLE :: res,symm
  INTEGER, DIMENSION(:), ALLOCATABLE :: trace_type 
  REAL(DP), DIMENSION(:), ALLOCATABLE :: t0_in,r_in,z_in,p0_in,p1_in,rho_in
  REAL(SP):: tSHOT=1.,edgecorr=1.
  LOGICAL :: dbswitch=.FALSE.,vmec=.true.,force_period=.TRUE.
  CHARACTER(10) :: config,code_name,expaug='aaa',diag='aaa'

!for FINDIF code only
  ! taken from input
  INTEGER :: out_gik=1 ! 0, 1, 2
  INTEGER :: in_type=2 ! 0 or 2
	! nb_tor is the array of toroidal neighbors of primary point i
	INTEGER, DIMENSION(:,:), ALLOCATABLE :: nb_tor

  !For output (filled in load mode)
  INTEGER :: num_startp

!read_boris

  REAL(DP), DIMENSION(:), ALLOCATABLE :: iota_boris 

!prep_cuts

  INTEGER, DIMENSION(:), ALLOCATABLE :: cut
  INTEGER :: resol

!matrix_d

  REAL(DP) :: detc
  REAL(DP), DIMENSION(3,3) :: c,d

!gij_contravar,gij_covar

  REAL(DP), DIMENSION(3,3) :: gij,g_ij

!quantities

  REAL(DP) :: Bfield,Br,Bp,Bz,Bt,divB,jac,inv
  REAL(DP) :: k1,k2,k_norm,k_geo,kappa_sq
  REAL(DP) :: zero_kappa,zero_cl1,zero_cl2,zero_cl3
  REAL(DP) :: b3,dBdz,dBdr,dBdp,dBdv1,dBdv2,dBdtau
  REAL(DP) :: dBrdr,dBrdp,dBrdz,dBpdr,dBpdp,dBpdz,dBzdr,dBzdp,dBzdz
  REAL(DP) :: jac_chain,jac_contravar,jac_covar
  REAL(DP) :: fac1,fac2,fac3,kappa_sq_imp,cart_x,cart_y,theta

!tracer

  REAL(DP) :: r,phi,dphi,z,dist,r_norm,rmin,rmax,zmin,zmax
  INTEGER :: count=0,genecnt=0

!prep_out

  REAL(DP), DIMENSION(:), ALLOCATABLE :: r_out,p_out,z_out,cart_x_out,cart_y_out,theta_out
  REAL(DP), DIMENSION(:), ALLOCATABLE :: Bfield_out,Br_out,b3_out,Bt_out
  REAL(DP), DIMENSION(:), ALLOCATABLE :: Bp_out,divb_out,zero_kappa_out,zero_cl1_out
  REAL(DP), DIMENSION(:), ALLOCATABLE :: zero_cl2_out,zero_cl3_out,dBdr_out,dBdp_out
  REAL(DP), DIMENSION(:), ALLOCATABLE :: dBdz_out,dBdtau_out,dBdv1_out,dBdv2_out
  REAL(DP), DIMENSION(:), ALLOCATABLE :: Bz_out,jac_out,k1_out,k2_out
  REAL(DP), DIMENSION(:), ALLOCATABLE :: k_norm_out,k_geo_out,kappa_sq_out,inv_out
  REAL(DP), DIMENSION(:), ALLOCATABLE :: g12_b_out,g22_b_out,g23_b_out,jac_b_out
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: gij_out,g_ij_out,c_out,d_out
  INTEGER :: resol_f,resol_b

!load_mod

  REAL(DP) :: temp1,omti,dens1,omni,temp2,omte,dens2,omne,alpha,Btor,beta,r_gs2,b0,ftrap
  REAL(DP), DIMENSION(:), ALLOCATABLE :: safetyf,shatf,c11_rho
  REAL(DP), DIMENSION(3) :: grads, refs

!gs2
  REAL(DP) :: aref_bel,bref_bel

!vmec_iota

  INTEGER :: npol

!vmec_mod

  INTEGER :: sign_of_bp
  REAL :: qpr

!output

  REAL(DP) :: safety,shat,drtotdrho,shatgij
  INTEGER :: res_out
!Final output
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: gij_outf
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: det_g_outf,det_g_2d_outf
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Bfield_outf,b3_outf
  

!eqintrp

  LOGICAL :: vmec_6_90=.FALSE.

!shear

  REAL(DP) :: slope
  REAL(DP), DIMENSION(:), ALLOCATABLE:: sloc

!d3d_iterdb and data_aug

  INTEGER :: shot=1,ishot
  REAL(DP):: time

!aug
  REAL(SP) :: Phi_Sep
  REAL(DP) :: Psi_Sep,Psi_Ax,zsym
  LOGICAL :: swsymm
END MODULE specifications_mod

