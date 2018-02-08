
MODULE benchmark_mod
USE type_mod
implicit none

real(dp),dimension(:),allocatable::tht_booz,zet_booz 
real(dp),dimension(:),allocatable::g11_vv,g12_vv,g22_vv,bdgrad_vv,curvdr_vv
real(dp),dimension(:),allocatable::dBdv1_vv,dBdv2_vv,jac_vv,s_vv
real(dp),dimension(:),allocatable::g11_gene,g12_gene,g22_gene,dBdx_vv,dBdy_vv,bmod_vv
real(dp),dimension(:),allocatable::dBdx_tr,dBdy_tr

END MODULE benchmark_mod