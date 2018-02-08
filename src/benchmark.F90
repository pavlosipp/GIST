!!=================================================================================!!
!! This Module deals with the comparison tests between the MHD and tracing methods !!
!!=================================================================================!!

MODULE benchmark
USE benchmark_mod
USE specifications_mod
implicit none

!=======================================================!
                       CONTAINS
!=======================================================!

SUBROUTINE bench
integer::m,n
real(dp)::s0,q0

s0=s_in
q0=safety

CALL curvs_tracer(dBdx_tr,dBdy_tr)

write(36,"(A)") "#******************  Benchmarking file  ******************"
write(36,"(A)") "#theta g11_T g^ss g12_T g^sa g22_T g^aa B_T B_B dBdx_T dBdx_B dBdy_T dBdy_B "//&
"bad_T bad_B geo_T geo_B" 

do m=1,resol
write(36,"(17F12.7)") &
p_out(m)/q0,& !same with s_vv for alpha0=0
gij_out(1,1,m),alpha**2/4./s0*g11_vv(m), &
gij_out(1,2,m),alpha**2/2./safety*g12_vv(m), &
gij_out(2,2,m),alpha**2*s0/q0**2*g22_vv(m),&
Bfield_out(m)/ABS(b0),Bmod_vv(m)/ABS(b0),&
alpha*k1_out(m)*ABS(Bfield_out(m)/b0),-alpha**2*sqrt(s0)/q0*dBdv1_vv(m)*parity,&
alpha*k2_out(m)*ABS(Bfield_out(m)/b0),alpha**2/2./sqrt(s0)*dBdv2_vv(m)*parity,&
alpha*dBdx_tr(m)*ABS(Bfield_out(m)/b0),-alpha**2*sqrt(s0)/q0*curvdr_vv(m)*parity,&
alpha*dBdy_tr(m)*ABS(Bfield_out(m)/b0),alpha**2/2./sqrt(s0)*dBdv2_vv(m)*parity
enddo

END SUBROUTINE bench
!==========================================================!
SUBROUTINE curvs_tracer(k_1,k_2)

integer :: k,m
real(dp) :: q,h 
real(dp), dimension(:),allocatable ::  jac,r,dr,dz,rp,zp,dr2,dz2,rpp,zpp,aux,tr,tp,tz
real(dp), dimension(:),allocatable ::  bf,g11,g12,kgeo,knorm,kr,kp,kz,c11,c12,c22,c21
real(dp), dimension(:),allocatable ::  k_1vv,k_2vv,curvdr,curvgeo
real(dp), dimension(:),allocatable,intent(out) :: k_1,k_2

allocate(dr(resol),dz(resol),rp(resol),zp(resol),dr2(resol),dz2(resol),rpp(resol),zpp(resol))
allocate(aux(resol),tr(resol),tp(resol),tz(resol))
allocate(r(resol),jac(resol),bf(resol),g11(resol),g12(resol),kgeo(resol))
allocate(kr(resol),kp(resol),kz(resol),c11(resol),c12(resol),c22(resol))
allocate(knorm(resol),k_1(resol),k_2(resol),k_1vv(resol),k_2vv(resol),curvdr(resol+1),curvgeo(resol),c21(resol))

!local definitions
q=safety
r=r_out(1:resol)
jac=jac_out(1:resol)*q !jac_out=jacobian wrt phi, jac=jacobian wrt chi=phi/q
bf=bfield_out(1:resol)
g11=gij_out(1,1,1:resol)
g12=gij_out(1,2,1:resol)
c11=c_out(1,1,1:resol)
c12=c_out(1,2,1:resol)
c21=c_out(2,1,1:resol)
c22=c_out(2,2,1:resol)

h=2*pi/real(resol-1)

!first derivative of cylindrical coords
do k=1,resol
dr(k)=r_out(k+1)-r_out(k-1)
dz(k)=z_out(k+1)-z_out(k-1)
enddo

 DO k=1,resol
       rp(k) = dr(k)/(2.*h)/q !!NB: Tracer gives r(chi), so dr_dphi=dr_dchi/q. Here rp=dr_dphi 
       zp(k) = dz(k)/(2.*h)/q 
    END DO

!second derivative
do k=1,resol
dr2(k)=r_out(k+1)+r_out(k-1)-2*r_out(k)
dz2(k)=z_out(k+1)+z_out(k-1)-2*z_out(k)
enddo

DO k=1,resol
       rpp(k) = dr2(k)/h**2/q**2 
       zpp(k) = dz2(k)/h**2/q**2 
    END DO

!physical cylindrical components of tangent vector

aux=sqrt(rp**2+r**2+zp**2)
tr=rp/aux ; tp=r/aux ; tz=zp/aux

!physical cylindrical components of curvature vector

kr=(rpp*zp**2+rpp*r**2-2*r*rp**2-r*zp**2-r**3-zp*zpp*rp)/aux**4
kp=(2*rp**3+2*rp*zp**2+rp*r**2-r*rp*rpp-zp*zpp*r)/aux**4
kz=(zpp*rp**2+r**2*zpp-rp*zp*rpp-r*rp*zp)/aux**4

!geodesic curvature
kgeo=jac/q/r*bf/abs(B0)*(c12*kr-c11*kz)/sqrt(g11)

!normal curvature
knorm=jac/q/r/sqrt(g11)*(kr*(g11*c22-g12*c12)+kz*(c11*g12-c21*g11))

!covariant components of curvature vector wrt Tracer system
k_1=jac/q/r*(c22*kr-c21*kz)
k_2=jac/q/r*(c11*kz-c12*kr)

END SUBROUTINE curvs_tracer
!=====================================================================!


END MODULE benchmark

