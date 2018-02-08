  MODULE fslib_mod
  USE type_mod
  IMPLICIT none

!-----------------------------------------------!
                  CONTAINS
!-----------------------------------------------!

SUBROUTINE viewer_lib
USE specifications_mod, only : alpha0,nz0,pol_turns,vmec_dir,vmec_file, &
                             & alpha0_start,alpha0_end,nalpha0,gviu,pest,boozer, &
                             & sloc_unif,x0min,x0max,x0ref,nsurf,global

  real(dp),dimension(3) :: sflCrd     
  real(dp),dimension(3) :: gradS,gradThetaStar,gradPhi,mag,gradAlpha
  real(dp),dimension(3) :: Bfield,gradB,gradF,gradTh,gradPh
  real(dp),dimension(3) :: es,ea,et,wrk
  real(dp) :: Fa,a,Ba,iot,q,qprim,dpdx
  real(dp) :: th,thetaStar,dtheta,shat
  real(dp) :: gss,g11,g12,g22,Bhat,jac,dBdt,L1,L2,L2_sloc,knorm 
  real(dp) :: g13,g23,c,dBds,dBda
  real(dp) :: gsa, gsa_sloc, gst, gaa, gat, jac1, test, dalpha, dsurf
  real(dp) :: alpha0_start_,xxx
  real, dimension(nsurf) :: x_ind,q_ind,shat_ind,dpdx_ind
  real(dp) :: xaux,iotaux,qprimaux,dxval
  character(10) :: ns_string
  character(40) :: form
  integer  :: nalpha0_
  integer  :: pol_turns_gene
  integer, parameter :: ncols=6 
  logical  :: first=.true.,first1=.true.

  EXTERNAL mcload,mcfree,mcsetb0,mctorflux,mciota,mcreff
  real(dp) ::  mctorflux,mciota,mcreff,mcIotaprime,mcPressurePrime
  integer  :: mcload,mcsavereduced
  character*80 :: fname

  integer(8) mcObj/0/   
  integer(4) FALSE/0/,TRUE/1/
  integer :: i,n,ialpha,is  
  
  real(dp),parameter ::  mu0=1.2566370614e-6
  real(dp) :: phi0,alph0
  integer  :: maxPnt

!Data from gist.inp
  maxPnt=nz0
  pol_turns_gene=pol_turns
  if (nsurf == 1) then
  x0min=x0ref ; x0max=x0ref ; dxval=0.
  else
  dxval = (x0max-x0min)/(nsurf-1)
  endif
  
!Preparation
  fname = trim(vmec_dir)//'/'//trim(vmec_file)  
  if  (mcload (mcObj,fname)==FALSE ) then 
      WRITE(*,"(/,2X,A,/)")   "--------------------- VMEC file not found. GIST will stop. -----------------------"
  stop 
  endif 

  Fa = mctorflux(mcObj,1d0) !Toroidal flux at the lcfs
    a  = mcreff(mcObj,1d0)    !Minor radius
  Ba = ABS(Fa/(pi*a**2))          !Normalization for modB
  dtheta = 2*pi*pol_turns_gene/maxPnt
  dalpha = (alpha0_end - alpha0_start) / real(nalpha0)


!Treat everything as global and make local a special case
if (global) then
nalpha0_=nalpha0
alpha0_start_=-alpha0_start
else
nalpha0_=1
alpha0_start_=alpha0
dalpha=0.
endif

!Keep log for output
do is=1,nsurf
 xaux=x0min+dxval*(is-1)
x_ind(is)=xaux
iotaux = mciota(mcObj,xaux**2)
q_ind(is)=1/iotaux
qprimaux = -mciotaprime(mcObj,xaux**2)*q_ind(is)**2
shat_ind(is)=2*xaux**2/q_ind(is)*qprimaux
dpdx_ind(is)=-4.*xaux/Ba**2 * mcPressurePrime(mcObj,xaux**2) * mu0
enddo

!------------------------------------------------------------------------------!
!-----------------------------------     Main loop    -------------------------!
!------------------------------------------------------------------------------!

global_loop: do is=1,nsurf
xxx=x0min+dxval*(is-1)
iot = mciota(mcObj,xxx**2)
q   = 1/iot
qprim = -mciotaprime(mcObj,xxx**2)*q**2
shat=2*xxx**2/q*qprim
dpdx = -4.*xxx/Ba**2 * mcPressurePrime(mcObj,xxx**2) * mu0


binormal_loop:  do ialpha=1,nalpha0_
  alph0=alpha0_start_ - (ialpha-1) * dalpha
  phi0=alph0

follow_line:  DO i=1, maxPnt
              th = -pi*pol_turns_gene + (i-1)*dtheta

    sflCrd(1) = xxx**2  
    sflCrd(2) = th                              
    sflCrd(3) = phi0 + q*th
    
  if (pest) then
  CALL mcSFLcontraBasis(mcObj,sflCrd,gradS,gradThetaStar,gradPhi,mag)    
  else
  CALL mcBoozerContraBasis(mcObj,sflCrd,gradS,gradThetaStar,gradPhi,mag)
  endif

    thetaStar = sflCrd(2)
    gradAlpha=qprim*thetaStar*gradS+q*gradThetaStar-gradPhi 
!--------------------------------------------------------------------------------!

! Metrics and jacobian   
    gss =  dot_product(gradS,gradS)         
    gsa =  dot_product(gradS,gradAlpha)     
    gst =  dot_product(gradS,gradThetaStar)     
    gaa =  dot_product(gradAlpha,gradAlpha) 
    gat =  dot_product(gradAlpha,gradThetaStar) 
    !for uniform local shear if requested
    gsa_sloc=shat*thetaStar*gss/(2.*xxx**2*iot)

    CALL cross_product(gradS,gradAlpha, wrk) 
    jac1 =  1/dot_product(wrk,gradThetaStar)

! Gradient of B wrt cylindrical coordinates
   CALL mcgetBandGradients(mcObj,mag,Bfield,gradB,gradF,gradTh,gradPh)
    gradB(:)  = gradB(:)/Ba
    Bfield(:) = Bfield(:)/Ba
    Bhat = sqrt(dot_product(Bfield,Bfield))

! GENE metrics
    g11 = gss*a**2/(4*xxx**2)
    if (sloc_unif) then ! Apply uniform local shear if requested
    g12 =  gsa_sloc*a**2*iot/2    
    else
    g12 =  gsa*a**2*iot/2    
    endif
    !g22 =  gaa*a**2*s*iot**2 ! standard
    g22=(Bhat**2+g12**2)/g11 ! using field-alignment property
    jac =  jac1*2*q/a**3

!Covariant basis wrt cylindrical coordinates
    CALL cross_product(gradAlpha,gradThetaStar, es) 
    CALL cross_product(gradThetaStar,gradS,     ea) 
    CALL cross_product(gradS,gradAlpha,         et)  

    es(:) = es(:)*jac1
    ea(:) = ea(:)*jac1
    et(:) = et(:)*jac1
        
    dBds = dot_product(gradB,es)
    dBda = dot_product(gradB,ea)
    dBdt = dot_product(gradB,et)

!Curvatures
    c = iot*iot*a**4
    L1 = q/xxx*(dBda + c*(gss*gat-gsa*gst)*dBdt/(4*Bhat**2))
    L2 = 2*xxx*(dBds + c*(gaa*gst-gsa*gat)*dBdt/(4*Bhat**2))
    knorm=(sqrt(g11)*L2-gsa*a**2*iot/2*L1)/Bhat
    L2_sloc=(Bhat*knorm+g12*L1)/sqrt(g11)
    if (sloc_unif) L2=L2_sloc

!Output file
       if (first1) then
              WRITE(58,"(A)") "&parameters"
       if (pest) then
              WRITE(58,"(A,3I7,3X,A,3X,A)") "(#lines, #points, #cols, coordinates, code) =", &
              &       nalpha0, MaxPnt, ncols, "PEST", "GIST" 
       else if (boozer) then
              WRITE(58,"(A,3I7,3X,A,3X,A)") "(#lines, #points, #cols, coordinates, code) =", &
              &              nalpha0, MaxPnt, ncols, "BOOZER", "GIST"
       endif 

                   WRITE(58,"(A,I4)") "#surfaces = ",nsurf
                   write(ns_string, '(I4)') nsurf ! for format
                   form="(A,"//trim(ns_string)//"F12.7)" ! for format
                   WRITE(58,form) "s0 = ", x_ind**2
                   WRITE(58,form) "q0 = ", q_ind
              WRITE(58,"(A)") "/"
              WRITE(58,"(8(A,X))") "theta","y","modB","kbad","kgeo","g11","g12","g22"       
              first1=.false.
       endif

            WRITE(58,"(8ES20.8)") th,phi0,Bhat*Ba,L2-dpdx/2/Bhat,L1,g11,g12,g22

ENDDO follow_line

ENDDO binormal_loop

ENDDO global_loop

!Output screen
write (*,"(/,2X,A,/)") "------------------------- GViU file successfully generated ----------------------"

!---------------------------------    
CALL mcfree(mcObj)  
       
END SUBROUTINE viewer_lib

!********************************************!
  SUBROUTINE cross_product(A,B,AxB)
  IMPLICIT none 
  REAL(dp),INTENT(in) :: A(3),B(3)
  REAL(dp),INTENT(out) :: AxB(3)

  AxB = (/ (A(2)*B(3)-A(3)*B(2)), & 
           (A(3)*B(1)-A(1)*B(3)), & 
          (A(1)*B(2)-A(2)*B(1)) /)
 
 END SUBROUTINE cross_product
  

!====================================================!

END MODULE fslib_mod

