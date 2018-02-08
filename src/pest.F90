MODULE pest_mod
  USE vmec_mod
  USE type_mod
  IMPLICIT none

  integer(8),PRIVATE::mcObj=0 

!-----------------------------------------------!
                  CONTAINS
!-----------------------------------------------!

SUBROUTINE gene_pest
USE specifications_mod, only :alpha0,nz0,pol_turns,vmec_dir,vmec_file, &
                            & alpha0_start,alpha0_end,nalpha0,gviu,pest,boozer,nodrift,notrap, &
                            & nsurf,global,x0min,x0max,fake3d

  real(dp),dimension(3) :: sflCrd,cylCrd
  real(dp),dimension(3) :: gradS,gradThetaStar,gradPhi,mag,gradAlpha
  real(dp),dimension(3) :: Bfield,gradB,gradF,gradTh,gradPh
  real(dp),dimension(3) :: es,ea,et,wrk
  real(dp) :: Fa,a,Ba,iot,q,qprim,dpdx,alpha,beta
  real(dp) :: th,thetaStar,dtheta,shat
  real(dp) :: gss,g11,g12,g22,g33,Bhat,jac,dBdt,L1,L2,L2_sloc,knorm 
  real(dp) :: g13,g23,c,dBds,dBda
  real(dp) :: gsa, gst, gtt,gsa_sloc,gaa,gat,jac1,jac2,temp,test,dalpha,dxval
  real(dp) :: alpha0_start_, zero=0.
  real(dp) :: xxx,xaux,iotaux,qprimaux,x0ref
  integer  :: nalpha0_
  integer  :: pol_turns_gene
  integer, parameter :: ncols=5 
  integer  :: nfp
  character(4) :: firstchar
  character(20) :: dumchar
  real, dimension(7) :: dumreal
  real, dimension(nsurf) :: x_ind,q_ind,shat_ind,dpdx_ind
  real :: q0_ref, s0_ref
  integer :: midsurf
  character(10) :: ns_string
  character(40) :: form
  integer :: dumint
  logical  :: first=.true., arc=.false.
  real(dp) :: dtharc, tharc, larc, arcth0
  real(dp) :: arc1, arc2, arc3, arc1_0, arc2_0, arc3_0
  real(dp), dimension(3) :: cylarc, cylarc0
  integer :: iarc
  integer, parameter :: Narc=100 

  EXTERNAL mcload,mcfree,mcsetb0,mctorflux,mciota,mcreff,mcR0,mcBoozerBmn
  real(dp) ::  mctorflux,mciota,mcreff,mcIotaprime,mcPressurePrime,mcPressure,mcR0,mcBoozerBmn
  integer  :: mcload,mcsavereduced
  character*80 :: fname

  integer(8) mcObj/0/   
  integer(4) FALSE/0/,TRUE/1/
  integer :: i,n,ialpha,is
  
  real(dp),parameter :: mu0=1.2566370614e-6
  real(dp) :: phi0,R0
  real(dp) :: alph0
  integer  :: maxPnt


!Set number of periods from equilibrium file (VMEC or EFIT)
open(8,file=trim(vmec_dir)//'/'//trim(vmec_file),action="read", err=100)
read(8,*) firstchar
close(8) 

!Differentiate between VMEC, Boozer file and EFIT
open(8,file=trim(vmec_dir)//'/'//trim(vmec_file),action="read", err=100)
if (firstchar == 'VMEC') then
        read (8,*) dumchar, dumchar, dumchar, dumchar
        read (8,*) (dumreal(i),i=1,7) 
        read (8,*) nfp
elseif (firstchar == 'EFIT') then
nfp=1 !tokamak case
elseif (firstchar == 'CC') then !Boozer file
        read (8,*)
        read (8,*)
        read (8,*)
        read (8,*)
        read (8,*)
        read (8,*) dumint, dumint, dumint, nfp
    else
WRITE(*,"(/,2X,A,/)")   "-----------------  Equilibrium format unknown. GIST will stop.  ------------------"
STOP
endif
close(8)    

!Define alpha0_start and alpha0_end based on nfp 
alpha0_start=-pi/nfp ; alpha0_end=pi/nfp

!Data from gist.inp
  maxPnt=nz0
  pol_turns_gene=pol_turns
  if ((nsurf == 1) .or. fake3d) then
   x0max=x0min ; dxval=0.
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
  R0=mcR0(mcObj)            !Major radius
  Ba = ABS(Fa/(pi*a**2))    !Normalization for modB

  dtheta = 2*pi*pol_turns_gene/maxPnt
  dalpha = 2*pi / nfp / real(nalpha0)


!Treat default global and make local a special case
if (global) then
nalpha0_=nalpha0
alpha0_start_=-alpha0_start !GENE convention: we go from plus to minus
else
nalpha0_=1
alpha0_start_=alpha0
dalpha=0.
endif



!Keep log for output
do is=1,nsurf
xaux=x0min +dxval*(is-1)
x_ind(is)=xaux
iotaux = mciota(mcObj,xaux**2) !mciota is a fuction of s=x^2
q_ind(is)=1/iotaux
qprimaux = -mciotaprime(mcObj,xaux**2)*q_ind(is)**2 !q'(s) 
shat_ind(is)=2*xaux**2/q_ind(is)*qprimaux
dpdx_ind(is)=-4.*xaux/Ba**2 * mcPressurePrime(mcObj,xaux**2) * mu0  
enddo

!Reference values: s0_ref q0_ref
!If nsurf is odd then s0 is the mid surface, if nsurf is even, then s0 the surface closest to (inexisting) middle
if (MOD(nsurf,2) == 0) then !even no. of surfaces
midsurf=nsurf/2
elseif (MOD(nsurf,2) == 1) then !odd no. of surfaces
midsurf=nsurf/2+1
endif
!s0_ref=s_ind(midsurf)
x0ref=x_ind(midsurf)
q0_ref=1/mciota(mcObj,x0ref**2)

!------------------------------------------------------------------------------!
!-----------------------------------     Main loop    -------------------------!
!------------------------------------------------------------------------------!

global_loop: do is=1,nsurf
   xxx=x0min + dxval*(is-1)
iot = mciota(mcObj,xxx**2)
q   = 1/iot
qprim = -mciotaprime(mcObj,xxx**2)*q**2 
shat=2*xxx**2/q*qprim
dpdx = -4.*xxx/Ba**2 * mcPressurePrime(mcObj,xxx**2) * mu0  
beta=1/Ba**2*mcPressure(mcObj,xxx**2)*mu0 !!Total beta is two times this but we need local for GENE

binormal_loop: DO ialpha=1,nalpha0_
  alph0=alpha0_start_ - (ialpha-1) * dalpha !Fix line
  phi0=alph0 !alph=phi-q*(theta-theta0) (Yuriy) ; then, alph0=phi0 

follow_line:  DO i=1, maxPnt
              th = -pi*pol_turns_gene + (i-1)*dtheta !theta = parallel coordinate

    sflCrd(1) = xxx**2  
    sflCrd(2) = th                               
    sflCrd(3) = phi0 + q*th !equation of line
  
!-------------------------------------------------------------------------------!
!Create Boozer or PEST coordinates

  if (pest) then
  CALL mcSFLcontraBasis(mcObj,sflCrd,gradS,gradThetaStar,gradPhi,mag)    
  else !Boozer
  CALL mcBoozerContraBasis(mcObj,sflCrd,gradS,gradThetaStar,gradPhi,mag)
  endif
  
  CALL mcmag2cyl(mcObj,sflCrd,cylCrd) !extract cylindrical coordinates
!--------------------------------------------------------------------------------!

    thetastar = sflCrd(2)
    gradAlpha=qprim*thetaStar*gradS+q*gradThetaStar-gradPhi !alpha=q*theta - phi

!Metrics
    gss =  dot_product(gradS,gradS)         
    gsa =  dot_product(gradS,gradAlpha)     
    gst =  dot_product(gradS,gradThetaStar)     
    gaa =  dot_product(gradAlpha,gradAlpha) 
    gat =  dot_product(gradAlpha,gradThetaStar)
    gtt =  dot_product(gradThetaStar,gradThetaStar)

!Jacobian
     CALL cross_product(gradS,gradAlpha, wrk) 
     jac1 =  1/dot_product(wrk,gradThetaStar)     

! Gradient of B wrt cylindrical coordinates
   CALL mcgetBandGradients(mcObj,mag,Bfield,gradB,gradF,gradTh,gradPh)
    gradB(:)  = gradB(:)/Ba
    Bfield(:) = Bfield(:)/Ba
    Bhat = sqrt(dot_product(Bfield,Bfield))

! GENE metrics
    g11 = gss*a**2/(4*xxx**2) 
    g12 =  gsa*a**2*iot/2
    g22 =  gaa*a**2*xxx**2*iot**2 
    g13 = 0.5*a**2/xxx*gst
    g23 = a**2*xxx*iot*gat
    g33 = a**2*gtt    
    jac =  jac1*2*q/a**3
!alternative calculation of jacobian via metrics; gives same answer with jac
!    temp=g11*g22*g33+2.*g12*g23*g13-(g11*g23**2+g22*g13**2+g33*g12**2)
!    jac=sqrt(1./abs(temp))
                
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

  if (notrap) dBdt=0.

! Curvatures
    c = iot*iot*a**4
    L1 = q/xxx*(dBda + c*(gss*gat-gsa*gst)*dBdt/(4*Bhat**2))
    L2 = 2*xxx*(dBds + c*(gaa*gst-gsa*gat)*dBdt/(4*Bhat**2))   
    !knorm=(g11*L2+gsa*a**2*iot/2*L1)/Bhat !L1=-kgeo
    !knorm=(g11*L2+g12*L1)/Bhat !L1=-kgeo
    !L2_sloc=(Bhat*knorm-g12*L1)/g11

! Remove drift      
 if (nodrift)  then
     L2=0. 
    dpdx=0. 
 endif


! Output file
 if (first) then
                   WRITE(56,"(A)") "&parameters"
                   if (pest .and. (nfp .ne. 1)) then 
                    WRITE(56,"(A)") "!PEST coordinates"
                   else if (boozer .and. (nfp .ne. 1)) then
                    WRITE(56,"(A)") "!BOOZER coordinates"
                   endif
                   WRITE(56,"(A,2F12.7)") "!major radius[m], minor radius[m]= ",R0, a
             if (global .and. (nsurf >=2)) then
                   WRITE(56,"(A)") "!----------------------    GLOBAL mode   ----------------------"
                   WRITE(56,"(A,2F12.7)") "!alpha0_start, alpha0_end = ", &
                                            &alpha0_start, alpha0_end
                   WRITE(56,"(A,I3)") "n0_global = ",nfp
                   WRITE(56,"(A,I4)") "nalpha0 = ",nalpha0
                   WRITE(56,"(A,F12.7)") "x0 = ", x0ref
                   WRITE(56,"(A,F12.7)") "q0 = ", q0_ref
                   WRITE(56,"(A,I4)") "ns = ",nsurf
                   write(ns_string, '(I4)') nsurf ! for format
                   form="(A,"//trim(ns_string)//"F12.7)" ! for format
                   WRITE(56,form) "xval_a = ", x_ind
                   WRITE(56,form) "q_prof = ", q_ind
                   WRITE(56,form) "dqdx_prof = ", shat_ind
                   WRITE(56,form) "dpdx_pm_arr = ", dpdx_ind
                elseif (global .and. (nsurf == 1)) then
                   WRITE(56,"(A)") "!----------------------    SURFACE mode   ----------------------"
                   WRITE(56,"(A,2F12.7)") "!alpha0_start, alpha0_end = ", &
                                            &alpha0_start, alpha0_end
                   WRITE(56,"(A,I3)") "n0_global = ",nfp
                   WRITE(56,"(A,I4)") "nalpha0 = ",nalpha0
                   WRITE(56,"(A,F12.7)") "s0 = ", x0min**2
                   WRITE(56,"(A,F12.7)") "q0 = ", q0_ref
                   write(ns_string, '(I4)') nsurf ! for format
                   form="(A,"//trim(ns_string)//"F12.7)" ! for format
                   WRITE(56,form) "shat = ", shat_ind
                   WRITE(56,form) "my_dpdx = ", dpdx_ind
                else
                   WRITE(56,"(A,3F12.7,I5)") "!---------   FLUX-TUBE   ------------"
                   WRITE(56,"(A,5F12.7)") "!Bref, B_00, B_01, B_10, B_11 = ",Ba, mcBoozerBmn(mcObj,0,0,xxx**2), &
                   & mcBoozerBmn(mcObj,0,1,xxx**2),mcBoozerBmn(mcObj,1,0,xxx**2),mcBoozerBmn(mcObj,1,1,xxx**2)
                   WRITE(56,"(A,F6.3)") "s0 = ",xxx**2
                   WRITE(56,"(A,F6.3)") "!alpha0 = ",alpha0
                   WRITE(56,"(A,F12.7)") "!beta = ",beta
                   WRITE(56,"(A,F12.7)") "my_dpdx = ",dpdx
                   WRITE(56,"(A,F12.7)") "q0 = ",ABS(q)                   
                   WRITE(56,"(A,F12.7)") "shat = ",shat            
                endif
                   WRITE(56,"(A,I5,A)") "gridpoints = ",maxPnt," !Parallel resolution"
                   WRITE(56,"(A,I3,A)") "n_pol = ",pol_turns_gene," !Poloidal turns"
                   WRITE(56,"(A)") "/"

first=.false. 
endif
                   
                   WRITE(56,"(8ES20.10)") g11,g12,g22,Bhat,ABS(jac),&
                                          &L2,L1,dBdt !,cylCrd(1)
 
                   !! The bad curvature is given by L2-dpdx/2./Bhat
ENDDO follow_line

ENDDO binormal_loop

ENDDO global_loop

! Circumference of a cut
if (arc) then
   arcth0=-pi
   dtharc=2.*pi/Narc
      !initialize
      larc=0.
      call arclength(.5,arcth0,pi/3.,cylarc0)

  do iarc=1,Narc
  arc1_0=cylarc0(1) ; arc2_0=cylarc0(2) ; arc3_0=cylarc0(3)
  tharc=arcth0+iarc*dtharc
  call arclength(.5,tharc,pi/3.,cylarc)
  arc1=cylarc(1) ; arc2=cylarc(2) ; arc3=cylarc(3)
  larc=larc+sqrt((arc1-arc1_0)**2+(arc2-arc2_0)**2+(arc3-arc3_0)**2)
  cylarc0=cylarc
  enddo
  print *,"larc=",larc
  endif

! Output on screen

if (.not. global) then
 WRITE(*,"(22X,A,2F9.4,//)") '    iota, shear = ',ABS(iot),shat
endif
              
! WRITE(*,"(2X,A,/)") '--------------------------------- GENE --------------------------------' 
! WRITE(*,"(18X,A,2F9.4)") ' Bref[T], Lref(minor radius)[m] =',Ba,a

write (*,"(/,2X,A,/)") "------------------------- GIST file successfully generated -----------------------"

!---------------------------------    
CALL mcfree(mcObj)         
return
100    WRITE(*,"(/,2X,A,/)")   "--------------------- VMEC file not found. GIST will stop. -----------------------" 
stop
END SUBROUTINE gene_pest

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

END MODULE pest_mod

