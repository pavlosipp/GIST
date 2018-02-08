
      program gist
      USE WMPI
      USE pest_mod
      USE type_mod
      USE tracer_mod
      USE field_mod
      USE rk5_mod
      USE load_mod
      USE output_mod
      USE specifications_mod
      use tprb_param
      use cidata
      use crdata
      use ceigen
      use cresis
      use geneyout_mod 
      use gs2_mod
      use fslib_mod

      implicit none

  INTERFACE
     SUBROUTINE rhs(phi,y,yp)             
       USE type_mod
       IMPLICIT NONE
       REAL(DP), dimension(:), INTENT(OUT) :: yp
       REAL(DP), dimension(:), INTENT(IN)  :: y
       REAL(DP), INTENT(IN) :: phi
     END SUBROUTINE rhs
  
     SUBROUTINE rhs1(phi,y,yp)             
       USE type_mod
       IMPLICIT NONE
       REAL(DP), dimension(:), INTENT(OUT) :: yp
       REAL(DP), dimension(:), INTENT(IN)  :: y
       REAL(DP), INTENT(IN) :: phi
     END SUBROUTINE rhs1

  END INTERFACE

  REAL(DP) :: pend_b,pend_f,hn1,r0_aug,c11,err,Bp0,phi_period
  INTEGER :: index,ierr=0,geneend=1
  INTEGER :: ij,kk,ll,nEDIT=0,j,i
  INTEGER, PARAMETER :: neqn=8 !Number of ODE's !do not change!
  REAL(DP), DIMENSION(neqn):: y 
  REAL(DP), DIMENSION(2):: w !This solves only the field line equation 
  REAL(DP), PARAMETER :: hn1min=0.
  REAL(SP) :: rhoPF,rhoTF,fPF,fTF,rrr,zzz
  REAL(DP) :: dalpha,alpha0_start_ 
  INTRINSIC TINY, MOD
   real :: tim (20)  

      external EQINVM, mtaskb, second, veqrec
      external datain, driver 

!=======================================================================!

!!Read parameters

    NAMELIST /coordinates/ &
     &pest,boozer

    NAMELIST /in_out/ &
     &  itrmax ,ci ,mercie &                
     &  ,xmax   ,abserr ,relerr ,gmin   ,gmax   ,ginit  &                
     &  ,dg1    ,dg2    ,gscal1 ,gscal2 ,jstart ,jstop  &                
     &  ,lresis ,sresis ,dlta   ,sk     ,x0    &                
     &  ,modes  ,lmodel ,ldiag  ,nstart   ,nstop, benchmarks  &
     &  ,initialize, name, gs2 , gs2_press, vvbal, vmec_6_90, tracer &
     &  ,vmec_file, table_tag, out_tag, vmec_dir, global, gviu, stellopt &
     &  ,terpsichore, nodrift,notrap, sloc_unif, sloc_curv, jorge, wilcox, fake3d
   
                   
    NAMELIST /setup/ &
    &   pol_turns,theta0,surf,rho_tr,z_tr,t_tr,alpha0,nz0,profw7, &
    &   nalpha0,alpha0_start,alpha0_end,x0ref,x0min,x0max,nsurf
         
    open(1,file='../inp/gist.inp',action="read")
    READ(1, NML=coordinates)        
    READ(1, NML=in_out)
    READ(1, NML=setup)                
    close(1)

!=======================================================================!
write (*,"(2X,A)") "================================================================================="
write (*,"(2X,A)")       "=========        Geometry Interface for Stellarators and Tokamaks       ========="
write (*,"(2X,A)")       "=========          (P.Xanthopoulos, W.A.Cooper, and Yu.Turkin)          ========="
write (*,"(2X,A)")       "=========          pax@ipp.mpg.de      www.pavlosipp.com                ========="
write (*,"(2X,A)") "================================================================================="


!==============================================================================================!
!------------------------------------------  LIBRARY ------------------------------------------!
!                                          ONLY for GENE                                       !
!==============================================================================================!



!--------------------------------------   Starting SANITY checks   -------------------------------------!

IF ((.not. pest) .and. (.not. boozer)) then
write (*,"(/,2X,A,/)") "--------------- No output. Please select Boozer or PEST coordinates. ------------"
stop
ENDIF


IF (pest .and. boozer) then
write (*,"(/,2X,A,/)") "--------------- No output. Please select Boozer or PEST coordinates. ------------"
stop
ENDIF

!for all modes
if (x0min < 0) then
   write (*,"(/,2X,A,/)") "--------------- Stop. Prescribe x0min. ------------"
   stop
endif

if (nz0 < 0) then
   write (*,"(/,2X,A,/)") "--------------- Stop. Prescribe nz0. ------------"
   stop
endif

IF (.not. global) nsurf=1

IF (global) THEN
   if ((x0max < 0) .and. (nsurf > 1)) then
     write (*,"(/,2X,A,/)") "-------------------     Stop. Prescribe x0max.     ----------------"
     stop
  endif

  if ((x0max <= x0min) .and. (nsurf > 1)) then
     write (*,"(/,2X,A,/)") "-------------------     Stop. x0max must be larger than x0min.     ----------------"
     stop
  endif
       
  if (nsurf < 0) then
     write (*,"(/,2X,A,/)") "--------------- Stop. Prescribe nsurf. ------------"
     stop
  endif

  if (nalpha0 == 0) then
     write (*,"(/,2X,A,/)") "--------------- Stop. Prescribe nalpha0. ------------"
     stop
  endif
ENDIF

IF ((.not. gs2) .and. (.not. terpsichore)) THEN
   if (.not. global) then
   WRITE(*,"(2X,A,/)")   "============================== FLUX-TUBE mode ==================================="
   WRITE(*,"(8X,A,2F8.4,//)") '                     s0, alpha0= ',x0min**2,alpha0
elseif (global .and. (nsurf==1)) then
   write (*,"(//,2X,A,/)")  "================================= SURFACE mode ===================================="
else
   WRITE(*,"(2X,A,/)")   "============================== GLOBAL mode ==================================="
endif

!------------------------------------------   Finishing SANITY checks   -----------------------------!
 
   if (global) then
      if (nsurf ==1) then
         open(56,file='../out/gist_surf_'//trim(vmec_file)//'_' &
                       &  //trim(out_tag),action='write')
    else
         open(56,file='../out/gist_glob_'//trim(vmec_file)//'_' &
                       &  //trim(out_tag),action='write')
      endif
   endif

  if (.not. global) then
         open(56,file='../out/gist_tube_'//trim(vmec_file)//'_' &
                       &  //trim(out_tag),action='write')     
      endif
    
       CALL gene_pest
       close(56)
       close(77) 
       close(78)
!Viewer
       if (global .and. gviu) then
       open(58,file='../out/gviu_'//trim(vmec_file)//'_' &
     &  //trim(out_tag),action='write')
       CALL viewer_lib
       close(58) 
    endif
stop
ENDIF





!==============================================================================================!
!------------------------------------------ TERPSICHORE ---------------------------------------!
!                                            GENE & GS2                                        !
!==============================================================================================!

!----------------------------------------------------------------------------------------------!
!------------------------------------------- GS2 ----------------------------------------------!
!                 ONLY flux-tube using Terpsichore (ie. no full-surface, no library)           !
!----------------------------------------------------------------------------------------------!

if (gs2) then 
write (*,"(/,2X,A,/)") "======================== Preparing geometry for GS2 code ========================"
   terpsichore=.true. ; vvbal=.true. ; tracer=.false. 

!Check
    if (pest) then
write (*,"(/,2X,A,/)") "--- PEST coordinates not available for GS2. Using Boozer coordinates instead. ---"
    endif

    if (global) then
write (*,"(/,2X,A,/)") "------- Full surface setup not available for GS2. Please set global=.f. --------"
    stop
   endif

   endif

!--------------------------------!
!This is not really necessary!
   if (stellopt) then
   terpsichore=.true. ; vvbal=.false. ; tracer=.true. ; initialize=.false.
   vmec_6_90=.true. ; global=.false. ; gviu=.false.
   pol_turns=1. ; nz0=128
   endif
!--------------------------------!
   
   s_tr=rho_tr**2

!----------------------------------------- TERPSICHORE begins -------------------------------------!
terp: IF (terpsichore) THEN

if (.not. global) then
write (*,"(/,2X,A,/)") "===================== BOOZER coordinates in flux tube mode ======================"
else
write (*,"(/,2X,A,/)") "===================== BOOZER coordinates in global mode ======================="
endif

   if (surf .eq. 0) then
write (*,"(/,2X,A,/)") "--- Please select a surface (the number of surfaces is found in the wout file) ---"
    stop
   endif


      IF (initialize) THEN 
write (*,"(/,2X,A,/)")  "========================== Initializing equilibrium ============================="
    
      CALL eqintrp_08_i
     
      do i=1,20
         tim(i) = 0.0
      end do

      call TPRALL_BAL
      call eqinvm  
      call second (tim (1) ) 

      call veqrec (tim)  
      call second (tim (4) )
      call mtaskb (tim)  
      call second (tim (9) )  

  100 continue  

  ELSE
write (*,"(2X,A,/)")    "======================== Using initialized equilibrium =========================="

        ENDIF 

!!========================================================================!!
!!================================ GLOBAL ==============================!!
!!========================================================================!!

GLOBAL_loop : if (global) then

call combal_alloc
allocate(alpha0_vec(0:nalpha0-1))

open(35,file='../out/gist_gene_global_'//trim(vmec_file)//'_' &
     &        //trim(table_tag)//'_'//trim(out_tag),action='write')

if (gviu) then
open(95,file='../out/gist_gviu_'//trim(vmec_file)//'_' &
     &        //trim(table_tag)//'_'//trim(out_tag),action='write')
endif


   dalpha = (alpha0_end - alpha0_start) / real(nalpha0)
   alpha0_start_=-alpha0_start

   DO nalpha=0,nalpha0-1
      alpha0_vec(nalpha)=alpha0_start_ - nalpha*dalpha
   ENDDO 

WRITE(*,"(/,4X,A,I4,2F10.4,I5,/)") '  surf, alpha0_start, alpha0_end, nalpha0 = ' &
,surf, alpha0_start, alpha0_end, nalpha0 

   DO nalpha=0,nalpha0-1
      alpha0=alpha0_vec(nalpha)
      alpha0_out=alpha0

   call datain
   call driver 
   call geneyout
   
   ENDDO

close(35)

else 


!!================================= VVBAL ==================================!!

  	IF (vvbal) THEN 

        call combal_alloc
        call datain
        call driver 
 	ENDIF	

!!================================= TRACER ===================================!! 


tracing:    	IF (tracer) THEN
!write (*,"(//,A,//)") ,"================== Running TRACER =================="       
  CALL Init_MPI
  CALL Comm_rank_MPI(COMM_WORLD_MPI,g_rank)
  CALL Comm_size_MPI(COMM_WORLD_MPI,g_num_proc)

ALLOCATE(cut(cuts+1))

!Read input file

    CALL read_input

!Open output files

  IF (out_mode==0) CALL open_out

!=========================================================================!
!                            LOOP OVER LINES                              !
!=========================================================================!


lines: DO ij=1,num_lines

geneloop: DO genecnt=1,geneend

     r=r_in(ij) ; z=z_in(ij) ; 
     phi=p0_in(ij)
     dist=p1_in(ij) ; resol=res(ij) 
     dist=dist*Pi  
     
      CALL field(r,phi,z,Br,Bp,Bz,dBrdr,dBrdp,dBrdz,dBpdr,dBpdp,&
               &dBpdz,dBzdr,dBzdp,dBzdz,device)

!Safety
     IF (ABS(Bp) < TINY(1.0)) STOP 'STOP(tracer): !Check integrity of field file!'   


#if 0
 IF (symm(ij) == 1) THEN        
        CALL symm_plane(r,phi,zmin,zmax,zsym,ierr)

        IF (ierr == 0)  THEN
           z=zsym
           CALL field(r,phi,zsym,br,bp,bz,dbrdr,dbrdp,dbrdz,dbpdr,dbpdp,&
                dbpdz,dbzdr,dbzdp,dbzdz,device)   
          WRITE(*,"(1X,A,E16.9,1X,A,E16.9)") &
                      &'Symmetry plane at z0 =',zsym,'with Br =', Br
        ENDIF
     ENDIF
#endif

!Define starting z

     z0=z

!Allocate output arrays

     !necessary due to changing resol
     IF (genecnt==2.AND.out_mode==0)  CALL dealloc_out
     IF (out_mode==0) CALL alloc_out 

!Resolutions for forward/backward tracing
     CALL set_resolution(ij)
     
!--------------------------------------------------
!Set dphi (initially)
     CALL set_dphi(code,ij,0)



!--------------------- WRITE OUT PARAMETERS ---------------!

WRITE(*,"(/,2X,A,/)")      "------------------------- Midpoint position of flux tube ------------------------"
 WRITE(*,"(15X,A,2F8.4,/)") '            rho/a, alpha = ',sqrt(s_in),alpha0

 ! WRITE(*,"(12X,A,3F8.3,/)") 'CYLINDRICAL (r0,phi0,zeta0) = ',r,phi,z
! WRITE(*,"(4X,A,1F8.3,/)") 'Fraction of trapped particles = ',ftrap

     hn1=0.1*dphi !guessed step for integrator

!=========================================================================!
!                          On the cut phi = pstart                        !
!=========================================================================!

  index=resol_b + 1

!Prescribe values for C^i_j
        Bp0=Bp
        c=0.
        c(1,1)=c11_rho(ij)  
        c(1,3)=-r*Br/Bp
        c(3,3)=1.      
        c(2,2)=-1./c(1,1)*Bp0/b0
        c(2,3)=+r*Bz/Bp/c(1,1)*Bp0/b0
!        c(2,2)=Bp/c(1,1)*Bp0/b0
!        c(2,3)=-r*Bz/c(1,1)*Bp0/b0
        
!Matrix D=INV(C)

     CALL matrix_d

!Contravariant Clebsch metrics 
     
     CALL gij_contravar

!Covariant Clebsch metmrics 

     CALL gij_covar    

!Calculate quantities

     CALL quantities

!Prepare output
 
     CALL prep_out(index,ij)


!=========================================================================!
!                      TRACING ALONG FIELD LINE                           !
!=========================================================================!


!Set accuracy for odeint

     err=1.e-11

!--------------------------  Tracing BACKWARDS  --------------------------!

! Skip backward tracing if required

IF(trace_type(ij) == 1) GOTO 101 

!Prepare cuts

     CALL prep_cuts(resol_b)    

     cuts_b : DO ll=1,cuts

!Prescribe initial conditions

           y(1) = r
           y(2) = z
           y(3)=c11_rho(ij)
           y(4) = 0.
           y(5)=-r*Br/Bp
           y(6) = 0.
           y(7) = -1./y(3)*bp0/b0
           y(8) = +r*Bz/Bp/y(3)*bp0/b0 
  !         y(7) = -Bp/y(3)*Bp0/b0
  !         y(8) = r*Bz/y(3)*bp0/b0
           
!Initial conditions for w

           w(1) = r
           w(2) = z

resolution_b : DO kk=cut(ll),cut(ll+1)

    index=resol_b-kk+1 !for output
    CALL set_dphi(code,ij,-kk)
    
    phi=phi-dphi
    
	if (midplane) then !create metrics 
    CALL odeint(y,phi+dphi,phi,err,hn1,hn1min,rhs)
        endif 
    CALL odeint(w,phi+dphi,phi,err,hn1,hn1min,rhs1) !!!!kludge 

!Update
        
    r = y(1)
    z = y(2)
    r = w(1)
    z = w(2)
    c(1,1) = y(3)
    c(1,2) = y(4)
    c(1,3) = y(5)
    c(2,1) = y(6)
    c(2,2) = y(7)
    c(2,3) = y(8)

!!$    !Implicit rule
!!$    CALL field(r,phi,z,br,bp,bz,dbrdr,dbrdp,dbrdz,dbpdr,dbpdp,&
!!$                dbpdz,dbzdr,dbzdp,dbzdz,device) 
!!$    c(1,3) = -(Br*c(1,1)+Bz*c(1,2))*r/Bp
!!$    c(2,3) = -(Br*c(2,1)+Bz*c(2,2))*r/Bp

!Matrix D=INV(C)

    CALL matrix_d

!Contravariant Clebsch metrics 
            
    CALL gij_contravar

!Covariant Clebsch metrics 

    CALL gij_covar     

!Calculate field in new position

    CALL field(r,phi,z,Br,Bp,Bz,dBrdr,dBrdp,dBrdz,dBpdr,dBpdp,&
              &dBpdz,dBzdr,dBzdp,dBzdz,device)

!Calculate quantities

    CALL quantities

!Prepare output

    CALL prep_out(index,ij)


 ENDDO resolution_b

ENDDO cuts_b

!--------------------------  Tracing FORWARDS  ---------------------------!

!Skip forward tracing if required

IF(trace_type(ij) == -1) GOTO 200 

!Restart
101 r=r_in(ij) ; phi=p0_in(ij)
IF (symm(ij) == 1 .AND. ierr == 0) THEN
   z=zsym
ELSE 
   z=z_in(ij)
ENDIF

CALL field(r,phi,z,Br,Bp,Bz,dBrdr,dBrdp,dBrdz,dBpdr,dBpdp,&
          &dBpdz,dBzdr,dBzdp,dBzdz,device) 

!Prepare cuts

CALL prep_cuts(resol_f)    


cuts_f : DO ll=1,cuts

!Prescribe initial conditions

           y(1) = r
           y(2) = z
           y(3)=c11_rho(ij)
           y(4) = 0.
           y(5)=-r*Br/Bp
           y(6) = 0.
           y(7) = -1./y(3)*bp0/b0
           y(8) = +r*Bz/Bp/y(3)*bp0/b0   
  !         y(7) = -Bp/y(3)*Bp0/b0
  !         y(8) = r*Bz/y(3)*Bp0/b0
          
!Initial conditions for w

           w(1) = r
           w(2) = z
           
resolution_f : DO kk=cut(ll),cut(ll+1)         

   index=resol_b+kk+1 !for output             

   CALL set_dphi(code,ij,kk)
   phi=phi+dphi


	if (midplane) then !create metrics 
    CALL odeint(y,phi-dphi,phi,err,hn1,hn1min,rhs) 
        endif
    CALL odeint(w,phi-dphi,phi,err,hn1,hn1min,rhs1) !!!kludge

!Update
   r = y(1)
   z = y(2)
   r = w(1)
   z = w(2)
   c(1,1) = y(3)
   c(1,2) = y(4)
   c(1,3) = y(5)
   c(2,1) = y(6)
   c(2,2) = y(7)
   c(2,3) = y(8)


!!$   !Implicit rule
!!$    CALL field(r,phi,z,br,bp,bz,dbrdr,dbrdp,dbrdz,dbpdr,dbpdp,&
!!$                dbpdz,dbzdr,dbzdp,dbzdz,device) 
!!$    c(1,3) = -(Br*c(1,1)+Bz*c(1,2))*r/Bp
!!$    c(2,3) = -(Br*c(2,1)+Bz*c(2,2))*r/Bp

!Matrix D=INV(C)

   CALL matrix_d

!Contravariant Clebsch metrics 
            
   CALL gij_contravar

!Covariant Clebsch metrics 

   CALL gij_covar     

!Calculate field in new position

   CALL field(r,phi,Z,Br,Bp,Bz,dBrdr,dBrdp,dBrdz,dBpdr,dBpdp,&
             &dBpdz,dBzdr,dBzdp,dBzdz,device)
 
!Calculate quantities

   CALL quantities

!Prepare output

   CALL prep_out(index,ij)

!For tokamaks: 
   IF (geneend/=1.AND.force_period) CALL Check_periodicity(index,phi_period)

ENDDO resolution_f

ENDDO cuts_f

!=========================================================================!
!                                OUTPUT                                   !
!=========================================================================!

!Calculate safety factor

200 IF (genecnt==1.AND..NOT.force_period) CALL safety_factor

!Calculate shat

IF (genecnt==1.AND..NOT.(code == 2).AND..NOT.force_period) CALL shear

!Produce output

ENDDO geneloop

IF(.NOT.(code == 2))CALL local_shear(p_out,sloc)

IF (out_mode==0) CALL output(ij)
!Deallocate output arrays

IF (out_mode==0) CALL dealloc_out

ENDDO lines

IF (out_mode==0) CALL close_out

CALL Finalize_MPI

	ENDIF  tracing
!=========================================================================!

IF (gs2) CALL gs2_out

endif GLOBAL_LOOP

ENDIF terp

        END PROGRAM gist
