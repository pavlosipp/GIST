
MODULE load_mod
  USE UtilCDF
  USE specifications_mod
  USE vmec_mod, ONLY : vmec_data


  IMPLICIT NONE


!=========================================================================!
CONTAINS               !                MODULE SUBPROGRAMS                !
!=========================================================================!

!===================== READ INPUT FILE tracer.inp ========================!

  SUBROUTINE read_input

    INTEGER :: i,j, scan
    INTEGER, DIMENSION(:), ALLOCATABLE:: pts,sym,tt
    REAL(DP) :: r0,c11,shat_loc
    REAL(DP), DIMENSION(:), ALLOCATABLE:: rstart,rstep,rend,z,p0,p1,lines
   
!    OPEN(65, FILE='../inp/tracer.inp', ACTION='READ', ERR=100)
!    READ(65,*) 
!    READ(65,*) num_lines,scan

     scan=0 !hardwired

    !If 'number of scans'=0, rstep and rend are read but ignored
    IF (scan==0) THEN
       ALLOCATE(rho_in(num_lines),r_in(num_lines),rstep(num_lines),t0_in(num_lines),&
       &rend(num_lines),z_in(num_lines),res(num_lines),symm(num_lines))
       ALLOCATE(trace_type(num_lines))
       ALLOCATE(p0_in(num_lines),p1_in(num_lines))
   !    READ(65,*)
       DO i=1,num_lines
   !       READ(65,*) r_in(i),rstep(i),rend(i),z_in(i),p0_in(i),&
   !            &p1_in(i),res(i),symm(i),trace_type(i)
       ENDDO
    ELSE
       !If 'number of scans'/=0, do a tracing for every field line between rstart and
       !rend using steps of rstep.
       ALLOCATE(rstart(scan),rstep(scan),rend(scan),z(scan),p0(scan),p1(scan),&
            &lines(0:scan),pts(scan),sym(scan),tt(scan))
       lines(0)=0
       READ(65,*)
       DO i=1,scan
          READ(65,*) rstart(i),rstep(i),rend(i),z(i),p0(i),p1(i),pts(i),sym(i),tt(i)
       ENDDO
       !line numbering
       DO i=1,scan
          lines(i)=lines(i-1)+int((rend(i)-rstart(i))/rstep(i))+1
       ENDDO
       num_lines=lines(scan)
       ALLOCATE(rho_in(num_lines),r_in(num_lines),z_in(num_lines))
       ALLOCATE(res(num_lines),symm(num_lines))
       ALLOCATE(trace_type(num_lines),p0_in(num_lines),p1_in(num_lines))
       DO j=1,scan
          DO i=lines(j-1)+1,lines(j)
             r_in(i)=rstart(j)+(i-lines(j-1)-1)*rstep(j)
             z_in=z(j)
             p0_in=p0(j)
             p1_in=p1(j)
             res=pts(j)
             symm=sym(j)
             trace_type=tt(j)
          ENDDO
       ENDDO
    ENDIF
 !   CLOSE(65)


!================================= VMEC ======================================!

       ALLOCATE(safetyf(num_lines),shatf(num_lines),c11_rho(num_lines))

!from prepar.f 
       if (vvbal) then
       r_in(1)=s_tracer
       res(1)=nz0
       z_in(1)=theta0
       p1_in(1)=2.*pol_turns*q0_vvbal
       p0_in(1)=alpha0+q0_vvbal*theta0
 
     else

       r_in(1)=s_tr 
       z_in(1)=t_tr
       p0_in(1)=z_tr
       res(1)=nz0
       endif

!decide whether we run also tracer for metrics
       if ((p0_in(1) == 0) .and. (z_in(1) == 0)) midplane=.true.
if (stellopt) midplane=.true.


       s_in=r_in(1)

!set by default
       symm(1)=0
       trace_type(1)=0

       DO j=1,num_lines  
         CALL vmec_data(r_in(j),z_in(j),p0_in(j),c11,safety,alpha,shat_loc,ftrap,grads,refs)
         c11_rho(j)=c11
         safetyf(j)=safety
         shatf(j)=shat_loc
         q0_tracer=safety
         shat_tracer=shat_loc
       if(tracer .and. (.not. vvbal)) p1_in(j)=2*pol_turns*safety 

       ENDDO
  
!=============================================================================!

    RETURN

  END SUBROUTINE read_input

!-------------------------------------------------------------------------!

!===================== READ MESH FILE for findif code ====================!

  SUBROUTINE read_mesh

    INTEGER :: ii,i,idummy,c,n,nl,pnl,ci,max_n,pplot,nplots
    REAL(DP) :: dphi,pl_dist,x,y,ddym(2)
    LOGICAL*1 :: ptype !If false - end or beginning of open line
    LOGICAL :: bfirst !Start point for open line
    !-------------------------------------------------
    INTEGER*1 :: i1dym(1)
    INTEGER :: dimids(2),idym(1),npoints,np,startl,numl,fid
    INTEGER, DIMENSION(:), ALLOCATABLE :: fl_size,closeind
#ifdef NETCDF
    in_type_if: If(in_type==2) THEN
    CALL open_cdf('../input/'//TRIM(config)//'/mesh.inp',READ_ACC_CDF,fid)
    IF(exist_cdf(fid,ATT_CDF,'part')) THEN
      WRITE(*,*) 'Partial points list!!!'
      STOP
    ENDIF
    CALL get_dim_cdf(fid,'npoints',dimids(1),npoints)
    CALL get_dim_cdf(fid,'nlines',dimids(2),nl)
    ALLOCATE(fl_size(nl))
    ALLOCATE(closeind(nl))
    CALL get_var_cdf(fid,'line_size',fl_size)
    CALL get_var_cdf(fid,'closeind',closeind)
    !--------Get start position and number of points----------------------
    startl = (nl/g_num_proc)*g_rank+1
    IF (g_rank==g_num_proc-1) THEN
      numl= nl - startl + 1
    ELSE
      numl= nl/g_num_proc
    ENDIF
    num_startp=1
    num_lines=0
    DO i=1,startl+numl-1
      IF(i<startl) num_startp=num_startp+fl_size(i)
      IF(i>=startl) num_lines=num_lines+fl_size(i)
    ENDDO
    !---------------------------------------------------------------------

    ! Allocating all needed resources
    ALLOCATE(r_in(num_lines),z_in(num_lines))
    ALLOCATE(res(num_lines),symm(num_lines))
    ALLOCATE(trace_type(num_lines))
    ALLOCATE(p0_in(num_lines),p1_in(num_lines))
    !---------------------------------------------------------------------
    IF(.NOT.get_att_cdf(fid,get_varid_cdf(fid,'plot'),'valid_range',ddym)) THEN
      WRITE(*,*) 'Error in cdf file!!! (load-MOD.F90)'
      STOP
    ENDIF
    pl_dist= 2./ddym(2)
    CALL seek_cdf(get_varid_cdf(fid,'plot'),num_startp,SEEK_SET_CDF)
    ! CALL seek_cdf(get_varid_cdf(fid,'phi'),num_startp,SEEK_SET_CDF)
    CALL seek_cdf(get_varid_cdf(fid,'R'),num_startp,SEEK_SET_CDF)
    CALL seek_cdf(get_varid_cdf(fid,'z'),num_startp,SEEK_SET_CDF)
    CALL seek_cdf(get_varid_cdf(fid,'line'),num_startp,SEEK_SET_CDF)
    CALL seek_cdf(get_varid_cdf(fid,'type'),num_startp,SEEK_SET_CDF)
    CALL seek_cdf(get_varid_cdf(fid,'closeind'),num_startp,SEEK_SET_CDF)
    CALL seek_cdf(get_varid_cdf(fid,'dxi'),num_startp,SEEK_SET_CDF)
    pnl=0
    bfirst=.false.
    trace_type=0
    DO ii = 1, num_lines
      CALL get_var_cdf(fid,'x',ddym,1)
      x=ddym(1)
      CALL get_var_cdf(fid,'y',ddym,1)
      y=ddym(1)
      CALL get_var_cdf(fid,'R',ddym,1)
      r_in(ii)=ddym(1)
      CALL get_var_cdf(fid,'z',ddym,1)
      z_in(ii)=ddym(1)
      CALL get_var_cdf(fid,'line',idym,1)
      nl=idym(1)
      CALL get_var_cdf(fid,'type',i1dym,1)
      ptype=i1dym(1)
      CALL get_var_cdf(fid,'plot',idym,1)
      CALL get_var_cdf(fid,'dxi',ddym,1)
      ci=closeind(nl)
      IF(bfirst) THEN ! ii is second point
        bfirst = .false.
        p0_in(ii)=pl_dist*REAL(idym(1)-1)
        p1_in(ii)=pl_dist*2.
        p1_in(ii-1) = ddym(1)/Pi !Fill the distance for the start point
        res(ii)=3
        trace_type(ii) = -2
      ENDIF
      IF((pnl/=nl).AND.(ci==1).AND.(.NOT.ptype)) THEN
        p0_in(ii) = Phi_xy(x,y)/Pi
        trace_type(ii)=1 !Beginning of the line (trace forwards)
        bfirst = .true.
        res(ii)=2
      ELSEIF((pnl==nl).AND.(ci==1).AND.(.NOT.ptype)) THEN
        p0_in(ii) = Phi_xy(x,y)/Pi
        p1_in(ii) = ddym(1)/Pi
        res(ii)=2
        trace_type(ii)=-1 !End of the line (trace backwards)
        trace_type(ii-1)=2 !We mark previous point
      ELSEIF(trace_type(ii)/=-2) THEN
        p0_in(ii)=pl_dist*REAL(idym(1)-1)
        p1_in(ii)=pl_dist*2.
        res(ii)=3
        trace_type(ii)=0
      ENDIF
      symm(ii)=0
      pnl=nl
    ENDDO
    CALL close_cdf(fid)
    DEALLOCATE(fl_size)
    ELSEIF(in_type==1) THEN
#endif
    OPEN(65, FILE='../input/'//TRIM(config)//'/mesh.inp', &
             ACTION='READ', ERR=200)
    READ(65,*)
    READ(65,*) num_lines
    READ(65,*)
    READ(65,*) nplots
    READ(65,*)
    READ(65,*) 
    READ(65,*)
    ALLOCATE(r_in(num_lines),z_in(num_lines))
    ALLOCATE(res(num_lines),symm(num_lines))
    ALLOCATE(trace_type(num_lines))
    ALLOCATE(p0_in(num_lines),p1_in(num_lines))
    c = 0
    max_n = 0
    pnl=0
    pl_dist=2/REAL(nplots)
    bfirst=.false.
    trace_type=0
    DO ii=1,num_lines
      READ(65, *) i, pplot, x, y, & 
         & r_in(i+1), z_in(i+1), dphi,ptype, idummy,nl,n,ci
      IF(bfirst) THEN ! ii is second point
        bfirst = .false.
        p0_in(ii)=pl_dist*REAL(pplot-1)
        p1_in(ii)=pl_dist*2.
        p1_in(ii-1) = dphi/Pi !Fill the distance for the start point
        res(ii)=3
        trace_type(ii) = -2
      ENDIF
      IF((pnl/=nl).AND.(ci==1).AND.(.NOT.ptype)) THEN
        p0_in(ii) = Phi_xy(x,y)/Pi
        trace_type(ii)=1 !Beginning of the line (trace forwards)
        bfirst = .true.
        res(ii)=2
      ELSEIF((pnl==nl).AND.(ci==1).AND.(.NOT.ptype)) THEN
        p0_in(ii) = Phi_xy(x,y)/Pi
        p1_in(ii) = dphi/Pi
        res(ii)=2
        trace_type(ii)=-1 !End of the line (trace backwards)
        trace_type(ii-1)=2
      ELSEIF(trace_type(ii)/=-2) THEN
        p0_in(ii)=pl_dist*REAL(pplot-1)
        p1_in(ii)=pl_dist*2.
        res(ii)=3
        trace_type(ii)=0
      ENDIF
      symm(ii)=0
      pnl=nl
    ENDDO
    CLOSE (65)

    WRITE(*,*) 'Findif mesh was read.'

#ifdef NETCDF
    ENDIF in_type_if
#endif
    resol = res(1)
    RETURN
200 STOP '!!Mesh file not found!!'

  END SUBROUTINE read_mesh

!===================== READ NEIGHBORS FILE for findif code ===============!

  SUBROUTINE read_neighbors

		! Counters
    INTEGER :: ii,i,fid,idym(1)
    CHARACTER(56) :: str_
#ifdef NETCDF
    in_type_if: If(in_type==2) THEN
      CALL open_cdf('../input/'//TRIM(config)//'/torn.inp',READ_ACC_CDF,fid)
      IF(exist_cdf(fid,ATT_CDF,'part')) THEN
        WRITE(*,*) 'Partial points list!!!'
        STOP
      ENDIF
      ALLOCATE(nb_tor(num_lines,1:2))
      CALL seek_cdf(get_varid_cdf(fid,'torn_behind'),num_startp,SEEK_SET_CDF)
      CALL seek_cdf(get_varid_cdf(fid,'torn_front'),num_startp,SEEK_SET_CDF)
      DO ii=1,num_lines
        CALL get_var_cdf(fid,'torn_behind',idym,1)
        nb_tor(ii,1)=idym(1)
        CALL get_var_cdf(fid,'torn_front',idym,1)
        nb_tor(ii,2)=idym(1)
      ENDDO
      CALL close_cdf(fid)
    ELSEIF(in_type==1) THEN
#endif
      OPEN(1, FILE='../input/'//TRIM(config)//'/torn.inp',&
              ACTION='READ', ERR=300)
      
      ! Reading of toroidal neighbors
      ALLOCATE(nb_tor(num_lines,1:2))
      READ(1,'(A25)') str_
      WRITE(*,*) str_
      DO ii = 1, num_lines
        READ(1,*) i, nb_tor(i+1,1:2) 
      END DO
      !=================
      ! The rest we don't need here :)
      !=================
      ! we shift the arrays by one so that the indices held in the arrays 
      ! range from 1 to k, not from 0 to k-1
      nb_tor   = nb_tor + 1
      WRITE(*,*) 'Findif neighbors file was read. '
#ifdef NETCDF
    ENDIF in_type_if
#endif
    RETURN
300 STOP '!!Neighbors file not found!!'

  END SUBROUTINE read_neighbors

!-------------------------------------------------------------------------!
  REAL(DP) FUNCTION Phi_xy(x_,y_)
    REAL(DP), INTENT(IN) :: x_,y_
    IF (ABS(x_).GT.ABS(y_)) THEN
      Phi_xy = ATAN(y_/x_)
      IF(x_ .LT. 0.) Phi_xy = Phi_xy + Pi
    ELSE
        Phi_xy = ATAN(x_/y_)
      IF (y_ .LT. 0.) Phi_xy = Phi_xy + Pi
        Phi_xy = Pi/2. - Phi_xy
    END IF
    IF (Phi_xy.LT.0) Phi_xy = Phi_xy+2.*Pi
  END FUNCTION Phi_xy


!================================ BORIS ===================================!

SUBROUTINE read_boris

  INTEGER :: i

  OPEN(65, FILE='../input/'//TRIM(config)//'/tracer.inp_boris', ACTION='READ', ERR=100)
  READ(65,*) 
  READ(65, '(I3)') num_lines
  ALLOCATE(r_in(num_lines),z_in(num_lines),res(num_lines),symm(num_lines))
  ALLOCATE(trace_type(num_lines),iota_boris(num_lines))
  ALLOCATE(p0_in(num_lines),p1_in(num_lines))
  READ(65,*)
  DO i=1,num_lines
     READ(65,*) r_in(i),iota_boris(i),z_in(i),p0_in(i),p1_in(i),res(i)
  ENDDO
  CLOSE(65)

!Hardwire parameters
  trace_type=0
  symm=1

  RETURN
100 STOP 'STOP(read_input): !Input file tracer.inp_boris not found!'

END SUBROUTINE read_boris


!----------------------------------------------------------------------------!

END MODULE load_mod
