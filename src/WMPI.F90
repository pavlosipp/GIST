!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	purpose: 		Module for wrapping functions over MPI
!
!	author:			Comrad
!	date of creation: 	25.01.2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE WMPI
	USE type_mod
  USE specifications_mod
  ! USE tag_params_MPI
#ifdef USE_MPI
  USE MPI
#endif
	IMPLICIT NONE

  PUBLIC :: Finalize_MPI,Init_MPI,Comm_rank_MPI,Comm_size_MPI,COMM_WORLD_MPI,&
           &INT_MPI,DOUBLE_MPI,LOGICAL_MPI,CHARACTER_MPI,REAL_MPI,&
           &STATUS_SIZE_MPI,SOURCE_MPI,SUCCESS_MPI,&
           &Iprobe_MPI,Notify_all,Recv_MPI,Allgather_alive_MPI,&
           &Allgatherv_alive_MPI,BarrierAlive_MPI,Barrier_MPI,&
           &BARRIER_TAG,WORK_FINISHED_TAG,&
           &Register_type_MPI,Abort_MPI,point_3d_mpi,Isend_MPI

INTERFACE Allgatherv_alive_MPI
  MODULE PROCEDURE Allgatherv_alive_MPI_d
  MODULE PROCEDURE Allgatherv_alive_MPI_i
END INTERFACE

INTERFACE Allgather_alive_MPI
  MODULE PROCEDURE Allgather_alive_MPI_d
  MODULE PROCEDURE Allgather_alive_MPI_i
  MODULE PROCEDURE Allgather_alive_MPI_i1
END INTERFACE

INTERFACE Recv_MPI
  MODULE PROCEDURE Recv_MPI_d
  MODULE PROCEDURE Recv_MPI_i
  MODULE PROCEDURE Recv_MPI_i1
END INTERFACE

INTERFACE Isend_MPI
  MODULE PROCEDURE Isend_MPI_d
  MODULE PROCEDURE Isend_MPI_i
  MODULE PROCEDURE Isend_MPI_i1
END INTERFACE

#ifdef USE_MPI
  INTEGER(ni), PARAMETER :: COMM_WORLD_MPI = MPI_COMM_WORLD
  INTEGER(ni), PARAMETER :: STATUS_SIZE_MPI = MPI_STATUS_SIZE
  INTEGER(ni), PARAMETER :: SOURCE_MPI = MPI_SOURCE
  INTEGER(ni), PARAMETER :: SUCCESS_MPI = MPI_SUCCESS

  INTEGER(ni), PARAMETER :: INT_MPI = MPI_INTEGER
  INTEGER(ni), PARAMETER :: REAL_MPI = MPI_REAL
  INTEGER(ni), PARAMETER :: DOUBLE_MPI = MPI_DOUBLE_PRECISION
  INTEGER(ni), PARAMETER :: CHARACTER_MPI = MPI_CHARACTER
  INTEGER(ni), PARAMETER :: LOGICAL_MPI = MPI_LOGICAL

#else
  INTEGER(ni), PARAMETER :: COMM_WORLD_MPI =  1
  INTEGER(ni), PARAMETER :: STATUS_SIZE_MPI =  1
  INTEGER(ni), PARAMETER :: SOURCE_MPI =  1
  INTEGER(ni), PARAMETER :: SUCCESS_MPI =  1

  INTEGER(ni), PARAMETER :: INT_MPI =  1
  INTEGER(ni), PARAMETER :: REAL_MPI =  1
  INTEGER(ni), PARAMETER :: DOUBLE_MPI =  1
  INTEGER(ni), PARAMETER :: CHARACTER_MPI =  1
  INTEGER(ni), PARAMETER :: LOGICAL_MPI =  1
#endif

  INTEGER(ni), PARAMETER :: BARRIER_TAG = 10
  INTEGER(ni), PARAMETER :: WORK_FINISHED_TAG = 11
  INTEGER(ni) :: point_3d_mpi

  PRIVATE

  INTEGER(ni) :: g_ierr
  INTEGER(ni) :: g_num_alive

  CONTAINS

  SUBROUTINE Init_MPI()
#ifdef USE_MPI
    CALL MPI_Init(g_ierr)
    CALL Check_Err_MPI()
#else
    ! OPEN(5)
    ! READ(5,*) g_rank,g_size
    ! CLOSE(5)
    ! IF(g_rank<0.OR.g_rank>=g_size.OR.g_size<1) THEN
    !   WRITE(*,*) 'Error! Wrong rank and size!!!'
    !   STOP
    ! ENDIF
    g_rank=0
    g_num_proc=1
#endif
  END SUBROUTINE Init_MPI

  SUBROUTINE Finalize_MPI()
#ifdef USE_MPI
    CALL MPI_Finalize(g_ierr)
    CALL Check_Err_MPI()
#endif
  END SUBROUTINE Finalize_MPI

  SUBROUTINE Comm_rank_MPI(mpi_comm_, rank_)
    INTEGER(ni), INTENT(IN) :: mpi_comm_
    INTEGER(ni), INTENT(OUT) :: rank_
#ifdef USE_MPI
    CALL MPI_Comm_rank(mpi_comm_,rank_,g_ierr)
    CALL Check_Err_MPI()
    g_rank=rank_
#else
    rank_ = g_rank
#endif
  END SUBROUTINE Comm_rank_MPI

  SUBROUTINE Comm_size_MPI(mpi_comm_, size_)
    INTEGER(ni), INTENT(IN) :: mpi_comm_
    INTEGER(ni), INTENT(OUT) :: size_
#ifdef USE_MPI
    CALL MPI_Comm_size(mpi_comm_,size_,g_ierr)
    CALL Check_Err_MPI()
    g_num_proc=size_
    g_num_alive=size_
#else
    size_=g_num_proc
#endif
  END SUBROUTINE Comm_size_MPI

  SUBROUTINE Barrier_MPI
#ifdef USE_MPI
    CALL MPI_Barrier(COMM_WORLD_MPI,g_ierr)
    CALL Check_Err_MPI
#endif
  END SUBROUTINE

  SUBROUTINE Iprobe_MPI(src_,tag_,comm_,b_,st_)
    !IN
    INTEGER(ni) :: src_,tag_,comm_
    !OUT
    INTEGER(ni) :: st_(STATUS_SIZE_MPI)
    LOGICAL :: b_

#ifdef USE_MPI
    CALL MPI_iprobe(src_,tag_,comm_,b_,&
                   &st_,g_ierr)
    CALL Check_Err_MPI
#endif
  END SUBROUTINE

  SUBROUTINE Notify_all(sval_,tagval_,size_,is_alive_procs_)
    !IN
    INTEGER(ni), INTENT(IN) :: sval_,tagval_,size_
    LOGICAL, INTENT(IN), OPTIONAL :: is_alive_procs_(size_)
    !Internal
    INTEGER(ni) :: i,rq
#ifdef USE_MPI
    IF(PRESENT(is_alive_procs_))THEN
      DO i=1,size_
        IF(is_alive_procs_(i)) THEN
          CALL Mpi_isend(sval_,1,INT_MPI,i-1,tagval_,&
                        &COMM_WORLD_MPI,rq,g_ierr)
          CALL Check_Err_MPI
        ENDIF
      ENDDO
    ELSE
      DO i=1,size_
        CALL Mpi_isend(sval_,1,INT_MPI,i-1,tagval_,&
                      &COMM_WORLD_MPI,rq,g_ierr)
        CALL Check_Err_MPI
      ENDDO
    ENDIF
#endif
  END SUBROUTINE Notify_all

  SUBROUTINE BarrierAlive_MPI(is_alive_procs_,size_)
    !IN
    INTEGER(ni), INTENT(IN) :: size_
    !OUT
    LOGICAL :: is_alive_procs_(size_)
    !Internal
    INTEGER :: i,itmp,n,rq,st(STATUS_SIZE_MPI),st1(STATUS_SIZE_MPI)
    LOGICAL :: b,b1
    !Initial check if somebody have finished the work
    b1=.TRUE.
#ifdef USE_MPI
    DO WHILE (b1)
      CALL Iprobe_MPI(MPI_ANY_SOURCE,WORK_FINISHED_TAG,COMM_WORLD_MPI,b1,st1)
      IF(b1) THEN
        CALL Recv_MPI(itmp,1,INT_MPI,st1(SOURCE_MPI),WORK_FINISHED_TAG,&
                     &COMM_WORLD_MPI)
        is_alive_procs_(st1(SOURCE_MPI)+1)=.FALSE.
        g_num_alive=g_num_alive-1
      ENDIF
    ENDDO
    DO i=1,size_
      IF(is_alive_procs_(i)) THEN
        CALL Mpi_isend(0,1,INT_MPI,i-1,BARRIER_TAG,&
                      &COMM_WORLD_MPI,rq,g_ierr)
        CALL Check_Err_MPI
      ENDIF
    ENDDO
    !Start waiting until either BARRIER_TAG or WORK_FINISHED_TAG will come
    ! from all alive 
    n=0
    DO WHILE (n<g_num_alive)
      !Check if somebody have finished the work
      b1=.TRUE.
      DO WHILE (b1)
        CALL Iprobe_MPI(MPI_ANY_SOURCE,WORK_FINISHED_TAG,COMM_WORLD_MPI,b1,st1)
        IF(b1) THEN
          CALL Recv_MPI(itmp,1,INT_MPI,st1(SOURCE_MPI),WORK_FINISHED_TAG,&
                       &COMM_WORLD_MPI)
          is_alive_procs_(st1(SOURCE_MPI)+1)=.FALSE.
          g_num_alive=g_num_alive-1
        ENDIF
      ENDDO
      !Check barrier
      CALL Iprobe_MPI(MPI_ANY_SOURCE,BARRIER_TAG,COMM_WORLD_MPI,b,st)
      IF(b) THEN
        CALL Recv_MPI(itmp,1,INT_MPI,st(SOURCE_MPI),BARRIER_TAG,&
                     &COMM_WORLD_MPI)
        n=n+1
      ENDIF
    ENDDO
#endif
  END SUBROUTINE BarrierAlive_MPI

  SUBROUTINE Isend_MPI_i1(ival_,inum_,type_,rank_,tag_,comm_)
    !IN
    INTEGER(ni), INTENT(IN) :: ival_,inum_,type_,rank_,tag_,comm_
    !Internal
    INTEGER(ni) :: rq
#ifdef USE_MPI
    CALL Mpi_isend(ival_,1,type_,rank_,tag_,&
                  &comm_,rq,g_ierr)
    CALL Check_Err_MPI
#endif
  END SUBROUTINE Isend_MPI_i1

  SUBROUTINE Isend_MPI_i(ival_,inum_,type_,rank_,tag_,comm_)
    !IN
    INTEGER(ni), INTENT(IN) :: inum_,type_,rank_,tag_,comm_
    INTEGER(ni), INTENT(IN) :: ival_(inum_)
    !Internal
    INTEGER(ni) :: rq
#ifdef USE_MPI
    CALL Mpi_isend(ival_,inum_,type_,rank_,tag_,&
                  &comm_,rq,g_ierr)
    CALL Check_Err_MPI
#endif
  END SUBROUTINE Isend_MPI_i

  SUBROUTINE Isend_MPI_d(dval_,dnum_,type_,rank_,tag_,comm_)
    !IN
    INTEGER(ni), INTENT(IN) :: dnum_,type_,rank_,tag_,comm_
    REAL(DP), INTENT(IN) :: dval_(dnum_)
    !Internal
    INTEGER(ni) :: rq
#ifdef USE_MPI
    CALL Mpi_isend(dval_,dnum_,type_,rank_,tag_,&
                  &comm_,rq,g_ierr)
    CALL Check_Err_MPI
#endif
  END SUBROUTINE Isend_MPI_d

  SUBROUTINE Recv_MPI_i1(ival_,inum_,type_,rank_,tag_,comm_)
    !IN
    INTEGER(ni), INTENT(IN) :: inum_,type_,rank_,tag_,comm_
    !OUT
    INTEGER(ni) :: ival_
    !Internal
    INTEGER(ni) :: st(STATUS_SIZE_MPI)
#ifdef USE_MPI
    CALL Mpi_recv(ival_,1,type_,rank_,tag_,&
                  &comm_,st,g_ierr)
    CALL Check_Err_MPI
#endif
  END SUBROUTINE Recv_MPI_i1

  SUBROUTINE Recv_MPI_d(dval_,dnum_,type_,rank_,tag_,comm_)
    !IN
    INTEGER(ni), INTENT(IN) :: dnum_,type_,rank_,tag_,comm_
    !OUT
    REAL(DP) :: dval_(dnum_)
    !Internal
    INTEGER(ni) :: st(STATUS_SIZE_MPI)
#ifdef USE_MPI
    CALL Mpi_recv(dval_,dnum_,type_,rank_,tag_,&
                  &comm_,st,g_ierr)
    CALL Check_Err_MPI
#endif
  END SUBROUTINE Recv_MPI_d

  SUBROUTINE Recv_MPI_i(ival_,inum_,type_,rank_,tag_,comm_)
    !IN
    INTEGER(ni), INTENT(IN) :: inum_,type_,rank_,tag_,comm_
    !OUT
    INTEGER(ni) :: ival_(inum_)
    !Internal
    INTEGER(ni) :: st(STATUS_SIZE_MPI)
#ifdef USE_MPI
    CALL Mpi_recv(ival_,inum_,type_,rank_,tag_,&
                  &comm_,st,g_ierr)
    CALL Check_Err_MPI
#endif
  END SUBROUTINE Recv_MPI_i

  SUBROUTINE  Allgatherv_alive_MPI_d(dvals_,dnums_,&
                                    &dvalr_,dnumr_,dnumrt_,displ_,&
                                    &type_,size_,tag_,&
                                    &comm_,is_alive_procs_)
    !IN
    INTEGER(ni), INTENT(IN) :: dnums_,dnumrt_,type_,tag_,size_,comm_
    INTEGER(ni), INTENT(IN) :: dnumr_(size_)
    INTEGER(ni), INTENT(IN) :: displ_(size_+1)
    REAL(DP), INTENT(IN) :: dvals_(dnums_)
    LOGICAL, INTENT(IN), OPTIONAL :: is_alive_procs_(size_)
    !OUT
    REAL(DP) :: dvalr_(dnumrt_)
    !Internal
    INTEGER(ni) :: i
#ifdef USE_MPI
    DO i=1,size_
      IF(is_alive_procs_(i)) THEN
        CALL Isend_MPI(dvals_,dnums_,type_,i-1,tag_,&
                      &comm_)
      ENDIF
    ENDDO
    DO i=1,size_
      IF(is_alive_procs_(i)) THEN
        CALL Recv_MPI(dvalr_(displ_(i):displ_(i+1)),dnumr_(i),type_,i-1,tag_,&
                      &comm_)
      ELSE
        dvalr_(i)=0
      ENDIF
    ENDDO
#endif
  END SUBROUTINE Allgatherv_alive_MPI_d

  SUBROUTINE Allgatherv_alive_MPI_i(ivals_,inums_,&
                                   &ivalr_,inumr_,inumrt_,displ_,&
                                   &type_,size_,tag_,&
                                   &comm_,is_alive_procs_)
    !IN
    INTEGER(ni), INTENT(IN) :: inums_,inumrt_,type_,tag_,size_,comm_
    INTEGER(ni), INTENT(IN) :: inumr_(size_)
    INTEGER(ni), INTENT(IN) :: displ_(size_+1)
    INTEGER(ni), INTENT(IN) :: ivals_(inums_)
    LOGICAL, INTENT(IN), OPTIONAL :: is_alive_procs_(size_)
    INTEGER(ni) :: ivalr_(inumrt_)
    !Internal
    INTEGER(ni) :: i
#ifdef USE_MPI
    DO i=1,size_
      IF(is_alive_procs_(i)) THEN
        CALL Isend_MPI(ivals_,inums_,type_,i-1,tag_,&
                      &comm_)
      ENDIF
    ENDDO
    DO i=1,size_
      IF(is_alive_procs_(i)) THEN
        CALL Recv_MPI(ivalr_(displ_(i):displ_(i+1)),inumr_(i),type_,i-1,tag_,&
                      &comm_)
      ELSE
        ivalr_(i)=0
      ENDIF
    ENDDO
#endif
  END SUBROUTINE Allgatherv_alive_MPI_i
!====================================================
  SUBROUTINE Allgather_alive_MPI_d(dvals_,dnums_,dvalr_,dnumr_,&
                                   &type_,size_,tag_,&
                                   &comm_,is_alive_procs_)
    !IN
    INTEGER(ni), INTENT(IN) :: dnums_,dnumr_,type_,tag_,size_,comm_
    REAL(DP), INTENT(IN) :: dvals_(dnums_)
    LOGICAL, INTENT(IN), OPTIONAL :: is_alive_procs_(size_)
    !OUT
    REAL(DP):: dvalr_(dnumr_)
    !Internal
    INTEGER(ni) :: i
#ifdef USE_MPI
    DO i=1,size_
      IF(is_alive_procs_(i)) THEN
        CALL Isend_MPI(dvals_,dnums_,type_,i-1,tag_,&
                      &comm_)
      ENDIF
    ENDDO
    DO i=1,size_
      IF(is_alive_procs_(i)) THEN
        CALL Recv_MPI(dvalr_((i-1)*dnums_+1:i*dnums_),dnums_,type_,i-1,tag_,&
                      &comm_)
      ELSE
        dvalr_(i)=0
      ENDIF
    ENDDO
#endif
  END SUBROUTINE Allgather_alive_MPI_d

  SUBROUTINE Allgather_alive_MPI_i(ivals_,inums_,ivalr_,inumr_,&
                                   &type_,size_,tag_,&
                                   &comm_,is_alive_procs_)
    !IN
    INTEGER(ni), INTENT(IN) :: inums_,inumr_,type_,tag_,size_,comm_
    INTEGER(ni), INTENT(IN) :: ivals_(inums_)
    LOGICAL, INTENT(IN), OPTIONAL :: is_alive_procs_(size_)
    !OUT
    INTEGER(ni), INTENT(OUT) :: ivalr_(inumr_)
    !Internal
    INTEGER(ni) :: i
    DO i=1,size_
      IF(is_alive_procs_(i)) THEN
        CALL Isend_MPI(ivals_,inums_,type_,i-1,tag_,&
                      &comm_)
      ENDIF
    ENDDO
    DO i=1,size_
      IF(is_alive_procs_(i)) THEN
        CALL Recv_MPI(ivalr_((i-1)*inums_+1:i*inums_),inums_,type_,i-1,tag_,&
                      &comm_)
      ELSE
        ivalr_(i)=0
      ENDIF
    ENDDO
  END SUBROUTINE Allgather_alive_MPI_i

  SUBROUTINE Allgather_alive_MPI_i1(ivals_,inums_,ivalr_,inumr_,&
                                   &type_,size_,tag_,&
                                   &comm_,is_alive_procs_)
    !IN
    INTEGER(ni), INTENT(IN) :: ivals_,inums_,inumr_,type_,tag_,size_,comm_
    LOGICAL, INTENT(IN), OPTIONAL :: is_alive_procs_(size_)
    !OUT
    INTEGER(ni), INTENT(OUT) :: ivalr_(inumr_)
    !Internal
    INTEGER(ni) :: i
    DO i=1,size_
      IF(is_alive_procs_(i)) THEN
        CALL Isend_MPI(ivals_,1,type_,i-1,tag_,&
                      &comm_)
      ENDIF
    ENDDO
    DO i=1,size_
      IF(is_alive_procs_(i)) THEN
        CALL Recv_MPI(ivalr_(i),1,type_,i-1,tag_,&
                      &comm_)
      ELSE
        ivalr_(i)=0
      ENDIF
    ENDDO
  END SUBROUTINE Allgather_alive_MPI_i1

  SUBROUTINE Register_type_MPI
    !Internal
    TYPE(POINT_3D) :: p
    INTEGER(ni) :: ierr,blocklen(3), typep(3)
#ifdef USE_MPI
    INTEGER(MPI_ADDRESS_KIND) :: disp(3), base

    CALL MPI_get_address(p%x,disp(1),ierr)
    CALL MPI_get_address(p%y,disp(2),ierr)
    CALL MPI_get_address(p%z,disp(3),ierr)

    base=disp(1)
    disp(1)=disp(1)-base
    disp(2)=disp(2)-base
    disp(3)=disp(3)-base

    blocklen(:)=1
    
    typep(:)=DOUBLE_MPI

    CALL MPI_type_create_struct(3, blocklen, disp, typep,point_3d_mpi, g_ierr)
    CALL Check_Err_MPI()
    CALL MPI_type_commit(point_3d_mpi, g_ierr)
    CALL Check_Err_MPI()
#endif
    WRITE(*,*) 'Type was registered.'

  END SUBROUTINE Register_type_MPI

  SUBROUTINE Abort_MPI
#ifdef USE_MPI
    CALL MPI_abort(COMM_WORLD_MPI,911,g_ierr)
    CALL Check_Err_MPI()
#else
    STOP
#endif
  END SUBROUTINE
  ! SUBROUTINE Isend_iMPI(buf_,count_,datatype_,dest_,&
  !                     &tag_,comm_,request_,ierror_)
  !   INTEGER(ni) :: buf_
  !   INTEGER(ni) :: count_
  !   INTEGER(ni) :: datatype_
  !   INTEGER(ni) :: dest_
  !   INTEGER(ni) :: tag_
  !   INTEGER(ni) :: comm_
  !   INTEGER(ni) :: request_
  !   INTEGER(ni) :: ierror_
  ! END SUBROUTINE

#ifdef USE_MPI
  SUBROUTINE Check_Err_MPI()
    INTEGER(ni) :: ierr, reslen
    CHARACTER (MPI_MAX_ERROR_STRING) :: str
    IF(g_ierr/=MPI_SUCCESS) THEN
      CALL MPI_Error_string(g_ierr,str,reslen,ierr)
      WRITE(*,*) 'MPI error.'
      WRITE(*,'(A)') str(1:reslen)
      IF(ierr/=MPI_SUCCESS) THEN
        WRITE(*,*) "Total Kaput!!!"
      ENDIF
      STOP
    ENDIF
  END SUBROUTINE Check_Err_MPI
#endif

END MODULE WMPI
