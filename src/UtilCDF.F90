!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! purpose:    Module for a CDF/netCDF wrapping
!
! author:     Comrad
! date of creation:   3.03.2006
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE UtilCDF
USE type_mod
! USE netcdf
IMPLICIT NONE
#ifdef NETCDF

#define ASSERT_CDF(x)  CALL check_cdf( (x), __FILE__, __LINE__ )

! One should use this file from sources if the
!   library was compiled with gnu fortran
! Otherwise one should use the original include file and change all
!   calls and folowing defines in order to substruct '_' at the end

!This include must be exactly here!!!
INCLUDE 'netcdf.inc'
INTERFACE get_att_cdf
  MODULE PROCEDURE get_att_int_cdf
  MODULE PROCEDURE get_att_double_cdf
END INTERFACE
INTERFACE put_att_cdf
  MODULE PROCEDURE put_att_text_cdf
  MODULE PROCEDURE put_att_int_cdf
  MODULE PROCEDURE put_att_double_cdf
END INTERFACE
INTERFACE get_var_cdf
  MODULE PROCEDURE get_m2ivar_cdf
  MODULE PROCEDURE get_m2dvar_cdf
  MODULE PROCEDURE get_mdvar_cdf
  MODULE PROCEDURE get_m1dvar_cdf
  MODULE PROCEDURE get_dvar_cdf
  MODULE PROCEDURE get_i1var_cdf
  MODULE PROCEDURE get_ivar_cdf
END INTERFACE
INTERFACE put_var_cdf
  MODULE PROCEDURE put_m4dvar_cdf
  MODULE PROCEDURE put_m2dvar_cdf
  MODULE PROCEDURE put_dvar_cdf
  MODULE PROCEDURE put_i1var_cdf
  MODULE PROCEDURE put_ivar_cdf
END INTERFACE
#endif

INTEGER(ni), PARAMETER :: WRITE_ACC_CDF = 1
INTEGER(ni), PARAMETER :: READ_ACC_CDF = 0

INTEGER(ni), PARAMETER :: ATT_CDF = 1
INTEGER(ni), PARAMETER :: DIM_CDF = 2
INTEGER(ni), PARAMETER :: VAR_CDF = 3

INTEGER(ni), PARAMETER :: SEEK_SET_CDF= 1
INTEGER(ni), PARAMETER :: SEEK_CUR_CDF = 2
! PUBLIC :: def_cdf_dim, enddef_cdf,open_cdf,close_cdf,create_cdf_var
! PRIVATE
INTEGER(ni), ALLOCATABLE, DIMENSION(:) :: p_cdf_vars

#ifdef NETCDF
INTEGER(ni), PARAMETER :: CDF_DOUBLE=nf_double
INTEGER(ni), PARAMETER :: CDF_BYTE=nf_byte
INTEGER(ni), PARAMETER :: CDF_INT=nf_int
INTEGER(ni), PARAMETER :: CDF_SHORT=nf_short
INTEGER(ni), PARAMETER :: CDF_LONG=nf_int
INTEGER(ni), PARAMETER :: CDF_GLOBAL=nf_global
#else
INTEGER(ni), PARAMETER :: CDF_DOUBLE=0
INTEGER(ni), PARAMETER :: CDF_BYTE=0
INTEGER(ni), PARAMETER :: CDF_INT=0
INTEGER(ni), PARAMETER :: CDF_SHORT=0
INTEGER(ni), PARAMETER :: CDF_LONG=0
INTEGER(ni), PARAMETER :: CDF_GLOBAL=0
#endif

#ifdef NETCDF
CONTAINS
!-----------------------------------------------------------------------------
! Sets the position for variable
!   n_ - number of cdf variable
!   offset_ - the offset
!   whence_ - where to count from
!-----------------------------------------------------------------------------
  SUBROUTINE seek_cdf(n_,offset_,whence_)
    INTEGER(ni), INTENT(IN):: n_,offset_
    INTEGER(ni), INTENT(IN), OPTIONAL :: whence_
    ! Internal
    INTEGER(ni) :: whence
    IF(PRESENT(whence_)) THEN
      whence = whence_
    ELSE
      whence = SEEK_SET_CDF
    ENDIF
    IF(whence==SEEK_SET_CDF) THEN
      p_cdf_vars(n_) = offset_
    ELSEIF(whence==SEEK_CUR_CDF) THEN
      p_cdf_vars(n_) = p_cdf_vars(n_) + offset_
    ENDIF
  END SUBROUTINE seek_cdf
!-----------------------------------------------------------------------------
! Sets number of cdf vars
!   n_ - number of cdf variables
!-----------------------------------------------------------------------------
  SUBROUTINE set_num_vars_cdf(n_)
    INTEGER(ni), INTENT(IN):: n_
    ALLOCATE(p_cdf_vars(n_))
    p_cdf_vars = 1
  END SUBROUTINE set_num_vars_cdf
!-----------------------------------------------------------------------------
! Adds a new dimension to an open cdf dataset
!   fid_,name_,len_ - in values
!     fid_ - file id
!     name_ - dimension name
!     len_ - Length of dimension; that is, number of values for this dimension
!       as an index to variables that use it. This should be either a positive 
!       integer or the predefined constant NF90_UNLIMITED. 
!   dimid_ - out value
!     dimid_ - Returned dimension ID
!-----------------------------------------------------------------------------
  SUBROUTINE def_cdf_dim(fid_,name_,len_,dimid_)
    INTEGER(ni), INTENT(IN) :: fid_, len_
    CHARACTER(LEN=*), INTENT(IN) :: name_
    INTEGER(ni), INTENT(OUT) :: dimid_
    ASSERT_CDF( nf_def_dim(fid_,name_,len_,dimid_) )
  END SUBROUTINE def_cdf_dim
!-----------------------------------------------------------------------------
! Open a CDF file
!   fid_ - file id
!-----------------------------------------------------------------------------
  SUBROUTINE enddef_cdf(fid_)
    INTEGER(ni), INTENT(IN) :: fid_
    ASSERT_CDF( nf_enddef(fid_) )
  END SUBROUTINE enddef_cdf
!-----------------------------------------------------------------------------
! Open a CDF file
!   filename_ - just a file name :-)
!   acc_ - access type (read or write)
!   fid_ - file id
!-----------------------------------------------------------------------------
  SUBROUTINE open_cdf(filename_,acc_,fid_)
    CHARACTER(LEN=*), INTENT(IN) :: filename_
    INTEGER(ni), INTENT(IN) :: acc_
    INTEGER(ni), INTENT(OUT) :: fid_
    ! Internal
    INTEGER(ni) :: nvars
    IF(acc_==WRITE_ACC_CDF) THEN
      ASSERT_CDF( nf_create(TRIM(filename_),nf_clobber,fid_) )
    ELSEIF(acc_==READ_ACC_CDF) THEN
      ASSERT_CDF( nf_open(TRIM(filename_),nf_nowrite,fid_) )
      ASSERT_CDF( nf_inq_nvars(fid_,nvars) )
      CALL set_num_vars_cdf(nvars)
    ELSE
      WRITE(*,*) 'Wrong passed access type !!!'
      STOP
    ENDIF
  END SUBROUTINE open_cdf
!-----------------------------------------------------------------------------
! Close a CDF file
!   fid_ - file id
!-----------------------------------------------------------------------------
  SUBROUTINE close_cdf(fid_)
    INTEGER(ni), INTENT(IN) :: fid_
    ASSERT_CDF( nf_close(fid_) )
    IF(ALLOCATED(p_cdf_vars)) DEALLOCATE(p_cdf_vars)
  END SUBROUTINE close_cdf
!-----------------------------------------------------------------------------
! Return the dimension id and length
!   fid_ - file id
!   dimname_ - name of the attribute
!   dimid_ - returned dimension id
!   len_ - returned length
!-----------------------------------------------------------------------------
 SUBROUTINE get_dim_cdf(fid_,dimname_,dimid_,len_)
    INTEGER(ni), INTENT(IN) :: fid_
    CHARACTER(LEN=*), INTENT(IN) :: dimname_
    INTEGER(ni), INTENT(OUT) :: dimid_,len_
    ! Internal
    INTEGER(ni) :: ierr
    ASSERT_CDF(nf_inq_dimid(fid_,dimname_,dimid_))
    ASSERT_CDF(nf_inq_dimlen(fid_,dimid_,len_))
  END SUBROUTINE get_dim_cdf
!-----------------------------------------------------------------------------
! Return the integer attribute either for global or for variable with varid_ id
!   fid_ - file id
!   varid_ - variable id
!   name_ - name of the attribute
!   attint_ - returned attribute value
!-----------------------------------------------------------------------------
  LOGICAL FUNCTION get_att_double_cdf(fid_,varid_,name_,attd_)
    INTEGER(ni), INTENT(IN) :: fid_,varid_
    CHARACTER(LEN=*), INTENT(IN) :: name_
    REAL(DP), INTENT(OUT) :: attd_(*)
    ! Internal
    INTEGER(ni) :: ierr
    get_att_double_cdf=.TRUE.
    ierr = nf_get_att_double(fid_,varid_,name_,attd_)
    IF(ierr/=nf_noerr) THEN
      get_att_double_cdf = .FALSE.
      RETURN
    ENDIF
  END FUNCTION get_att_double_cdf
!-----------------------------------------------------------------------------
! Return the integer attribute either for global or for variable with varid_ id
!   fid_ - file id
!   varid_ - variable id
!   name_ - name of the attribute
!   attint_ - returned attribute value
!-----------------------------------------------------------------------------
  LOGICAL FUNCTION get_att_int_cdf(fid_,varid_,name_,attint_)
    INTEGER(ni), INTENT(IN) :: fid_,varid_
    CHARACTER(LEN=*), INTENT(IN) :: name_
    INTEGER(ni), INTENT(OUT) :: attint_
    ! Internal
    INTEGER(ni) :: ierr
    get_att_int_cdf=.TRUE.
    ierr = nf_get_att_int(fid_,varid_,name_,attint_)
    IF(ierr/=nf_noerr) THEN
      get_att_int_cdf = .FALSE.
      RETURN
    ENDIF
  END FUNCTION get_att_int_cdf
!-----------------------------------------------------------------------------
! Checks existance of global (ONLY!) attributes, dimensions, variables, etc...
!   fid_ - file id
!   type_ - type of different CDF substances
!   name_ - name of the substance
!-----------------------------------------------------------------------------
  LOGICAL FUNCTION exist_cdf(fid_,type_,name_)
    INTEGER(ni), INTENT(IN) :: fid_,type_
    CHARACTER(LEN=*), INTENT(IN) :: name_
    ! Internal
    INTEGER(ni) :: ierr,id
    exist_cdf=.FALSE.
    IF(type_==ATT_CDF) THEN
      ierr = nf_inq_attid(fid_,nf_global,name_,id)
      IF(ierr==nf_noerr) exist_cdf = .TRUE.
    ELSEIF(type_==DIM_CDF) THEN
    ELSEIF(type_==VAR_CDF) THEN
    ELSE
    ENDIF
  END FUNCTION exist_cdf

!-----------------------------------------------------------------------------
! Adds or changes a variable attribute or global attribute
!   fid_,varid_,name_,len_,values_
!     fid_ - cdf file id
!     varid_ - variable id
!     name_ - name of the attribute
!     len_ - number of values provided for the attribute
!     values_ - an array of len_ attribute values
!-----------------------------------------------------------------------------
  SUBROUTINE put_att_text_cdf(fid_,varid_,name_,len_,values_)
    INTEGER(ni), INTENT(IN) :: fid_,varid_,len_
    CHARACTER(LEN=*), INTENT(IN) :: name_,values_
    ASSERT_CDF( nf_put_att_text(fid_,varid_,name_,len_,values_) )
  END SUBROUTINE put_att_text_cdf

  SUBROUTINE put_att_int_cdf(fid_,varid_,name_,len_,values_)
    INTEGER(ni), INTENT(IN) :: fid_,varid_,len_
    CHARACTER(LEN=*), INTENT(IN) :: name_
    INTEGER(ni) :: values_(:)
    ASSERT_CDF( nf_put_att_int(fid_,varid_,name_,nf_int,len_,values_) )
  END SUBROUTINE put_att_int_cdf

  SUBROUTINE put_att_double_cdf(fid_,varid_,name_,len_,values_)
    INTEGER(ni), INTENT(IN) :: fid_,varid_,len_
    CHARACTER(LEN=*), INTENT(IN) :: name_
    REAL(DP) :: values_(:)
    ASSERT_CDF( nf_put_att_double(fid_,varid_,name_,&
                  & nf_double,len_,values_) )
  END SUBROUTINE put_att_double_cdf

!-----------------------------------------------------------------------------
! Create cdf variable with different atributes
!   ncid_,varname_,vartype_,ndims_,vdims_,units_,long_name_,minmax_ - in values 
!     ncid_ - cdf file id
!     varname_ - name of the variable
!     vartype_ - type of the variable
!     ndims_ - number of the variable dimensions
!     vdims_ - vector of the valid dimensions ids
!     units_ - units of the variable
!     long_name_ - long name of the variable (description)
!   varid_ - out value
!     varid_ - variable id
!   minmax_ - valid range of the variable (optional)
!-----------------------------------------------------------------------------
  SUBROUTINE create_cdf_var(ncid_,varname_,vartype_,ndims_,vdims_,units_,&
        &long_name_,varid_,minmax_)
    INTEGER(ni), INTENT(IN) :: ncid_,vartype_,ndims_,vdims_(ndims_)
    CHARACTER(LEN=*), INTENT(IN) :: varname_,units_,long_name_
    INTEGER(ni), INTENT(OUT) :: varid_
    REAL,OPTIONAL, INTENT(IN) :: minmax_(2)
    ASSERT_CDF( nf_def_var(ncid_,varname_,vartype_,ndims_,vdims_,&
      &varid_) )
    ASSERT_CDF( nf_put_att_text(ncid_,varid_,'units',LEN(units_),units_) )
    ASSERT_CDF( nf_put_att_text(ncid_,varid_,'long_name',LEN(long_name_),&
                  & long_name_) )
    IF(PRESENT(minmax_) ) THEN
      ! WRITE(*,*) 'minmax',minmax_(1),minmax_(2)
      ASSERT_CDF( nf_put_att_double(ncid_,varid_,'valid_range',&
                    & nf_double,2,minmax_) )
    ENDIF
  END SUBROUTINE create_cdf_var

!-----------------------------------------------------------------------------
! Return the variable id 
!-----------------------------------------------------------------------------
  INTEGER(ni) FUNCTION get_varid_cdf(fid_,varname_)
    INTEGER(ni), INTENT(IN) :: fid_
    CHARACTER(*), INTENT(IN) :: varname_
    !Internal
    INTEGER(ni) :: id
    ASSERT_CDF( nf_inq_varid(fid_,varname_,id) )
    get_varid_cdf = id
  END FUNCTION get_varid_cdf
!-----------------------------------------------------------------------------
! Return the variable according to the name
!-----------------------------------------------------------------------------
  SUBROUTINE get_dvar_cdf(fid_,varname_,dvals_,count_)
    INTEGER(ni), INTENT(IN) :: fid_
    CHARACTER(*), INTENT(IN) :: varname_
    INTEGER(ni), OPTIONAL, INTENT(IN) :: count_
    REAL(DP), INTENT(OUT) :: dvals_(*)
    !Internal
    INTEGER(ni) :: icount(1),istart(1),varid
    ASSERT_CDF( nf_inq_varid(fid_,varname_,varid) )
    IF(PRESENT(count_)) THEN
      icount(1) = count_
      istart(1) = p_cdf_vars(varid)
      ASSERT_CDF( nf_get_vara_double(fid_,varid,istart,&
                    & icount,dvals_) )
      p_cdf_vars(varid) = p_cdf_vars(varid)+count_
    ELSE
      ASSERT_CDF( nf_get_var_double(fid_,varid,dvals_) )
    ENDIF
  END SUBROUTINE get_dvar_cdf
!-----------------------------------------------------------------------------
  SUBROUTINE get_m2dvar_cdf(fid_,varname_,dvals_)
    INTEGER(ni), INTENT(IN) :: fid_
    CHARACTER(*), INTENT(IN) :: varname_
    REAL(DP), INTENT(OUT) :: dvals_(:,:)
    !Internal
    INTEGER(ni) :: varid
    ASSERT_CDF( nf_inq_varid(fid_,varname_,varid) )
    ASSERT_CDF( nf_get_var_double(fid_,varid,dvals_) )
  END SUBROUTINE get_m2dvar_cdf
!-----------------------------------------------------------------------------
  SUBROUTINE get_m2ivar_cdf(fid_,varname_,ivals_)
    INTEGER(ni), INTENT(IN) :: fid_
    CHARACTER(*), INTENT(IN) :: varname_
    INTEGER(ni),DIMENSION(:,:), INTENT(OUT) :: ivals_
    !Internal
    INTEGER(ni) :: varid
    ASSERT_CDF( nf_inq_varid(fid_,varname_,varid) )
    ASSERT_CDF( nf_get_var_double(fid_,varid,ivals_) )
  END SUBROUTINE get_m2ivar_cdf
!-----------------------------------------------------------------------------
  SUBROUTINE get_m1dvar_cdf(fid_,varname_,dval_,start_,count_)
    INTEGER(ni), INTENT(IN) :: fid_
    CHARACTER(*), INTENT(IN) :: varname_
    INTEGER(ni), INTENT(IN) :: start_(*),count_(*)
    REAL(DP), INTENT(OUT) :: dval_
    !Internal
    INTEGER(ni) :: varid
    ASSERT_CDF( nf_inq_varid(fid_,varname_,varid) )
    ASSERT_CDF( nf_get_vara_double(fid_,varid,start_,&
                  & count_,dval_) )
  END SUBROUTINE get_m1dvar_cdf
!-----------------------------------------------------------------------------
  SUBROUTINE get_mdvar_cdf(fid_,varname_,dvals_,start_,count_)
    INTEGER(ni), INTENT(IN) :: fid_
    CHARACTER(*), INTENT(IN) :: varname_
    INTEGER(ni), INTENT(IN) :: start_(*),count_(*)
    REAL(DP),DIMENSION(*), INTENT(OUT) :: dvals_
    !Internal
    INTEGER(ni) :: varid
    ASSERT_CDF( nf_inq_varid(fid_,varname_,varid) )
    ASSERT_CDF( nf_get_vara_double(fid_,varid,start_,&
                  & count_,dvals_) )
  END SUBROUTINE get_mdvar_cdf
!-----------------------------------------------------------------------------
  SUBROUTINE get_i1var_cdf(fid_,varname_,ivals_,count_)
    INTEGER(ni), INTENT(IN) :: fid_
    CHARACTER(*), INTENT(IN) :: varname_
    INTEGER*1, INTENT(IN) :: ivals_(*)
    INTEGER(ni), OPTIONAL, INTENT(IN) :: count_
    !Internal
    INTEGER(ni) :: icount(1),istart(1),varid
    ASSERT_CDF( nf_inq_varid(fid_,varname_,varid) )
    IF(PRESENT(count_)) THEN
      icount(1) = count_
      istart(1) = p_cdf_vars(varid)
      ASSERT_CDF( nf_get_vara_int1(fid_,varid,istart,&
                    & icount,ivals_) )
      p_cdf_vars(varid) = p_cdf_vars(varid)+count_
    ELSE
      ASSERT_CDF( nf_get_var_int1(fid_,varid,ivals_) )
    ENDIF
  END SUBROUTINE get_i1var_cdf
!!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
  SUBROUTINE get_ivar_cdf(fid_,varname_,ivals_,count_)
    INTEGER(ni), INTENT(IN) :: fid_
    CHARACTER(*), INTENT(IN) :: varname_
    INTEGER(ni), INTENT(IN) :: ivals_(*)
    INTEGER(ni), OPTIONAL, INTENT(IN) :: count_
    !Internal
    INTEGER(ni) :: icount(1),istart(1),varid
    ASSERT_CDF( nf_inq_varid(fid_,varname_,varid) )
    IF(PRESENT(count_)) THEN
      icount(1) = count_
      istart(1) = p_cdf_vars(varid)
      ASSERT_CDF( nf_get_vara_int(fid_,varid,istart,&
                    & icount,ivals_) )
      p_cdf_vars(varid) = p_cdf_vars(varid)+count_
    ELSE
      ASSERT_CDF( nf_get_var_int(fid_,varid,ivals_) )
    ENDIF
  END SUBROUTINE get_ivar_cdf
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
  SUBROUTINE put_m4dvar_cdf(fid_,varid_,dvals_)
    INTEGER(ni), INTENT(IN) :: fid_
    INTEGER(ni), INTENT(IN) :: varid_
    REAL(DP), DIMENSION(:,:,:,:), INTENT(IN) :: dvals_
    ASSERT_CDF( nf_put_var_double(fid_,varid_,dvals_))
  END SUBROUTINE put_m4dvar_cdf
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
  SUBROUTINE put_m2dvar_cdf(fid_,varid_,dvals_)
    INTEGER(ni), INTENT(IN) :: fid_
    INTEGER(ni), INTENT(IN) :: varid_
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: dvals_
    ASSERT_CDF( nf_put_var_double(fid_,varid_,dvals_))
  END SUBROUTINE put_m2dvar_cdf
!-----------------------------------------------------------------------------
  SUBROUTINE put_dvar_cdf(fid_,varid_,dvals_,count_)
    INTEGER(ni), INTENT(IN) :: fid_
    INTEGER(ni), INTENT(IN) :: varid_
    REAL(DP), INTENT(IN) :: dvals_(*)
    INTEGER(ni), INTENT(IN) :: count_
    !Internal
    INTEGER(ni) :: icount(1),istart(1)
    icount(1) = count_
    istart(1) = p_cdf_vars(varid_)
    ASSERT_CDF( nf_put_vara_double(fid_,varid_,istart,&
                  & icount,dvals_) )
    p_cdf_vars(varid_) = p_cdf_vars(varid_)+count_
  END SUBROUTINE put_dvar_cdf
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
  SUBROUTINE put_i1var_cdf(fid_,varid_,ivals_,count_)
    INTEGER(ni), INTENT(IN) :: fid_
    INTEGER(ni), INTENT(IN) :: varid_
    INTEGER*1, INTENT(IN) :: ivals_(*)
    INTEGER(ni), INTENT(IN) :: count_
    !Internal
    INTEGER(ni) :: icount(1),istart(1)
    icount(1) = count_
    istart(1) = p_cdf_vars(varid_)
    ASSERT_CDF( nf_put_vara_int1(fid_,varid_,istart,&
                  & icount,ivals_) )
    p_cdf_vars(varid_) = p_cdf_vars(varid_)+count_
  END SUBROUTINE put_i1var_cdf
!!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
  SUBROUTINE put_ivar_cdf(fid_,varid_,ivals_,count_)
    INTEGER(ni), INTENT(IN) :: fid_
    INTEGER(ni), INTENT(IN) :: varid_
    INTEGER(ni), INTENT(IN) :: ivals_(*)
    INTEGER(ni), INTENT(IN) :: count_
    !Internal
    INTEGER(ni) :: icount(1),istart(1)
    icount(1) = count_
    istart(1) = p_cdf_vars(varid_)
    ASSERT_CDF( nf_put_vara_int(fid_,varid_,istart,&
                  & icount,ivals_) )
    p_cdf_vars(varid_) = p_cdf_vars(varid_)+count_
  END SUBROUTINE put_ivar_cdf
!-----------------------------------------------------------------------------
! Check cdf operation status
!   status - id of the cdf operation status
!-----------------------------------------------------------------------------
  SUBROUTINE check_cdf(status_,filen_,linen_)
    INTEGER(ni), INTENT (IN) :: status_,linen_
    CHARACTER(*), INTENT(IN) :: filen_

    IF(status_ /= nf_noerr) THEN
      WRITE(*,*) 'Error in netcdf.' 
      WRITE(*,*) 'See file:        ',filen_
      WRITE(*,*) 'line number:',linen_
      WRITE(*,*) 'Reason: '
      WRITE(*,*)  trim(nf_strerror(status_))
      STOP "Stopped"
    ENDIF
  END SUBROUTINE check_cdf
#endif
#ifdef HDF5
#endif
END MODULE UtilCDF
