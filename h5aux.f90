!> @file h5aux.f90
!!
! © 2023. Triad National Security, LLC. All rights reserved.
! This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National
! Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of
! Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad
! National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration.
! The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up,
! irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute
! copies to the public, perform publicly and display publicly, and to permit others to do so.
!!
!! Functions and subroutines to simplify HDF5 interface
!! OK 2025.03.27
!! ----------------------
module h5aux
  use hdf5
  implicit none
  !!!include 'mpif.h'

  integer :: number_of_open_files= 0
  private number_of_open_files

  interface h5_read
    module procedure &
      h5aux_get_1darr, h5aux_get_1darr_i4, h5aux_get_1darr_i8, &
      h5aux_get_2darr, h5aux_get_2darr_i4, h5aux_get_2darr_i8, &
      h5aux_get_3darr, h5aux_get_3darr_i4, h5aux_get_3darr_i8, &
      h5aux_get_4darr, h5aux_get_4darr_i4, h5aux_get_4darr_i8
  end interface

  interface h5_write
     module procedure &
      h5aux_write_1darr, h5aux_write_1darr_i4, h5aux_write_1darr_i8, &
      h5aux_write_2darr, h5aux_write_2darr_i4, h5aux_write_2darr_i8, &
      h5aux_write_3darr, h5aux_write_3darr_i4, h5aux_write_3darr_i8, &
      h5aux_write_4darr, h5aux_write_4darr_i4, h5aux_write_4darr_i8
  end interface

  interface h5_write_attr
     module procedure h5aux_add_attr, h5aux_add_attr_i4, h5aux_add_attr_i8
  end interface

  interface h5_read_attr
     module procedure h5aux_get_attr, h5aux_get_attr_i4, h5aux_get_attr_i8
  end interface

contains

!> creates an HDF5 file
function h5_create (fname) result(retval)
character(*), intent(IN) :: fname  ! HDF5 file name
integer(HID_T) :: plist_id, retval
!
integer(4) :: err

   if (number_of_open_files.eq.0) then
      call h5open_f(err)
      if(err.NE.0) THEN
         print "(A,I3,A)", "ERROR(",err,"): cannot initialize HDF5 interface"
         error stop
      endif
   endif

   ! open the file
   call h5fcreate_f(trim(fname), H5F_ACC_TRUNC_F, retval, err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)","ERROR(",err,"): cannot create file ",trim(fname)
      error stop
   endif

   ! increment files counter
   number_of_open_files= number_of_open_files + 1

end function h5_create


!!! !> creates an HDF5 file with parallel access
!!! function h5_create_mpi (fname) result(retval)
!!! character(*), intent(IN) :: fname  ! HDF5 file name
!!! integer(HID_T) :: plist_id, retval
!!! !
!!! integer(4) :: err, comm, info
!!!
!!!    comm= MPI_COMM_WORLD
!!!    info= MPI_INFO_NULL
!!!
!!!    if (number_of_open_files.eq.0) then
!!!       call h5open_f(err)
!!!       if(err.NE.0) THEN
!!!          print "(A,I3,A)", "ERROR(",err,"): cannot initialize HDF5 interface"
!!!          error stop
!!!       endif
!!!    endif
!!!
!!!    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, err)
!!!    if(err.NE.0) THEN
!!!       print "(A,I3,A)", "ERROR(",err,"): cannot set up file access plist"
!!!       error stop
!!!    endif
!!!
!!!    call h5pset_fapl_mpio_f(plist_id, comm, info, err)
!!!    if(err.NE.0) THEN
!!!       print "(A,I3,A)", "ERROR(",err, &
!!!                         "): cannot store MPI comm in file access plist"
!!!       error stop
!!!    endif
!!!
!!!    ! open the file
!!!    call h5fcreate_f(trim(fname), H5F_ACC_TRUNC_F, retval, err, &
!!!                     access_prp=plist_id)
!!!    if(err.NE.0) THEN
!!!       print "(A,I3,A,A)","ERROR(",err,"): cannot create file ",trim(fname)
!!!       error stop
!!!    endif
!!!
!!!    call h5pclose_f(plist_id, err)
!!!
!!!    ! increment files counter
!!!    number_of_open_files= number_of_open_files + 1
!!!
!!! end function h5_create_mpi


!> open existing HDF5 file
function h5_open (fname, readonly_) result(retval)
character(*), intent(IN) :: fname          ! HDF5 file name
logical, intent(IN), optional :: readonly_ ! open readonly or r&w?
integer(HID_T) :: retval
!
integer(4) :: err
logical :: readonly

   readonly= .false.; if (present(readonly_)) readonly= readonly_

   if (number_of_open_files.eq.0) then
      call h5open_f(err)
      if(err.NE.0) THEN
         print "(A,I3,A)", "ERROR(",err,"): cannot initialize HDF5 interface"
         error stop
      endif
   endif

   ! open the file
   if (readonly) THEN
      call h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, retval, err)
   else
      call h5fopen_f(trim(fname), H5F_ACC_RDWR_F, retval, err)
   endif
   if(err.NE.0) THEN
      print "(A,I3,A,A)","ERROR(",err,"): cannot open your file ",trim(fname)
      error stop
   endif

   ! increment files counter
   number_of_open_files= number_of_open_files + 1

end function h5_open


!!! !> open existing HDF5 file in parallel mode
!!! function h5_open_mpi (fname, readonly_) result(retval)
!!! character(*), intent(IN) :: fname          ! HDF5 file name
!!! logical, intent(IN), optional :: readonly_ ! open readonly or r&w?
!!! integer(HID_T) :: plist_id, retval
!!! !
!!! integer(4) :: err, comm, info
!!! logical :: readonly
!!!
!!!    comm= MPI_COMM_WORLD
!!!    info= MPI_INFO_NULL
!!!    readonly= .false.; if (present(readonly_)) readonly= readonly_
!!!
!!!    if (number_of_open_files.eq.0) then
!!!       call h5open_f(err)
!!!       if(err.NE.0) THEN
!!!          print "(A,I3,A)", "ERROR(",err,"): cannot initialize HDF5 interface"
!!!          error stop
!!!       endif
!!!    endif
!!!
!!!    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, err)
!!!    if(err.NE.0) THEN
!!!       print "(A,I3,A)", "ERROR(",err,"): cannot set up file access plist"
!!!       error stop
!!!    endif
!!!
!!!    call h5pset_fapl_mpio_f(plist_id, comm, info, err)
!!!    if(err.NE.0) THEN
!!!       print "(A,I3,A)", "ERROR(",err, &
!!!                         "): cannot store MPI comm in file access plist"
!!!       error stop
!!!    endif
!!!
!!!    ! open the file
!!!    if (readonly) THEN
!!!       call h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, retval, err, &
!!!                     access_prp=plist_id)
!!!    else
!!!       call h5fopen_f(trim(fname), H5F_ACC_RDWR_F, retval, err, &
!!!                     access_prp=plist_id)
!!!    endif
!!!    if(err.NE.0) THEN
!!!       print "(A,I3,A,A)","ERROR(",err,"): cannot open your file ",trim(fname)
!!!       error stop
!!!    endif
!!!
!!!    call h5pclose_f(plist_id, err)
!!!
!!!    ! increment files counter
!!!    number_of_open_files= number_of_open_files + 1
!!!
!!! end function h5_open_mpi


subroutine h5_close (file_id)
integer(HID_T) :: file_id
!
integer(4) :: err

   ! close resources
   call h5fclose_f(file_id, err)
   if(err.NE.0) THEN
      print "(A,I3,A)", "ERROR(",err,"): failed closing HDF5 file"
      error stop
   endif

   ! close hdf5 interface
   if (number_of_open_files.lt.0) then
      call h5close_f(err)
      if(err.NE.0) THEN
         print "(A,I3,A)", "ERROR(",err,"): failed shutting down HDF5 iface"
         error stop
      endif
   endif

   ! increment files counter
   number_of_open_files= number_of_open_files - 1
   if (number_of_open_files.lt.0) &
      error stop "ERROR: internal error with the number of open files"

end subroutine h5_close


!> Adds a double attribute to the root group of the HDF5 file
subroutine h5aux_add_attr (obj_id, attr_name, attr_val)
integer(HID_T), intent(IN)   :: obj_id      !< object identifier
character(*), intent(IN)     :: attr_name   !< name of the attribute
double precision, intent(IN) :: attr_val    !< value to assign to attribute
!
integer(4) :: err
integer(HSIZE_T) :: attr_dims(1)
integer(HID_T) :: asp_id, attr_id

   attr_dims(1)= 1
   call h5screate_simple_f(1_4, attr_dims, asp_id, err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)", &
            "ERROR(",err,"): cannot create attr. space: ",attr_name
      error stop
   endif

   call h5acreate_f(obj_id, attr_name, H5T_NATIVE_DOUBLE, asp_id, attr_id, err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)", &
            "ERROR(",err,"): cannot create attribute: ",attr_name
      error stop
   endif

   call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_val, attr_dims, err)
   if(err.NE.0) THEN
      print "(A,I3,A,I0.1)", &
            "ERROR(",err,"): failed to write attribute value: ",attr_val
      error stop
   endif

   call h5aclose_f(attr_id, err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)", &
            "ERROR(",err,"): failed to close attribute: ",attr_name
      error stop
   endif

end subroutine h5aux_add_attr


!> Read double attribute
subroutine h5aux_get_attr (retval, obj_id, attr_name)
integer(HID_T), intent(IN) :: obj_id     ! object id (object = file or group)
character(*), intent(IN)   :: attr_name  ! attribute name
double precision, intent(OUT) :: retval
!
integer(HSIZE_T) :: adimsf(1)   ! attribute set dimensions
integer(HID_T) :: aid,atid,asid ! attribute, type and space identifiers
integer(4) :: err

   ! attribute wrappers
   adimsf(1)= 1
   call h5screate_simple_f(1_4, adimsf, asid, err)
   if(err.NE.0) THEN
      print "(A,I3,A)", "ERROR(",err,"): cannot create attribute dataspace"
      error stop
   endif

   call h5aopen_f(obj_id, trim(attr_name), aid, err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)", "ERROR(",err,"): open fail on attribute ",attr_name
      error stop
   endif

   call h5aread_f(aid, H5T_NATIVE_DOUBLE, retval, adimsf, err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)", "ERROR(",err,"): read fail on attribute ",attr_name
      error stop
   endif

   ! close attribute
   call h5aclose_f(aid, err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)", "ERROR(",err,"): close fail on attribute",attr_name
      error stop
   endif

   ! close dataspace
   call h5sclose_f(asid, err)
   if(err.NE.0) THEN
      print "(A,I3,A)", "ERROR(",err,"): close fail on dataspace"
      error stop
   endif

end subroutine h5aux_get_attr


!> Adds an integer attribute to the root group of the HDF5 file
subroutine h5aux_add_attr_i4 (obj_id, attr_name, attr_val)
integer(HID_T), intent(IN)  :: obj_id      !< object identifier
character(*), intent(IN)    :: attr_name   !< name of the attribute
integer(kind=4), intent(IN) :: attr_val    !< value to assign to attribute
!
integer(kind=4) :: err
integer(HSIZE_T) :: attr_dims(1)
integer(HID_T) :: asp_id, attr_id

   attr_dims(1)= 1
   call h5screate_simple_f(1_4, attr_dims, asp_id, err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)", &
            "ERROR(",err,"): cannot create attr. space: ",attr_name
      error stop
   endif

   call h5acreate_f(obj_id, attr_name, H5T_STD_I32LE, asp_id, attr_id, err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)", &
            "ERROR(",err,"): cannot create attribute: ",attr_name
      error stop
   endif

   call h5awrite_f(attr_id, H5T_STD_I32LE, attr_val, attr_dims, err)
   if(err.NE.0) THEN
      print "(A,I3,A,I0.1)", &
            "ERROR(",err,"): failed to write attribute value: ",attr_val
      error stop
   endif

   call h5aclose_f(attr_id, err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)", &
            "ERROR(",err,"): failed to close attribute: ",attr_name
      error stop
   endif

end subroutine h5aux_add_attr_i4


!> Adds an integer attribute to the root group of the HDF5 file
subroutine h5aux_add_attr_i8 (obj_id, attr_name, attr_val)
integer(HID_T), intent(IN)  :: obj_id      !< object identifier
character(*), intent(IN)    :: attr_name   !< name of the attribute
integer(kind=8), intent(IN) :: attr_val    !< value to assign to attribute
!
integer(4) :: err
integer(HSIZE_T) :: attr_dims(1)
integer(HID_T) :: asp_id, attr_id

   attr_dims(1)= 1
   call h5screate_simple_f(1_4, attr_dims, asp_id, err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)", &
            "ERROR(",err,"): cannot create attr. space: ",attr_name
      error stop
   endif

   call h5acreate_f(obj_id, attr_name, H5T_STD_I64LE, asp_id, attr_id, err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)", &
            "ERROR(",err,"): cannot create attribute: ",attr_name
      error stop
   endif

   call h5awrite_f(attr_id, H5T_STD_I64LE, attr_val, attr_dims, err)
   if(err.NE.0) THEN
      print "(A,I3,A,I0.1)", &
            "ERROR(",err,"): failed to write attribute value: ",attr_val
      error stop
   endif

   call h5aclose_f(attr_id, err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)", &
            "ERROR(",err,"): failed to close attribute: ",attr_name
      error stop
   endif

end subroutine h5aux_add_attr_i8


!> Return integer attribute
subroutine h5aux_get_attr_i4 (retval, file_id, attr_name)
integer(kind=4), intent(OUT) :: retval
integer(HID_T), intent(IN)   :: file_id    ! file identifier
character(*), intent(IN)     :: attr_name  ! attribute name
!
integer(HSIZE_T) :: adimsf(1)   ! attribute set dimensions
integer(HID_T) :: aid,atid,asid ! attribute, type and space identifiers
integer(4) :: err

   ! attribute wrappers
   adimsf(1)= 1
   call h5screate_simple_f(int(1,kind=4), adimsf, asid, err)
   if(err.NE.0) THEN
      print "(A,I3,A)", "ERROR(",err,"): cannot create attribute dataspace"
      error stop
   endif

   call h5aopen_f(file_id, trim(attr_name), aid, err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)", "ERROR(",err,"): open fail on attribute ",attr_name
      error stop
   endif

   call h5aread_f(aid, H5T_STD_I32LE, retval, adimsf, err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)", "ERROR(",err,"): read fail on attribute ",attr_name
      error stop
   endif

   ! close attribute
   call h5aclose_f(aid, err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)", "ERROR(",err,"): close fail on attribute",attr_name
      error stop
   endif

   ! close dataspace
   call h5sclose_f(asid, err)
   if(err.NE.0) THEN
      print "(A,I3,A)", "ERROR(",err,"): close fail on dataspace"
      error stop
   endif

end subroutine h5aux_get_attr_i4


subroutine h5aux_get_attr_i8 (retval, file_id, attr_name)
integer(kind=8), intent(OUT) :: retval
integer(HID_T), intent(IN)   :: file_id    ! file identifier
character(*), intent(IN)     :: attr_name  ! attribute name
!
integer(HSIZE_T) :: adimsf(1)   ! attribute set dimensions
integer(HID_T) :: aid,atid,asid ! attribute, type and space identifiers
integer(4) :: err
   ! attribute wrappers
   adimsf(1)= 1
   call h5screate_simple_f(int(1,kind=4), adimsf, asid, err)
   if(err.NE.0) THEN
      print "(A,I3,A)", "ERROR(",err,"): cannot create attribute dataspace"
      error stop
   endif

   call h5aopen_f(file_id, trim(attr_name), aid, err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)", "ERROR(",err,"): open fail on attribute ",attr_name
      error stop
   endif

   call h5aread_f(aid, H5T_STD_I64LE, retval, adimsf, err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)", "ERROR(",err,"): read fail on attribute ",attr_name
      error stop
   endif

   ! close attribute
   call h5aclose_f(aid, err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)", "ERROR(",err,"): close fail on attribute",attr_name
      error stop
   endif

   ! close dataspace
   call h5sclose_f(asid, err)
   if(err.NE.0) THEN
      print "(A,I3,A)", "ERROR(",err,"): close fail on dataspace"
      error stop
   endif

end subroutine h5aux_get_attr_i8


!> Writes arbitrary double array to the HDF5 file
subroutine h5aux_write_data (obj_id, dset_name, dset_data, dims)
integer(HID_T), intent(IN)   :: obj_id       !< object identifier
character(*), intent(IN)     :: dset_name    !< name of the dataset
double precision, intent(IN) :: dset_data(*) !< data to write
integer, intent(IN)          :: dims(:)      !< data dimensions
!
integer(4) :: err, ds_rank
integer(HID_T) :: dataspace, dset_id    ! dataspace & dataset handles
integer(HSIZE_T), allocatable :: dims_(:)

   ! create the dataspace
   dims_= int(dims, HSIZE_T)
   call h5screate_simple_f(size(dims_,kind=4), dims_, dataspace, err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)", &
            "ERROR(",err,"): cannot create dataspace for ",dset_name
      error stop
   endif

   ! create the dataset
   call h5dcreate_f(obj_id, dset_name, H5T_NATIVE_DOUBLE, dataspace,dset_id,err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)", "ERROR(",err,"): cannot create dataset: ",dset_name
      error stop
   endif

   ! write the dataset
   call h5dwrite_f (dset_id, H5T_NATIVE_DOUBLE, dset_data, dims_, err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)", "ERROR(",err,"): cannot write dataset: ",dset_name
      error stop
   endif

   ! close dataset and dataspace
   call h5dclose_f (dset_id, err)
   call h5sclose_f (dataspace, err)
   deallocate(dims_)

end subroutine h5aux_write_data


!> Writes 1D dataset to the HDF5 file
subroutine h5aux_write_1darr (obj_id, dset_name, dset_data)
integer(HID_T), intent(IN)   :: obj_id       !< object identifier
character(*), intent(IN)     :: dset_name    !< name of the dataset
double precision, intent(IN) :: dset_data(:) !< data to write
   call h5aux_write_data (obj_id, dset_name, dset_data, shape(dset_data))
end subroutine h5aux_write_1darr


!> Writes 2D dataset to the HDF5 file
subroutine h5aux_write_2darr (obj_id, dset_name, dset_data)
integer(HID_T), intent(IN)   :: obj_id       !< object identifier
character(*), intent(IN)     :: dset_name    !< name of the dataset
double precision, intent(IN) :: dset_data(:,:) !< data to write
   call h5aux_write_data (obj_id, dset_name, dset_data, shape(dset_data))
end subroutine h5aux_write_2darr


!> Writes 3D dataset to the HDF5 file
subroutine h5aux_write_3darr (obj_id, dset_name, dset_data)
integer(HID_T), intent(IN)   :: obj_id       !< object identifier
character(*), intent(IN)     :: dset_name    !< name of the dataset
double precision, intent(IN) :: dset_data(:,:,:) !< data to write
   call h5aux_write_data (obj_id, dset_name, dset_data, shape(dset_data))
end subroutine h5aux_write_3darr


!> Writes 4D dataset to the HDF5 file
subroutine h5aux_write_4darr (obj_id, dset_name, dset_data)
integer(HID_T), intent(IN)   :: obj_id       !< object identifier
character(*), intent(IN)     :: dset_name    !< name of the dataset
double precision, intent(IN) :: dset_data(:,:,:,:) !< data to write
   call h5aux_write_data (obj_id, dset_name, dset_data, shape(dset_data))
end subroutine h5aux_write_4darr


!> Return dataset shape
function h5_get_data_shape (file_id, dset_name) result(dims)
integer(HID_T), intent(IN) :: file_id    ! file identifier
character(*), intent(IN)   :: dset_name  ! dataset name
integer, allocatable, dimension(:) :: dims
!
integer(4) :: err, ds_rank
integer(HID_T) :: dset_id       ! dataset identifier
integer(HID_T) :: dataspace     ! dataspace handle
integer(HSIZE_T), allocatable :: dims_(:), maxdims(:)

   ! open the dataset
   call h5dopen_f(file_id, dset_name, dset_id, err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)", "ERROR(",err,"): cannot open dataset: ",dset_name
      error stop
   endif

   ! get dataset's dataspace handle
   call h5dget_space_f(dset_id, dataspace, err)
   if(err.NE.0) THEN
      print "(A,I3,A)", "ERROR(",err,"): cannot get dataspace handle"
      error stop
   endif

   ! get dataspace's rank
   call h5sget_simple_extent_ndims_f(dataspace, ds_rank, err)
   if(err.NE.0) THEN
      print "(A,I3,A)", "ERROR(",err,"): cannot determine dataspace rank"
      error stop
   endif

   ! allocate dims
   allocate(dims_(ds_rank), dims(ds_rank), maxdims(ds_rank))
   dims_(1:ds_rank)= 0 ! initialize to avoid warnings

   ! query dataspace's dimensinons
   call h5sget_simple_extent_dims_f(dataspace, dims_, maxdims, err)
   if(err.NE.ds_rank) THEN
      print "(A,I3,A)", "ERROR(",err,"): failed to get dataspace dimensions"
      error stop
   endif
   dims= int(dims_)

   ! close dataspace
   call h5sclose_f(dataspace, err)
   if(err.NE.0) THEN
      print "(A,I3,A)", "ERROR(",err,"): close fail on dataspace"
      error stop
   endif

   ! close dataset
   call h5dclose_f(dset_id, err)
   if(err.NE.0) THEN
      print "(A,I3,A)", "ERROR(",err,"): close fail on dataset"
      error stop
   endif

end function h5_get_data_shape


!> Read an array of doubles into 1D array
subroutine h5aux_get_1darr (data1, obj_id, dset_name)
double precision, allocatable, intent(OUT) :: data1(:)
integer(HID_T), intent(IN) :: obj_id    !< file identifier
character(*), intent(IN)   :: dset_name  !< dataset name
!
integer(4) :: err, ds_rank
integer(HID_T) :: dset_id       ! dataset identifier
integer(HID_T) :: dataspace     ! dataspace handle
integer(HSIZE_T), allocatable, dimension(:) :: dims, maxdims ! max dimensions

   ! open the dataset
   call h5dopen_f(obj_id, dset_name, dset_id, err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)", "ERROR(",err,"): cannot open dataset: ",dset_name
      error stop
   endif

   ! get dataset's dataspace handle
   call h5dget_space_f(dset_id, dataspace, err)
   if(err.NE.0) THEN
      print "(A,I3,A)", "ERROR(",err,"): cannot get dataspace handle"
      error stop
   endif

   ! get dataspace's rank
   call h5sget_simple_extent_ndims_f(dataspace, ds_rank, err)
   if(err.NE.0) THEN
      print "(A,I3,A)", "ERROR(",err,"): cannot determine dataspace rank"
      error stop
   endif

   ! allocate dims
   allocate(dims(ds_rank), maxdims(ds_rank))
   dims(1:ds_rank)= 0 ! initialize to avoid warnings

   ! query dataspace's dimensinons
   call h5sget_simple_extent_dims_f(dataspace, dims, maxdims, err)
   if(err.NE.ds_rank) THEN
      print "(A,I3,A)", "ERROR(",err,"): failed to get dataspace dimensions"
      error stop
   endif

   ! read the data into data1
   if (allocated(data1)) deallocate(data1)
   allocate(data1(product(dims(1:ds_rank))))
   call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data1, dims, err)
   if(err.NE.0) THEN
      print "(A,I3,A)", "ERROR(",err,"): in h5dread_f"
      error stop
   endif

   ! close dataspace
   call h5sclose_f(dataspace, err)
   if(err.NE.0) THEN
      print "(A,I3,A)", "ERROR(",err,"): close fail on dataspace"
      error stop
   endif

   ! close dataset
   call h5dclose_f(dset_id, err)
   if(err.NE.0) THEN
      print "(A,I3,A)", "ERROR(",err,"): close fail on dataset"
      error stop
   endif

end subroutine h5aux_get_1darr


!> Read and return integer data as a 1D array
subroutine h5aux_get_1darr_i4 (data1, file_id, dset_name)
integer(kind=4), allocatable, intent(OUT) :: data1(:)
integer(HID_T), intent(IN) :: file_id    ! file identifier
character(*), intent(IN)   :: dset_name  ! dataset name
!
integer(kind=4) :: err, ds_rank
integer(HID_T) :: dset_id       ! dataset identifier
integer(HID_T) :: dataspace     ! dataspace handle
integer(HSIZE_T), allocatable, dimension(:) :: dims, maxdims ! max dimensions

   ! open the dataset
   call h5dopen_f(file_id, dset_name, dset_id, err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)", "ERROR(",err,"): cannot open dataset: ",dset_name
      error stop
   endif

   ! get dataset's dataspace handle
   call h5dget_space_f(dset_id, dataspace, err)
   if(err.NE.0) THEN
      print "(A,I3,A)", "ERROR(",err,"): cannot get dataspace handle"
      error stop
   endif

   ! get dataspace's rank
   call h5sget_simple_extent_ndims_f(dataspace, ds_rank, err)
   if(err.NE.0) THEN
      print "(A,I3,A)", "ERROR(",err,"): cannot determine dataspace rank"
      error stop
   endif

   ! allocate dims
   allocate(dims(ds_rank), maxdims(ds_rank))
   dims(1:ds_rank)= 0 ! initialize to avoid warnings

   ! query dataspace's dimensinons
   call h5sget_simple_extent_dims_f(dataspace, dims, maxdims, err)
   if(err.NE.ds_rank) THEN
      print "(A,I3,A)", "ERROR(",err,"): failed to get dataspace dimensions"
      error stop
   endif

   ! read the data into data1
   if (allocated(data1)) deallocate(data1)
   allocate(data1(product(dims(1:ds_rank))))
   call h5dread_f(dset_id, H5T_STD_I32LE, data1, dims, err)
   if(err.NE.0) THEN
      print "(A,I3,A)", "ERROR(",err,"): in h5dread_f"
      error stop
   endif

   ! close dataspace
   call h5sclose_f(dataspace, err)
   if(err.NE.0) THEN
      print "(A,I3,A)", "ERROR(",err,"): close fail on dataspace"
      error stop
   endif

   ! close dataset
   call h5dclose_f(dset_id, err)
   if(err.NE.0) THEN
      print "(A,I3,A)", "ERROR(",err,"): close fail on dataset"
      error stop
   endif

end subroutine h5aux_get_1darr_i4


!> Read and return integer data as a 1D array
subroutine h5aux_get_1darr_i8 (data1, file_id, dset_name)
integer(kind=8), allocatable, intent(OUT) :: data1(:)
integer(HID_T), intent(IN) :: file_id    ! file identifier
character(*), intent(IN)   :: dset_name  ! dataset name
!
integer(kind=4) :: err, ds_rank
integer(HID_T) :: dset_id       ! dataset identifier
integer(HID_T) :: dataspace     ! dataspace handle
integer(HSIZE_T), allocatable, dimension(:) :: dims, maxdims ! max dimensions

   ! open the dataset
   call h5dopen_f(file_id, dset_name, dset_id, err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)", "ERROR(",err,"): cannot open dataset: ",dset_name
      error stop
   endif

   ! get dataset's dataspace handle
   call h5dget_space_f(dset_id, dataspace, err)
   if(err.NE.0) THEN
      print "(A,I3,A)", "ERROR(",err,"): cannot get dataspace handle"
      error stop
   endif

   ! get dataspace's rank
   call h5sget_simple_extent_ndims_f(dataspace, ds_rank, err)
   if(err.NE.0) THEN
      print "(A,I3,A)", "ERROR(",err,"): cannot determine dataspace rank"
      error stop
   endif

   ! allocate dims
   allocate(dims(ds_rank), maxdims(ds_rank))
   dims(1:ds_rank)= 0 ! initialize to avoid warnings

   ! query dataspace's dimensinons
   call h5sget_simple_extent_dims_f(dataspace, dims, maxdims, err)
   if(err.NE.ds_rank) THEN
      print "(A,I3,A)", "ERROR(",err,"): failed to get dataspace dimensions"
      error stop
   endif

   ! read the data into data1
   if (allocated(data1)) deallocate(data1)
   allocate(data1(product(dims(1:ds_rank))))
   call h5dread_f(dset_id, H5T_STD_I64LE, data1, dims, err)
   if(err.NE.0) THEN
      print "(A,I3,A)", "ERROR(",err,"): in h5dread_f"
      error stop
   endif

   ! close dataspace
   call h5sclose_f(dataspace, err)
   if(err.NE.0) THEN
      print "(A,I3,A)", "ERROR(",err,"): close fail on dataspace"
      error stop
   endif

   ! close dataset
   call h5dclose_f(dset_id, err)
   if(err.NE.0) THEN
      print "(A,I3,A)", "ERROR(",err,"): close fail on dataset"
      error stop
   endif

end subroutine h5aux_get_1darr_i8


!> Read a 2D array
subroutine h5aux_get_2darr (data2, obj_id, dset_name)
double precision, allocatable, intent(OUT) :: data2(:,:)
integer(HID_T), intent(IN) :: obj_id    !< file identifier
character(*), intent(IN)   :: dset_name  !< dataset name
!
integer :: dims(2)
double precision, allocatable :: data1(:)

   call h5aux_get_1darr(data1, obj_id, dset_name)
   dims= h5_get_data_shape(obj_id, dset_name)
   if (size(dims).ne.2) error stop "ERROR: dataset "//dset_name//" is not 2D"
   data2= reshape(data1, dims)

end subroutine h5aux_get_2darr


!> Read a 2D int array
subroutine h5aux_get_2darr_i4 (data2, obj_id, dset_name)
integer(kind=4), allocatable, intent(OUT) :: data2(:,:)
integer(HID_T), intent(IN) :: obj_id    !< file identifier
character(*), intent(IN)   :: dset_name  !< dataset name
!
integer :: dims(2)
integer(kind=4), allocatable :: data1(:)

   call h5aux_get_1darr_i4(data1, obj_id, dset_name)
   dims= h5_get_data_shape(obj_id, dset_name)
   if (size(dims).ne.2) error stop "ERROR: dataset "//dset_name//" is not 2D"
   data2= reshape(data1, dims)

end subroutine h5aux_get_2darr_i4


subroutine h5aux_get_2darr_i8 (data2, obj_id, dset_name)
integer(kind=8), allocatable, intent(OUT) :: data2(:,:)
integer(HID_T), intent(IN) :: obj_id    !< file identifier
character(*), intent(IN)   :: dset_name  !< dataset name
!
integer :: dims(2)
integer(kind=8), allocatable :: data1(:)

   call h5aux_get_1darr_i8(data1, obj_id, dset_name)
   dims= h5_get_data_shape(obj_id, dset_name)
   if (size(dims).ne.2) error stop "ERROR: dataset "//dset_name//" is not 2D"
   data2= reshape(data1, dims)

end subroutine h5aux_get_2darr_i8


!> Read a 3D array
subroutine h5aux_get_3darr (data3, obj_id, dset_name)
double precision, allocatable, intent(OUT) :: data3(:,:,:)
integer(HID_T), intent(IN) :: obj_id    !< file identifier
character(*), intent(IN)   :: dset_name  !< dataset name
!
integer :: dims(3)
double precision, allocatable :: data1(:)

   call h5aux_get_1darr(data1, obj_id, dset_name)
   dims= h5_get_data_shape(obj_id, dset_name)
   if (size(dims).ne.3) error stop "ERROR: dataset "//dset_name//" is not 3D"
   data3= reshape(data1, dims)

end subroutine h5aux_get_3darr


!> Read a 3D int array
subroutine h5aux_get_3darr_i4 (data3, obj_id, dset_name)
integer(kind=4), allocatable, intent(OUT) :: data3(:,:,:)
integer(HID_T), intent(IN) :: obj_id    !< file identifier
character(*), intent(IN)   :: dset_name  !< dataset name
!
integer :: dims(3)
integer(kind=4), allocatable :: data1(:)

   call h5aux_get_1darr_i4(data1, obj_id, dset_name)
   dims= h5_get_data_shape(obj_id, dset_name)
   if (size(dims).ne.3) error stop "ERROR: dataset "//dset_name//" is not 3D"
   data3= reshape(data1, dims)

end subroutine h5aux_get_3darr_i4


subroutine h5aux_get_3darr_i8 (data3, obj_id, dset_name)
integer(kind=8), allocatable, intent(OUT) :: data3(:,:,:)
integer(HID_T), intent(IN) :: obj_id    !< file identifier
character(*), intent(IN)   :: dset_name  !< dataset name
!
integer :: dims(3)
integer(kind=8), allocatable :: data1(:)

   call h5aux_get_1darr_i8(data1, obj_id, dset_name)
   dims= h5_get_data_shape(obj_id, dset_name)
   if (size(dims).ne.3) error stop "ERROR: dataset "//dset_name//" is not 3D"
   data3= reshape(data1, dims)

end subroutine h5aux_get_3darr_i8


!> Read a 4D array
subroutine h5aux_get_4darr (data4, obj_id, dset_name)
double precision, allocatable, intent(OUT) :: data4(:,:,:,:)
integer(HID_T), intent(IN) :: obj_id    !< file identifier
character(*), intent(IN)   :: dset_name  !< dataset name
!
integer :: dims(4)
double precision, allocatable :: data1(:)

   call h5aux_get_1darr(data1, obj_id, dset_name)
   dims= h5_get_data_shape(obj_id, dset_name)
   if (size(dims).ne.4) error stop "ERROR: dataset "//dset_name//" is not 4D"
   data4= reshape(data1, dims)

end subroutine h5aux_get_4darr


!> Read a 4D int array
subroutine h5aux_get_4darr_i4 (data4, obj_id, dset_name)
integer(kind=4), allocatable, intent(OUT) :: data4(:,:,:,:)
integer(HID_T), intent(IN) :: obj_id    !< file identifier
character(*), intent(IN)   :: dset_name  !< dataset name
!
integer :: dims(4)
integer(kind=4), allocatable :: data1(:)

   call h5aux_get_1darr_i4(data1, obj_id, dset_name)
   dims= h5_get_data_shape(obj_id, dset_name)
   if (size(dims).ne.4) error stop "ERROR: dataset "//dset_name//" is not 4D"
   data4= reshape(data1, dims)

end subroutine h5aux_get_4darr_i4


subroutine h5aux_get_4darr_i8 (data4, obj_id, dset_name)
integer(kind=8), allocatable, intent(OUT) :: data4(:,:,:,:)
integer(HID_T), intent(IN) :: obj_id    !< file identifier
character(*), intent(IN)   :: dset_name  !< dataset name
!
integer :: dims(4)
integer(kind=8), allocatable :: data1(:)

   call h5aux_get_1darr_i8(data1, obj_id, dset_name)
   dims= h5_get_data_shape(obj_id, dset_name)
   if (size(dims).ne.4) error stop "ERROR: dataset "//dset_name//" is not 4D"
   data4= reshape(data1, dims)

end subroutine h5aux_get_4darr_i8


subroutine h5aux_write_data_i4 (obj_id, dset_name, dset_data, dims)
integer(HID_T), intent(IN)   :: obj_id       !< object identifier
character(*), intent(IN)     :: dset_name    !< name of the dataset
integer(kind=4), intent(IN)  :: dset_data(*) !< data to write
integer, intent(IN)          :: dims(:)      !< data dimensions
!
integer(4) :: err
integer :: ds_rank
integer(HID_T) :: dataspace, dset_id    ! dataspace handle
integer(HSIZE_T), allocatable :: dims_(:)

   ! create the dataspace
   dims_= int(dims, HSIZE_T)
   call h5screate_simple_f(size(dims_,kind=4), dims_, dataspace, err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)", &
            "ERROR(",err,"): cannot create dataspace for ",dset_name
      error stop
   endif

   ! create the dataset
   call h5dcreate_f(obj_id,dset_name,H5T_STD_I32LE,dataspace,dset_id,err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)", "ERROR(",err,"): cannot create dataset: ",dset_name
      error stop
   endif

   ! write the dataset
   call h5dwrite_f (dset_id, H5T_STD_I32LE, dset_data, dims_, err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)", "ERROR(",err,"): cannot write dataset: ",dset_name
      error stop
   endif

   ! close dataset and dataspace
   call h5dclose_f (dset_id, err)
   call h5sclose_f (dataspace, err)
   deallocate(dims_)

end subroutine h5aux_write_data_i4


subroutine h5aux_write_data_i8 (obj_id, dset_name, dset_data, dims)
integer(HID_T), intent(IN)   :: obj_id       !< object identifier
character(*), intent(IN)     :: dset_name    !< name of the dataset
integer(kind=8), intent(IN)  :: dset_data(*) !< data to write
integer, intent(IN)          :: dims(:)      !< data dimensions
!
integer(4) :: err
integer :: ds_rank
integer(HID_T) :: dataspace, dset_id    ! dataspace handle
integer(HSIZE_T), allocatable :: dims_(:)

   ! create the dataspace
   dims_= int(dims, HSIZE_T)
   call h5screate_simple_f(size(dims_,kind=4), dims_, dataspace, err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)", &
            "ERROR(",err,"): cannot create dataspace for ",dset_name
      error stop
   endif

   ! create the dataset
   call h5dcreate_f(obj_id,dset_name,H5T_STD_I64LE,dataspace,dset_id,err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)", "ERROR(",err,"): cannot create dataset: ",dset_name
      error stop
   endif

   ! write the dataset
   call h5dwrite_f (dset_id, H5T_STD_I64LE, dset_data, dims_, err)
   if(err.NE.0) THEN
      print "(A,I3,A,A)", "ERROR(",err,"): cannot write dataset: ",dset_name
      error stop
   endif

   ! close dataset and dataspace
   call h5dclose_f (dset_id, err)
   call h5sclose_f (dataspace, err)
   deallocate(dims_)

end subroutine h5aux_write_data_i8


subroutine h5aux_write_1darr_i4 (obj_id, dset_name, dset_data)
integer(HID_T), intent(IN)   :: obj_id       !< object identifier
character(*), intent(IN)     :: dset_name    !< name of the dataset
integer(kind=4), intent(IN)  :: dset_data(:) !< data to write
   call h5aux_write_data_i4 (obj_id, dset_name, dset_data, shape(dset_data))
end subroutine h5aux_write_1darr_i4


subroutine h5aux_write_1darr_i8 (obj_id, dset_name, dset_data)
integer(HID_T), intent(IN)   :: obj_id       !< object identifier
character(*), intent(IN)     :: dset_name    !< name of the dataset
integer(kind=8), intent(IN)  :: dset_data(:) !< data to write
   call h5aux_write_data_i8 (obj_id, dset_name, dset_data, shape(dset_data))
end subroutine h5aux_write_1darr_i8


subroutine h5aux_write_2darr_i4 (obj_id, dset_name, dset_data)
integer(HID_T), intent(IN)   :: obj_id         !< object identifier
character(*), intent(IN)     :: dset_name      !< name of the dataset
integer(kind=4), intent(IN)  :: dset_data(:,:) !< data to write
   call h5aux_write_data_i4 (obj_id, dset_name, dset_data, shape(dset_data))
end subroutine h5aux_write_2darr_i4


subroutine h5aux_write_2darr_i8 (obj_id, dset_name, dset_data)
integer(HID_T), intent(IN)   :: obj_id         !< object identifier
character(*), intent(IN)     :: dset_name      !< name of the dataset
integer(kind=8), intent(IN)  :: dset_data(:,:) !< data to write
   call h5aux_write_data_i8 (obj_id, dset_name, dset_data, shape(dset_data))
end subroutine h5aux_write_2darr_i8


subroutine h5aux_write_3darr_i4 (obj_id, dset_name, dset_data)
integer(HID_T), intent(IN)   :: obj_id         !< object identifier
character(*), intent(IN)     :: dset_name      !< name of the dataset
integer(kind=4), intent(IN)  :: dset_data(:,:,:) !< data to write
   call h5aux_write_data_i4 (obj_id, dset_name, dset_data, shape(dset_data))
end subroutine h5aux_write_3darr_i4


subroutine h5aux_write_3darr_i8 (obj_id, dset_name, dset_data)
integer(HID_T), intent(IN)   :: obj_id         !< object identifier
character(*), intent(IN)     :: dset_name      !< name of the dataset
integer(kind=8), intent(IN)  :: dset_data(:,:,:) !< data to write
   call h5aux_write_data_i8 (obj_id, dset_name, dset_data, shape(dset_data))
end subroutine h5aux_write_3darr_i8


subroutine h5aux_write_4darr_i4 (obj_id, dset_name, dset_data)
integer(HID_T), intent(IN)   :: obj_id         !< object identifier
character(*), intent(IN)     :: dset_name      !< name of the dataset
integer(kind=4), intent(IN)  :: dset_data(:,:,:,:) !< data to write
   call h5aux_write_data_i4 (obj_id, dset_name, dset_data, shape(dset_data))
end subroutine h5aux_write_4darr_i4


subroutine h5aux_write_4darr_i8 (obj_id, dset_name, dset_data)
integer(HID_T), intent(IN)   :: obj_id         !< object identifier
character(*), intent(IN)     :: dset_name      !< name of the dataset
integer(kind=8), intent(IN)  :: dset_data(:,:,:,:) !< data to write
   call h5aux_write_data_i8 (obj_id, dset_name, dset_data, shape(dset_data))
end subroutine h5aux_write_4darr_i8


!>  Create a group "group_name" in an HDF5 object
function h5_create_group (obj_id, group_path) result(group_id)
use stringmod, only: split
integer(HID_T), intent(IN) :: obj_id      !< file identifier
character(*), intent(IN)   :: group_path  !< name of the data group
integer(HID_T) :: group_id                !< group identifier
!
integer(4) :: err
integer :: i, num_groups
integer(HID_T) :: gid
character(len=:), allocatable :: groups(:)

   call split(groups, trim(group_path), "/")
   num_groups= size(groups)
   gid= obj_id
   do i=1,num_groups
      call h5gcreate_f(gid, trim(groups(i)), group_id, err)
      if(err.NE.0) THEN
         print "(A,I3,A,A)","ERROR(",err,"): cannot create group: ",groups(i)
         error stop
      endif
      gid= group_id
   enddo

end function h5_create_group


!> Open an existing group
function h5_open_group (obj_id, group_path) result(group_id)
use stringmod, only: split
integer(HID_T), intent(IN) :: obj_id      !< file identifier
character(*), intent(IN)   :: group_path  !< name of the data group
integer(HID_T) :: group_id                !< group identifier
!
integer(4) :: err
integer :: i, num_groups
integer(HID_T) :: gid
character(len=:), allocatable :: groups(:)

   call split(groups, trim(group_path), "/")
   num_groups= size(groups)
   gid= obj_id
   do i=1,num_groups
      call h5gopen_f(gid, trim(groups(i)), group_id, err)
      if(err.NE.0) THEN
         print "(A,I3,A,A)","ERROR(",err,"): cannot open group: ",groups(i)
         error stop
      endif
      gid= group_id
   enddo

end function h5_open_group


!> Close group by group_id (wrapper around h5gclose_f
subroutine h5_close_group (group_id)
integer(HID_T), intent(IN) :: group_id    !< group identifier
integer(4) :: err

   call h5gclose_f(group_id, err)
   if(err.NE.0) THEN
      print "(A,I3,A,I10)","ERROR(",err,"): cannot close group: ",group_id
      error stop
   endif

end subroutine h5_close_group

end module h5aux
