! © 2023. Triad National Security, LLC. All rights reserved.
! This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National
! Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of
! Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad
! National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration.
! The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up,
! irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute
! copies to the public, perform publicly and display publicly, and to permit others to do so.
!---------------------------------------------
!stub version of h5aux, if HDF5 is not enabled
!---------------------------------------------
module h5aux
  implicit none

  integer, parameter :: HID_T = selected_int_kind(2)

  interface h5_read
     module procedure &
          h5aux_get_1darr, h5aux_get_1darr_i4, h5aux_get_1darr_i8, &
          h5aux_get_2darr, h5aux_get_2darr_i4, h5aux_get_2darr_i8, &
          h5aux_get_3darr, h5aux_get_3darr_i4, h5aux_get_3darr_i8, &
          h5aux_get_4darr, h5aux_get_4darr_i4, h5aux_get_4darr_i8
  end interface h5_read

contains

  function h5_open (fname, readonly_) result(retval)
    character(*), intent(IN) :: fname          ! HDF5 file name
    logical, intent(IN), optional :: readonly_ ! open readonly or r&w?
    integer(HID_T) :: retval
    error stop 'Attempting to open HDF5 file without HDF5 enabled.'
  end function h5_open

  subroutine h5_close (file_id)
    integer(HID_T) :: file_id
    error stop 'Attempting to open HDF5 file without HDF5 enabled.'
  end subroutine h5_close

  subroutine h5aux_get_1darr (data1, obj_id, dset_name)
    double precision, allocatable, intent(OUT) :: data1(:)
    integer(HID_T), intent(IN) :: obj_id    !< file identifier
    character(*), intent(IN)   :: dset_name  !< dataset name
    error stop 'Attempting to read HDF5 file without HDF5 enabled.'
  end subroutine h5aux_get_1darr

  subroutine h5aux_get_1darr_i4 (data1, file_id, dset_name)
    integer(kind=4), allocatable, intent(OUT) :: data1(:)
    integer(HID_T), intent(IN) :: file_id    ! file identifier
    character(*), intent(IN)   :: dset_name  ! dataset name
    error stop 'Attempting to read HDF5 file without HDF5 enabled.'
  end subroutine h5aux_get_1darr_i4

  subroutine h5aux_get_1darr_i8 (data1, file_id, dset_name)
    integer(kind=8), allocatable, intent(OUT) :: data1(:)
    integer(HID_T), intent(IN) :: file_id    ! file identifier
    character(*), intent(IN)   :: dset_name  ! dataset name
    error stop 'Attempting to read HDF5 file without HDF5 enabled.'
  end subroutine h5aux_get_1darr_i8

  subroutine h5aux_get_2darr (data2, obj_id, dset_name)
    double precision, allocatable, intent(OUT) :: data2(:,:)
    integer(HID_T), intent(IN) :: obj_id    !< file identifier
    character(*), intent(IN)   :: dset_name  !< dataset name
    error stop 'Attempting to read HDF5 file without HDF5 enabled.'
  end subroutine h5aux_get_2darr

  subroutine h5aux_get_2darr_i4 (data2, obj_id, dset_name)
    integer(kind=4), allocatable, intent(OUT) :: data2(:,:)
    integer(HID_T), intent(IN) :: obj_id    !< file identifier
    character(*), intent(IN)   :: dset_name  !< dataset name
    error stop 'Attempting to read HDF5 file without HDF5 enabled.'
  end subroutine h5aux_get_2darr_i4

  subroutine h5aux_get_2darr_i8 (data2, obj_id, dset_name)
    integer(kind=8), allocatable, intent(OUT) :: data2(:,:)
    integer(HID_T), intent(IN) :: obj_id    !< file identifier
    character(*), intent(IN)   :: dset_name  !< dataset name
    error stop 'Attempting to read HDF5 file without HDF5 enabled.'
  end subroutine h5aux_get_2darr_i8

  subroutine h5aux_get_3darr (data3, obj_id, dset_name)
    double precision, allocatable, intent(OUT) :: data3(:,:,:)
    integer(HID_T), intent(IN) :: obj_id    !< file identifier
    character(*), intent(IN)   :: dset_name  !< dataset name
    error stop 'Attempting to read HDF5 file without HDF5 enabled.'
  end subroutine h5aux_get_3darr

  subroutine h5aux_get_3darr_i4 (data3, obj_id, dset_name)
    integer(kind=4), allocatable, intent(OUT) :: data3(:,:,:)
    integer(HID_T), intent(IN) :: obj_id    !< file identifier
    character(*), intent(IN)   :: dset_name  !< dataset name
    error stop 'Attempting to read HDF5 file without HDF5 enabled.'
  end subroutine h5aux_get_3darr_i4

  subroutine h5aux_get_3darr_i8 (data3, obj_id, dset_name)
    integer(kind=8), allocatable, intent(OUT) :: data3(:,:,:)
    integer(HID_T), intent(IN) :: obj_id    !< file identifier
    character(*), intent(IN)   :: dset_name  !< dataset name
    error stop 'Attempting to read HDF5 file without HDF5 enabled.'
  end subroutine h5aux_get_3darr_i8

  subroutine h5aux_get_4darr (data4, obj_id, dset_name)
    double precision, allocatable, intent(OUT) :: data4(:,:,:,:)
    integer(HID_T), intent(IN) :: obj_id    !< file identifier
    character(*), intent(IN)   :: dset_name  !< dataset name
    error stop 'Attempting to read HDF5 file without HDF5 enabled.'
  end subroutine h5aux_get_4darr

  subroutine h5aux_get_4darr_i4 (data4, obj_id, dset_name)
    integer(kind=4), allocatable, intent(OUT) :: data4(:,:,:,:)
    integer(HID_T), intent(IN) :: obj_id    !< file identifier
    character(*), intent(IN)   :: dset_name  !< dataset name
    error stop 'Attempting to read HDF5 file without HDF5 enabled.'
  end subroutine h5aux_get_4darr_i4

  subroutine h5aux_get_4darr_i8 (data4, obj_id, dset_name)
    integer(kind=8), allocatable, intent(OUT) :: data4(:,:,:,:)
    integer(HID_T), intent(IN) :: obj_id    !< file identifier
    character(*), intent(IN)   :: dset_name  !< dataset name
    error stop 'Attempting to read HDF5 file without HDF5 enabled.'
  end subroutine h5aux_get_4darr_i8

end module h5aux
