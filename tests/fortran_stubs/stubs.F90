! Minimal stand-ins for the CAM modules used by marine_inp_cam.F90 (compile test only).
module shr_kind_mod
  integer, parameter :: shr_kind_r8 = selected_real_kind(12)
end module
module spmd_utils
  logical :: masterproc = .true.
end module
module cam_logfile
  integer :: iulog = 6
end module
module cam_abortutils
contains
  subroutine endrun(msg)
    character(len=*), intent(in) :: msg
    print *, msg
    stop 1
  end subroutine
end module
module mpishorthand
  integer :: mpiint, mpir8, mpichar, mpicom
end module
