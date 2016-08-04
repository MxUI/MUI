program main
  use iso_c_binding, only : c_ptr,c_null_char, c_double
  use mui_3df, only : mui_create_uniface3d_f, &
    mui_push_f,&
    mui_destroy_uniface3d_f
  use iso_fortran_env, only : error_unit
  implicit none
  type(c_ptr),target :: uniface
  character(len=1024) :: arg

  real(c_double) :: zero=0.0_c_double
  if (command_argument_count()==1) then
    call get_command_argument(1,arg)
  else
    ! this shall be properly handled via MPI
    write(error_unit,*)"Wrong number of arguments passed"
    stop -101
  endif 

  call mui_create_uniface3d_f(uniface, trim(arg)//c_null_char)
  call mui_push_f(uniface, "position"//c_null_char, zero, zero, zero, zero)

  call mui_destroy_uniface3d_f(uniface)
end program main
