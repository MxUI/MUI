program main
  implicit none
  integer, pointer :: uniface
  character(len=1024) arg
  external mui_create_uniface3d_f, mui_push_f, mui_destroy_uniface3d_f

  call getarg(2,arg)
  call mui_create_uniface3d_f(uniface, arg)

  call mui_push_f(uniface, "position", 0., 0., 0., 0.)

  call mui_destroy_uniface3d_f(uniface)
end program main
