module mui_3df
  use iso_c_binding, only : c_ptr, c_char,c_double
  implicit none
  public 
  interface
    subroutine mui_create_uniface3d_f(uniface,domain) bind(c)
      import :: c_ptr,c_char
      type(c_ptr), intent(out), target :: uniface(*)
      character(kind=c_char), intent(in) :: domain(*)
    end subroutine mui_create_uniface3d_f

    subroutine mui_push_f(uniface,attr,x,y,z,val) bind(c)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in),value :: uniface
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double),intent(in) :: x,y,z,val
    end subroutine mui_push_f

    subroutine mui_destroy_uniface3d_f(uniface) bind(c)
      import :: c_ptr
      type(c_ptr), intent(in),value :: uniface
    end subroutine mui_destroy_uniface3d_f
  end interface 
end module mui_3df
