!Multiscale Universal Interface Code Coupling Library
!
!Copyright (C) 2017 Y. H. Tang, S. Kudo, X. Bian, Z. Li, G. E. Karniadakis
!
!This software is jointly licensed under the Apache License, Version 2.0
!and the GNU General Public License version 3, you may use it according
!to either.
!
!** Apache License, version 2.0 **
!
!Licensed under the Apache License, Version 2.0 (the "License");
!you may not use this file except in compliance with the License.
!You may obtain a copy of the License at
!
!http://www.apache.org/licenses/LICENSE-2.0
!
!Unless required by applicable law or agreed to in writing, software
!distributed under the License is distributed on an "AS IS" BASIS,
!WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!See the License for the specific language governing permissions and
!limitations under the License.
!
!** GNU General Public License, version 3 **
!
!This program is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!** File Details **
!
!Filename: mui.F90
!Created:
!Author: S. Kudo
!Description: Fortran wrapper to create 3D MUI uniface.

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
