!******************************************************************************
!* Multiscale Universal Interface Code Coupling Library                       *
!*                                                                            *
!* Copyright (C) 2021 Y. H. Tang, S. Kudo, X. Bian, Z. Li, G. E. Karniadakis, *
!*                    S. M. Longshaw                                          *
!*                                                                            *
!* This software is jointly licensed under the Apache License, Version 2.0    *
!* and the GNU General Public License version 3, you may use it according     *
!* to either.                                                                 *
!*                                                                            *
!* ** Apache License, version 2.0 **                                          *
!*                                                                            *
!* Licensed under the Apache License, Version 2.0 (the "License");            *
!* you may not use this file except in compliance with the License.           *
!* You may obtain a copy of the License at                                    *
!*                                                                            *
!* http://www.apache.org/licenses/LICENSE-2.0                                 *
!*                                                                            *
!* Unless required by applicable law or agreed to in writing, software        *
!* distributed under the License is distributed on an "AS IS" BASIS,          *
!* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
!* See the License for the specific language governing permissions and        *
!* limitations under the License.                                             *
!*                                                                            *
!* ** GNU General Public License, version 3 **                                *
!*                                                                            *
!* This program is free software: you can redistribute it and/or modify       *
!* it under the terms of the GNU General Public License as published by       *
!* the Free Software Foundation, either version 3 of the License, or          *
!* (at your option) any later version.                                        *
!*                                                                            *
!* This program is distributed in the hope that it will be useful,            *
!* but WITHOUT ANY WARRANTY; without even the implied warranty of             *
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
!* GNU General Public License for more details.                               *
!*                                                                            *
!* You should have received a copy of the GNU General Public License          *
!* along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
!******************************************************************************
!
!** File Details **
!
!Filename: unit_test.F90
!Created: 18 October 2021
!Author: S. M. Longshaw (derived from original 3D unit test by S. Kudo)
!Description: Unit test for Fortran wrapper to create and manage 1D MUI interfaces
!             and associated sampler objects

program main
  use iso_c_binding, only : c_ptr,c_null_char,c_double
  use iso_fortran_env, only : error_unit
  use mui_1d_f

  implicit none

  !Local variables
  character(len=1024) :: arg
  type(c_ptr), target :: uniface_1f
  real(c_double) :: zero=0.0_c_double

  !Read in URI from command line
  if (command_argument_count()==1) then
    call get_command_argument(1,arg)
  else
    ! this shall be properly handled via MPI
    write(error_unit,*)"Wrong number of arguments passed"
    stop 0
  endif

  call mui_create_uniface_1f_f(uniface_1f, trim(arg)//c_null_char)

  call mui_push_1f_f(uniface_1f, "position"//c_null_char, zero, zero)

  call mui_destroy_uniface_1f_f(uniface_1f)
end program main
