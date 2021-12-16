!******************************************************************************
!* Multiscale Universal Interface Code Coupling Library                       *
!*                                                                            *
!* Copyright (C) 2021 Y. H. Tang, S. Kudo, X. Bian, Z. Li, G. E. Karniadakis, *
!*                    S. M. Longshaw, W. Liu                                  *
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
!Filename: mui_f_wrapper_general.f90
!Created: 15 September 2021
!Author: S. M. Longshaw (derived from original 3D wrapper by S. Kudo)
!Description: Fortran wrapper to create and manage MUI functions without
!             dimensionality
!

module mui_general_f
  use iso_c_binding
  implicit none
  public

  interface
    !Function to split MPI communicator and return new, local communicator
    subroutine mui_mpi_split_by_app_f(communicator) bind(C)
      import :: c_ptr
      type(c_ptr), intent(out), target :: communicator(*)
    end subroutine mui_mpi_split_by_app_f

    !Function to split MPI communicator and return new, local communicator
    subroutine mui_mpi_split_by_app_threaded_f(communicator,argc,argv,threadType,thread_support) bind(C)
      import :: c_ptr,c_char,c_int
      type(c_ptr), intent(out), target :: thread_support(*)
      type(c_ptr), intent(out), target :: communicator(*)
      integer(kind=c_int), intent(in) :: argc, threadType
      character(kind=c_char), intent(in), dimension(argc) :: argv(*)
    end subroutine mui_mpi_split_by_app_threaded_f

  end interface 

end module mui_general_f
