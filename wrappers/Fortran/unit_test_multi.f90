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
!Filename: unit_test_multi.f90
!Created: 25 November 2021
!Author: S. M. Longshaw (derived from original 3D unit test by S. Kudo)
!Description: Unit test for Fortran wrapper to create and manage 1D MUI interfaces
!             and associated sampler objects using helper function to create multiple
!             for a single domain

program main
  use iso_c_binding
  use, intrinsic:: iso_fortran_env, only: stdout=>output_unit, stdin=>input_unit, stderr=>error_unit
  use mui_1d_f

  implicit none

  !Local variables
  character(len=1024) :: arg_domain
  character(len=1024) :: arg_interface
  character(:), allocatable :: domain
  integer(c_int) :: num_interfaces=0
  type(c_ptr), target :: uniface_1d=c_null_ptr
  type(c_ptr), target :: spatial_sampler_exact_1d=c_null_ptr
  type(c_ptr), target :: chrono_sampler_exact_1d=c_null_ptr
  real(c_double) :: tolerance=1e-37_c_double
  real(c_double) :: push_point_1=0.0_c_double
  real(c_double) :: fetch_point_1=0.0_c_double
  real(c_double) :: commit_time=0.0_c_double
  real(c_double) :: fetch_time=0.0_c_double
  real(c_double) :: send_value=1.0_c_double
  real(c_double) :: fetch_result=-1_c_double

  !Read in URI from command line
  if (command_argument_count()==3) then
    call get_command_argument(1,arg_domain)
    call get_command_argument(2,arg_interface)
    call get_command_argument(3,num_interfaces)
  else
    print *,"Wrong number of arguments passed: [domain] [interface] [interface_count]"
    stop 0
  endif

  !Initialise MPI (returns an MPI comm world we might use in a local MPI_Init() call)
  mui_mpi_split_by_app_f();

  !***********************
  !* Multi 1D interfaces *
  !***********************

  domain = trim(arg_domain) // "_1d"

  !Create MUI interface
  call mui_create_uniface_1d_f(uniface_1d, trim(uri)//c_null_char)

  !Push send_value=1 through the MUI interface with the tag "position" at location={push_point}={0}
  call mui_push_1d_f(uniface_1d, "position"//c_null_char, push_point_1, send_value)

  !Commit (transmit) the pushed value at commit_time=0
  call mui_commit_1d_f(uniface_1d, commit_time)

  !Create spatial and temporal samplers for fetch operation
  call mui_create_sampler_exact_1d_f(spatial_sampler_exact_1d, tolerance)
  call mui_create_chrono_sampler_exact_1d_f(chrono_sampler_exact_1d, tolerance)

  !Fetch the value for tag "position" at location={fetch_point}={0} at fetch_time=0
  call mui_fetch_exact_exact_1d_f(uniface_1d, "position"//c_null_char, fetch_point_1, fetch_time, &
       spatial_sampler_exact_1d, chrono_sampler_exact_1d, fetch_result)

  print *, "Fetched 1D interface value = ",fetch_result

  !Destroy created 1D MUI objects
  call mui_destroy_sampler_exact_1d_f(spatial_sampler_exact_1d)
  call mui_destroy_chrono_sampler_exact_1d_f(chrono_sampler_exact_1d)
  call mui_destroy_uniface_1d_f(uniface_1d)

end program main