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
!Filename: unit_test_single.f90
!Created: 18 October 2021
!Author: S. M. Longshaw (derived from original 3D unit test by S. Kudo)
!Description: Unit test for Fortran wrapper to create and manage 1D MUI interfaces
!             and associated sampler objects

program main
  use iso_c_binding
  use, intrinsic:: iso_fortran_env, only: stdout=>output_unit, stdin=>input_unit, stderr=>error_unit
  use mui_1d_f
  use mui_2d_f
  use mui_3d_f

  implicit none

  !Local variables
  character(len=1024) :: arg_domain
  character(len=1024) :: arg_interface
  character(:), allocatable :: uri
  type(c_ptr), target :: uniface_1d=c_null_ptr
  type(c_ptr), target :: uniface_2d=c_null_ptr
  type(c_ptr), target :: uniface_3d=c_null_ptr
  type(c_ptr), target :: spatial_sampler_exact_1d=c_null_ptr
  type(c_ptr), target :: spatial_sampler_exact_2d=c_null_ptr
  type(c_ptr), target :: spatial_sampler_exact_3d=c_null_ptr
  type(c_ptr), target :: chrono_sampler_exact_1d=c_null_ptr
  type(c_ptr), target :: chrono_sampler_exact_2d=c_null_ptr
  type(c_ptr), target :: chrono_sampler_exact_3d=c_null_ptr
  real(c_double) :: tolerance=1e-37_c_double
  real(c_double) :: push_point_1=0.0_c_double
  real(c_double) :: push_point_2=0.0_c_double
  real(c_double) :: push_point_3=0.0_c_double
  real(c_double) :: fetch_point_1=0.0_c_double
  real(c_double) :: fetch_point_2=0.0_c_double
  real(c_double) :: fetch_point_3=0.0_c_double
  real(c_double) :: commit_time=0.0_c_double
  real(c_double) :: fetch_time=0.0_c_double
  real(c_double) :: send_value_1d=1.0_c_double
  real(c_double) :: send_value_2d=2.0_c_double
  real(c_double) :: send_value_3d=3.0_c_double
  real(c_double) :: fetch_result_1d=-1_c_double
  real(c_double) :: fetch_result_2d=-2_c_double
  real(c_double) :: fetch_result_3d=-3_c_double

  !Read in URI from command line
  if (command_argument_count()==2) then
    call get_command_argument(1,arg_domain)
    call get_command_argument(2,arg_interface)
  else
    print *,"Wrong number of arguments passed: [domain] [interface]"
    stop 0
  endif

  !No need to call mui_mpi_split_by_app() function first as with the _multi example
  !as direct uniface creation also called MPI_Init()

  !***********************
  !* Single 1D interface *
  !***********************

  uri = "mpi://" // trim(arg_domain) // "/" // trim(arg_interface) // "_1d"

  !Create MUI interface
  call mui_create_uniface_1d_f(uniface_1d, trim(uri)//c_null_char)

  !Push send_value_1d=1 through the MUI interface with the tag "position" at location={push_point}={0}
  call mui_push_1d_f(uniface_1d, "position"//c_null_char, push_point_1, send_value_1d)

  !Commit (transmit) the pushed value at commit_time=0
  call mui_commit_1d_f(uniface_1d, commit_time)

  !Create spatial and temporal samplers for fetch operation
  call mui_create_sampler_exact_1d_f(spatial_sampler_exact_1d, tolerance)
  call mui_create_chrono_sampler_exact_1d_f(chrono_sampler_exact_1d, tolerance)

  !Fetch the value for tag "position" at location={fetch_point}={0} at fetch_time=0
  call mui_fetch_exact_exact_1d_f(uniface_1d, "position"//c_null_char, fetch_point_1, fetch_time, &
       spatial_sampler_exact_1d, chrono_sampler_exact_1d, fetch_result_1d)

  print *, "Fetched 1D interface value = ",fetch_result_1d

  !Destroy created 1D MUI objects
  call mui_destroy_sampler_exact_1d_f(spatial_sampler_exact_1d)
  call mui_destroy_chrono_sampler_exact_1d_f(chrono_sampler_exact_1d)

  !***********************
  !* Single 2D interface *
  !***********************

  uri = "mpi://" // trim(arg_domain) // "/" // trim(arg_interface) // "_2d"

  !Create MUI interface
  call mui_create_uniface_2d_f(uniface_2d, trim(uri)//c_null_char)

  !Push send_value_2d=2 through the MUI interface with the tag "position" at location={push_point}={0,0}
  call mui_push_2d_f(uniface_2d, "position"//c_null_char, push_point_1, push_point_2, send_value_2d)

  !Commit (transmit) the pushed value at commit_time=0
  call mui_commit_2d_f(uniface_2d, commit_time)

  !Create spatial and temporal samplers for fetch operation
  call mui_create_sampler_exact_2d_f(spatial_sampler_exact_2d, tolerance)
  call mui_create_chrono_sampler_exact_2d_f(chrono_sampler_exact_2d, tolerance)

  !Fetch the value for tag "position" at location={fetch_point}={0,0} at fetch_time=0
  call mui_fetch_exact_exact_2d_f(uniface_2d, "position"//c_null_char, fetch_point_1, fetch_point_2, fetch_time, &
       spatial_sampler_exact_2d, chrono_sampler_exact_2d, fetch_result_2d)

  print *, "Fetched 2D interface value = ",fetch_result_2d

  !Destroy created 2D MUI objects
  call mui_destroy_sampler_exact_2d_f(spatial_sampler_exact_2d)
  call mui_destroy_chrono_sampler_exact_2d_f(chrono_sampler_exact_2d)

  !***********************
  !* Single 3D interface *
  !***********************

  uri = "mpi://" // trim(arg_domain) // "/" // trim(arg_interface) // "_3d"

  !Create MUI interface
  call mui_create_uniface_3d_f(uniface_3d, trim(uri)//c_null_char)

  !Push send_value_3d=3 through the MUI interface with the tag "position" at location={push_point}={0,0,0}
  call mui_push_3d_f(uniface_3d, "position"//c_null_char, push_point_1, push_point_2, push_point_3, send_value_3d)

  !Commit (transmit) the pushed value at commit_time=0
  call mui_commit_3d_f(uniface_3d, commit_time)

  !Create spatial and temporal samplers for fetch operation
  call mui_create_sampler_exact_3d_f(spatial_sampler_exact_3d, tolerance)
  call mui_create_chrono_sampler_exact_3d_f(chrono_sampler_exact_3d, tolerance)

  !Fetch the value for tag "position" at location={fetch_point}={0,0,0} at fetch_time=0
  call mui_fetch_exact_exact_3d_f(uniface_3d, "position"//c_null_char, fetch_point_1, fetch_point_2, fetch_point_3, fetch_time, &
       spatial_sampler_exact_3d, chrono_sampler_exact_3d, fetch_result_3d)

  print *, "Fetched 3D interface value = ",fetch_result_3d

  !Destroy created 3D MUI objects
  call mui_destroy_sampler_exact_3d_f(spatial_sampler_exact_3d)
  call mui_destroy_chrono_sampler_exact_3d_f(chrono_sampler_exact_3d)
  
  !Destroy created MUI interfaces note: calls MPI_Finalize(), so need to do last
  call mui_destroy_uniface_1d_f(uniface_1d)
  call mui_destroy_uniface_2d_f(uniface_2d)
  call mui_destroy_uniface_3d_f(uniface_3d)

end program main
