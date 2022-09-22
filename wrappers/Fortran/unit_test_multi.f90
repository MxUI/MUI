!******************************************************************************
!* Multiscale Universal Interface Code Coupling Library                       *
!*                                                                            *
!* Copyright (C) 2021 Y. H. Tang, S. Kudo, X. Bian, Z. Li, G. E. Karniadakis, *
!*                    W. Liu, S. M. Longshaw                                  *
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
!Created: 05 September 2022
!Author:  W, Liu (derived from original 3D unit test by S. Kudo and unit test 
!             single by S. M. Longshaw)
!Description: Unit test for Fortran wrapper to create and manage multiple MUI interfaces
!             and associated sampler objects

program main
  use iso_c_binding
  use, intrinsic:: iso_fortran_env, only: stdout=>output_unit, stdin=>input_unit, stderr=>error_unit
  use mui_1d_f
  use mui_2d_f
  use mui_3d_f
  use mui_general_f

  implicit none

  !Local variables
  character(len=1024) :: arg_domain
  character(len=1024) :: arg_interface
  character(len=1024) :: arg_interface_count
  integer(c_int) :: interface_count
  character(:), allocatable, target :: interfaces1d(:)
  character(:), allocatable, target :: interfaces2d(:)
  character(:), allocatable, target :: interfaces3d(:)
  character(:), allocatable :: domain1d
  character(:), allocatable :: domain2d
  character(:), allocatable :: domain3d
  character(len=1024) :: numberSuffix
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
  real(c_double) :: send_value_1d=3.0_c_double
  real(c_double) :: send_value_2d=3.0_c_double
  real(c_double) :: send_value_3d=3.0_c_double
  real(c_double) :: fetch_result_1d=-3_c_double
  real(c_double) :: fetch_result_2d=-3_c_double
  real(c_double) :: fetch_result_3d=-3_c_double
  type(c_ptr) :: MUI_COMM_WORLD
  integer(kind=8) :: i

  !Read in URI from command line
  if (command_argument_count()==3) then
    call get_command_argument(1,arg_domain)
    call get_command_argument(2,arg_interface)
    call get_command_argument(3,arg_interface_count)
  else
    print *,"Wrong number of arguments passed: [domain] [interface] [interface_count]"
    stop 0
  endif

  !Convert the third command argument into interger as the interface count
  read(arg_interface_count,"(I8)") interface_count

  !Call mui_mpi_split_by_app_f() function to init MPI
  call mui_mpi_split_by_app_f(MUI_COMM_WORLD)

  !*************************
  !* multiple 1D interface *
  !*************************

  !Allociate memory based on number of interfaces
  allocate(character(len_trim(arg_interface)+5) :: interfaces1d(interface_count))
  !For multi-domain function, "uniface_pointers_1d" should be used to collect the array of
  ! MUI uniface pointers. It is decleared in the MUI FORTRAN wrapper.
  allocate(uniface_pointers_1d(interface_count))

  ! Obtain the domain name
  domain1d = trim(arg_domain)

  ! Create interface names
  do i = 1, interface_count
    !Generate character type of number suffix
    if (i < 10) then
        write (numberSuffix, "(I1)") i
    else if ((i < 100) .and. (i > 9)) then
        write (numberSuffix, "(I2)") i
    else if ((i < 1000) .and. (i > 99)) then
        write (numberSuffix, "(I3)") i
    else
        write (numberSuffix, "(I4)") i
    endif

    !Create and collect interface names
    interfaces1d(i) = trim(arg_interface) // "_" // trim(numberSuffix)
  end do 

  !Create MUI interfaces. MUI interfaces will be collected by the "uniface_pointers_1d" after this subroutine
  call create_and_get_uniface_multi_1d_f(uniface_pointers_1d, domain1d, interfaces1d, interface_count)

  !Push send_value_1d=3 through all the MUI interfaces with the tag "position" at location={push_point}={0}
  do i = 1, interface_count
    call mui_push_1d_f(uniface_pointers_1d(i)%ptr, "position"//c_null_char, push_point_1, &
    send_value_1d)
  end do

  !Commit (transmit) the pushed value at commit_time=0 for all MUI interfaces
  do i = 1, interface_count
    call mui_commit_1d_f(uniface_pointers_1d(i)%ptr, commit_time)
  end do

  !Create spatial and temporal samplers for fetch operation
  call mui_create_sampler_exact_1d_f(spatial_sampler_exact_1d, tolerance)
  call mui_create_chrono_sampler_exact_1d_f(chrono_sampler_exact_1d, tolerance)

  !Fetch the value for tag "position" at location={fetch_point}={0} at fetch_time=0 for all MUI interfaces
  do i = 1, interface_count
     call mui_fetch_exact_exact_1d_f(uniface_pointers_1d(i)%ptr, "position"//c_null_char, fetch_point_1,  &
           fetch_time, spatial_sampler_exact_1d, chrono_sampler_exact_1d, fetch_result_1d)
     print *, "Fetched 1D interface value = ",fetch_result_1d, " at the interface of ", interfaces1d(i),             &
           " domian of ", domain1d
  end do

  !Destroy created 1D MUI objects
  call mui_destroy_sampler_exact_1d_f(spatial_sampler_exact_1d)
  call mui_destroy_chrono_sampler_exact_1d_f(chrono_sampler_exact_1d)

  ! Destroy created MUI interfaces note: calls MPI_Finalize(), so need to do last
  do i = 1, interface_count
     call mui_destroy_uniface_1d_f(uniface_pointers_1d(i)%ptr)
  end do

  !Release the memory on unifaces and interface names
  deallocate(uniface_pointers_1d)
  deallocate(interfaces1d)

  !*************************
  !* multiple 2D interface *
  !*************************

  !Allociate memory based on number of interfaces
  allocate(character(len_trim(arg_interface)+5) :: interfaces2d(interface_count))
  !For multi-domain function, "uniface_pointers_2d" should be used to collect the array of
  ! MUI uniface pointers. It is decleared in the MUI FORTRAN wrapper.
  allocate(uniface_pointers_2d(interface_count))

  ! Obtain the domain name
  domain2d = trim(arg_domain)

  ! Create interface names
  do i = 1, interface_count
    !Generate character type of number suffix
    if (i < 10) then
        write (numberSuffix, "(I1)") i
    else if ((i < 100) .and. (i > 9)) then
        write (numberSuffix, "(I2)") i
    else if ((i < 1000) .and. (i > 99)) then
        write (numberSuffix, "(I3)") i
    else
        write (numberSuffix, "(I4)") i
    endif

    !Create and collect interface names
    interfaces2d(i) = trim(arg_interface) // "_" // trim(numberSuffix)
  end do 

  !Create MUI interfaces. MUI interfaces will be collected by the "uniface_pointers_2d" after this subroutine
  call create_and_get_uniface_multi_2d_f(uniface_pointers_2d, domain2d, interfaces2d, interface_count)

  !Push send_value_2d=3 through all the MUI interfaces with the tag "position" at location={push_point}={0,0}
  do i = 1, interface_count
    call mui_push_2d_f(uniface_pointers_2d(i)%ptr, "position"//c_null_char, push_point_1, &
    push_point_2, send_value_2d)
  end do

  !Commit (transmit) the pushed value at commit_time=0 for all MUI interfaces
  do i = 1, interface_count
    call mui_commit_2d_f(uniface_pointers_2d(i)%ptr, commit_time)
  end do

  !Create spatial and temporal samplers for fetch operation
  call mui_create_sampler_exact_2d_f(spatial_sampler_exact_2d, tolerance)
  call mui_create_chrono_sampler_exact_2d_f(chrono_sampler_exact_2d, tolerance)

  !Fetch the value for tag "position" at location={fetch_point}={0,0} at fetch_time=0 for all MUI interfaces
  do i = 1, interface_count
     call mui_fetch_exact_exact_2d_f(uniface_pointers_2d(i)%ptr, "position"//c_null_char, fetch_point_1, fetch_point_2, &
           fetch_time, spatial_sampler_exact_2d, chrono_sampler_exact_2d, fetch_result_2d)
     print *, "Fetched 2D interface value = ",fetch_result_2d, " at the interface of ", interfaces2d(i),             &
           " domian of ", domain2d
  end do

  !Destroy created 2D MUI objects
  call mui_destroy_sampler_exact_2d_f(spatial_sampler_exact_2d)
  call mui_destroy_chrono_sampler_exact_2d_f(chrono_sampler_exact_2d)

  ! Destroy created MUI interfaces note: calls MPI_Finalize(), so need to do last
  do i = 1, interface_count
     call mui_destroy_uniface_2d_f(uniface_pointers_2d(i)%ptr)
  end do

  !Release the memory on unifaces and interface names
  deallocate(uniface_pointers_2d)
  deallocate(interfaces2d)


  !*************************
  !* multiple 3D interface *
  !*************************

  !Allociate memory based on number of interfaces
  allocate(character(len_trim(arg_interface)+5) :: interfaces3d(interface_count))
  !For multi-domain function, "uniface_pointers_3d" should be used to collect the array of
  ! MUI uniface pointers. It is decleared in the MUI FORTRAN wrapper.
  allocate(uniface_pointers_3d(interface_count))

  ! Obtain the domain name
  domain3d = trim(arg_domain)

  ! Create interface names
  do i = 1, interface_count
    !Generate character type of number suffix
    if (i < 10) then
        write (numberSuffix, "(I1)") i
    else if ((i < 100) .and. (i > 9)) then
        write (numberSuffix, "(I2)") i
    else if ((i < 1000) .and. (i > 99)) then
        write (numberSuffix, "(I3)") i
    else
        write (numberSuffix, "(I4)") i
    endif

    !Create and collect interface names
    interfaces3d(i) = trim(arg_interface) // "_" // trim(numberSuffix)
  end do 

  !Create MUI interfaces. MUI interfaces will be collected by the "uniface_pointers_3d" after this subroutine
  call create_and_get_uniface_multi_3d_f(uniface_pointers_3d, domain3d, interfaces3d, interface_count)

  !Push send_value_3d=3 through all the MUI interfaces with the tag "position" at location={push_point}={0,0,0}
  do i = 1, interface_count
    call mui_push_3d_f(uniface_pointers_3d(i)%ptr, "position"//c_null_char, push_point_1, &
    push_point_2, push_point_3, send_value_3d)
  end do

  !Commit (transmit) the pushed value at commit_time=0 for all MUI interfaces
  do i = 1, interface_count
    call mui_commit_3d_f(uniface_pointers_3d(i)%ptr, commit_time)
  end do

  !Create spatial and temporal samplers for fetch operation
  call mui_create_sampler_exact_3d_f(spatial_sampler_exact_3d, tolerance)
  call mui_create_chrono_sampler_exact_3d_f(chrono_sampler_exact_3d, tolerance)

  !Fetch the value for tag "position" at location={fetch_point}={0,0,0} at fetch_time=0 for all MUI interfaces
  do i = 1, interface_count
     call mui_fetch_exact_exact_3d_f(uniface_pointers_3d(i)%ptr, "position"//c_null_char, fetch_point_1, fetch_point_2, &
           fetch_point_3, fetch_time, spatial_sampler_exact_3d, chrono_sampler_exact_3d, fetch_result_3d)
     print *, "Fetched 3D interface value = ",fetch_result_3d, " at the interface of ", interfaces3d(i),             &
           " domian of ", domain3d
  end do

  !Destroy created 3D MUI objects
  call mui_destroy_sampler_exact_3d_f(spatial_sampler_exact_3d)
  call mui_destroy_chrono_sampler_exact_3d_f(chrono_sampler_exact_3d)

  ! Destroy created MUI interfaces note: calls MPI_Finalize(), so need to do last
  do i = 1, interface_count
     call mui_destroy_uniface_3d_f(uniface_pointers_3d(i)%ptr)
  end do

  !Release the memory on unifaces and interface names
  deallocate(uniface_pointers_3d)
  deallocate(interfaces3d)

end program main
