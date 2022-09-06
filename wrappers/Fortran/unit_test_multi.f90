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
  character(:), allocatable :: uri
  type(c_ptr), target :: uniface_1d=c_null_ptr
  type(c_ptr), target :: uniface_2d=c_null_ptr
 ! type(c_ptr), target :: uniface_3d=c_null_ptr
  ! type uniface_3d_p_array
   ! type(c_ptr) :: uniface_p=c_null_ptr
  ! end type uniface_3d_p_array
!  type(uniface_p), target, allocatable :: uniface_3d(:)
  type(c_ptr), target, allocatable :: interfaces3d_ptr(:)
  character(:), allocatable, target :: interfaces3d(:)
  character(1), allocatable, target :: interfaces3dArray(:,:)
  character(:), allocatable :: domain3d
  character(len=1024) :: numbrtSuffix
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
  integer(kind=8) :: i,j

  type(c_ptr) :: MUI_COMM_WORLD

  !Read in URI from command line
  if (command_argument_count()==3) then
    call get_command_argument(1,arg_domain)
    call get_command_argument(2,arg_interface)
    call get_command_argument(3,arg_interface_count)
  else
    print *,"Wrong number of arguments passed: [domain] [interface] [interface_count]"
    stop 0
  endif

  read(arg_interface_count,"(I8)") interface_count

  allocate(uniface_p_array(interface_count))
  allocate(interfaces3d_ptr(interface_count))
  allocate(character(len_trim(arg_interface)+5) :: interfaces3d(interface_count))
  allocate(interfaces3dArray(interface_count,(len_trim(arg_interface)+5)))

  !No need to call mui_mpi_split_by_app() function first as with the _multi example
  !as direct uniface creation also called MPI_Init()

  call mui_mpi_split_by_app_f(MUI_COMM_WORLD)
  !***********************
  !* Single 1D interface *
  !***********************


  !***********************
  !* Single 2D interface *
  !***********************



  !***********************
  !* Single 3D interface *
  !***********************

  uri = "mpi://" // trim(arg_domain) // "/" // trim(arg_interface) // "_3d"

  domain3d = trim(arg_domain)

  do i = 1, interface_count

    if (i < 10) then
        write (numbrtSuffix, "(I1)") i
    else if ((i < 100) .and. (i > 9)) then
        write (numbrtSuffix, "(I2)") i
    else if ((i < 1000) .and. (i > 99)) then
        write (numbrtSuffix, "(I3)") i
    else
        write (numbrtSuffix, "(I4)") i
    endif

    interfaces3d(i) = trim(arg_interface) // "_" // trim(numbrtSuffix)
    ! interfaces3d_ptr(i)= c_loc(interfaces3d(i))

    ! do j=1,len_trim(interfaces3d(i))
	    ! interfaces3dArray(i, j)=interfaces3d(i)(j:j)
		! print *, "interfaces3dArray: ", interfaces3dArray(i, j)
    ! end do
    if (i==1) then
    interfaces3d(i) = "abc"
    else if (i==2) then
    interfaces3d(i) = "222223333305"
    endif
    
	print*,  interfaces3d(i)

  end do 
  
 ! interfaces3d_ptr = c_loc(interfaces3d(1))
  
print *, "interface_count: ", interface_count
print *, "interfaces3d(1): ", interfaces3d(1)
  !Create MUI interface
!  call mui_create_uniface_3d_f(uniface_3d(1), trim(uri)//c_null_char)
   call create_uniface_multi_3d(uniface_p_array, domain3d, interfaces3d, interface_count)
 
! do i = 1, interface_count
    ! print *, "uniface_3d: ", uniface_3d(i), " at ", i
! end do
  do i = 1, interface_count
    !Push send_value_3d=3 through the MUI interface with the tag "position" at location={push_point}={0,0,0}
    call mui_push_3d_f(uniface_p_array(i)%uniface_ptr, "position"//c_null_char, push_point_1, &
    push_point_2, push_point_3, send_value_3d)
  end do


  do i = 1, interface_count
    !Commit (transmit) the pushed value at commit_time=0
    call mui_commit_3d_f(uniface_p_array(i)%uniface_ptr, commit_time)
  end do


  !Create spatial and temporal samplers for fetch operation
  call mui_create_sampler_exact_3d_f(spatial_sampler_exact_3d, tolerance)
  call mui_create_chrono_sampler_exact_3d_f(chrono_sampler_exact_3d, tolerance)


  do i = 1, interface_count
      !Fetch the value for tag "position" at location={fetch_point}={0,0,0} at fetch_time=0
     call mui_fetch_exact_exact_3d_f(uniface_p_array(i)%uniface_ptr, "position"//c_null_char, fetch_point_1, fetch_point_2, &
           fetch_point_3, fetch_time, spatial_sampler_exact_3d, chrono_sampler_exact_3d, fetch_result_3d)

     print *, "Fetched 3D interface value = ",fetch_result_3d
  end do

  ! !Destroy created 3D MUI objects
  ! call mui_destroy_sampler_exact_3d_f(spatial_sampler_exact_3d)
  ! call mui_destroy_chrono_sampler_exact_3d_f(chrono_sampler_exact_3d)
  
! !  Destroy created MUI interfaces note: calls MPI_Finalize(), so need to do last
 ! deallocate(uniface_3d)
 ! deallocate(interfaces3d)
 ! call mui_destroy_uniface_3d_f(uniface_3d)

end program main
