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
!Filename: mui_f_wrapper_2d.f90
!Created: 08 December 2021
!Author: S. M. Longshaw (derived from original 3D wrapper by S. Kudo)
!Description: Fortran wrapper to create and manage 2D MUI interfaces
!             and associated sampler objects
!
!             NOTE: Any point co-ordinates are enumerated rather than assuming
!                   Cartesian form, i.e. {1, 2, 3} rather than {x, y, z}.

module mui_2d_f
  use iso_c_binding
  implicit none
  public

  !Define pointer type to collcet uniface pointers for multi-domian function
  type ptr_typ_2d
    type(c_ptr) :: ptr
  end type ptr_typ_2d

  !Create an allocatable array to collect MUI uniface pointers with the type of
  ! ptr_typ_2d for multi-domian function
  type(ptr_typ_2d), target, save, allocatable :: uniface_pointers_2d(:)

  interface
    !****************************************
    !* Create MUI interfaces                *
    !****************************************

    !2D interface with float=single and int=int32
    subroutine mui_create_uniface_2f_f(uniface,domain) bind(C)
      import :: c_ptr,c_char
      type(c_ptr), intent(out), target :: uniface(*)
      character(kind=c_char), intent(in) :: domain(*)
    end subroutine mui_create_uniface_2f_f

    !2D interface with float=single and int=int64
    subroutine mui_create_uniface_2fx_f(uniface,domain) bind(C)
      import :: c_ptr,c_char
      type(c_ptr), intent(out), target :: uniface(*)
      character(kind=c_char), intent(in) :: domain(*)
    end subroutine mui_create_uniface_2fx_f

    !2D interface with float=double and int=int32
    subroutine mui_create_uniface_2d_f(uniface,domain) bind(C)
      import :: c_ptr,c_char
      type(c_ptr), intent(out), target :: uniface(*)
      character(kind=c_char), intent(in) :: domain(*)
    end subroutine mui_create_uniface_2d_f

    !2D interface with float=double and int=int64
    subroutine mui_create_uniface_2dx_f(uniface,domain) bind(C)
      import :: c_ptr,c_char
      type(c_ptr), intent(out), target :: uniface(*)
      character(kind=c_char), intent(in) :: domain(*)
    end subroutine mui_create_uniface_2dx_f

    !2D interface using config from config_f_wrapper.h
    subroutine mui_create_uniface_2t_f(uniface,domain) bind(C)
      import :: c_ptr,c_char
      type(c_ptr), intent(out), target :: uniface(*)
      character(kind=c_char), intent(in) :: domain(*)
    end subroutine mui_create_uniface_2t_f

    !Set of 2D interfaces with float=single and int=int32
    !Recomend to use the create_and_get_uniface_multi_2f_f(*) subroutine instead of use
    ! this subroutine directly
    subroutine mui_create_uniface_multi_2f_f(domain, interfaces, interface_count) bind(C)
      import :: c_char,c_int
      character(kind=c_char), intent(in) :: domain(*)
      character(kind=c_char,len=*), intent(in) :: interfaces(*)
      integer(kind=c_int), VALUE :: interface_count
    end subroutine mui_create_uniface_multi_2f_f

    !Set of 2D interfaces with float=single and int=int64
    !Recomend to use the create_and_get_uniface_multi_2fx_f(*) subroutine instead of use
    ! this subroutine directly
    subroutine mui_create_uniface_multi_2fx_f(domain, interfaces, interface_count) bind(C)
      import :: c_char,c_int
      character(kind=c_char), intent(in) :: domain(*)
      character(kind=c_char,len=*), intent(in) :: interfaces(*)
      integer(kind=c_int), VALUE :: interface_count
    end subroutine mui_create_uniface_multi_2fx_f

    !Set of 2D interfaces with float=double and int=int32
    !Recomend to use the create_and_get_uniface_multi_2d_f(*) subroutine instead of use
    ! this subroutine directly
    subroutine mui_create_uniface_multi_2d_f(domain, interfaces, interface_count) bind(C)
      import :: c_char,c_int
      character(kind=c_char), intent(in) :: domain(*)
      character(kind=c_char,len=*), intent(in) :: interfaces(*)
      integer(kind=c_int), VALUE :: interface_count
    end subroutine mui_create_uniface_multi_2d_f

    !Set of 2D interfaces with float=double and int=int64
    !Recomend to use the create_and_get_uniface_multi_2dx_f(*) subroutine instead of use
    ! this subroutine directly
    subroutine mui_create_uniface_multi_2dx_f(domain, interfaces, interface_count) bind(C)
      import :: c_char,c_int
      character(kind=c_char), intent(in) :: domain(*)
      character(kind=c_char,len=*), intent(in) :: interfaces(*)
      integer(kind=c_int), VALUE :: interface_count
    end subroutine mui_create_uniface_multi_2dx_f

    !Set of 2D interfaces using config from config_c_wrapper.h
    !Recomend to use the create_and_get_uniface_multi_2t_f(*) subroutine instead of use
    ! this subroutine directly
    subroutine mui_create_uniface_multi_2t_f(domain, interfaces, interface_count) bind(C)
      import :: c_char,c_int
      character(kind=c_char), intent(in) :: domain(*)
      character(kind=c_char,len=*), intent(in) :: interfaces(*)
      integer(kind=c_int), VALUE :: interface_count
    end subroutine mui_create_uniface_multi_2t_f

    !Access to MUI set of 2D interfaces with float=single and int=int32
    !Recomend to use the create_and_get_uniface_multi_2f_f(*) subroutine instead of use
    ! this subroutine directly
    type(c_ptr) function get_mui_uniface_multi_2f_f(interface_count) bind(C)
      import :: c_ptr,c_int
      integer(kind=c_int), VALUE :: interface_count
    end function get_mui_uniface_multi_2f_f

    !Access to MUI set of 2D interfaces with float=single and int=int64
    !Recomend to use the create_and_get_uniface_multi_2fx_f(*) subroutine instead of use
    ! this subroutine directly
    type(c_ptr) function get_mui_uniface_multi_2fx_f(interface_count) bind(C)
      import :: c_ptr,c_int
      integer(kind=c_int), VALUE :: interface_count
    end function get_mui_uniface_multi_2fx_f

    !Access to MUI set of 2D interfaces with float=double and int=int32
    !Recomend to use the create_and_get_uniface_multi_2d_f(*) subroutine instead of use
    ! this subroutine directly
    type(c_ptr) function get_mui_uniface_multi_2d_f(interface_count) bind(C)
      import :: c_ptr,c_int
      integer(kind=c_int), VALUE :: interface_count
    end function get_mui_uniface_multi_2d_f

    !Access to MUI set of 2D interfaces with float=double and int=int64
    !Recomend to use the create_and_get_uniface_multi_2dx_f(*) subroutine instead of use
    ! this subroutine directly
    type(c_ptr) function get_mui_uniface_multi_2dx_f(interface_count) bind(C)
      import :: c_ptr,c_int
      integer(kind=c_int), VALUE :: interface_count
    end function get_mui_uniface_multi_2dx_f

    !Access to MUI set of 2D interfaces using config from config_f_wrapper.h
    !Recomend to use the create_and_get_uniface_multi_2t_f(*) subroutine instead of use
    ! this subroutine directly
    type(c_ptr) function get_mui_uniface_multi_2t_f(interface_count) bind(C)
      import :: c_ptr,c_int
      integer(kind=c_int), VALUE :: interface_count
    end function get_mui_uniface_multi_2t_f

    !****************************************
    !* Destroy MUI interface                *
    !****************************************

    subroutine mui_destroy_uniface_2f_f(uniface) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: uniface
    end subroutine mui_destroy_uniface_2f_f

    subroutine mui_destroy_uniface_2fx_f(uniface) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: uniface
    end subroutine mui_destroy_uniface_2fx_f

    subroutine mui_destroy_uniface_2d_f(uniface) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: uniface
    end subroutine mui_destroy_uniface_2d_f

    subroutine mui_destroy_uniface_2dx_f(uniface) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: uniface
    end subroutine mui_destroy_uniface_2dx_f

    subroutine mui_destroy_uniface_2t_f(uniface) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: uniface
    end subroutine mui_destroy_uniface_2t_f

    !******************************************
    !* Create 2D spatial samplers             *
    !******************************************

    !Exact sampler
    subroutine mui_create_sampler_exact_2f_f(sampler,tolerance) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_float), intent(in) :: tolerance
    end subroutine mui_create_sampler_exact_2f_f

    subroutine mui_create_sampler_exact_2fx_f(sampler,tolerance) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_float), intent(in) :: tolerance
    end subroutine mui_create_sampler_exact_2fx_f

    subroutine mui_create_sampler_exact_2d_f(sampler,tolerance) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double), intent(in) :: tolerance
    end subroutine mui_create_sampler_exact_2d_f

    subroutine mui_create_sampler_exact_2dx_f(sampler,tolerance) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double), intent(in) :: tolerance
    end subroutine mui_create_sampler_exact_2dx_f

    subroutine mui_create_sampler_exact_2t_f(sampler,tolerance) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double), intent(in) :: tolerance
    end subroutine mui_create_sampler_exact_2t_f

    !Gauss sampler
    subroutine mui_create_sampler_gauss_2f_f(sampler,r,h) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_float), intent(in) :: r,h
    end subroutine mui_create_sampler_gauss_2f_f

    subroutine mui_create_sampler_gauss_2fx_f(sampler,r,h) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_float), intent(in) :: r,h
    end subroutine mui_create_sampler_gauss_2fx_f

    subroutine mui_create_sampler_gauss_2d_f(sampler,r,h) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double), intent(in) :: r,h
    end subroutine mui_create_sampler_gauss_2d_f

    subroutine mui_create_sampler_gauss_2dx_f(sampler,r,h) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double), intent(in) :: r,h
    end subroutine mui_create_sampler_gauss_2dx_f

    subroutine mui_create_sampler_gauss_2t_f(sampler,r,h) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double), intent(in) :: r,h
    end subroutine mui_create_sampler_gauss_2t_f

    !Moving average sampler
    subroutine mui_create_sampler_moving_average_2f_f(sampler,bbox_1,bbox_2) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_float), intent(in) :: bbox_1,bbox_2
    end subroutine mui_create_sampler_moving_average_2f_f

    subroutine mui_create_sampler_moving_average_2fx_f(sampler,bbox_1,bbox_2) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_float), intent(in) :: bbox_1,bbox_2
    end subroutine mui_create_sampler_moving_average_2fx_f

    subroutine mui_create_sampler_moving_average_2d_f(sampler,bbox_1,bbox_2) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double), intent(in) :: bbox_1,bbox_2
    end subroutine mui_create_sampler_moving_average_2d_f

    subroutine mui_create_sampler_moving_average_2dx_f(sampler,bbox_1,bbox_2) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double), intent(in) :: bbox_1,bbox_2
    end subroutine mui_create_sampler_moving_average_2dx_f

    subroutine mui_create_sampler_moving_average_2t_f(sampler,bbox_1,bbox_2) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double), intent(in) :: bbox_1,bbox_2
    end subroutine mui_create_sampler_moving_average_2t_f

    !Nearest neighbour sampler
    subroutine mui_create_sampler_nearest_neighbor_2f_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(out), target :: sampler(*)
    end subroutine mui_create_sampler_nearest_neighbor_2f_f

    subroutine mui_create_sampler_nearest_neighbor_2fx_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(out), target :: sampler(*)
    end subroutine mui_create_sampler_nearest_neighbor_2fx_f

    subroutine mui_create_sampler_nearest_neighbor_2d_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(out), target :: sampler(*)
    end subroutine mui_create_sampler_nearest_neighbor_2d_f

    subroutine mui_create_sampler_nearest_neighbor_2dx_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(out), target :: sampler(*)
    end subroutine mui_create_sampler_nearest_neighbor_2dx_f

    subroutine mui_create_sampler_nearest_neighbor_2t_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(out), target :: sampler(*)
    end subroutine mui_create_sampler_nearest_neighbor_2t_f

    !Pseudo-linear n^2 interpolation sampler
    subroutine mui_create_sampler_pseudo_n2_linear_2f_f(sampler,r) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_float), intent(in), target :: r
    end subroutine mui_create_sampler_pseudo_n2_linear_2f_f

    subroutine mui_create_sampler_pseudo_n2_linear_2fx_f(sampler,r) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_float), intent(in), target :: r
    end subroutine mui_create_sampler_pseudo_n2_linear_2fx_f

    subroutine mui_create_sampler_pseudo_n2_linear_2d_f(sampler,r) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double), intent(in), target :: r
    end subroutine mui_create_sampler_pseudo_n2_linear_2d_f

    subroutine mui_create_sampler_pseudo_n2_linear_2dx_f(sampler,r) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double), intent(in), target :: r
    end subroutine mui_create_sampler_pseudo_n2_linear_2dx_f

    subroutine mui_create_sampler_pseudo_n2_linear_2t_f(sampler,r) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double), intent(in), target :: r
    end subroutine mui_create_sampler_pseudo_n2_linear_2t_f

    !Pseudo-linear nearest neighbour interpolation sampler
    subroutine mui_create_sampler_pseudo_nearest_neighbor_2f_f(sampler,h) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_float), intent(in), target :: h
    end subroutine mui_create_sampler_pseudo_nearest_neighbor_2f_f

    subroutine mui_create_sampler_pseudo_nearest_neighbor_2fx_f(sampler,h) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_float), intent(in), target :: h
    end subroutine mui_create_sampler_pseudo_nearest_neighbor_2fx_f

    subroutine mui_create_sampler_pseudo_nearest_neighbor_2d_f(sampler,h) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double), intent(in), target :: h
    end subroutine mui_create_sampler_pseudo_nearest_neighbor_2d_f

    subroutine mui_create_sampler_pseudo_nearest_neighbor_2dx_f(sampler,h) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double), intent(in), target :: h
    end subroutine mui_create_sampler_pseudo_nearest_neighbor_2dx_f

    subroutine mui_create_sampler_pseudo_nearest_neighbor_2t_f(sampler,h) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double), intent(in), target :: h
    end subroutine mui_create_sampler_pseudo_nearest_neighbor_2t_f

    !Shepard interpolation with a quintic kernel sampler
    subroutine mui_create_sampler_shepard_quintic_2f_f(sampler,r) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_float), intent(in), target :: r
    end subroutine mui_create_sampler_shepard_quintic_2f_f

    subroutine mui_create_sampler_shepard_quintic_2fx_f(sampler,r) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_float), intent(in), target :: r
    end subroutine mui_create_sampler_shepard_quintic_2fx_f

    subroutine mui_create_sampler_shepard_quintic_2d_f(sampler,r) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double), intent(in), target :: r
    end subroutine mui_create_sampler_shepard_quintic_2d_f

    subroutine mui_create_sampler_shepard_quintic_2dx_f(sampler,r) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double), intent(in), target :: r
    end subroutine mui_create_sampler_shepard_quintic_2dx_f

    subroutine mui_create_sampler_shepard_quintic_2t_f(sampler,r) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double), intent(in), target :: r
    end subroutine mui_create_sampler_shepard_quintic_2t_f

    !SPH derived interpolation method with a quintic spline kernel sampler
    subroutine mui_create_sampler_sph_quintic_2f_f(sampler,r) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_float), intent(in), target :: r
    end subroutine mui_create_sampler_sph_quintic_2f_f

    subroutine mui_create_sampler_sph_quintic_2fx_f(sampler,r) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_float), intent(in), target :: r
    end subroutine mui_create_sampler_sph_quintic_2fx_f

    subroutine mui_create_sampler_sph_quintic_2d_f(sampler,r) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double), intent(in), target :: r
    end subroutine mui_create_sampler_sph_quintic_2d_f

    subroutine mui_create_sampler_sph_quintic_2dx_f(sampler,r) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double), intent(in), target :: r
    end subroutine mui_create_sampler_sph_quintic_2dx_f

    subroutine mui_create_sampler_sph_quintic_2t_f(sampler,r) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double), intent(in), target :: r
    end subroutine mui_create_sampler_sph_quintic_2t_f

    !Summation with a quintic kernel sampler
    subroutine mui_create_sampler_sum_quintic_2f_f(sampler,r) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_float), intent(in), target :: r
    end subroutine mui_create_sampler_sum_quintic_2f_f

    subroutine mui_create_sampler_sum_quintic_2fx_f(sampler,r) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_float), intent(in), target :: r
    end subroutine mui_create_sampler_sum_quintic_2fx_f

    subroutine mui_create_sampler_sum_quintic_2d_f(sampler,r) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double), intent(in), target :: r
    end subroutine mui_create_sampler_sum_quintic_2d_f

    subroutine mui_create_sampler_sum_quintic_2dx_f(sampler,r) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double), intent(in), target :: r
    end subroutine mui_create_sampler_sum_quintic_2dx_f

    subroutine mui_create_sampler_sum_quintic_2t_f(sampler,r) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double), intent(in), target :: r
    end subroutine mui_create_sampler_sum_quintic_2t_f

#ifdef USE_RBF
    !Radial Basis Function sampler
    subroutine mui_create_sampler_rbf_2f_f(sampler,r,points_1,points_2,points_count, &
               basis_func,conservative,polynomial,smoothFunc,readMatrix, &
               file_address,cutoff) bind(C)
      import :: c_ptr,c_int,c_float
      type(c_ptr), intent(out), target :: sampler(*)
      character(kind=c_char), intent(in) :: file_address(*)
      type(c_float), intent(in), dimension(points_count), target :: points_1,points_2
      integer(kind=c_int), intent(in), target :: points_count,basis_func,conservative,polynomial,smoothFunc,readMatrix,writeMatrix,cgSolve,pouSize
      real(kind=c_float), intent(in), target :: r,cutoff,cgSolveTol
    end subroutine mui_create_sampler_rbf_2f_f

    subroutine mui_create_sampler_rbf_2fx_f(sampler,r,points_1,points_2,points_count, &
               basis_func,conservative,polynomial,smoothFunc,readMatrix, &
               file_address,cutoff) bind(C)
      import :: c_ptr,c_int,c_float
      type(c_ptr), intent(out), target :: sampler(*)
      character(kind=c_char), intent(in) :: file_address(*)
      type(c_float), intent(in), dimension(points_count), target :: points_1,points_2
      integer(kind=c_int), intent(in), target :: points_count,basis_func,conservative,polynomial,smoothFunc,readMatrix,writeMatrix,cgSolve,pouSize
      real(kind=c_float), intent(in), target :: r,cutoff,cgSolveTol
    end subroutine mui_create_sampler_rbf_2fx_f

    subroutine mui_create_sampler_rbf_2d_f(sampler,r,points_1,points_2,points_count, &
               basis_func,conservative,polynomial,smoothFunc,readMatrix, &
               file_address,cutoff) bind(C)
      import :: c_ptr,c_int,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      character(kind=c_char), intent(in) :: file_address(*)
      type(c_double), intent(in), dimension(points_count), target :: points_1,points_2
      integer(kind=c_int), intent(in), target :: points_count,basis_func,conservative,polynomial,smoothFunc,readMatrix,writeMatrix,cgSolve,pouSize
      real(kind=c_double), intent(in), target :: r,cutoff,cgSolveTol
    end subroutine mui_create_sampler_rbf_2d_f

    subroutine mui_create_sampler_rbf_2dx_f(sampler,r,points_1,points_2,points_count, &
               basis_func,conservative,polynomial,smoothFunc,readMatrix, &
               file_address,cutoff) bind(C)
      import :: c_ptr,c_int,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      character(kind=c_char), intent(in) :: file_address(*)
      type(c_double), intent(in), dimension(points_count), target :: points_1,points_2
      integer(kind=c_int), intent(in), target :: points_count,basis_func,conservative,polynomial,smoothFunc,readMatrix,writeMatrix,cgSolve,pouSize
      real(kind=c_double), intent(in), target :: r,cutoff,cgSolveTol
    end subroutine mui_create_sampler_rbf_2dx_f

    subroutine mui_create_sampler_rbf_2t_f(sampler,r,points_1,points_2,points_count, &
               basis_func,conservative,polynomial,smoothFunc,readMatrix, &
               file_address,cutoff) bind(C)
      import :: c_ptr,c_int,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      character(kind=c_char), intent(in) :: file_address(*)
      type(c_double), intent(in), dimension(points_count), target :: points_1,points_2
      integer(kind=c_int), intent(in), target :: points_count,basis_func,conservative,polynomial,smoothFunc,readMatrix,writeMatrix,cgSolve,pouSize
      real(kind=c_double), intent(in), target :: r,cutoff,cgSolveTol
    end subroutine mui_create_sampler_rbf_2t_f
#endif

    !******************************************
    !* Destroy 2D spatial samplers            *
    !******************************************

    !Exact sampler
    subroutine mui_destroy_sampler_exact_2f_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_exact_2f_f

    subroutine mui_destroy_sampler_exact_2fx_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_exact_2fx_f

    subroutine mui_destroy_sampler_exact_2d_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_exact_2d_f

    subroutine mui_destroy_sampler_exact_2dx_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_exact_2dx_f

    subroutine mui_destroy_sampler_exact_2t_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_exact_2t_f

    !Gaussian sampler
    subroutine mui_destroy_sampler_gauss_2f_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_gauss_2f_f

    subroutine mui_destroy_sampler_gauss_2fx_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_gauss_2fx_f

    subroutine mui_destroy_sampler_gauss_2d_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_gauss_2d_f

    subroutine mui_destroy_sampler_gauss_2dx_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_gauss_2dx_f

    subroutine mui_destroy_sampler_gauss_2t_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_gauss_2t_f

    !Moving average sampler
    subroutine mui_destroy_sampler_moving_average_2f_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_moving_average_2f_f

    subroutine mui_destroy_sampler_moving_average_2fx_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_moving_average_2fx_f

    subroutine mui_destroy_sampler_moving_average_2d_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_moving_average_2d_f

    subroutine mui_destroy_sampler_moving_average_2dx_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_moving_average_2dx_f

    subroutine mui_destroy_sampler_moving_average_2t_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_moving_average_2t_f

    !Nearest neighbour sampler
    subroutine mui_destroy_sampler_nearest_neighbor_2f_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_nearest_neighbor_2f_f

    subroutine mui_destroy_sampler_nearest_neighbor_2fx_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_nearest_neighbor_2fx_f

    subroutine mui_destroy_sampler_nearest_neighbor_2d_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_nearest_neighbor_2d_f

    subroutine mui_destroy_sampler_nearest_neighbor_2dx_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_nearest_neighbor_2dx_f

    subroutine mui_destroy_sampler_nearest_neighbor_2t_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_nearest_neighbor_2t_f

    !Pseudo-linear n^2 interpolation sampler
    subroutine mui_destroy_sampler_pseudo_nearest2_linear_2f_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_pseudo_nearest2_linear_2f_f

    subroutine mui_destroy_sampler_pseudo_nearest2_linear_2fx_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_pseudo_nearest2_linear_2fx_f

    subroutine mui_destroy_sampler_pseudo_nearest2_linear_2d_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_pseudo_nearest2_linear_2d_f

    subroutine mui_destroy_sampler_pseudo_nearest2_linear_2dx_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_pseudo_nearest2_linear_2dx_f

    subroutine mui_destroy_sampler_pseudo_nearest2_linear_2t_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_pseudo_nearest2_linear_2t_f

    !Pseudo-linear nearest neighbour interpolation sampler
    subroutine mui_destroy_sampler_pseudo_nearest_neighbor_2f_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_pseudo_nearest_neighbor_2f_f

    subroutine mui_destroy_sampler_pseudo_nearest_neighbor_2fx_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_pseudo_nearest_neighbor_2fx_f

    subroutine mui_destroy_sampler_pseudo_nearest_neighbor_2d_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_pseudo_nearest_neighbor_2d_f

    subroutine mui_destroy_sampler_pseudo_nearest_neighbor_2dx_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_pseudo_nearest_neighbor_2dx_f

    subroutine mui_destroy_sampler_pseudo_nearest_neighbor_2t_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_pseudo_nearest_neighbor_2t_f

    !Shepard interpolation with a quintic kernel sampler
    subroutine mui_destroy_sampler_shepard_quintic_2f_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_shepard_quintic_2f_f

    subroutine mui_destroy_sampler_shepard_quintic_2fx_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_shepard_quintic_2fx_f

    subroutine mui_destroy_sampler_shepard_quintic_2d_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_shepard_quintic_2d_f

    subroutine mui_destroy_sampler_shepard_quintic_2dx_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_shepard_quintic_2dx_f

    subroutine mui_destroy_sampler_shepard_quintic_2t_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_shepard_quintic_2t_f

    !SPH derived interpolation method with a quintic spline kernel sampler
    subroutine mui_destroy_sampler_sph_quintic_2f_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_sph_quintic_2f_f

    subroutine mui_destroy_sampler_sph_quintic_2fx_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_sph_quintic_2fx_f

    subroutine mui_destroy_sampler_sph_quintic_2d_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_sph_quintic_2d_f

    subroutine mui_destroy_sampler_sph_quintic_2dx_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_sph_quintic_2dx_f

    subroutine mui_destroy_sampler_sph_quintic_2t_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_sph_quintic_2t_f

    !Summation with a quintic kernel sampler
    subroutine mui_destroy_sampler_sum_quintic_2f_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_sum_quintic_2f_f

    subroutine mui_destroy_sampler_sum_quintic_2fx_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_sum_quintic_2fx_f

    subroutine mui_destroy_sampler_sum_quintic_2d_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_sum_quintic_2d_f

    subroutine mui_destroy_sampler_sum_quintic_2dx_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_sum_quintic_2dx_f

    subroutine mui_destroy_sampler_sum_quintic_2t_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_sum_quintic_2t_f

#ifdef USE_RBF
    !Radial Basis Function sampler
    subroutine mui_destroy_sampler_rbf_2f_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_rbf_2f_f

    subroutine mui_destroy_sampler_rbf_2fx_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_rbf_2fx_f

    subroutine mui_destroy_sampler_rbf_2d_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_rbf_2d_f

    subroutine mui_destroy_sampler_rbf_2dx_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_rbf_2dx_f

    subroutine mui_destroy_sampler_rbf_2t_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_sampler_rbf_2t_f
#endif

    !******************************************
    !* Create temporal samplers               *
    !******************************************

    !Exact temporal sampler
    subroutine mui_create_chrono_sampler_exact_2f_f(sampler,tolerance) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_float),intent(in) :: tolerance
    end subroutine mui_create_chrono_sampler_exact_2f_f

    subroutine mui_create_chrono_sampler_exact_2fx_f(sampler,tolerance) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_float),intent(in) :: tolerance
    end subroutine mui_create_chrono_sampler_exact_2fx_f

    subroutine mui_create_chrono_sampler_exact_2d_f(sampler,tolerance) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double),intent(in) :: tolerance
    end subroutine mui_create_chrono_sampler_exact_2d_f

    subroutine mui_create_chrono_sampler_exact_2dx_f(sampler,tolerance) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double),intent(in) :: tolerance
    end subroutine mui_create_chrono_sampler_exact_2dx_f

    subroutine mui_create_chrono_sampler_exact_2t_f(sampler,tolerance) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double),intent(in) :: tolerance
    end subroutine mui_create_chrono_sampler_exact_2t_f

    !Gauss temporal sampler
    subroutine mui_create_chrono_sampler_gauss_2f_f(sampler,cutoff,sigma) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_float),intent(in) :: cutoff,sigma
    end subroutine mui_create_chrono_sampler_gauss_2f_f

    subroutine mui_create_chrono_sampler_gauss_2fx_f(sampler,cutoff,sigma) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_float),intent(in) :: cutoff,sigma
    end subroutine mui_create_chrono_sampler_gauss_2fx_f

    subroutine mui_create_chrono_sampler_gauss_2d_f(sampler,cutoff,sigma) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double),intent(in) :: cutoff,sigma
    end subroutine mui_create_chrono_sampler_gauss_2d_f

    subroutine mui_create_chrono_sampler_gauss_2dx_f(sampler,cutoff,sigma) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double),intent(in) :: cutoff,sigma
    end subroutine mui_create_chrono_sampler_gauss_2dx_f

    subroutine mui_create_chrono_sampler_gauss_2t_f(sampler,cutoff,sigma) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double),intent(in) :: cutoff,sigma
    end subroutine mui_create_chrono_sampler_gauss_2t_f

    !Mean temporal sampler
    subroutine mui_create_chrono_sampler_mean_2f_f(sampler,lower,upper) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_float),intent(in) :: lower,upper
    end subroutine mui_create_chrono_sampler_mean_2f_f

    subroutine mui_create_chrono_sampler_mean_2fx_f(sampler,lower,upper) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_float),intent(in) :: lower,upper
    end subroutine mui_create_chrono_sampler_mean_2fx_f

    subroutine mui_create_chrono_sampler_mean_2d_f(sampler,lower,upper) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double),intent(in) :: lower,upper
    end subroutine mui_create_chrono_sampler_mean_2d_f

    subroutine mui_create_chrono_sampler_mean_2dx_f(sampler,lower,upper) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double),intent(in) :: lower,upper
    end subroutine mui_create_chrono_sampler_mean_2dx_f

    subroutine mui_create_chrono_sampler_mean_2t_f(sampler,lower,upper) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double),intent(in) :: lower,upper
    end subroutine mui_create_chrono_sampler_mean_2t_f

    !Sum temporal sampler
    subroutine mui_create_chrono_sampler_sum_2f_f(sampler,lower,upper) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_float),intent(in) :: lower,upper
    end subroutine mui_create_chrono_sampler_sum_2f_f

    subroutine mui_create_chrono_sampler_sum_2fx_f(sampler,lower,upper) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_float),intent(in) :: lower,upper
    end subroutine mui_create_chrono_sampler_sum_2fx_f

    subroutine mui_create_chrono_sampler_sum_2d_f(sampler,lower,upper) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double),intent(in) :: lower,upper
    end subroutine mui_create_chrono_sampler_sum_2d_f

    subroutine mui_create_chrono_sampler_sum_2dx_f(sampler,lower,upper) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double),intent(in) :: lower,upper
    end subroutine mui_create_chrono_sampler_sum_2dx_f

    subroutine mui_create_chrono_sampler_sum_2t_f(sampler,lower,upper) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: sampler(*)
      real(kind=c_double),intent(in) :: lower,upper
    end subroutine mui_create_chrono_sampler_sum_2t_f

    !******************************************
    !* Destroy temporal samplers              *
    !******************************************

    !Exact temporal sampler
    subroutine mui_destroy_chrono_sampler_exact_2f_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_chrono_sampler_exact_2f_f

    subroutine mui_destroy_chrono_sampler_exact_2fx_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_chrono_sampler_exact_2fx_f

    subroutine mui_destroy_chrono_sampler_exact_2d_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_chrono_sampler_exact_2d_f

    subroutine mui_destroy_chrono_sampler_exact_2dx_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_chrono_sampler_exact_2dx_f

    subroutine mui_destroy_chrono_sampler_exact_2t_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_chrono_sampler_exact_2t_f

    !Gauss temporal sampler
    subroutine mui_destroy_chrono_sampler_gauss_2f_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_chrono_sampler_gauss_2f_f

    subroutine mui_destroy_chrono_sampler_gauss_2fx_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_chrono_sampler_gauss_2fx_f

    subroutine mui_destroy_chrono_sampler_gauss_2d_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_chrono_sampler_gauss_2d_f

    subroutine mui_destroy_chrono_sampler_gauss_2dx_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_chrono_sampler_gauss_2dx_f

    subroutine mui_destroy_chrono_sampler_gauss_2t_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_chrono_sampler_gauss_2t_f

    !Mean temporal sampler
    subroutine mui_destroy_chrono_sampler_mean_2f_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_chrono_sampler_mean_2f_f

    subroutine mui_destroy_chrono_sampler_mean_2fx_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_chrono_sampler_mean_2fx_f

    subroutine mui_destroy_chrono_sampler_mean_2d_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_chrono_sampler_mean_2d_f

    subroutine mui_destroy_chrono_sampler_mean_2dx_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_chrono_sampler_mean_2dx_f

    subroutine mui_destroy_chrono_sampler_mean_2t_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_chrono_sampler_mean_2t_f

    !Sum temporal sampler
    subroutine mui_destroy_chrono_sampler_sum_2f_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_chrono_sampler_sum_2f_f

    subroutine mui_destroy_chrono_sampler_sum_2fx_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_chrono_sampler_sum_2fx_f

    subroutine mui_destroy_chrono_sampler_sum_2d_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_chrono_sampler_sum_2d_f

    subroutine mui_destroy_chrono_sampler_sum_2dx_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_chrono_sampler_sum_2dx_f

    subroutine mui_destroy_chrono_sampler_sum_2t_f(sampler) bind(C)
      import :: c_ptr
      type(c_ptr), intent(in), value :: sampler
    end subroutine mui_destroy_chrono_sampler_sum_2t_f

    !******************************************
    !* MUI functions for data push            *
    !******************************************

    !Standard push functions
    subroutine mui_push_2f_f(uniface,attr,point_1,point_2,value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,value
    end subroutine mui_push_2f_f

    subroutine mui_push_2fx_f(uniface,attr,point_1,point_2,value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,value
    end subroutine mui_push_2fx_f

    subroutine mui_push_2d_f(uniface,attr,point_1,point_2,value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,value
    end subroutine mui_push_2d_f

    subroutine mui_push_2dx_f(uniface,attr,point_1,point_2,value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,value
    end subroutine mui_push_2dx_f

    subroutine mui_push_2t_f(uniface,attr,point_1,point_2,value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,value
    end subroutine mui_push_2t_f

    !Single parameter push functions
    subroutine mui_push_2f_param_f(uniface,attr,value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: value
    end subroutine mui_push_2f_param_f

    subroutine mui_push_2fx_param_f(uniface,attr,value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: value
    end subroutine mui_push_2fx_param_f

    subroutine mui_push_2d_param_f(uniface,attr,value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: value
    end subroutine mui_push_2d_param_f

    subroutine mui_push_2dx_param_f(uniface,attr,value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: value
    end subroutine mui_push_2dx_param_f

    subroutine mui_push_2t_param_f(uniface,attr,value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: value
    end subroutine mui_push_2t_param_f

    !******************************************
    !* MUI functions for data commit          *
    !******************************************

    !Commit using one time value
    subroutine mui_commit_2f_f(uniface,t) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in) :: t
    end subroutine mui_commit_2f_f

    subroutine mui_commit_2fx_f(uniface,t) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in) :: t
    end subroutine mui_commit_2fx_f

    subroutine mui_commit_2d_f(uniface,t) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in) :: t
    end subroutine mui_commit_2d_f

    subroutine mui_commit_2dx_f(uniface,t) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in) :: t
    end subroutine mui_commit_2dx_f

    subroutine mui_commit_2t_f(uniface,t) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in) :: t
    end subroutine mui_commit_2t_f

    !Commit using two time values
    subroutine mui_commit_2f_pair_f(uniface,t1,t2) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in) :: t1,t2
    end subroutine mui_commit_2f_pair_f

    subroutine mui_commit_2fx_pair_f(uniface,t1,t2) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in) :: t1,t2
    end subroutine mui_commit_2fx_pair_f

    subroutine mui_commit_2d_pair_f(uniface,t1,t2) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in) :: t1,t2
    end subroutine mui_commit_2d_pair_f

    subroutine mui_commit_2dx_pair_f(uniface,t1,t2) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in) :: t1,t2
    end subroutine mui_commit_2dx_pair_f

    subroutine mui_commit_2t_pair_f(uniface,t1,t2) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in) :: t1,t2
    end subroutine mui_commit_2t_pair_f

    !******************************************
    !* MUI functions for data forecast        *
    !******************************************

    !Forecast using one time value
    subroutine mui_forecast_2f_f(uniface,t) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in) :: t
    end subroutine mui_forecast_2f_f

    subroutine mui_forecast_2fx_f(uniface,t) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in) :: t
    end subroutine mui_forecast_2fx_f

    subroutine mui_forecast_2d_f(uniface,t) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in) :: t
    end subroutine mui_forecast_2d_f

    subroutine mui_forecast_2dx_f(uniface,t) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in) :: t
    end subroutine mui_forecast_2dx_f

    subroutine mui_forecast_2t_f(uniface,t) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in) :: t
    end subroutine mui_forecast_2t_f

    !Commit using two time values
    subroutine mui_forecast_2f_pair_f(uniface,t1,t2) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in) :: t1,t2
    end subroutine mui_forecast_2f_pair_f

    subroutine mui_forecast_2fx_pair_f(uniface,t1,t2) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in) :: t1,t2
    end subroutine mui_forecast_2fx_pair_f

    subroutine mui_forecast_2d_pair_f(uniface,t1,t2) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in) :: t1,t2
    end subroutine mui_forecast_2d_pair_f

    subroutine mui_forecast_2dx_pair_f(uniface,t1,t2) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in) :: t1,t2
    end subroutine mui_forecast_2dx_pair_f

    subroutine mui_forecast_2t_pair_f(uniface,t1,t2) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in) :: t1,t2
    end subroutine mui_forecast_2t_pair_f

    !*********************************************************
    !* MUI functions for 2D data fetch using one time value  *
    !*********************************************************

    !Spatial sampler: exact; temporal sampler: exact
    subroutine mui_fetch_exact_exact_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_exact_exact_2f_f

    subroutine mui_fetch_exact_exact_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_exact_exact_2fx_f

    subroutine mui_fetch_exact_exact_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_exact_exact_2d_f

    subroutine mui_fetch_exact_exact_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_exact_exact_2dx_f

    subroutine mui_fetch_exact_exact_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_exact_exact_2t_f

    !Spatial sampler: exact; temporal sampler: gauss
    subroutine mui_fetch_exact_gauss_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_exact_gauss_2f_f

    subroutine mui_fetch_exact_gauss_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_exact_gauss_2fx_f

    subroutine mui_fetch_exact_gauss_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_exact_gauss_2d_f

    subroutine mui_fetch_exact_gauss_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_exact_gauss_2dx_f

    subroutine mui_fetch_exact_gauss_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_exact_gauss_2t_f

    !Spatial sampler: exact; temporal sampler: mean
    subroutine mui_fetch_exact_mean_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_exact_mean_2f_f

    subroutine mui_fetch_exact_mean_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_exact_mean_2fx_f

    subroutine mui_fetch_exact_mean_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_exact_mean_2d_f

    subroutine mui_fetch_exact_mean_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_exact_mean_2dx_f

    subroutine mui_fetch_exact_mean_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_exact_mean_2t_f

    !Spatial sampler: exact; temporal sampler: sum
    subroutine mui_fetch_exact_sum_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_exact_sum_2f_f

    subroutine mui_fetch_exact_sum_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_exact_sum_2fx_f

    subroutine mui_fetch_exact_sum_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_exact_sum_2d_f

    subroutine mui_fetch_exact_sum_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_exact_sum_2dx_f

    subroutine mui_fetch_exact_sum_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_exact_sum_2t_f

    !Spatial sampler: gauss; temporal sampler: exact
    subroutine mui_fetch_gauss_exact_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_gauss_exact_2f_f

    subroutine mui_fetch_gauss_exact_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_gauss_exact_2fx_f

    subroutine mui_fetch_gauss_exact_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_gauss_exact_2d_f

    subroutine mui_fetch_gauss_exact_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_gauss_exact_2dx_f

    subroutine mui_fetch_gauss_exact_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_gauss_exact_2t_f

    !Spatial sampler: gauss; temporal sampler: gauss
    subroutine mui_fetch_gauss_gauss_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_gauss_gauss_2f_f

    subroutine mui_fetch_gauss_gauss_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_gauss_gauss_2fx_f

    subroutine mui_fetch_gauss_gauss_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_gauss_gauss_2d_f

    subroutine mui_fetch_gauss_gauss_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_gauss_gauss_2dx_f

    subroutine mui_fetch_gauss_gauss_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_gauss_gauss_2t_f

    !Spatial sampler: gauss; temporal sampler: mean
    subroutine mui_fetch_gauss_mean_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_gauss_mean_2f_f

    subroutine mui_fetch_gauss_mean_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_gauss_mean_2fx_f

    subroutine mui_fetch_gauss_mean_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_gauss_mean_2d_f

    subroutine mui_fetch_gauss_mean_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_gauss_mean_2dx_f

    subroutine mui_fetch_gauss_mean_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_gauss_mean_2t_f

    !Spatial sampler: gauss; temporal sampler: sum
    subroutine mui_fetch_gauss_sum_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_gauss_sum_2f_f

    subroutine mui_fetch_gauss_sum_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_gauss_sum_2fx_f

    subroutine mui_fetch_gauss_sum_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_gauss_sum_2d_f

    subroutine mui_fetch_gauss_sum_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_gauss_sum_2dx_f

    subroutine mui_fetch_gauss_sum_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_gauss_sum_2t_f

    !Spatial sampler: moving average; temporal sampler: exact
    subroutine mui_fetch_moving_average_exact_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_exact_2f_f

    subroutine mui_fetch_moving_average_exact_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_exact_2fx_f

    subroutine mui_fetch_moving_average_exact_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_exact_2d_f

    subroutine mui_fetch_moving_average_exact_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_exact_2dx_f

    subroutine mui_fetch_moving_average_exact_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_exact_2t_f

    !Spatial sampler: moving average; temporal sampler: gauss
    subroutine mui_fetch_moving_average_gauss_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_gauss_2f_f

    subroutine mui_fetch_moving_average_gauss_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_gauss_2fx_f

    subroutine mui_fetch_moving_average_gauss_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_gauss_2d_f

    subroutine mui_fetch_moving_average_gauss_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_gauss_2dx_f

    subroutine mui_fetch_moving_average_gauss_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_gauss_2t_f

    !Spatial sampler: moving average; temporal sampler: mean
    subroutine mui_fetch_moving_average_mean_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_mean_2f_f

    subroutine mui_fetch_moving_average_mean_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_mean_2fx_f

    subroutine mui_fetch_moving_average_mean_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_mean_2d_f

    subroutine mui_fetch_moving_average_mean_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_mean_2dx_f

    subroutine mui_fetch_moving_average_mean_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_mean_2t_f

    !Spatial sampler: moving average; temporal sampler: sum
    subroutine mui_fetch_moving_average_sum_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_sum_2f_f

    subroutine mui_fetch_moving_average_sum_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_sum_2fx_f

    subroutine mui_fetch_moving_average_sum_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_sum_2d_f

    subroutine mui_fetch_moving_average_sum_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_sum_2dx_f

    subroutine mui_fetch_moving_average_sum_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_sum_2t_f

    !Spatial sampler: nearest neighbor; temporal sampler: exact
    subroutine mui_fetch_nearest_neighbor_exact_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_exact_2f_f

    subroutine mui_fetch_nearest_neighbor_exact_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_exact_2fx_f

    subroutine mui_fetch_nearest_neighbor_exact_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_exact_2d_f

    subroutine mui_fetch_nearest_neighbor_exact_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_exact_2dx_f

    subroutine mui_fetch_nearest_neighbor_exact_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_exact_2t_f

    !Spatial sampler: nearest neighbor; temporal sampler: gauss
    subroutine mui_fetch_nearest_neighbor_gauss_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_gauss_2f_f

    subroutine mui_fetch_nearest_neighbor_gauss_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_gauss_2fx_f

    subroutine mui_fetch_nearest_neighbor_gauss_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_gauss_2d_f

    subroutine mui_fetch_nearest_neighbor_gauss_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_gauss_2dx_f

    subroutine mui_fetch_nearest_neighbor_gauss_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_gauss_2t_f

    !Spatial sampler: nearest neighbor; temporal sampler: mean
    subroutine mui_fetch_nearest_neighbor_mean_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_mean_2f_f

    subroutine mui_fetch_nearest_neighbor_mean_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_mean_2fx_f

    subroutine mui_fetch_nearest_neighbor_mean_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_mean_2d_f

    subroutine mui_fetch_nearest_neighbor_mean_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_mean_2dx_f

    subroutine mui_fetch_nearest_neighbor_mean_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_mean_2t_f

    !Spatial sampler: nearest neighbor; temporal sampler: sum
    subroutine mui_fetch_nearest_neighbor_sum_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_sum_2f_f

    subroutine mui_fetch_nearest_neighbor_sum_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_sum_2fx_f

    subroutine mui_fetch_nearest_neighbor_sum_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_sum_2d_f

    subroutine mui_fetch_nearest_neighbor_sum_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_sum_2dx_f

    subroutine mui_fetch_nearest_neighbor_sum_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_sum_2t_f

    !Spatial sampler: pseudo nearest neighbor; temporal sampler: exact
    subroutine mui_fetch_pseudo_nearest_neighbor_exact_2f_f(uniface,attr,point_1,point_2,t,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_exact_2f_f

    subroutine mui_fetch_pseudo_nearest_neighbor_exact_2fx_f(uniface,attr,point_1,point_2,t,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_exact_2fx_f

    subroutine mui_fetch_pseudo_nearest_neighbor_exact_2d_f(uniface,attr,point_1,point_2,t,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_exact_2d_f

    subroutine mui_fetch_pseudo_nearest_neighbor_exact_2dx_f(uniface,attr,point_1,point_2,t,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_exact_2dx_f

    subroutine mui_fetch_pseudo_nearest_neighbor_exact_2t_f(uniface,attr,point_1,point_2,t,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_exact_2t_f

    !Spatial sampler: pseudo nearest neighbor; temporal sampler: gauss
    subroutine mui_fetch_pseudo_nearest_neighbor_gauss_2f_f(uniface,attr,point_1,point_2,t,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_gauss_2f_f

    subroutine mui_fetch_pseudo_nearest_neighbor_gauss_2fx_f(uniface,attr,point_1,point_2,t,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_gauss_2fx_f

    subroutine mui_fetch_pseudo_nearest_neighbor_gauss_2d_f(uniface,attr,point_1,point_2,t,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_gauss_2d_f

    subroutine mui_fetch_pseudo_nearest_neighbor_gauss_2dx_f(uniface,attr,point_1,point_2,t,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_gauss_2dx_f

    subroutine mui_fetch_pseudo_nearest_neighbor_gauss_2t_f(uniface,attr,point_1,point_2,t,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_gauss_2t_f

    !Spatial sampler: pseudo nearest neighbor; temporal sampler: mean
    subroutine mui_fetch_pseudo_nearest_neighbor_mean_2f_f(uniface,attr,point_1,point_2,t,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_mean_2f_f

    subroutine mui_fetch_pseudo_nearest_neighbor_mean_2fx_f(uniface,attr,point_1,point_2,t,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_mean_2fx_f

    subroutine mui_fetch_pseudo_nearest_neighbor_mean_2d_f(uniface,attr,point_1,point_2,t,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_mean_2d_f

    subroutine mui_fetch_pseudo_nearest_neighbor_mean_2dx_f(uniface,attr,point_1,point_2,t,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_mean_2dx_f

    subroutine mui_fetch_pseudo_nearest_neighbor_mean_2t_f(uniface,attr,point_1,point_2,t,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_mean_2t_f

    !Spatial sampler: pseudo nearest neighbor; temporal sampler: sum
    subroutine mui_fetch_pseudo_nearest_neighbor_sum_2f_f(uniface,attr,point_1,point_2,t,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_sum_2f_f

    subroutine mui_fetch_pseudo_nearest_neighbor_sum_2fx_f(uniface,attr,point_1,point_2,t,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_sum_2fx_f

    subroutine mui_fetch_pseudo_nearest_neighbor_sum_2d_f(uniface,attr,point_1,point_2,t,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_sum_2d_f

    subroutine mui_fetch_pseudo_nearest_neighbor_sum_2dx_f(uniface,attr,point_1,point_2,t,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_sum_2dx_f

    subroutine mui_fetch_pseudo_nearest_neighbor_sum_2t_f(uniface,attr,point_1,point_2,t,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_sum_2t_f

    !Spatial sampler: shepard quintic; temporal sampler: exact
    subroutine mui_fetch_shepard_quintic_exact_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_exact_2f_f

    subroutine mui_fetch_shepard_quintic_exact_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_exact_2fx_f

    subroutine mui_fetch_shepard_quintic_exact_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_exact_2d_f

    subroutine mui_fetch_shepard_quintic_exact_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_exact_2dx_f

    subroutine mui_fetch_shepard_quintic_exact_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_exact_2t_f

    !Spatial sampler: shepard quintic; temporal sampler: gauss
    subroutine mui_fetch_shepard_quintic_gauss_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_gauss_2f_f

    subroutine mui_fetch_shepard_quintic_gauss_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_gauss_2fx_f

    subroutine mui_fetch_shepard_quintic_gauss_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_gauss_2d_f

    subroutine mui_fetch_shepard_quintic_gauss_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_gauss_2dx_f

    subroutine mui_fetch_shepard_quintic_gauss_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_gauss_2t_f

    !Spatial sampler: shepard quintic; temporal sampler: mean
    subroutine mui_fetch_shepard_quintic_mean_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_mean_2f_f

    subroutine mui_fetch_shepard_quintic_mean_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_mean_2fx_f

    subroutine mui_fetch_shepard_quintic_mean_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_mean_2d_f

    subroutine mui_fetch_shepard_quintic_mean_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_mean_2dx_f

    subroutine mui_fetch_shepard_quintic_mean_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_mean_2t_f

    !Spatial sampler: shepard quintic; temporal sampler: sum
    subroutine mui_fetch_shepard_quintic_sum_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_sum_2f_f

    subroutine mui_fetch_shepard_quintic_sum_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_sum_2fx_f

    subroutine mui_fetch_shepard_quintic_sum_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_sum_2d_f

    subroutine mui_fetch_shepard_quintic_sum_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_sum_2dx_f

    subroutine mui_fetch_shepard_quintic_sum_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_sum_2t_f

    !Spatial sampler: sph-derived quintic; temporal sampler: exact
    subroutine mui_fetch_sph_quintic_exact_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_exact_2f_f

    subroutine mui_fetch_sph_quintic_exact_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_exact_2fx_f

    subroutine mui_fetch_sph_quintic_exact_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_exact_2d_f

    subroutine mui_fetch_sph_quintic_exact_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_exact_2dx_f

    subroutine mui_fetch_sph_quintic_exact_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_exact_2t_f

    !Spatial sampler: sph-derived quintic; temporal sampler: gauss
    subroutine mui_fetch_sph_quintic_gauss_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_gauss_2f_f

    subroutine mui_fetch_sph_quintic_gauss_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_gauss_2fx_f

    subroutine mui_fetch_sph_quintic_gauss_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_gauss_2d_f

    subroutine mui_fetch_sph_quintic_gauss_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_gauss_2dx_f

    subroutine mui_fetch_sph_quintic_gauss_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_gauss_2t_f

    !Spatial sampler: sph-derived quintic; temporal sampler: mean
    subroutine mui_fetch_sph_quintic_mean_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_mean_2f_f

    subroutine mui_fetch_sph_quintic_mean_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_mean_2fx_f

    subroutine mui_fetch_sph_quintic_mean_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_mean_2d_f

    subroutine mui_fetch_sph_quintic_mean_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_mean_2dx_f

    subroutine mui_fetch_sph_quintic_mean_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_mean_2t_f

    !Spatial sampler: sph-derived quintic; temporal sampler: sum
    subroutine mui_fetch_sph_quintic_sum_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_sum_2f_f

    subroutine mui_fetch_sph_quintic_sum_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_sum_2fx_f

    subroutine mui_fetch_sph_quintic_sum_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_sum_2d_f

    subroutine mui_fetch_sph_quintic_sum_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_sum_2dx_f

    subroutine mui_fetch_sph_quintic_sum_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_sum_2t_f

    !Spatial sampler: summation quintic; temporal sampler: exact
    subroutine mui_fetch_sum_quintic_exact_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_exact_2f_f

    subroutine mui_fetch_sum_quintic_exact_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_exact_2fx_f

    subroutine mui_fetch_sum_quintic_exact_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_exact_2d_f

    subroutine mui_fetch_sum_quintic_exact_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_exact_2dx_f

    subroutine mui_fetch_sum_quintic_exact_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_exact_2t_f

    !Spatial sampler: summation quintic; temporal sampler: gauss
    subroutine mui_fetch_sum_quintic_gauss_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_gauss_2f_f

    subroutine mui_fetch_sum_quintic_gauss_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_gauss_2fx_f

    subroutine mui_fetch_sum_quintic_gauss_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_gauss_2d_f

    subroutine mui_fetch_sum_quintic_gauss_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_gauss_2dx_f

    subroutine mui_fetch_sum_quintic_gauss_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_gauss_2t_f

    !Spatial sampler: summation quintic; temporal sampler: mean
    subroutine mui_fetch_sum_quintic_mean_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_mean_2f_f

    subroutine mui_fetch_sum_quintic_mean_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_mean_2fx_f

    subroutine mui_fetch_sum_quintic_mean_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_mean_2d_f

    subroutine mui_fetch_sum_quintic_mean_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_mean_2dx_f

    subroutine mui_fetch_sum_quintic_mean_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_mean_2t_f

    !Spatial sampler: summation quintic; temporal sampler: sum
    subroutine mui_fetch_sum_quintic_sum_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_sum_2f_f

    subroutine mui_fetch_sum_quintic_sum_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_sum_2fx_f

    subroutine mui_fetch_sum_quintic_sum_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_sum_2d_f

    subroutine mui_fetch_sum_quintic_sum_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_sum_2dx_f

    subroutine mui_fetch_sum_quintic_sum_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_sum_2t_f

#ifdef USE_RBF
    !Spatial sampler: radial basis function; temporal sampler: exact
    subroutine mui_fetch_rbf_exact_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_rbf_exact_2f_f

    subroutine mui_fetch_rbf_exact_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_rbf_exact_2fx_f

    subroutine mui_fetch_rbf_exact_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_rbf_exact_2d_f

    subroutine mui_fetch_rbf_exact_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_rbf_exact_2dx_f

    subroutine mui_fetch_rbf_exact_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_rbf_exact_2t_f

    !Spatial sampler: radial basis function; temporal sampler: gauss
    subroutine mui_fetch_rbf_gauss_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_rbf_gauss_2f_f

    subroutine mui_fetch_rbf_gauss_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_rbf_gauss_2fx_f

    subroutine mui_fetch_rbf_gauss_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_rbf_gauss_2d_f

    subroutine mui_fetch_rbf_gauss_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_rbf_gauss_2dx_f

    subroutine mui_fetch_rbf_gauss_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_rbf_gauss_2t_f

    !Spatial sampler: radial basis function; temporal sampler: mean
    subroutine mui_fetch_rbf_mean_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_rbf_mean_2f_f

    subroutine mui_fetch_rbf_mean_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_rbf_mean_2fx_f

    subroutine mui_fetch_rbf_mean_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_rbf_mean_2d_f

    subroutine mui_fetch_rbf_mean_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_rbf_mean_2dx_f

    subroutine mui_fetch_rbf_mean_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_rbf_mean_2t_f

    !Spatial sampler: radial basis function; temporal sampler: sum
    subroutine mui_fetch_rbf_sum_2f_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_rbf_sum_2f_f

    subroutine mui_fetch_rbf_sum_2fx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_rbf_sum_2fx_f

    subroutine mui_fetch_rbf_sum_2d_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_rbf_sum_2d_f

    subroutine mui_fetch_rbf_sum_2dx_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_rbf_sum_2dx_f

    subroutine mui_fetch_rbf_sum_2t_f(uniface,attr,point_1,point_2,t,spatial_sampler,temporal_sampler, &
      return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_rbf_sum_2t_f
#endif

    !*********************************************************
    !* MUI functions for 2D data fetch using two time values *
    !*********************************************************

    !Spatial sampler: exact; temporal sampler: exact
    subroutine mui_fetch_exact_exact_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_exact_exact_2f_pair_f

    subroutine mui_fetch_exact_exact_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_exact_exact_2fx_pair_f

    subroutine mui_fetch_exact_exact_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_exact_exact_2d_pair_f

    subroutine mui_fetch_exact_exact_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_exact_exact_2dx_pair_f

    subroutine mui_fetch_exact_exact_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_exact_exact_2t_pair_f

    !Spatial sampler: exact; temporal sampler: gauss
    subroutine mui_fetch_exact_gauss_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_exact_gauss_2f_pair_f

    subroutine mui_fetch_exact_gauss_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_exact_gauss_2fx_pair_f

    subroutine mui_fetch_exact_gauss_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_exact_gauss_2d_pair_f

    subroutine mui_fetch_exact_gauss_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_exact_gauss_2dx_pair_f

    subroutine mui_fetch_exact_gauss_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_exact_gauss_2t_pair_f

    !Spatial sampler: exact; temporal sampler: mean
    subroutine mui_fetch_exact_mean_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_exact_mean_2f_pair_f

    subroutine mui_fetch_exact_mean_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_exact_mean_2fx_pair_f

    subroutine mui_fetch_exact_mean_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_exact_mean_2d_pair_f

    subroutine mui_fetch_exact_mean_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_exact_mean_2dx_pair_f

    subroutine mui_fetch_exact_mean_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_exact_mean_2t_pair_f

    !Spatial sampler: exact; temporal sampler: sum
    subroutine mui_fetch_exact_sum_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_exact_sum_2f_pair_f

    subroutine mui_fetch_exact_sum_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_exact_sum_2fx_pair_f

    subroutine mui_fetch_exact_sum_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_exact_sum_2d_pair_f

    subroutine mui_fetch_exact_sum_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_exact_sum_2dx_pair_f

    subroutine mui_fetch_exact_sum_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_exact_sum_2t_pair_f

    !Spatial sampler: gauss; temporal sampler: exact
    subroutine mui_fetch_gauss_exact_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_gauss_exact_2f_pair_f

    subroutine mui_fetch_gauss_exact_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_gauss_exact_2fx_pair_f

    subroutine mui_fetch_gauss_exact_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_gauss_exact_2d_pair_f

    subroutine mui_fetch_gauss_exact_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_gauss_exact_2dx_pair_f

    subroutine mui_fetch_gauss_exact_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_gauss_exact_2t_pair_f

    !Spatial sampler: gauss; temporal sampler: gauss
    subroutine mui_fetch_gauss_gauss_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_gauss_gauss_2f_pair_f

    subroutine mui_fetch_gauss_gauss_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_gauss_gauss_2fx_pair_f

    subroutine mui_fetch_gauss_gauss_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_gauss_gauss_2d_pair_f

    subroutine mui_fetch_gauss_gauss_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_gauss_gauss_2dx_pair_f

    subroutine mui_fetch_gauss_gauss_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_gauss_gauss_2t_pair_f

    !Spatial sampler: gauss; temporal sampler: mean
    subroutine mui_fetch_gauss_mean_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_gauss_mean_2f_pair_f

    subroutine mui_fetch_gauss_mean_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_gauss_mean_2fx_pair_f

    subroutine mui_fetch_gauss_mean_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_gauss_mean_2d_pair_f

    subroutine mui_fetch_gauss_mean_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_gauss_mean_2dx_pair_f

    subroutine mui_fetch_gauss_mean_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_gauss_mean_2t_pair_f

    !Spatial sampler: gauss; temporal sampler: sum
    subroutine mui_fetch_gauss_sum_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_gauss_sum_2f_pair_f

    subroutine mui_fetch_gauss_sum_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_gauss_sum_2fx_pair_f

    subroutine mui_fetch_gauss_sum_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_gauss_sum_2d_pair_f

    subroutine mui_fetch_gauss_sum_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_gauss_sum_2dx_pair_f

    subroutine mui_fetch_gauss_sum_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_gauss_sum_2t_pair_f

    !Spatial sampler: moving average; temporal sampler: exact
    subroutine mui_fetch_moving_average_exact_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_exact_2f_pair_f

    subroutine mui_fetch_moving_average_exact_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_exact_2fx_pair_f

    subroutine mui_fetch_moving_average_exact_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_exact_2d_pair_f

    subroutine mui_fetch_moving_average_exact_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_exact_2dx_pair_f

    subroutine mui_fetch_moving_average_exact_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_exact_2t_pair_f

    !Spatial sampler: moving average; temporal sampler: gauss
    subroutine mui_fetch_moving_average_gauss_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_gauss_2f_pair_f

    subroutine mui_fetch_moving_average_gauss_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_gauss_2fx_pair_f

    subroutine mui_fetch_moving_average_gauss_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_gauss_2d_pair_f

    subroutine mui_fetch_moving_average_gauss_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_gauss_2dx_pair_f

    subroutine mui_fetch_moving_average_gauss_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_gauss_2t_pair_f

    !Spatial sampler: moving average; temporal sampler: mean
    subroutine mui_fetch_moving_average_mean_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_mean_2f_pair_f

    subroutine mui_fetch_moving_average_mean_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_mean_2fx_pair_f

    subroutine mui_fetch_moving_average_mean_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_mean_2d_pair_f

    subroutine mui_fetch_moving_average_mean_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_mean_2dx_pair_f

    subroutine mui_fetch_moving_average_mean_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_mean_2t_pair_f

    !Spatial sampler: moving average; temporal sampler: sum
    subroutine mui_fetch_moving_average_sum_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_sum_2f_pair_f

    subroutine mui_fetch_moving_average_sum_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_sum_2fx_pair_f

    subroutine mui_fetch_moving_average_sum_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_sum_2d_pair_f

    subroutine mui_fetch_moving_average_sum_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_sum_2dx_pair_f

    subroutine mui_fetch_moving_average_sum_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_moving_average_sum_2t_pair_f

    !Spatial sampler: nearest neighbor; temporal sampler: exact
    subroutine mui_fetch_nearest_neighbor_exact_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_exact_2f_pair_f

    subroutine mui_fetch_nearest_neighbor_exact_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_exact_2fx_pair_f

    subroutine mui_fetch_nearest_neighbor_exact_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_exact_2d_pair_f

    subroutine mui_fetch_nearest_neighbor_exact_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_exact_2dx_pair_f

    subroutine mui_fetch_nearest_neighbor_exact_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_exact_2t_pair_f

    !Spatial sampler: nearest neighbor; temporal sampler: gauss
    subroutine mui_fetch_nearest_neighbor_gauss_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_gauss_2f_pair_f

    subroutine mui_fetch_nearest_neighbor_gauss_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_gauss_2fx_pair_f

    subroutine mui_fetch_nearest_neighbor_gauss_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_gauss_2d_pair_f

    subroutine mui_fetch_nearest_neighbor_gauss_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_gauss_2dx_pair_f

    subroutine mui_fetch_nearest_neighbor_gauss_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_gauss_2t_pair_f

    !Spatial sampler: nearest neighbor; temporal sampler: mean
    subroutine mui_fetch_nearest_neighbor_mean_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_mean_2f_pair_f

    subroutine mui_fetch_nearest_neighbor_mean_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_mean_2fx_pair_f

    subroutine mui_fetch_nearest_neighbor_mean_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_mean_2d_pair_f

    subroutine mui_fetch_nearest_neighbor_mean_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_mean_2dx_pair_f

    subroutine mui_fetch_nearest_neighbor_mean_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_mean_2t_pair_f

    !Spatial sampler: nearest neighbor; temporal sampler: sum
    subroutine mui_fetch_nearest_neighbor_sum_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_sum_2f_pair_f

    subroutine mui_fetch_nearest_neighbor_sum_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_sum_2fx_pair_f

    subroutine mui_fetch_nearest_neighbor_sum_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_sum_2d_pair_f

    subroutine mui_fetch_nearest_neighbor_sum_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_sum_2dx_pair_f

    subroutine mui_fetch_nearest_neighbor_sum_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_nearest_neighbor_sum_2t_pair_f

    !Spatial sampler: pseudo nearest neighbor; temporal sampler: exact
    subroutine mui_fetch_pseudo_nearest_neighbor_exact_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_exact_2f_pair_f

    subroutine mui_fetch_pseudo_nearest_neighbor_exact_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_exact_2fx_pair_f

    subroutine mui_fetch_pseudo_nearest_neighbor_exact_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_exact_2d_pair_f

    subroutine mui_fetch_pseudo_nearest_neighbor_exact_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_exact_2dx_pair_f

    subroutine mui_fetch_pseudo_nearest_neighbor_exact_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_exact_2t_pair_f

    !Spatial sampler: pseudo nearest neighbor; temporal sampler: gauss
    subroutine mui_fetch_pseudo_nearest_neighbor_gauss_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_gauss_2f_pair_f

    subroutine mui_fetch_pseudo_nearest_neighbor_gauss_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_gauss_2fx_pair_f

    subroutine mui_fetch_pseudo_nearest_neighbor_gauss_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_gauss_2d_pair_f

    subroutine mui_fetch_pseudo_nearest_neighbor_gauss_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_gauss_2dx_pair_f

    subroutine mui_fetch_pseudo_nearest_neighbor_gauss_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_gauss_2t_pair_f

    !Spatial sampler: pseudo nearest neighbor; temporal sampler: mean
    subroutine mui_fetch_pseudo_nearest_neighbor_mean_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_mean_2f_pair_f

    subroutine mui_fetch_pseudo_nearest_neighbor_mean_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_mean_2fx_pair_f

    subroutine mui_fetch_pseudo_nearest_neighbor_mean_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_mean_2d_pair_f

    subroutine mui_fetch_pseudo_nearest_neighbor_mean_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_mean_2dx_pair_f

    subroutine mui_fetch_pseudo_nearest_neighbor_mean_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_mean_2t_pair_f

    !Spatial sampler: pseudo nearest neighbor; temporal sampler: sum
    subroutine mui_fetch_pseudo_nearest_neighbor_sum_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_sum_2f_pair_f

    subroutine mui_fetch_pseudo_nearest_neighbor_sum_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_sum_2fx_pair_f

    subroutine mui_fetch_pseudo_nearest_neighbor_sum_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_sum_2d_pair_f

    subroutine mui_fetch_pseudo_nearest_neighbor_sum_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_sum_2dx_pair_f

    subroutine mui_fetch_pseudo_nearest_neighbor_sum_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_pseudo_nearest_neighbor_sum_2t_pair_f

    !Spatial sampler: shepard quintic; temporal sampler: exact
    subroutine mui_fetch_shepard_quintic_exact_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_exact_2f_pair_f

    subroutine mui_fetch_shepard_quintic_exact_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_exact_2fx_pair_f

    subroutine mui_fetch_shepard_quintic_exact_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_exact_2d_pair_f

    subroutine mui_fetch_shepard_quintic_exact_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_exact_2dx_pair_f

    subroutine mui_fetch_shepard_quintic_exact_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_exact_2t_pair_f

    !Spatial sampler: shepard quintic; temporal sampler: gauss
    subroutine mui_fetch_shepard_quintic_gauss_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_gauss_2f_pair_f

    subroutine mui_fetch_shepard_quintic_gauss_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_gauss_2fx_pair_f

    subroutine mui_fetch_shepard_quintic_gauss_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_gauss_2d_pair_f

    subroutine mui_fetch_shepard_quintic_gauss_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_gauss_2dx_pair_f

    subroutine mui_fetch_shepard_quintic_gauss_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_gauss_2t_pair_f

    !Spatial sampler: shepard quintic; temporal sampler: mean
    subroutine mui_fetch_shepard_quintic_mean_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_mean_2f_pair_f

    subroutine mui_fetch_shepard_quintic_mean_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_mean_2fx_pair_f

    subroutine mui_fetch_shepard_quintic_mean_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_mean_2d_pair_f

    subroutine mui_fetch_shepard_quintic_mean_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_mean_2dx_pair_f

    subroutine mui_fetch_shepard_quintic_mean_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_mean_2t_pair_f

    !Spatial sampler: shepard quintic; temporal sampler: sum
    subroutine mui_fetch_shepard_quintic_sum_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_sum_2f_pair_f

    subroutine mui_fetch_shepard_quintic_sum_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_sum_2fx_pair_f

    subroutine mui_fetch_shepard_quintic_sum_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_sum_2d_pair_f

    subroutine mui_fetch_shepard_quintic_sum_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_sum_2dx_pair_f

    subroutine mui_fetch_shepard_quintic_sum_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_shepard_quintic_sum_2t_pair_f

    !Spatial sampler: sph-derived quintic; temporal sampler: exact
    subroutine mui_fetch_sph_quintic_exact_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_exact_2f_pair_f

    subroutine mui_fetch_sph_quintic_exact_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_exact_2fx_pair_f

    subroutine mui_fetch_sph_quintic_exact_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_exact_2d_pair_f

    subroutine mui_fetch_sph_quintic_exact_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_exact_2dx_pair_f

    subroutine mui_fetch_sph_quintic_exact_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_exact_2t_pair_f

    !Spatial sampler: sph-derived quintic; temporal sampler: gauss
    subroutine mui_fetch_sph_quintic_gauss_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_gauss_2f_pair_f

    subroutine mui_fetch_sph_quintic_gauss_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_gauss_2fx_pair_f

    subroutine mui_fetch_sph_quintic_gauss_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_gauss_2d_pair_f

    subroutine mui_fetch_sph_quintic_gauss_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_gauss_2dx_pair_f

    subroutine mui_fetch_sph_quintic_gauss_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_gauss_2t_pair_f

    !Spatial sampler: sph-derived quintic; temporal sampler: mean
    subroutine mui_fetch_sph_quintic_mean_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_mean_2f_pair_f

    subroutine mui_fetch_sph_quintic_mean_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_mean_2fx_pair_f

    subroutine mui_fetch_sph_quintic_mean_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_mean_2d_pair_f

    subroutine mui_fetch_sph_quintic_mean_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_mean_2dx_pair_f

    subroutine mui_fetch_sph_quintic_mean_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_mean_2t_pair_f

    !Spatial sampler: sph-derived quintic; temporal sampler: sum
    subroutine mui_fetch_sph_quintic_sum_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_sum_2f_pair_f

    subroutine mui_fetch_sph_quintic_sum_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_sum_2fx_pair_f

    subroutine mui_fetch_sph_quintic_sum_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_sum_2d_pair_f

    subroutine mui_fetch_sph_quintic_sum_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_sum_2dx_pair_f

    subroutine mui_fetch_sph_quintic_sum_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sph_quintic_sum_2t_pair_f

    !Spatial sampler: summation quintic; temporal sampler: exact
    subroutine mui_fetch_sum_quintic_exact_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_exact_2f_pair_f

    subroutine mui_fetch_sum_quintic_exact_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_exact_2fx_pair_f

    subroutine mui_fetch_sum_quintic_exact_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_exact_2d_pair_f

    subroutine mui_fetch_sum_quintic_exact_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_exact_2dx_pair_f

    subroutine mui_fetch_sum_quintic_exact_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_exact_2t_pair_f

    !Spatial sampler: summation quintic; temporal sampler: gauss
    subroutine mui_fetch_sum_quintic_gauss_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_gauss_2f_pair_f

    subroutine mui_fetch_sum_quintic_gauss_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_gauss_2fx_pair_f

    subroutine mui_fetch_sum_quintic_gauss_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_gauss_2d_pair_f

    subroutine mui_fetch_sum_quintic_gauss_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_gauss_2dx_pair_f

    subroutine mui_fetch_sum_quintic_gauss_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_gauss_2t_pair_f

    !Spatial sampler: summation quintic; temporal sampler: mean
    subroutine mui_fetch_sum_quintic_mean_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_mean_2f_pair_f

    subroutine mui_fetch_sum_quintic_mean_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_mean_2fx_pair_f

    subroutine mui_fetch_sum_quintic_mean_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_mean_2d_pair_f

    subroutine mui_fetch_sum_quintic_mean_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_mean_2dx_pair_f

    subroutine mui_fetch_sum_quintic_mean_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_mean_2t_pair_f

    !Spatial sampler: summation quintic; temporal sampler: sum
    subroutine mui_fetch_sum_quintic_sum_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_sum_2f_pair_f

    subroutine mui_fetch_sum_quintic_sum_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_sum_2fx_pair_f

    subroutine mui_fetch_sum_quintic_sum_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_sum_2d_pair_f

    subroutine mui_fetch_sum_quintic_sum_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_sum_2dx_pair_f

    subroutine mui_fetch_sum_quintic_sum_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_sum_quintic_sum_2t_pair_f

#ifdef USE_RBF
    !Spatial sampler: radial basis function; temporal sampler: exact
    subroutine mui_fetch_rbf_exact_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_rbf_exact_2f_pair_f

    subroutine mui_fetch_rbf_exact_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_rbf_exact_2fx_pair_f

    subroutine mui_fetch_rbf_exact_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_rbf_exact_2d_pair_f

    subroutine mui_fetch_rbf_exact_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_rbf_exact_2dx_pair_f

    subroutine mui_fetch_rbf_exact_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_rbf_exact_2t_pair_f

    !Spatial sampler: radial basis function; temporal sampler: gauss
    subroutine mui_fetch_rbf_gauss_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_rbf_gauss_2f_pair_f

    subroutine mui_fetch_rbf_gauss_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_rbf_gauss_2fx_pair_f

    subroutine mui_fetch_rbf_gauss_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_rbf_gauss_2d_pair_f

    subroutine mui_fetch_rbf_gauss_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_rbf_gauss_2dx_pair_f

    subroutine mui_fetch_rbf_gauss_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_rbf_gauss_2t_pair_f

    !Spatial sampler: radial basis function; temporal sampler: mean
    subroutine mui_fetch_rbf_mean_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_rbf_mean_2f_pair_f

    subroutine mui_fetch_rbf_mean_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_rbf_mean_2fx_pair_f

    subroutine mui_fetch_rbf_mean_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_rbf_mean_2d_pair_f

    subroutine mui_fetch_rbf_mean_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_rbf_mean_2dx_pair_f

    subroutine mui_fetch_rbf_mean_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_rbf_mean_2t_pair_f

    !Spatial sampler: radial basis function; temporal sampler: sum
    subroutine mui_fetch_rbf_sum_2f_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_rbf_sum_2f_pair_f

    subroutine mui_fetch_rbf_sum_2fx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_rbf_sum_2fx_pair_f

    subroutine mui_fetch_rbf_sum_2d_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_rbf_sum_2d_pair_f

    subroutine mui_fetch_rbf_sum_2dx_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_rbf_sum_2dx_pair_f

    subroutine mui_fetch_rbf_sum_2t_pair_f(uniface,attr,point_1,point_2,t_1,t_2,&
    spatial_sampler,temporal_sampler,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: point_1,point_2,t_1,t_2
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_rbf_sum_2t_pair_f
#endif

    !*******************************************************************
    !* MUI functions for 2D data point only fetch using one time value *
    !*******************************************************************

    !Temporal sampler: exact
    subroutine mui_fetch_points_exact_2f_f(uniface,attr,t,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_exact_2f_f

    subroutine mui_fetch_points_exact_2fx_f(uniface,attr,t,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_exact_2fx_f

    subroutine mui_fetch_points_exact_2d_f(uniface,attr,t,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_exact_2d_f

    subroutine mui_fetch_points_exact_2dx_f(uniface,attr,t,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_exact_2dx_f

    subroutine mui_fetch_points_exact_2t_f(uniface,attr,t,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_exact_2t_f

    !Temporal sampler: gauss
    subroutine mui_fetch_points_gauss_2f_f(uniface,attr,t,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_gauss_2f_f

    subroutine mui_fetch_points_gauss_2fx_f(uniface,attr,t,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_gauss_2fx_f

    subroutine mui_fetch_points_gauss_2d_f(uniface,attr,t,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_gauss_2d_f

    subroutine mui_fetch_points_gauss_2dx_f(uniface,attr,t,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_gauss_2dx_f

    subroutine mui_fetch_points_gauss_2t_f(uniface,attr,t,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_gauss_2t_f

    !Temporal sampler: mean
    subroutine mui_fetch_points_mean_2f_f(uniface,attr,t,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_mean_2f_f

    subroutine mui_fetch_points_mean_2fx_f(uniface,attr,t,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_mean_2fx_f

    subroutine mui_fetch_points_mean_2d_f(uniface,attr,t,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_mean_2d_f

    subroutine mui_fetch_points_mean_2dx_f(uniface,attr,t,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_mean_2dx_f

    subroutine mui_fetch_points_mean_2t_f(uniface,attr,t,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_mean_2t_f

    !Temporal sampler: sum
    subroutine mui_fetch_points_sum_2f_f(uniface,attr,t,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_sum_2f_f

    subroutine mui_fetch_points_sum_2fx_f(uniface,attr,t,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_sum_2fx_f

    subroutine mui_fetch_points_sum_2d_f(uniface,attr,t,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_sum_2d_f

    subroutine mui_fetch_points_sum_2dx_f(uniface,attr,t,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_sum_2dx_f

    subroutine mui_fetch_points_sum_2t_f(uniface,attr,t,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_sum_2t_f

    !********************************************************************
    !* MUI functions for 2D data point only fetch using two time values *
    !********************************************************************

    !Temporal sampler: exact
    subroutine mui_fetch_points_exact_2f_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_exact_2f_pair_f

    subroutine mui_fetch_points_exact_2fx_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_exact_2fx_pair_f

    subroutine mui_fetch_points_exact_2d_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_exact_2d_pair_f

    subroutine mui_fetch_points_exact_2dx_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_exact_2dx_pair_f

    subroutine mui_fetch_points_exact_2t_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_exact_2t_pair_f

    !Temporal sampler: gauss
    subroutine mui_fetch_points_gauss_2f_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_gauss_2f_pair_f

    subroutine mui_fetch_points_gauss_2fx_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_gauss_2fx_pair_f

    subroutine mui_fetch_points_gauss_2d_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_gauss_2d_pair_f

    subroutine mui_fetch_points_gauss_2dx_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_gauss_2dx_pair_f

    subroutine mui_fetch_points_gauss_2t_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_gauss_2t_pair_f

    !Temporal sampler: mean
    subroutine mui_fetch_points_mean_2f_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_mean_2f_pair_f

    subroutine mui_fetch_points_mean_2fx_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_mean_2fx_pair_f

    subroutine mui_fetch_points_mean_2d_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_mean_2d_pair_f

    subroutine mui_fetch_points_mean_2dx_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_mean_2dx_pair_f

    subroutine mui_fetch_points_mean_2t_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_mean_2t_pair_f

    !Temporal sampler: sum
    subroutine mui_fetch_points_sum_2f_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_sum_2f_pair_f

    subroutine mui_fetch_points_sum_2fx_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_sum_2fx_pair_f

    subroutine mui_fetch_points_sum_2d_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_sum_2d_pair_f

    subroutine mui_fetch_points_sum_2dx_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_sum_2dx_pair_f

    subroutine mui_fetch_points_sum_2t_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_points_1,ret_points_2,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_points_1(*),ret_points_2(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_points_sum_2t_pair_f

    !********************************************************************
    !* MUI functions for 2D data values only fetch using one time value *
    !********************************************************************

    !Temporal sampler: exact
    subroutine mui_fetch_values_exact_2f_f(uniface,attr,t,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_exact_2f_f

    subroutine mui_fetch_values_exact_2fx_f(uniface,attr,t,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_exact_2fx_f

    subroutine mui_fetch_values_exact_2d_f(uniface,attr,t,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_exact_2d_f

    subroutine mui_fetch_values_exact_2dx_f(uniface,attr,t,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_exact_2dx_f

    subroutine mui_fetch_values_exact_2t_f(uniface,attr,t,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_exact_2t_f

    !Temporal sampler: gauss
    subroutine mui_fetch_values_gauss_2f_f(uniface,attr,t,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_gauss_2f_f

    subroutine mui_fetch_values_gauss_2fx_f(uniface,attr,t,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_gauss_2fx_f

    subroutine mui_fetch_values_gauss_2d_f(uniface,attr,t,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_gauss_2d_f

    subroutine mui_fetch_values_gauss_2dx_f(uniface,attr,t,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_gauss_2dx_f

    subroutine mui_fetch_values_gauss_2t_f(uniface,attr,t,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_gauss_2t_f

    !Temporal sampler: mean
    subroutine mui_fetch_values_mean_2f_f(uniface,attr,t,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_mean_2f_f

    subroutine mui_fetch_values_mean_2fx_f(uniface,attr,t,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_mean_2fx_f

    subroutine mui_fetch_values_mean_2d_f(uniface,attr,t,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_mean_2d_f

    subroutine mui_fetch_values_mean_2dx_f(uniface,attr,t,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_mean_2dx_f

    subroutine mui_fetch_values_mean_2t_f(uniface,attr,t,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_mean_2t_f

    !Temporal sampler: sum
    subroutine mui_fetch_values_sum_2f_f(uniface,attr,t,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_sum_2f_f

    subroutine mui_fetch_values_sum_2fx_f(uniface,attr,t,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_sum_2fx_f

    subroutine mui_fetch_values_sum_2d_f(uniface,attr,t,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_sum_2d_f

    subroutine mui_fetch_values_sum_2dx_f(uniface,attr,t,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_sum_2dx_f

    subroutine mui_fetch_values_sum_2t_f(uniface,attr,t,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_sum_2t_f

    !*********************************************************************
    !* MUI functions for 2D data values only fetch using two time values *
    !*********************************************************************

    !Temporal sampler: exact
    subroutine mui_fetch_values_exact_2f_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_exact_2f_pair_f

    subroutine mui_fetch_values_exact_2fx_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_exact_2fx_pair_f

    subroutine mui_fetch_values_exact_2d_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_exact_2d_pair_f

    subroutine mui_fetch_values_exact_2dx_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_exact_2dx_pair_f

    subroutine mui_fetch_values_exact_2t_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_exact_2t_pair_f

    !Temporal sampler: gauss
    subroutine mui_fetch_values_gauss_2f_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_gauss_2f_pair_f

    subroutine mui_fetch_values_gauss_2fx_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_gauss_2fx_pair_f

    subroutine mui_fetch_values_gauss_2d_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_gauss_2d_pair_f

    subroutine mui_fetch_values_gauss_2dx_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_gauss_2dx_pair_f

    subroutine mui_fetch_values_gauss_2t_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_gauss_2t_pair_f

    !Temporal sampler: mean
    subroutine mui_fetch_values_mean_2f_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_mean_2f_pair_f

    subroutine mui_fetch_values_mean_2fx_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_mean_2fx_pair_f

    subroutine mui_fetch_values_mean_2d_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_mean_2d_pair_f

    subroutine mui_fetch_values_mean_2dx_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_mean_2dx_pair_f

    subroutine mui_fetch_values_mean_2t_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_mean_2t_pair_f

    !Temporal sampler: sum
    subroutine mui_fetch_values_sum_2f_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_sum_2f_pair_f

    subroutine mui_fetch_values_sum_2fx_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_sum_2fx_pair_f

    subroutine mui_fetch_values_sum_2d_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_sum_2d_pair_f

    subroutine mui_fetch_values_sum_2dx_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_sum_2dx_pair_f

    subroutine mui_fetch_values_sum_2t_pair_f(uniface,attr,t_1,t_2,&
        temporal_sampler,ret_values_1,num_points) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface,temporal_sampler
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t_1,t_2
      type(c_ptr), intent(out) :: ret_values_1(*)
      integer(kind=c_int), intent(in) :: num_points
    end subroutine mui_fetch_values_sum_2t_pair_f

    !********************************************
    !* MUI functions for single parameter fetch *
    !********************************************

    subroutine mui_fetch_2f_param_f(uniface,attr,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_2f_param_f

    subroutine mui_fetch_2fx_param_f(uniface,attr,return_value) bind(C)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(out) :: return_value
    end subroutine mui_fetch_2fx_param_f

    subroutine mui_fetch_2d_param_f(uniface,attr,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_2d_param_f

    subroutine mui_fetch_2dx_param_f(uniface,attr,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_2dx_param_f

    subroutine mui_fetch_2t_param_f(uniface,attr,return_value) bind(C)
      import :: c_ptr,c_char,c_double
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(out) :: return_value
    end subroutine mui_fetch_2t_param_f

    !******************************************
    !* MUI data receive test functions        *
    !******************************************

    !Data ready test using single time value
    subroutine mui_is_ready_2f_f(uniface,attr,t,return_value) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t
      integer(kind=c_int), intent(out) :: return_value
    end subroutine mui_is_ready_2f_f

    subroutine mui_is_ready_2fx_f(uniface,attr,t,return_value) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t
      integer(kind=c_int), intent(out) :: return_value
    end subroutine mui_is_ready_2fx_f

    subroutine mui_is_ready_2d_f(uniface,attr,t,return_value) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t
      integer(kind=c_int), intent(out) :: return_value
    end subroutine mui_is_ready_2d_f

    subroutine mui_is_ready_2dx_f(uniface,attr,t,return_value) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t
      integer(kind=c_int), intent(out) :: return_value
    end subroutine mui_is_ready_2dx_f

    subroutine mui_is_ready_2t_f(uniface,attr,t,return_value) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t
      integer(kind=c_int), intent(out) :: return_value
    end subroutine mui_is_ready_2t_f

    !Data ready test using two time values
    subroutine mui_is_ready_2f_pair_f(uniface,attr,t_1,t_2,return_value) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t_1,t_2
      integer(kind=c_int), intent(out) :: return_value
    end subroutine mui_is_ready_2f_pair_f

    subroutine mui_is_ready_2fx_pair_f(uniface,attr,t_1,t_2,return_value) bind(C)
      import :: c_ptr,c_char,c_int,c_float
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_float), intent(in) :: t_1,t_2
      integer(kind=c_int), intent(out) :: return_value
    end subroutine mui_is_ready_2fx_pair_f

    subroutine mui_is_ready_2d_pair_f(uniface,attr,t_1,t_2,return_value) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t_1,t_2
      integer(kind=c_int), intent(out) :: return_value
    end subroutine mui_is_ready_2d_pair_f

    subroutine mui_is_ready_2dx_pair_f(uniface,attr,t_1,t_2,return_value) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t_1,t_2
      integer(kind=c_int), intent(out) :: return_value
    end subroutine mui_is_ready_2dx_pair_f

    subroutine mui_is_ready_2t_pair_f(uniface,attr,t_1,t_2,return_value) bind(C)
      import :: c_ptr,c_char,c_int,c_double
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(in) :: attr(*)
      real(kind=c_double), intent(in) :: t_1,t_2
      integer(kind=c_int), intent(out) :: return_value
    end subroutine mui_is_ready_2t_pair_f

    !******************************************
    !* MUI Smart Send functions               *
    !******************************************

    !Send span announce using 2D box geometry
    subroutine mui_announce_send_span_2f_box_f(uniface,box_1_1,box_1_2,box_2_1,box_2_2,t_start, &
      t_timeout,synchronised) bind(C)
      import :: c_ptr,c_float,c_int
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in) :: box_1_1,box_1_2,box_2_1,box_2_2,t_start,t_timeout
      integer(kind=c_int), intent(in) :: synchronised
    end subroutine mui_announce_send_span_2f_box_f

    subroutine mui_announce_send_span_2fx_box_f(uniface,box_1_1,box_1_2,box_2_1,box_2_2,t_start, &
      t_timeout,synchronised) bind(C)
      import :: c_ptr,c_float,c_int
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in) :: box_1_1,box_1_2,box_2_1,box_2_2,t_start,t_timeout
      integer(kind=c_int), intent(in) :: synchronised
    end subroutine mui_announce_send_span_2fx_box_f

    subroutine mui_announce_send_span_2d_box_f(uniface,box_1_1,box_1_2,box_2_1,box_2_2,t_start, &
      t_timeout,synchronised) bind(C)
      import :: c_ptr,c_double,c_int
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in) :: box_1_1,box_1_2,box_2_1,box_2_2,t_start,t_timeout
      integer(kind=c_int), intent(in) :: synchronised
    end subroutine mui_announce_send_span_2d_box_f

    subroutine mui_announce_send_span_2dx_box_f(uniface,box_1_1,box_1_2,box_2_1,box_2_2,t_start, &
      t_timeout,synchronised) bind(C)
      import :: c_ptr,c_double,c_int
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in) :: box_1_1,box_1_2,box_2_1,box_2_2,t_start,t_timeout
      integer(kind=c_int), intent(in) :: synchronised
    end subroutine mui_announce_send_span_2dx_box_f

    subroutine mui_announce_send_span_2t_box_f(uniface,box_1_1,box_1_2,box_2_1,box_2_2,t_start, &
      t_timeout,synchronised) bind(C)
      import :: c_ptr,c_double,c_int
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in) :: box_1_1,box_1_2,box_2_1,box_2_2,t_start,t_timeout
      integer(kind=c_int), intent(in) :: synchronised
    end subroutine mui_announce_send_span_2t_box_f

    !Send span announce using 2D sphere geometry
    subroutine mui_announce_send_span_2f_sphere_f(uniface,centre_1,centre_2,radius,t_start, &
      t_timeout,synchronised) bind(C)
      import :: c_ptr,c_float,c_int
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in) :: centre_1,centre_2,radius,t_start,t_timeout
      integer(kind=c_int), intent(in) :: synchronised
    end subroutine mui_announce_send_span_2f_sphere_f

    subroutine mui_announce_send_span_2fx_sphere_f(uniface,centre_1,centre_2,radius,t_start, &
      t_timeout,synchronised) bind(C)
      import :: c_ptr,c_float,c_int
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in) :: centre_1,centre_2,radius,t_start,t_timeout
      integer(kind=c_int), intent(in) :: synchronised
    end subroutine mui_announce_send_span_2fx_sphere_f

    subroutine mui_announce_send_span_2d_sphere_f(uniface,centre_1,centre_2,radius,t_start, &
      t_timeout,synchronised) bind(C)
      import :: c_ptr,c_double,c_int
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in) :: centre_1,centre_2,radius,t_start,t_timeout
      integer(kind=c_int), intent(in) :: synchronised
    end subroutine mui_announce_send_span_2d_sphere_f

    subroutine mui_announce_send_span_2dx_sphere_f(uniface,centre_1,centre_2,radius,t_start, &
      t_timeout,synchronised) bind(C)
      import :: c_ptr,c_double,c_int
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in) :: centre_1,centre_2,radius,t_start,t_timeout
      integer(kind=c_int), intent(in) :: synchronised
    end subroutine mui_announce_send_span_2dx_sphere_f

    subroutine mui_announce_send_span_2t_sphere_f(uniface,centre_1,centre_2,radius,t_start, &
      t_timeout,synchronised) bind(C)
      import :: c_ptr,c_double,c_int
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in) :: centre_1,centre_2,radius,t_start,t_timeout
      integer(kind=c_int), intent(in) :: synchronised
    end subroutine mui_announce_send_span_2t_sphere_f

    !Receive span announce using 2D box geometry
    subroutine mui_announce_recv_span_2f_box_f(uniface,box_1_1,box_1_2,box_2_1,box_2_2,t_start, &
      t_timeout,synchronised) bind(C)
      import :: c_ptr,c_float,c_int
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in) :: box_1_1,box_1_2,box_2_1,box_2_2,t_start,t_timeout
      integer(kind=c_int), intent(in) :: synchronised
    end subroutine mui_announce_recv_span_2f_box_f

    subroutine mui_announce_recv_span_2fx_box_f(uniface,box_1_1,box_1_2,box_2_1,box_2_2,t_start, &
      t_timeout,synchronised) bind(C)
      import :: c_ptr,c_float,c_int
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in) :: box_1_1,box_1_2,box_2_1,box_2_2,t_start,t_timeout
      integer(kind=c_int), intent(in) :: synchronised
    end subroutine mui_announce_recv_span_2fx_box_f

    subroutine mui_announce_recv_span_2d_box_f(uniface,box_1_1,box_1_2,box_2_1,box_2_2,t_start, &
      t_timeout,synchronised) bind(C)
      import :: c_ptr,c_double,c_int
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in) :: box_1_1,box_1_2,box_2_1,box_2_2,t_start,t_timeout
      integer(kind=c_int), intent(in) :: synchronised
    end subroutine mui_announce_recv_span_2d_box_f

    subroutine mui_announce_recv_span_2dx_box_f(uniface,box_1_1,box_1_2,box_2_1,box_2_2,t_start, &
      t_timeout,synchronised) bind(C)
      import :: c_ptr,c_double,c_int
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in) :: box_1_1,box_1_2,box_2_1,box_2_2,t_start,t_timeout
      integer(kind=c_int), intent(in) :: synchronised
    end subroutine mui_announce_recv_span_2dx_box_f

    subroutine mui_announce_recv_span_2t_box_f(uniface,box_1_1,box_1_2,box_2_1,box_2_2,t_start, &
      t_timeout,synchronised) bind(C)
      import :: c_ptr,c_double,c_int
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in) :: box_1_1,box_1_2,box_2_1,box_2_2,t_start,t_timeout
      integer(kind=c_int), intent(in) :: synchronised
    end subroutine mui_announce_recv_span_2t_box_f

    !Receive span announce using 2D sphere geometry
    subroutine mui_announce_recv_span_2f_sphere_f(uniface,centre_1,centre_2,radius,t_start, &
      t_timeout,synchronised) bind(C)
      import :: c_ptr,c_float,c_int
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in) :: centre_1,centre_2,radius,t_start,t_timeout
      integer(kind=c_int), intent(in) :: synchronised
    end subroutine mui_announce_recv_span_2f_sphere_f

    subroutine mui_announce_recv_span_2fx_sphere_f(uniface,centre_1,centre_2,radius,t_start, &
      t_timeout,synchronised) bind(C)
      import :: c_ptr,c_float,c_int
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in) :: centre_1,centre_2,radius,t_start,t_timeout
      integer(kind=c_int), intent(in) :: synchronised
    end subroutine mui_announce_recv_span_2fx_sphere_f

    subroutine mui_announce_recv_span_2d_sphere_f(uniface,centre_1,centre_2,radius,t_start, &
      t_timeout,synchronised) bind(C)
      import :: c_ptr,c_double,c_int
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in) :: centre_1,centre_2,radius,t_start,t_timeout
      integer(kind=c_int), intent(in) :: synchronised
    end subroutine mui_announce_recv_span_2d_sphere_f

    subroutine mui_announce_recv_span_2dx_sphere_f(uniface,centre_1,centre_2,radius,t_start, &
      t_timeout,synchronised) bind(C)
      import :: c_ptr,c_double,c_int
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in) :: centre_1,centre_2,radius,t_start,t_timeout
      integer(kind=c_int), intent(in) :: synchronised
    end subroutine mui_announce_recv_span_2dx_sphere_f

    subroutine mui_announce_recv_span_2t_sphere_f(uniface,centre_1,centre_2,radius,t_start, &
      t_timeout,synchronised) bind(C)
      import :: c_ptr,c_double,c_int
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in) :: centre_1,centre_2,radius,t_start,t_timeout
      integer(kind=c_int), intent(in) :: synchronised
    end subroutine mui_announce_recv_span_2t_sphere_f

    !Send disable announce (local call per MPI rank)
    subroutine mui_announce_send_disable_2f_f(uniface,synchronised) bind(C)
      import :: c_ptr,c_int
      type(c_ptr), intent(in), value :: uniface
      integer(kind=c_int), intent(in) :: synchronised
    end subroutine mui_announce_send_disable_2f_f

    subroutine mui_announce_send_disable_2fx_f(uniface,synchronised) bind(C)
      import :: c_ptr,c_int
      type(c_ptr), intent(in), value :: uniface
      integer(kind=c_int), intent(in) :: synchronised
    end subroutine mui_announce_send_disable_2fx_f

    subroutine mui_announce_send_disable_2d_f(uniface,synchronised) bind(C)
      import :: c_ptr,c_int
      type(c_ptr), intent(in), value :: uniface
      integer(kind=c_int), intent(in) :: synchronised
    end subroutine mui_announce_send_disable_2d_f

    subroutine mui_announce_send_disable_2dx_f(uniface,synchronised) bind(C)
      import :: c_ptr,c_int
      type(c_ptr), intent(in), value :: uniface
      integer(kind=c_int), intent(in) :: synchronised
    end subroutine mui_announce_send_disable_2dx_f

    subroutine mui_announce_send_disable_2t_f(uniface,synchronised) bind(C)
      import :: c_ptr,c_int
      type(c_ptr), intent(in), value :: uniface
      integer(kind=c_int), intent(in) :: synchronised
    end subroutine mui_announce_send_disable_2t_f

    !Receive disable announce (local call per MPI rank)
    subroutine mui_announce_recv_disable_2f_f(uniface,synchronised) bind(C)
      import :: c_ptr,c_int
      type(c_ptr), intent(in), value :: uniface
      integer(kind=c_int), intent(in) :: synchronised
    end subroutine mui_announce_recv_disable_2f_f

    subroutine mui_announce_recv_disable_2fx_f(uniface,synchronised) bind(C)
      import :: c_ptr,c_int
      type(c_ptr), intent(in), value :: uniface
      integer(kind=c_int), intent(in) :: synchronised
    end subroutine mui_announce_recv_disable_2fx_f

    subroutine mui_announce_recv_disable_2d_f(uniface,synchronised) bind(C)
      import :: c_ptr,c_int
      type(c_ptr), intent(in), value :: uniface
      integer(kind=c_int), intent(in) :: synchronised
    end subroutine mui_announce_recv_disable_2d_f

    subroutine mui_announce_recv_disable_2dx_f(uniface,synchronised) bind(C)
      import :: c_ptr,c_int
      type(c_ptr), intent(in), value :: uniface
      integer(kind=c_int), intent(in) :: synchronised
    end subroutine mui_announce_recv_disable_2dx_f

    subroutine mui_announce_recv_disable_2t_f(uniface,synchronised) bind(C)
      import :: c_ptr,c_int
      type(c_ptr), intent(in), value :: uniface
      integer(kind=c_int), intent(in) :: synchronised
    end subroutine mui_announce_recv_disable_2t_f

    !******************************************
    !* MUI barrier functions                  *
    !******************************************

    !Barrier at single time value
    subroutine mui_barrier_2f_f(uniface,t) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in), value :: t
    end subroutine mui_barrier_2f_f

    subroutine mui_barrier_2fx_f(uniface,t) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in), value :: t
    end subroutine mui_barrier_2fx_f

    subroutine mui_barrier_2d_f(uniface,t) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in), value :: t
    end subroutine mui_barrier_2d_f

    subroutine mui_barrier_2dx_f(uniface,t) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in), value :: t
    end subroutine mui_barrier_2dx_f

    subroutine mui_barrier_2t_f(uniface,t) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in), value :: t
    end subroutine mui_barrier_2t_f

    !Barrier at two time values
    subroutine mui_barrier_2f_pair_f(uniface,t_1,t_2) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in), value :: t_1,t_2
    end subroutine mui_barrier_2f_pair_f

    subroutine mui_barrier_2fx_pair_f(uniface,t_1,t_2) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in), value :: t_1,t_2
    end subroutine mui_barrier_2fx_pair_f

    subroutine mui_barrier_2d_pair_f(uniface,t_1,t_2) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in), value :: t_1,t_2
    end subroutine mui_barrier_2d_pair_f

    subroutine mui_barrier_2dx_pair_f(uniface,t_1,t_2) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in), value :: t_1,t_2
    end subroutine mui_barrier_2dx_pair_f

    subroutine mui_barrier_2t_pair_f(uniface,t_1,t_2) bind(C)
      import :: c_ptr,c_double
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in), value :: t_1,t_2
    end subroutine mui_barrier_2t_pair_f

    !******************************************
    !* MUI forget functions                   *
    !******************************************

    !Forget log between [-inf, upper]
    subroutine mui_forget_upper_2f_f(uniface,upper,reset_log) bind(C)
      import :: c_ptr,c_int,c_float
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in), value :: upper
      integer(kind=c_int), intent(in), value :: reset_log
    end subroutine mui_forget_upper_2f_f

    subroutine mui_forget_upper_2fx_f(uniface,upper,reset_log) bind(C)
      import :: c_ptr,c_int,c_float
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in), value :: upper
      integer(kind=c_int), intent(in), value :: reset_log
    end subroutine mui_forget_upper_2fx_f

    subroutine mui_forget_upper_2d_f(uniface,upper,reset_log) bind(C)
      import :: c_ptr,c_int,c_double
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in), value :: upper
      integer(kind=c_int), intent(in), value :: reset_log
    end subroutine mui_forget_upper_2d_f

    subroutine mui_forget_upper_2dx_f(uniface,upper,reset_log) bind(C)
      import :: c_ptr,c_int,c_double
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in), value :: upper
      integer(kind=c_int), intent(in), value :: reset_log
    end subroutine mui_forget_upper_2dx_f

    subroutine mui_forget_upper_2t_f(uniface,upper,reset_log) bind(C)
      import :: c_ptr,c_int,c_double
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in), value :: upper
      integer(kind=c_int), intent(in), value :: reset_log
    end subroutine mui_forget_upper_2t_f

    !Forget log between [-inf, -inf], [upper_1, upper_2]
    subroutine mui_forget_upper_2f_pair_f(uniface,upper_1,upper_2,reset_log) bind(C)
      import :: c_ptr,c_int,c_float
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in), value :: upper_1,upper_2
      integer(kind=c_int), intent(in), value :: reset_log
    end subroutine mui_forget_upper_2f_pair_f

    subroutine mui_forget_upper_2fx_pair_f(uniface,upper_1,upper_2,reset_log) bind(C)
      import :: c_ptr,c_int,c_float
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in), value :: upper_1,upper_2
      integer(kind=c_int), intent(in), value :: reset_log
    end subroutine mui_forget_upper_2fx_pair_f

    subroutine mui_forget_upper_2d_pair_f(uniface,upper_1,upper_2,reset_log) bind(C)
      import :: c_ptr,c_int,c_double
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in), value :: upper_1,upper_2
      integer(kind=c_int), intent(in), value :: reset_log
    end subroutine mui_forget_upper_2d_pair_f

    subroutine mui_forget_upper_2dx_pair_f(uniface,upper_1,upper_2,reset_log) bind(C)
      import :: c_ptr,c_int,c_double
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in), value :: upper_1,upper_2
      integer(kind=c_int), intent(in), value :: reset_log
    end subroutine mui_forget_upper_2dx_pair_f

    subroutine mui_forget_upper_2t_pair_f(uniface,upper_1,upper_2,reset_log) bind(C)
      import :: c_ptr,c_int,c_double
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in), value :: upper_1,upper_2
      integer(kind=c_int), intent(in), value :: reset_log
    end subroutine mui_forget_upper_2t_pair_f

    !Forget log between [lower, upper]
    subroutine mui_forget_lower_upper_2f_f(uniface,lower,upper,reset_log) bind(C)
      import :: c_ptr,c_int,c_float
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in), value :: lower,upper
      integer(kind=c_int), intent(in), value :: reset_log
    end subroutine mui_forget_lower_upper_2f_f

    subroutine mui_forget_lower_upper_2fx_f(uniface,lower,upper,reset_log) bind(C)
      import :: c_ptr,c_int,c_float
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in), value :: lower,upper
      integer(kind=c_int), intent(in), value :: reset_log
    end subroutine mui_forget_lower_upper_2fx_f

    subroutine mui_forget_lower_upper_2d_f(uniface,lower,upper,reset_log) bind(C)
      import :: c_ptr,c_int,c_double
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in), value :: lower,upper
      integer(kind=c_int), intent(in), value :: reset_log
    end subroutine mui_forget_lower_upper_2d_f

    subroutine mui_forget_lower_upper_2dx_f(uniface,lower,upper,reset_log) bind(C)
      import :: c_ptr,c_int,c_double
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in), value :: lower,upper
      integer(kind=c_int), intent(in), value :: reset_log
    end subroutine mui_forget_lower_upper_2dx_f

    subroutine mui_forget_lower_upper_2t_f(uniface,lower,upper,reset_log) bind(C)
      import :: c_ptr,c_int,c_double
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in), value :: lower,upper
      integer(kind=c_int), intent(in), value :: reset_log
    end subroutine mui_forget_lower_upper_2t_f

    !Forget log between [lower_1, lower_2], [upper_1, upper_2]
    subroutine mui_forget_lower_upper_2f_pair_f(uniface,lower_1,lower_2,upper_1,upper_2, &
      reset_log) bind(C)
      import :: c_ptr,c_int,c_float
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in), value :: lower_1,lower_2,upper_1,upper_2
      integer(kind=c_int), intent(in), value :: reset_log
    end subroutine mui_forget_lower_upper_2f_pair_f

    subroutine mui_forget_lower_upper_2fx_pair_f(uniface,lower_1,lower_2,upper_1,upper_2, &
      reset_log) bind(C)
      import :: c_ptr,c_int,c_float
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in), value :: lower_1,lower_2,upper_1,upper_2
      integer(kind=c_int), intent(in), value :: reset_log
    end subroutine mui_forget_lower_upper_2fx_pair_f

    subroutine mui_forget_lower_upper_2d_pair_f(uniface,lower_1,lower_2,upper_1,upper_2, &
      reset_log) bind(C)
      import :: c_ptr,c_int,c_double
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in), value :: lower_1,lower_2,upper_1,upper_2
      integer(kind=c_int), intent(in), value :: reset_log
    end subroutine mui_forget_lower_upper_2d_pair_f

    subroutine mui_forget_lower_upper_2dx_pair_f(uniface,lower_1,lower_2,upper_1,upper_2, &
      reset_log) bind(C)
      import :: c_ptr,c_int,c_double
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in), value :: lower_1,lower_2,upper_1,upper_2
      integer(kind=c_int), intent(in), value :: reset_log
    end subroutine mui_forget_lower_upper_2dx_pair_f

    subroutine mui_forget_lower_upper_2t_pair_f(uniface,lower_1,lower_2,upper_1,upper_2, &
      reset_log) bind(C)
      import :: c_ptr,c_int,c_double
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_double), intent(in), value :: lower_1,lower_2,upper_1,upper_2
      integer(kind=c_int), intent(in), value :: reset_log
    end subroutine mui_forget_lower_upper_2t_pair_f

    !Set to forget log between [-inf, current-length] automatically
    subroutine mui_set_forget_length_2f_f(uniface,length) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in), value :: length
    end subroutine mui_set_forget_length_2f_f

    subroutine mui_set_forget_length_2fx_f(uniface,length) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in), value :: length
    end subroutine mui_set_forget_length_2fx_f

    subroutine mui_set_forget_length_2d_f(uniface,length) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in), value :: length
    end subroutine mui_set_forget_length_2d_f

    subroutine mui_set_forget_length_2dx_f(uniface,length) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in), value :: length
    end subroutine mui_set_forget_length_2dx_f

    subroutine mui_set_forget_length_2t_f(uniface,length) bind(C)
      import :: c_ptr,c_float
      type(c_ptr), intent(in), value :: uniface
      real(kind=c_float), intent(in), value :: length
    end subroutine mui_set_forget_length_2t_f

    !******************************************
    !* MUI URI functions                      *
    !******************************************

    !Obtain original URI host value from existing interface
    subroutine mui_uri_host_2f_f(uniface,return_val) bind(C)
      import :: c_ptr,c_char
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(out) :: return_val(*)
    end subroutine mui_uri_host_2f_f

    subroutine mui_uri_host_2fx_f(uniface,return_val) bind(C)
      import :: c_ptr,c_char
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(out) :: return_val(*)
    end subroutine mui_uri_host_2fx_f

    subroutine mui_uri_host_2d_f(uniface,return_val) bind(C)
      import :: c_ptr,c_char
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(out) :: return_val(*)
    end subroutine mui_uri_host_2d_f

    subroutine mui_uri_host_2dx_f(uniface,return_val) bind(C)
      import :: c_ptr,c_char
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(out) :: return_val(*)
    end subroutine mui_uri_host_2dx_f

    subroutine mui_uri_host_2t_f(uniface,return_val) bind(C)
      import :: c_ptr,c_char
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(out) :: return_val(*)
    end subroutine mui_uri_host_2t_f

    !Obtain original URI path value from existing interface
    subroutine mui_uri_path_2f_f(uniface,return_val) bind(C)
      import :: c_ptr,c_char
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(out) :: return_val(*)
    end subroutine mui_uri_path_2f_f

    subroutine mui_uri_path_2fx_f(uniface,return_val) bind(C)
      import :: c_ptr,c_char
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(out) :: return_val(*)
    end subroutine mui_uri_path_2fx_f

    subroutine mui_uri_path_2d_f(uniface,return_val) bind(C)
      import :: c_ptr,c_char
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(out) :: return_val(*)
    end subroutine mui_uri_path_2d_f

    subroutine mui_uri_path_2dx_f(uniface,return_val) bind(C)
      import :: c_ptr,c_char
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(out) :: return_val(*)
    end subroutine mui_uri_path_2dx_f

    subroutine mui_uri_path_2t_f(uniface,return_val) bind(C)
      import :: c_ptr,c_char
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(out) :: return_val(*)
    end subroutine mui_uri_path_2t_f

    ! Obtain original URI protocol value from existing interface
    subroutine mui_uri_protocol_2f_f(uniface,return_val) bind(C)
      import :: c_ptr,c_char
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(out) :: return_val(*)
    end subroutine mui_uri_protocol_2f_f

    subroutine mui_uri_protocol_2fx_f(uniface,return_val) bind(C)
      import :: c_ptr,c_char
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(out) :: return_val(*)
    end subroutine mui_uri_protocol_2fx_f

    subroutine mui_uri_protocol_2d_f(uniface,return_val) bind(C)
      import :: c_ptr,c_char
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(out) :: return_val(*)
    end subroutine mui_uri_protocol_2d_f

    subroutine mui_uri_protocol_2dx_f(uniface,return_val) bind(C)
      import :: c_ptr,c_char
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(out) :: return_val(*)
    end subroutine mui_uri_protocol_2dx_f

    subroutine mui_uri_protocol_2t_f(uniface,return_val) bind(C)
      import :: c_ptr,c_char
      type(c_ptr), intent(in), value :: uniface
      character(kind=c_char), intent(out) :: return_val(*)
    end subroutine mui_uri_protocol_2t_f

  end interface 

  contains

  !****************************************
  !* Create and get MUI interfaces for    *
  !* multi-domain function                *
  !****************************************

  !Create and access set of 2D interfaces with float=single and int=int32
  subroutine create_and_get_uniface_multi_2f_f(uniface_pointers_2d, domain, interfaces, &
    interface_count)
    use, intrinsic :: iso_c_binding
    implicit none

    type(ptr_typ_2d), target :: uniface_pointers_2d(:)
    character(kind=c_char), intent(in) :: domain(*)
    character(kind=c_char,len=*), intent(in) :: interfaces(*)
    integer(kind=c_int), VALUE :: interface_count
    integer :: i

    call mui_create_uniface_multi_2f_f(domain, interfaces, &
      interface_count)

    do i = 1, interface_count
      uniface_pointers_2d(i)%ptr = get_mui_uniface_multi_2f_f(i)
    end do
  end subroutine create_and_get_uniface_multi_2f_f

  !Create and access set of 2D interfaces with float=single and int=int64
  subroutine create_and_get_uniface_multi_2fx_f(uniface_pointers_2d, domain, interfaces, &
    interface_count)
    use, intrinsic :: iso_c_binding
    implicit none

    type(ptr_typ_2d), target :: uniface_pointers_2d(:)
    character(kind=c_char), intent(in) :: domain(*)
    character(kind=c_char,len=*), intent(in) :: interfaces(*)
    integer(kind=c_int), VALUE :: interface_count
    integer :: i

    call mui_create_uniface_multi_2fx_f(domain, interfaces, &
      interface_count)

    do i = 1, interface_count
      uniface_pointers_2d(i)%ptr = get_mui_uniface_multi_2fx_f(i)
    end do
  end subroutine create_and_get_uniface_multi_2fx_f

  !Create and access set of 2D interfaces with float=double and int=int32
  subroutine create_and_get_uniface_multi_2d_f(uniface_pointers_2d, domain, interfaces, &
    interface_count)
    use, intrinsic :: iso_c_binding
    implicit none

    type(ptr_typ_2d), target :: uniface_pointers_2d(:)
    character(kind=c_char), intent(in) :: domain(*)
    character(kind=c_char,len=*), intent(in) :: interfaces(*)
    integer(kind=c_int), VALUE :: interface_count
    integer :: i

    call mui_create_uniface_multi_2d_f(domain, interfaces, &
      interface_count)

    do i = 1, interface_count
      uniface_pointers_2d(i)%ptr = get_mui_uniface_multi_2d_f(i)
    end do
  end subroutine create_and_get_uniface_multi_2d_f

  !Create and access set of 2D interfaces with float=double and int=int64
  subroutine create_and_get_uniface_multi_2dx_f(uniface_pointers_2d, domain, interfaces, &
    interface_count)
    use, intrinsic :: iso_c_binding
    implicit none

    type(ptr_typ_2d), target :: uniface_pointers_2d(:)
    character(kind=c_char), intent(in) :: domain(*)
    character(kind=c_char,len=*), intent(in) :: interfaces(*)
    integer(kind=c_int), VALUE :: interface_count
    integer :: i

    call mui_create_uniface_multi_2dx_f(domain, interfaces, &
      interface_count)

    do i = 1, interface_count
      uniface_pointers_2d(i)%ptr = get_mui_uniface_multi_2dx_f(i)
    end do
  end subroutine create_and_get_uniface_multi_2dx_f

  !Create and access set of 2D interfaces using config from config_f_wrapper.h
  subroutine create_and_get_uniface_multi_2t_f(uniface_pointers_2d, domain, interfaces, &
    interface_count)
    use, intrinsic :: iso_c_binding
    implicit none

    type(ptr_typ_2d), target :: uniface_pointers_2d(:)
    character(kind=c_char), intent(in) :: domain(*)
    character(kind=c_char,len=*), intent(in) :: interfaces(*)
    integer(kind=c_int), VALUE :: interface_count
    integer :: i

    call mui_create_uniface_multi_2t_f(domain, interfaces, &
      interface_count)

    do i = 1, interface_count
      uniface_pointers_2d(i)%ptr = get_mui_uniface_multi_2t_f(i)
    end do
  end subroutine create_and_get_uniface_multi_2t_f

end module mui_2d_f
