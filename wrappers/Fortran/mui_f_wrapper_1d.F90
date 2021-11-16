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
!Filename: mui_f_wrapper_1d.f90
!Created: 15 September 2021
!Author: S. M. Longshaw (derived from original 3D wrapper by S. Kudo)
!Description: Fortran wrapper to create and manage 1D MUI interfaces
!             and associated sampler objects
!
!             NOTE: Any point co-ordinates are enumerated rather than assuming
!                   Cartesian form, i.e. {1, 2, 3} rather than {x, y, z}.

module mui_1d_f
  use iso_c_binding, only : c_ptr,c_char,c_int,c_float,c_double
  implicit none
  public

  type, bind(C) :: mui_point_1f_t
    real(c_float) :: point_1
  end type mui_point_1f_t

  type, bind(C) :: mui_point_1fx_t
    real(c_float) :: point_1
  end type mui_point_1fx_t

  type, bind(C) :: mui_point_1d_t
    real(c_double) :: point_1
  end type mui_point_1d_t

  type, bind(C) :: mui_point_1dx_t
    real(c_double) :: point_1
  end type mui_point_1dx_t

  type, bind(C) :: mui_point_1t_t
    real(c_double) :: point_1
  end type mui_point_1t_t

  interface
    !****************************************
    !* Create MUI interfaces                *
    !****************************************

    !1D interface with float=single and int=int32
    subroutine mui_create_uniface_1f_f(uniface,domain) bind(c)
      import :: c_ptr,c_char
      type(c_ptr), intent(out), target :: uniface(*)
      character(c_char), intent(in) :: domain(*)
    end subroutine mui_create_uniface_1f_f

    !1D interface with float=single and int=int64
    subroutine mui_create_uniface_1fx_f(uniface,domain) bind(c)
      import :: c_ptr,c_char
      type(c_ptr), intent(out), target :: uniface(*)
      character(c_char), intent(in) :: domain(*)
    end subroutine mui_create_uniface_1fx_f

    !1D interface with float=double and int=int32
    subroutine mui_create_uniface_1d_f(uniface,domain) bind(c)
      import :: c_ptr,c_char
      type(c_ptr), intent(out), target :: uniface(*)
      character(c_char), intent(in) :: domain(*)
    end subroutine mui_create_uniface_1d_f

    !1D interface with float=double and int=int64
    subroutine mui_create_uniface_1dx_f(uniface,domain) bind(c)
      import :: c_ptr,c_char
      type(c_ptr), intent(out), target :: uniface(*)
      character(c_char), intent(in) :: domain(*)
    end subroutine mui_create_uniface_1dx_f

    !1D interface using config from config_f_wrapper.h
    subroutine mui_create_uniface_1t_f(uniface,domain) bind(c)
      import :: c_ptr,c_char
      type(c_ptr), intent(out), target :: uniface(*)
      character(c_char), intent(in) :: domain(*)
    end subroutine mui_create_uniface_1t_f

    !****************************************
    !* Destroy MUI interface                *
    !****************************************

    subroutine mui_destroy_uniface_1f_f(uniface) bind(c)
      import :: c_ptr
      type(c_ptr), intent(in), value :: uniface
    end subroutine mui_destroy_uniface_1f_f

    subroutine mui_destroy_uniface_1fx_f(uniface) bind(c)
      import :: c_ptr
      type(c_ptr), intent(in), target :: uniface
    end subroutine mui_destroy_uniface_1fx_f

    subroutine mui_destroy_uniface_1d_f(uniface) bind(c)
      import :: c_ptr
      type(c_ptr), intent(in), target :: uniface
    end subroutine mui_destroy_uniface_1d_f

    subroutine mui_destroy_uniface_1dx_f(uniface) bind(c)
      import :: c_ptr
      type(c_ptr), intent(in), target :: uniface
    end subroutine mui_destroy_uniface_1dx_f

    subroutine mui_destroy_uniface_1t_f(uniface) bind(c)
      import :: c_ptr
      type(c_ptr), intent(in), target :: uniface
    end subroutine mui_destroy_uniface_1t_f

    !******************************************
    !* Create 1D spatial samplers             *
    !******************************************

    !Exact sampler

    subroutine mui_create_sampler_exact_1f_f(sampler,tolerance) bind(c)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: sampler(*)
      real(c_float), intent(in) :: tolerance
    end subroutine mui_create_sampler_exact_1f_f

    subroutine mui_create_sampler_exact_1fx_f(uniface,tolerance) bind(c)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_float), intent(in) :: tolerance
    end subroutine mui_create_sampler_exact_1fx_f

    subroutine mui_create_sampler_exact_1d_f(uniface,tolerance) bind(c)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_double), intent(in) :: tolerance
    end subroutine mui_create_sampler_exact_1d_f

    subroutine mui_create_sampler_exact_1dx_f(uniface,tolerance) bind(c)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_double), intent(in) :: tolerance
    end subroutine mui_create_sampler_exact_1dx_f

    subroutine mui_create_sampler_exact_1t_f(uniface,tolerance) bind(c)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_double), intent(in) :: tolerance
    end subroutine mui_create_sampler_exact_1t_f

    !Gauss sampler
    subroutine mui_create_sampler_gauss_1f_f(uniface,r,h) bind(c)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_float), intent(in) :: r,h
    end subroutine mui_create_sampler_gauss_1f_f

    subroutine mui_create_sampler_gauss_1fx_f(uniface,r,h) bind(c)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_float), intent(in) :: r,h
    end subroutine mui_create_sampler_gauss_1fx_f

    subroutine mui_create_sampler_gauss_1d_f(uniface,r,h) bind(c)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_double), intent(in) :: r,h
    end subroutine mui_create_sampler_gauss_1d_f

    subroutine mui_create_sampler_gauss_1dx_f(uniface,r,h) bind(c)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_double), intent(in) :: r,h
    end subroutine mui_create_sampler_gauss_1dx_f

    subroutine mui_create_sampler_gauss_1t_f(uniface,r,h) bind(c)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_double), intent(in) :: r,h
    end subroutine mui_create_sampler_gauss_1t_f

    !Moving average sampler
    subroutine mui_create_sampler_moving_average_1f_f(uniface,bbox_1) bind(c)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_float), intent(in) :: bbox_1
    end subroutine mui_create_sampler_moving_average_1f_f

    subroutine mui_create_sampler_moving_average_1fx_f(uniface,bbox_1) bind(c)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_float), intent(in) :: bbox_1
    end subroutine mui_create_sampler_moving_average_1fx_f

    subroutine mui_create_sampler_moving_average_1d_f(uniface,bbox_1) bind(c)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_double), intent(in) :: bbox_1
    end subroutine mui_create_sampler_moving_average_1d_f

    subroutine mui_create_sampler_moving_average_1dx_f(uniface,bbox_1) bind(c)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_double), intent(in) :: bbox_1
    end subroutine mui_create_sampler_moving_average_1dx_f

    subroutine mui_create_sampler_moving_average_1t_f(uniface,bbox_1) bind(c)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_double), intent(in) :: bbox_1
    end subroutine mui_create_sampler_moving_average_1t_f

    !Nearest neighbour sampler
    subroutine mui_create_sampler_nearest_neighbor_1f_f(uniface) bind(c)
      import :: c_ptr
      type(c_ptr), intent(out), target :: uniface(*)
    end subroutine mui_create_sampler_nearest_neighbor_1f_f

    subroutine mui_create_sampler_nearest_neighbor_1fx_f(uniface) bind(c)
      import :: c_ptr
      type(c_ptr), intent(out), target :: uniface(*)
    end subroutine mui_create_sampler_nearest_neighbor_1fx_f

    subroutine mui_create_sampler_nearest_neighbor_1d_f(uniface) bind(c)
      import :: c_ptr
      type(c_ptr), intent(out), target :: uniface(*)
    end subroutine mui_create_sampler_nearest_neighbor_1d_f

    subroutine mui_create_sampler_nearest_neighbor_1dx_f(uniface) bind(c)
      import :: c_ptr
      type(c_ptr), intent(out), target :: uniface(*)
    end subroutine mui_create_sampler_nearest_neighbor_1dx_f

    subroutine mui_create_sampler_nearest_neighbor_1t_f(uniface) bind(c)
      import :: c_ptr
      type(c_ptr), intent(out), target :: uniface(*)
    end subroutine mui_create_sampler_nearest_neighbor_1t_f

    !Pseudo-linear n^2 interpolation sampler
    subroutine mui_create_sampler_pseudo_n2_linear_1f_f(uniface,r) bind(c)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_float), intent(in), target :: r
    end subroutine mui_create_sampler_pseudo_n2_linear_1f_f

    subroutine mui_create_sampler_pseudo_n2_linear_1fx_f(uniface,r) bind(c)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_float), intent(in), target :: r
    end subroutine mui_create_sampler_pseudo_n2_linear_1fx_f

    subroutine mui_create_sampler_pseudo_n2_linear_1d_f(uniface,r) bind(c)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_double), intent(in), target :: r
    end subroutine mui_create_sampler_pseudo_n2_linear_1d_f

    subroutine mui_create_sampler_pseudo_n2_linear_1dx_f(uniface,r) bind(c)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_double), intent(in), target :: r
    end subroutine mui_create_sampler_pseudo_n2_linear_1dx_f

    subroutine mui_create_sampler_pseudo_n2_linear_1t_f(uniface,r) bind(c)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_double), intent(in), target :: r
    end subroutine mui_create_sampler_pseudo_n2_linear_1t_f

    !Pseudo-linear nearest neighbour interpolation sampler
    subroutine mui_create_sampler_pseudo_nearest_neighbor_1f_f(uniface,h) bind(c)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_float), intent(in), target :: h
    end subroutine mui_create_sampler_pseudo_nearest_neighbor_1f_f

    subroutine mui_create_sampler_pseudo_nearest_neighbor_1fx_f(uniface,h) bind(c)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_float), intent(in), target :: h
    end subroutine mui_create_sampler_pseudo_nearest_neighbor_1fx_f

    subroutine mui_create_sampler_pseudo_nearest_neighbor_1d_f(uniface,h) bind(c)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_double), intent(in), target :: h
    end subroutine mui_create_sampler_pseudo_nearest_neighbor_1d_f

    subroutine mui_create_sampler_pseudo_nearest_neighbor_1dx_f(uniface,h) bind(c)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_double), intent(in), target :: h
    end subroutine mui_create_sampler_pseudo_nearest_neighbor_1dx_f

    subroutine mui_create_sampler_pseudo_nearest_neighbor_1t_f(uniface,h) bind(c)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_double), intent(in), target :: h
    end subroutine mui_create_sampler_pseudo_nearest_neighbor_1t_f

    !Shepard interpolation with a quintic kernel sampler
    subroutine mui_create_sampler_shepard_quintic_1f_f(uniface,r) bind(c)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_float), intent(in), target :: r
    end subroutine mui_create_sampler_shepard_quintic_1f_f

    subroutine mui_create_sampler_shepard_quintic_1fx_f(uniface,r) bind(c)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_float), intent(in), target :: r
    end subroutine mui_create_sampler_shepard_quintic_1fx_f

    subroutine mui_create_sampler_shepard_quintic_1d_f(uniface,r) bind(c)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_double), intent(in), target :: r
    end subroutine mui_create_sampler_shepard_quintic_1d_f

    subroutine mui_create_sampler_shepard_quintic_1dx_f(uniface,r) bind(c)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_double), intent(in), target :: r
    end subroutine mui_create_sampler_shepard_quintic_1dx_f

    subroutine mui_create_sampler_shepard_quintic_1t_f(uniface,r) bind(c)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_double), intent(in), target :: r
    end subroutine mui_create_sampler_shepard_quintic_1t_f

    !SPH derived interpolation method with a quintic spline kernel sampler
    subroutine mui_create_sampler_sph_quintic_1f_f(uniface,r) bind(c)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_float), intent(in), target :: r
    end subroutine mui_create_sampler_sph_quintic_1f_f

    subroutine mui_create_sampler_sph_quintic_1fx_f(uniface,r) bind(c)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_float), intent(in), target :: r
    end subroutine mui_create_sampler_sph_quintic_1fx_f

    subroutine mui_create_sampler_sph_quintic_1d_f(uniface,r) bind(c)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_double), intent(in), target :: r
    end subroutine mui_create_sampler_sph_quintic_1d_f

    subroutine mui_create_sampler_sph_quintic_1dx_f(uniface,r) bind(c)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_double), intent(in), target :: r
    end subroutine mui_create_sampler_sph_quintic_1dx_f

    subroutine mui_create_sampler_sph_quintic_1t_f(uniface,r) bind(c)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_double), intent(in), target :: r
    end subroutine mui_create_sampler_sph_quintic_1t_f

    !Summation with a quintic kernel sampler
    subroutine mui_create_sampler_sum_quintic_1f_f(uniface,r) bind(c)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_float), intent(in), target :: r
    end subroutine mui_create_sampler_sum_quintic_1f_f

    subroutine mui_create_sampler_sum_quintic_1fx_f(uniface,r) bind(c)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_float), intent(in), target :: r
    end subroutine mui_create_sampler_sum_quintic_1fx_f

    subroutine mui_create_sampler_sum_quintic_1d_f(uniface,r) bind(c)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_double), intent(in), target :: r
    end subroutine mui_create_sampler_sum_quintic_1d_f

    subroutine mui_create_sampler_sum_quintic_1dx_f(uniface,r) bind(c)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_double), intent(in), target :: r
    end subroutine mui_create_sampler_sum_quintic_1dx_f

    subroutine mui_create_sampler_sum_quintic_1t_f(uniface,r) bind(c)
      import :: c_ptr,c_double
      type(c_ptr), intent(out), target :: uniface(*)
      real(c_double), intent(in), target :: r
    end subroutine mui_create_sampler_sum_quintic_1t_f

#ifdef USE_RBF
    !Radial Basis Function sampler
    subroutine mui_create_sampler_rbf_1f_f(uniface,r,points,points_count, &
               basis_func,conservative,polynomial,smoothFunc,readMatrix, &
               file_address,cutoff) bind(c)
      import :: c_ptr,c_int,c_float,mui_point_1f_t
      type(c_ptr), intent(out), target :: uniface(*)
      character(c_char), intent(in) :: file_address(*)
      type(mui_point_1f_t), intent(in), target :: points(*)
      integer(c_int), intent(in), target :: points_count,basis_func,conservative,polynomial,smoothFunc,readMatrix
      real(c_float), intent(in), target :: r,cutoff
    end subroutine mui_create_sampler_rbf_1f_f

    subroutine mui_create_sampler_rbf_1f_f(uniface,r,points,points_count, &
               basis_func,conservative,polynomial,smoothFunc,readMatrix, &
               file_address,cutoff) bind(c)
      import :: c_ptr,c_int,c_float,mui_point_1f_t
      type(c_ptr), intent(out), target :: uniface(*)
      character(c_char), intent(in) :: file_address(*)
      type(mui_point_1f_t), intent(in), target :: points(*)
      integer(c_int), intent(in), target :: points_count,basis_func,conservative,polynomial,smoothFunc,readMatrix
      real(c_float), intent(in), target :: r,cutoff
    end subroutine mui_create_sampler_rbf_1f_f

    subroutine mui_create_sampler_rbf_1fx_f(uniface,r,points,points_count, &
               basis_func,conservative,polynomial,smoothFunc,readMatrix, &
               file_address,cutoff) bind(c)
      import :: c_ptr,c_int,c_double,mui_point_1fx_t
      type(c_ptr), intent(out), target :: uniface(*)
      character(c_char), intent(in) :: file_address(*)
      type(mui_point_1fx_t), intent(in), target :: points(*)
      integer(c_int), intent(in), target :: points_count,basis_func,conservative,polynomial,smoothFunc,readMatrix
      real(c_double), intent(in), target :: r,cutoff
    end subroutine mui_create_sampler_rbf_1fx_f

    subroutine mui_create_sampler_rbf_1d_f(uniface,r,points,points_count, &
               basis_func,conservative,polynomial,smoothFunc,readMatrix, &
               file_address,cutoff) bind(c)
      import :: c_ptr,c_int,c_double,mui_point_1d_t
      type(c_ptr), intent(out), target :: uniface(*)
      character(c_char), intent(in) :: file_address(*)
      type(mui_point_1d_t), intent(in), target :: points(*)
      integer(c_int), intent(in), target :: points_count,basis_func,conservative,polynomial,smoothFunc,readMatrix
      real(c_double), intent(in), target :: r,cutoff
    end subroutine mui_create_sampler_rbf_1d_f

    subroutine mui_create_sampler_rbf_1dx_f(uniface,r,points,points_count, &
               basis_func,conservative,polynomial,smoothFunc,readMatrix, &
               file_address,cutoff) bind(c)
      import :: c_ptr,c_int,c_double,mui_point_1dx_t
      type(c_ptr), intent(out), target :: uniface(*)
      character(c_char), intent(in) :: file_address(*)
      type(mui_point_1dx_t), intent(in), target :: points(*)
      integer(c_int), intent(in), target :: points_count,basis_func,conservative,polynomial,smoothFunc,readMatrix
      real(c_double), intent(in), target :: r,cutoff
    end subroutine mui_create_sampler_rbf_1dx_f

    subroutine mui_create_sampler_rbf_1t_f(uniface,r,points,points_count, &
               basis_func,conservative,polynomial,smoothFunc,readMatrix, &
               file_address,cutoff) bind(c)
      import :: c_ptr,c_int,c_double,mui_point_1t_t
      type(c_ptr), intent(out), target :: uniface(*)
      character(c_char), intent(in) :: file_address(*)
      type(mui_point_1t_t), intent(in), target :: points(*)
      integer(c_int), intent(in), target :: points_count,basis_func,conservative,polynomial,smoothFunc,readMatrix
      real(c_double), intent(in), target :: r,cutoff
    end subroutine mui_create_sampler_rbf_1t_f
#endif

    subroutine mui_create_chrono_sampler_exact_1f_f(sampler,tolerance) bind(c)
      import :: c_ptr,c_float
      type(c_ptr), intent(out), target :: sampler(*)
      real(c_float),intent(in) :: tolerance
    end subroutine mui_create_chrono_sampler_exact_1f_f

    subroutine mui_push_1f_f(uniface,attr,point_1,value) bind(c)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface
      character(c_char), intent(in) :: attr(*)
      real(c_float), intent(in) :: point_1,value
    end subroutine mui_push_1f_f

    real(c_float) function mui_fetch_exact_exact_1f_f(uniface,attr,point_1,t,spatial_sampler,temporal_sampler) bind(c)
      import :: c_ptr,c_char,c_float
      type(c_ptr), intent(in), value :: uniface,spatial_sampler,temporal_sampler
      character(c_char), intent(in) :: attr(*)
      real(c_float), intent(in) :: point_1,t
    end function mui_fetch_exact_exact_1f_f

    subroutine mui_commit_1f(uniface,t) bind(c)
      import :: c_ptr,c_float
      type(c_ptr), intent(in), value :: uniface
      real(c_float), intent(in) :: t
    end subroutine mui_commit_1f

  end interface 

end module mui_1d_f
