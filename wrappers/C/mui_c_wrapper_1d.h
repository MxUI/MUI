/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library                       *
*                                                                            *
* Copyright (C) 2021 Y. H. Tang, S. Kudo, X. Bian, Z. Li, G. E. Karniadakis, *
*                    S. M. Longshaw                                          *
*                                                                            *
* This software is jointly licensed under the Apache License, Version 2.0    *
* and the GNU General Public License version 3, you may use it according     *
* to either.                                                                 *
*                                                                            *
* ** Apache License, version 2.0 **                                          *
*                                                                            *
* Licensed under the Apache License, Version 2.0 (the "License");            *
* you may not use this file except in compliance with the License.           *
* You may obtain a copy of the License at                                    *
*                                                                            *
* http://www.apache.org/licenses/LICENSE-2.0                                 *
*                                                                            *
* Unless required by applicable law or agreed to in writing, software        *
* distributed under the License is distributed on an "AS IS" BASIS,          *
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
* See the License for the specific language governing permissions and        *
* limitations under the License.                                             *
*                                                                            *
* ** GNU General Public License, version 3 **                                *
*                                                                            *
* This program is free software: you can redistribute it and/or modify       *
* it under the terms of the GNU General Public License as published by       *
* the Free Software Foundation, either version 3 of the License, or          *
* (at your option) any later version.                                        *
*                                                                            *
* This program is distributed in the hope that it will be useful,            *
* but WITHOUT ANY WARRANTY; without even the implied warranty of             *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
* GNU General Public License for more details.                               *
*                                                                            *
* You should have received a copy of the GNU General Public License          *
* along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
*****************************************************************************/

/**
 * @file mui_c_wrapper_1d.h
 * @author S. M. Longshaw (derived from original 3D wrapper by Y. H. Tang)
 * @date Jul 28, 2021
 * @brief Header for C wrapper to create and manage 1D MUI interfaces and
 *        associated sampler objects
 *
 *        NOTE: Any point co-ordinates are enumerated rather than assuming
 *              Cartesian form, i.e. {1, 2, 3} rather than {x, y, z}.
 */

#ifndef MUI_C_WRAPPER_1D_H_
#define MUI_C_WRAPPER_1D_H_

// C-defined simple struct versions of MUI point types
typedef struct mui_point_1f {
	float point_1;
} mui_point_1f;

typedef struct mui_point_1fx {
	float point_1;
} mui_point_1fx;

typedef struct mui_point_1d {
	double point_1;
} mui_point_1d;

typedef struct mui_point_1dx {
	double point_1;
} mui_point_1dx;

typedef struct mui_point_1t {
	double point_1;
} mui_point_1t;

// C access typedefs for uniface and sampler types
typedef struct mui_uniface_1f mui_uniface_1f;
typedef struct mui_uniface_1fx mui_uniface_1fx;
typedef struct mui_uniface_1d mui_uniface_1d;
typedef struct mui_uniface_1dx mui_uniface_1dx;
typedef struct mui_uniface_1t mui_uniface_1t;

typedef struct mui_sampler_exact_1f mui_sampler_exact_1f;
typedef struct mui_sampler_exact_1fx mui_sampler_exact_1fx;
typedef struct mui_sampler_exact_1d mui_sampler_exact_1d;
typedef struct mui_sampler_exact_1dx mui_sampler_exact_1dx;
typedef struct mui_sampler_exact_1t mui_sampler_exact_1t;

typedef struct mui_sampler_gauss_1f mui_sampler_gauss_1f;
typedef struct mui_sampler_gauss_1fx mui_sampler_gauss_1fx;
typedef struct mui_sampler_gauss_1d mui_sampler_gauss_1d;
typedef struct mui_sampler_gauss_1dx mui_sampler_gauss_1dx;
typedef struct mui_sampler_gauss_1t mui_sampler_gauss_1t;

typedef struct mui_sampler_moving_average_1f mui_sampler_moving_average_1f;
typedef struct mui_sampler_moving_average_1fx mui_sampler_moving_average_1fx;
typedef struct mui_sampler_moving_average_1d mui_sampler_moving_average_1d;
typedef struct mui_sampler_moving_average_1dx mui_sampler_moving_average_1dx;
typedef struct mui_sampler_moving_average_1t mui_sampler_moving_average_1t;

typedef struct mui_sampler_nearest_neighbor_1f mui_sampler_nearest_neighbor_1f;
typedef struct mui_sampler_nearest_neighbor_1fx mui_sampler_nearest_neighbor_1fx;
typedef struct mui_sampler_nearest_neighbor_1d mui_sampler_nearest_neighbor_1d;
typedef struct mui_sampler_nearest_neighbor_1dx mui_sampler_nearest_neighbor_1dx;
typedef struct mui_sampler_nearest_neighbor_1t mui_sampler_nearest_neighbor_1t;

typedef struct mui_sampler_pseudo_n2_linear_1f mui_sampler_pseudo_n2_linear_1f;
typedef struct mui_sampler_pseudo_n2_linear_1fx mui_sampler_pseudo_n2_linear_1fx;
typedef struct mui_sampler_pseudo_n2_linear_1d mui_sampler_pseudo_n2_linear_1d;
typedef struct mui_sampler_pseudo_n2_linear_1dx mui_sampler_pseudo_n2_linear_1dx;
typedef struct mui_sampler_pseudo_n2_linear_1t mui_sampler_pseudo_n2_linear_1t;

typedef struct mui_sampler_pseudo_nearest_neighbor_1f mui_sampler_pseudo_nearest_neighbor_1f;
typedef struct mui_sampler_pseudo_nearest_neighbor_1fx mui_sampler_pseudo_nearest_neighbor_1fx;
typedef struct mui_sampler_pseudo_nearest_neighbor_1d mui_sampler_pseudo_nearest_neighbor_1d;
typedef struct mui_sampler_pseudo_nearest_neighbor_1dx mui_sampler_pseudo_nearest_neighbor_1dx;
typedef struct mui_sampler_pseudo_nearest_neighbor_1t mui_sampler_pseudo_nearest_neighbor_1t;

typedef struct mui_sampler_shepard_quintic_1f mui_sampler_shepard_quintic_1f;
typedef struct mui_sampler_shepard_quintic_1fx mui_sampler_shepard_quintic_1fx;
typedef struct mui_sampler_shepard_quintic_1d mui_sampler_shepard_quintic_1d;
typedef struct mui_sampler_shepard_quintic_1dx mui_sampler_shepard_quintic_1dx;
typedef struct mui_sampler_shepard_quintic_1t mui_sampler_shepard_quintic_1t;

typedef struct mui_sampler_sph_quintic_1f mui_sampler_sph_quintic_1f;
typedef struct mui_sampler_sph_quintic_1fx mui_sampler_sph_quintic_1fx;
typedef struct mui_sampler_sph_quintic_1d mui_sampler_sph_quintic_1d;
typedef struct mui_sampler_sph_quintic_1dx mui_sampler_sph_quintic_1dx;
typedef struct mui_sampler_sph_quintic_1t mui_sampler_sph_quintic_1t;

typedef struct mui_sampler_sum_quintic_1f mui_sampler_sum_quintic_1f;
typedef struct mui_sampler_sum_quintic_1fx mui_sampler_sum_quintic_1fx;
typedef struct mui_sampler_sum_quintic_1d mui_sampler_sum_quintic_1d;
typedef struct mui_sampler_sum_quintic_1dx mui_sampler_sum_quintic_1dx;
typedef struct mui_sampler_sum_quintic_1t mui_sampler_sum_quintic_1t;

#ifdef USE_RBF
typedef struct mui_sampler_rbf_1f mui_sampler_rbf_1f;
typedef struct mui_sampler_rbf_1fx mui_sampler_rbf_1fx;
typedef struct mui_sampler_rbf_1d mui_sampler_rbf_1d;
typedef struct mui_sampler_rbf_1dx mui_sampler_rbf_1dx;
typedef struct mui_sampler_rbf_1t mui_sampler_rbf_1t;
#endif

typedef struct mui_chrono_sampler_exact_1f mui_chrono_sampler_exact_1f;
typedef struct mui_chrono_sampler_exact_1fx mui_chrono_sampler_exact_1fx;
typedef struct mui_chrono_sampler_exact_1d mui_chrono_sampler_exact_1d;
typedef struct mui_chrono_sampler_exact_1dx mui_chrono_sampler_exact_1dx;
typedef struct mui_chrono_sampler_exact_1t mui_chrono_sampler_exact_1t;

typedef struct mui_chrono_sampler_gauss_1f mui_chrono_sampler_gauss_1f;
typedef struct mui_chrono_sampler_gauss_1fx mui_chrono_sampler_gauss_1fx;
typedef struct mui_chrono_sampler_gauss_1d mui_chrono_sampler_gauss_1d;
typedef struct mui_chrono_sampler_gauss_1dx mui_chrono_sampler_gauss_1dx;
typedef struct mui_chrono_sampler_gauss_1t mui_chrono_sampler_gauss_1t;

typedef struct mui_chrono_sampler_mean_1f mui_chrono_sampler_mean_1f;
typedef struct mui_chrono_sampler_mean_1fx mui_chrono_sampler_mean_1fx;
typedef struct mui_chrono_sampler_mean_1d mui_chrono_sampler_mean_1d;
typedef struct mui_chrono_sampler_mean_1dx mui_chrono_sampler_mean_1dx;
typedef struct mui_chrono_sampler_mean_1t mui_chrono_sampler_mean_1t;

typedef struct mui_chrono_sampler_sum_1f mui_chrono_sampler_sum_1f;
typedef struct mui_chrono_sampler_sum_1fx mui_chrono_sampler_sum_1fx;
typedef struct mui_chrono_sampler_sum_1d mui_chrono_sampler_sum_1d;
typedef struct mui_chrono_sampler_sum_1dx mui_chrono_sampler_sum_1dx;
typedef struct mui_chrono_sampler_sum_1t mui_chrono_sampler_sum_1t;

// MUI single uniface creation
mui_uniface_1f* mui_create_uniface_1f(const char *URI);
mui_uniface_1fx* mui_create_uniface_1fx(const char *URI);
mui_uniface_1d* mui_create_uniface_1d(const char *URI);
mui_uniface_1dx* mui_create_uniface_1dx(const char *URI);
mui_uniface_1t* mui_create_uniface_1t(const char *URI);

// MUI multi uniface creation
mui_uniface_1f** mui_create_uniface_multi_1f( const char *domain, const char **interfaces, int interface_count );
mui_uniface_1fx** mui_create_uniface_multi_1fx( const char *domain, const char **interfaces, int interface_count );
mui_uniface_1d** mui_create_uniface_multi_1d( const char *domain, const char **interfaces, int interface_count );
mui_uniface_1dx** mui_create_uniface_multi_1dx( const char *domain, const char **interfaces, int interface_count );
mui_uniface_1t** mui_create_uniface_multi_1t( const char *domain, const char **interfaces, int interface_count );

// MUI uniface destruction
void mui_destroy_uniface_1f(mui_uniface_1f *uniface);
void mui_destroy_uniface_1fx(mui_uniface_1fx *uniface);
void mui_destroy_uniface_1d(mui_uniface_1d *uniface);
void mui_destroy_uniface_1dx(mui_uniface_1dx *uniface);
void mui_destroy_uniface_1t(mui_uniface_1t *uniface);

// MUI spatial samplers creation
mui_sampler_exact_1f* mui_create_sampler_exact_1f(float tolerance);
mui_sampler_exact_1fx* mui_create_sampler_exact_1fx(float tolerance);
mui_sampler_exact_1d* mui_create_sampler_exact_1d(double tolerance);
mui_sampler_exact_1dx* mui_create_sampler_exact_1dx(double tolerance);
mui_sampler_exact_1t* mui_create_sampler_exact_1t(double tolerance);
mui_sampler_gauss_1f* mui_create_sampler_gauss_1f(float r, float h);
mui_sampler_gauss_1fx* mui_create_sampler_gauss_1fx(float r, float h);
mui_sampler_gauss_1d* mui_create_sampler_gauss_1d(double r, double h);
mui_sampler_gauss_1dx* mui_create_sampler_gauss_1dx(double r, double h);
mui_sampler_gauss_1t* mui_create_sampler_gauss_1t(double r, double h);
mui_sampler_moving_average_1f* mui_create_sampler_moving_average_1f(float bbox_1);
mui_sampler_moving_average_1fx* mui_create_sampler_moving_average_1fx(float bbox_1);
mui_sampler_moving_average_1d* mui_create_sampler_moving_average_1d(double bbox_1);
mui_sampler_moving_average_1dx* mui_create_sampler_moving_average_1dx(double bbox_1);
mui_sampler_moving_average_1t* mui_create_sampler_moving_average_1t(double bbox_1);
mui_sampler_nearest_neighbor_1f* mui_create_sampler_nearest_neighbor_1f();
mui_sampler_nearest_neighbor_1fx* mui_create_sampler_nearest_neighbor_1fx();
mui_sampler_nearest_neighbor_1d* mui_create_sampler_nearest_neighbor_1d();
mui_sampler_nearest_neighbor_1dx* mui_create_sampler_nearest_neighbor_1dx();
mui_sampler_nearest_neighbor_1t* mui_create_sampler_nearest_neighbor_1t();
mui_sampler_pseudo_n2_linear_1f* mui_create_sampler_pseudo_n2_linear_1f(float r);
mui_sampler_pseudo_n2_linear_1fx* mui_create_sampler_pseudo_n2_linear_1fx(float r);
mui_sampler_pseudo_n2_linear_1d* mui_create_sampler_pseudo_n2_linear_1d(double r);
mui_sampler_pseudo_n2_linear_1dx* mui_create_sampler_pseudo_n2_linear_1dx(double r);
mui_sampler_pseudo_n2_linear_1t* mui_create_sampler_pseudo_n2_linear_t(double r);
mui_sampler_pseudo_nearest_neighbor_1f* mui_create_sampler_pseudo_nearest_neighbor_1f(float h);
mui_sampler_pseudo_nearest_neighbor_1fx* mui_create_sampler_pseudo_nearest_neighbor_1fx(float h);
mui_sampler_pseudo_nearest_neighbor_1d* mui_create_sampler_pseudo_nearest_neighbor_1d(double h);
mui_sampler_pseudo_nearest_neighbor_1dx* mui_create_sampler_pseudo_nearest_neighbor_1dx(double h);
mui_sampler_pseudo_nearest_neighbor_1t* mui_create_sampler_pseudo_nearest_neighbor_1t(double h);
mui_sampler_shepard_quintic_1f* mui_create_sampler_shepard_quintic_1f(float r);
mui_sampler_shepard_quintic_1fx* mui_create_sampler_shepard_quintic_1fx(float r);
mui_sampler_shepard_quintic_1d* mui_create_sampler_shepard_quintic_1d(double r);
mui_sampler_shepard_quintic_1dx* mui_create_sampler_shepard_quintic_1dx(double r);
mui_sampler_shepard_quintic_1t* mui_create_sampler_shepard_quintic_1t(double r);
mui_sampler_sph_quintic_1f* mui_create_sampler_sph_quintic_1f(float r);
mui_sampler_sph_quintic_1fx* mui_create_sampler_sph_quintic_1fx(float r);
mui_sampler_sph_quintic_1d* mui_create_sampler_sph_quintic_1d(double r);
mui_sampler_sph_quintic_1dx* mui_create_sampler_sph_quintic_1dx(double r);
mui_sampler_sph_quintic_1t* mui_create_sampler_sph_quintic_1t(double r);
mui_sampler_sum_quintic_1f* mui_create_sampler_sum_quintic_1f(float r);
mui_sampler_sum_quintic_1fx* mui_create_sampler_sum_quintic_1fx(float r);
mui_sampler_sum_quintic_1d* mui_create_sampler_sum_quintic_1d(double r);
mui_sampler_sum_quintic_1dx* mui_create_sampler_sum_quintic_1dx(double r);
mui_sampler_sum_quintic_1t* mui_create_sampler_sum_quintic_1t(double r);
#ifdef USE_RBF
mui_sampler_rbf_1f* mui_create_sampler_rbf_1f(float r, mui_point_1f *points, int points_count, int basis_func,
		int conservative, int polynomial, int smoothFunc, int readMatrix, int writeMatrix, const char *file_address,
		float cutoff, float cg_solve_tol, int cg_solve_it, int pou_size);
mui_sampler_rbf_1fx* mui_create_sampler_rbf_1fx(float r, mui_point_1fx *points, int points_count, int basis_func,
		int conservative, int polynomial, int smoothFunc, int readMatrix, int writeMatrix, const char *file_address,
		float cutoff, float cg_solve_tol, int cg_solve_it, int pou_size);
mui_sampler_rbf_1d* mui_create_sampler_rbf_1d(double r, mui_point_1d *points, int points_count, int basis_func,
		int conservative, int polynomial, int smoothFunc, int readMatrix, int writeMatrix, const char *file_address,
		double cutoff, double cg_solve_tol, int cg_solve_it, int pou_size);
mui_sampler_rbf_1dx* mui_create_sampler_rbf_1dx(double r, mui_point_1dx *points, int points_count, int basis_func,
		int conservative, int polynomial, int smoothFunc, int readMatrix, int writeMatrix, const char *file_address,
		double cutoff, double cg_solve_tol, int cg_solve_it, int pou_size);
mui_sampler_rbf_1t* mui_create_sampler_rbf_1t(double r, mui_point_1t *points, int points_count, int basis_func,
		int conservative, int polynomial, int smoothFunc, int readMatrix, int writeMatrix, const char *file_address,
		double cutoff, double cg_solve_tol, int cg_solve_it, int pou_size);
#endif

// MUI spatial samplers destruction
void mui_destroy_sampler_exact_1f(mui_sampler_exact_1f *sampler);
void mui_destroy_sampler_exact_1fx(mui_sampler_exact_1fx *sampler);
void mui_destroy_sampler_exact_1d(mui_sampler_exact_1d *sampler);
void mui_destroy_sampler_exact_1dx(mui_sampler_exact_1dx *sampler);
void mui_destroy_sampler_exact_1t(mui_sampler_exact_1t *sampler);
void mui_destroy_sampler_gauss_1f(mui_sampler_gauss_1f *sampler);
void mui_destroy_sampler_gauss_1fx(mui_sampler_gauss_1fx *sampler);
void mui_destroy_sampler_gauss_1d(mui_sampler_gauss_1d *sampler);
void mui_destroy_sampler_gauss_1dx(mui_sampler_gauss_1dx *sampler);
void mui_destroy_sampler_gauss_1t(mui_sampler_gauss_1t *sampler);
void mui_destroy_sampler_moving_average_1f(mui_sampler_moving_average_1f *sampler);
void mui_destroy_sampler_moving_average_1fx(mui_sampler_moving_average_1fx *sampler);
void mui_destroy_sampler_moving_average_1d(mui_sampler_moving_average_1d *sampler);
void mui_destroy_sampler_moving_average_1dx(mui_sampler_moving_average_1dx *sampler);
void mui_destroy_sampler_moving_average_1t(mui_sampler_moving_average_1t *sampler);
void mui_destroy_sampler_nearest_neighbor_1f(mui_sampler_nearest_neighbor_1f *sampler);
void mui_destroy_sampler_nearest_neighbor_1fx(mui_sampler_nearest_neighbor_1fx *sampler);
void mui_destroy_sampler_nearest_neighbor_1d(mui_sampler_nearest_neighbor_1d *sampler);
void mui_destroy_sampler_nearest_neighbor_1dx(mui_sampler_nearest_neighbor_1dx *sampler);
void mui_destroy_sampler_nearest_neighbor_1t(mui_sampler_nearest_neighbor_1t *sampler);
void mui_destroy_sampler_pseudo_nearest2_linear_1f(mui_sampler_pseudo_nearest_neighbor_1f *sampler);
void mui_destroy_sampler_pseudo_nearest2_linear_1fx(mui_sampler_pseudo_nearest_neighbor_1fx *sampler);
void mui_destroy_sampler_pseudo_nearest2_linear_1d(mui_sampler_pseudo_nearest_neighbor_1d *sampler);
void mui_destroy_sampler_pseudo_nearest2_linear_1dx(mui_sampler_pseudo_nearest_neighbor_1dx *sampler);
void mui_destroy_sampler_pseudo_nearest2_linear_1t(mui_sampler_pseudo_nearest_neighbor_1t *sampler);
void mui_destroy_sampler_pseudo_nearest_neighbor_1f(mui_sampler_pseudo_nearest_neighbor_1f *sampler);
void mui_destroy_sampler_pseudo_nearest_neighbor_1fx(mui_sampler_pseudo_nearest_neighbor_1fx *sampler);
void mui_destroy_sampler_pseudo_nearest_neighbor_1d(mui_sampler_pseudo_nearest_neighbor_1d *sampler);
void mui_destroy_sampler_pseudo_nearest_neighbor_1dx(mui_sampler_pseudo_nearest_neighbor_1dx *sampler);
void mui_destroy_sampler_pseudo_nearest_neighbor_1t(mui_sampler_pseudo_nearest_neighbor_1t *sampler);
void mui_destroy_sampler_shepard_quintic_1f(mui_sampler_shepard_quintic_1f *sampler);
void mui_destroy_sampler_shepard_quintic_1fx(mui_sampler_shepard_quintic_1fx *sampler);
void mui_destroy_sampler_shepard_quintic_1d(mui_sampler_shepard_quintic_1d *sampler);
void mui_destroy_sampler_shepard_quintic_1dx(mui_sampler_shepard_quintic_1dx *sampler);
void mui_destroy_sampler_shepard_quintic_1t(mui_sampler_shepard_quintic_1t *sampler);
void mui_destroy_sampler_sph_quintic_1f(mui_sampler_sph_quintic_1f *sampler);
void mui_destroy_sampler_sph_quintic_1fx(mui_sampler_sph_quintic_1fx *sampler);
void mui_destroy_sampler_sph_quintic_1d(mui_sampler_sph_quintic_1d *sampler);
void mui_destroy_sampler_sph_quintic_1dx(mui_sampler_sph_quintic_1dx *sampler);
void mui_destroy_sampler_sph_quintic_1t(mui_sampler_sph_quintic_1t *sampler);
void mui_destroy_sampler_sum_quintic_1f(mui_sampler_sum_quintic_1f *sampler);
void mui_destroy_sampler_sum_quintic_1fx(mui_sampler_sum_quintic_1fx *sampler);
void mui_destroy_sampler_sum_quintic_1d(mui_sampler_sum_quintic_1d *sampler);
void mui_destroy_sampler_sum_quintic_1dx(mui_sampler_sum_quintic_1dx *sampler);
void mui_destroy_sampler_sum_quintic_1t(mui_sampler_sum_quintic_1t *sampler);
#ifdef USE_RBF
void mui_destroy_sampler_rbf_1f(mui_sampler_rbf_1f *sampler);
void mui_destroy_sampler_rbf_1fx(mui_sampler_rbf_1fx *sampler);
void mui_destroy_sampler_rbf_1d(mui_sampler_rbf_1d *sampler);
void mui_destroy_sampler_rbf_1dx(mui_sampler_rbf_1dx *sampler);
void mui_destroy_sampler_rbf_1t(mui_sampler_rbf_1t *sampler);
#endif

// MUI temporal samplers creation
mui_chrono_sampler_exact_1f* mui_create_chrono_sampler_exact_1f(float tolerance);
mui_chrono_sampler_exact_1fx* mui_create_chrono_sampler_exact_1fx(float tolerance);
mui_chrono_sampler_exact_1d* mui_create_chrono_sampler_exact_1d(double tolerance);
mui_chrono_sampler_exact_1dx* mui_create_chrono_sampler_exact_1dx(double tolerance);
mui_chrono_sampler_exact_1t* mui_create_chrono_sampler_exact_1t(double tolerance);
mui_chrono_sampler_gauss_1f* mui_create_chrono_sampler_gauss_1f(float cutoff, float sigma);
mui_chrono_sampler_gauss_1fx* mui_create_chrono_sampler_gauss_1fx(float cutoff, float sigma);
mui_chrono_sampler_gauss_1d* mui_create_chrono_sampler_gauss_1d(double cutoff, double sigma);
mui_chrono_sampler_gauss_1dx* mui_create_chrono_sampler_gauss_1dx(double cutoff, double sigma);
mui_chrono_sampler_gauss_1t* mui_create_chrono_sampler_gauss_1t(double cutoff, double sigma);
mui_chrono_sampler_mean_1f* mui_create_chrono_sampler_mean_1f(float lower, float upper);
mui_chrono_sampler_mean_1fx* mui_create_chrono_sampler_mean_1fx(float lower, float upper);
mui_chrono_sampler_mean_1d* mui_create_chrono_sampler_mean_1d(double lower, double upper);
mui_chrono_sampler_mean_1dx* mui_create_chrono_sampler_mean_1dx(double lower, double upper);
mui_chrono_sampler_mean_1t* mui_create_chrono_sampler_mean_1t(double lower, double upper);
mui_chrono_sampler_sum_1f* mui_create_chrono_sampler_sum_1f(float lower, float upper);
mui_chrono_sampler_sum_1fx* mui_create_chrono_sampler_sum_1fx(float lower, float upper);
mui_chrono_sampler_sum_1d* mui_create_chrono_sampler_sum_1d(double lower, double upper);
mui_chrono_sampler_sum_1dx* mui_create_chrono_sampler_sum_1dx(double lower, double upper);
mui_chrono_sampler_sum_1t* mui_create_chrono_sampler_sum_1t(double lower, double upper);

// MUI temporal samplers destruction
void mui_destroy_chrono_sampler_exact_1f(mui_chrono_sampler_exact_1f *sampler);
void mui_destroy_chrono_sampler_exact_1fx(mui_chrono_sampler_exact_1fx *sampler);
void mui_destroy_chrono_sampler_exact_1d(mui_chrono_sampler_exact_1d *sampler);
void mui_destroy_chrono_sampler_exact_1dx(mui_chrono_sampler_exact_1dx *sampler);
void mui_destroy_chrono_sampler_exact_1t(mui_chrono_sampler_exact_1t *sampler);
void mui_destroy_chrono_sampler_gauss_1f(mui_chrono_sampler_gauss_1f *sampler);
void mui_destroy_chrono_sampler_gauss_1fx(mui_chrono_sampler_gauss_1fx *sampler);
void mui_destroy_chrono_sampler_gauss_1d(mui_chrono_sampler_gauss_1d *sampler);
void mui_destroy_chrono_sampler_gauss_1dx(mui_chrono_sampler_gauss_1dx *sampler);
void mui_destroy_chrono_sampler_gauss_1t(mui_chrono_sampler_gauss_1t *sampler);
void mui_destroy_chrono_sampler_mean_1f(mui_chrono_sampler_mean_1f *sampler);
void mui_destroy_chrono_sampler_mean_1fx(mui_chrono_sampler_mean_1fx *sampler);
void mui_destroy_chrono_sampler_mean_1d(mui_chrono_sampler_mean_1d *sampler);
void mui_destroy_chrono_sampler_mean_1dx(mui_chrono_sampler_mean_1dx *sampler);
void mui_destroy_chrono_sampler_mean_1t(mui_chrono_sampler_mean_1t *sampler);
void mui_destroy_chrono_sampler_sum_1f(mui_chrono_sampler_sum_1f *sampler);
void mui_destroy_chrono_sampler_sum_1fx(mui_chrono_sampler_sum_1fx *sampler);
void mui_destroy_chrono_sampler_sum_1d(mui_chrono_sampler_sum_1d *sampler);
void mui_destroy_chrono_sampler_sum_1dx(mui_chrono_sampler_sum_1dx *sampler);
void mui_destroy_chrono_sampler_sum_1t(mui_chrono_sampler_sum_1t *sampler);

// MUI functions for data push
void mui_push_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float value);
void mui_push_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float value);
void mui_push_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double value);
void mui_push_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double value);
void mui_push_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double value);
void mui_push_1f_param(mui_uniface_1f *uniface, const char *attr, float value);
void mui_push_1fx_param(mui_uniface_1fx *uniface, const char *attr, float value);
void mui_push_1d_param(mui_uniface_1d *uniface, const char *attr, double value);
void mui_push_1dx_param(mui_uniface_1dx *uniface, const char *attr, double value);
void mui_push_1t_param(mui_uniface_1t *uniface, const char *attr, double value);

// MUI functions for data commit
void mui_commit_1f(mui_uniface_1f *uniface, float t);
void mui_commit_1fx(mui_uniface_1fx *uniface, float t);
void mui_commit_1d(mui_uniface_1d *uniface, double t);
void mui_commit_1dx(mui_uniface_1dx *uniface, double t);
void mui_commit_1t(mui_uniface_1t *uniface, double t);
void mui_commit_1f_pair(mui_uniface_1f *uniface, float t_1, float t_2);
void mui_commit_1fx_pair(mui_uniface_1fx *uniface, float t_1, float t_2);
void mui_commit_1d_pair(mui_uniface_1d *uniface, double t_1, double t_2);
void mui_commit_1dx_pair(mui_uniface_1dx *uniface, double t_1, double t_2);
void mui_commit_1t_pair(mui_uniface_1t *uniface, double t_1, double t_2);

// MUI functions for data forecast
void mui_forecast_1f(mui_uniface_1f *uniface, float t);
void mui_forecast_1fx(mui_uniface_1fx *uniface, float t);
void mui_forecast_1d(mui_uniface_1d *uniface, double t);
void mui_forecast_1dx(mui_uniface_1dx *uniface, double t);
void mui_forecast_1t(mui_uniface_1t *uniface, double t);
void mui_forecast_1f_pair(mui_uniface_1f *uniface, float t_1, float t_2);
void mui_forecast_1fx_pair(mui_uniface_1fx *uniface, float t_1, float t_2);
void mui_forecast_1d_pair(mui_uniface_1d *uniface, double t_1, double t_2);
void mui_forecast_1dx_pair(mui_uniface_1dx *uniface, double t_1, double t_2);
void mui_forecast_1t_pair(mui_uniface_1t *uniface, double t_1, double t_2);

// MUI functions for data fetch
float mui_fetch_exact_exact_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_exact_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler);
float mui_fetch_exact_exact_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t,
		mui_sampler_exact_1fx *spatial_sampler, mui_chrono_sampler_exact_1fx *temporal_sampler);
double mui_fetch_exact_exact_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_exact_1d *spatial_sampler, mui_chrono_sampler_exact_1d *temporal_sampler);
double mui_fetch_exact_exact_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t,
		mui_sampler_exact_1dx *spatial_sampler, mui_chrono_sampler_exact_1dx *temporal_sampler);
double mui_fetch_exact_exact_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_exact_1t *spatial_sampler, mui_chrono_sampler_exact_1t *temporal_sampler);
float mui_fetch_exact_gauss_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_exact_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler);
float mui_fetch_exact_gauss_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t,
		mui_sampler_exact_1fx *spatial_sampler, mui_chrono_sampler_gauss_1fx *temporal_sampler);
double mui_fetch_exact_gauss_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_exact_1d *spatial_sampler, mui_chrono_sampler_gauss_1d *temporal_sampler);
double mui_fetch_exact_gauss_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t,
		mui_sampler_exact_1dx *spatial_sampler, mui_chrono_sampler_gauss_1dx *temporal_sampler);
double mui_fetch_exact_gauss_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_exact_1t *spatial_sampler, mui_chrono_sampler_gauss_1t *temporal_sampler);
float mui_fetch_exact_mean_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_exact_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler);
float mui_fetch_exact_mean_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t,
		mui_sampler_exact_1fx *spatial_sampler, mui_chrono_sampler_mean_1fx *temporal_sampler);
double mui_fetch_exact_mean_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_exact_1d *spatial_sampler, mui_chrono_sampler_mean_1d *temporal_sampler);
double mui_fetch_exact_mean_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t,
		mui_sampler_exact_1dx *spatial_sampler, mui_chrono_sampler_mean_1dx *temporal_sampler);
double mui_fetch_exact_mean_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_exact_1t *spatial_sampler, mui_chrono_sampler_mean_1t *temporal_sampler);
float mui_fetch_exact_sum_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_exact_1f *spatial_sampler, mui_chrono_sampler_sum_1f *temporal_sampler);
float mui_fetch_exact_sum_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t,
		mui_sampler_exact_1fx *spatial_sampler, mui_chrono_sampler_sum_1fx *temporal_sampler);
double mui_fetch_exact_sum_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_exact_1d *spatial_sampler, mui_chrono_sampler_sum_1d *temporal_sampler);
double mui_fetch_exact_sum_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t,
		mui_sampler_exact_1dx *spatial_sampler, mui_chrono_sampler_sum_1dx *temporal_sampler);
double mui_fetch_exact_sum_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_exact_1t *spatial_sampler, mui_chrono_sampler_sum_1t *temporal_sampler);
float mui_fetch_gauss_exact_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_gauss_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler);
float mui_fetch_gauss_exact_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t,
		mui_sampler_gauss_1fx *spatial_sampler, mui_chrono_sampler_exact_1fx *temporal_sampler);
double mui_fetch_gauss_exact_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_gauss_1d *spatial_sampler, mui_chrono_sampler_exact_1d *temporal_sampler);
double mui_fetch_gauss_exact_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t,
		mui_sampler_gauss_1dx *spatial_sampler, mui_chrono_sampler_exact_1dx *temporal_sampler);
double mui_fetch_gauss_exact_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_gauss_1t *spatial_sampler, mui_chrono_sampler_exact_1t *temporal_sampler);
float mui_fetch_gauss_gauss_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_gauss_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler);
float mui_fetch_gauss_gauss_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t,
		mui_sampler_gauss_1fx *spatial_sampler, mui_chrono_sampler_gauss_1fx *temporal_sampler);
double mui_fetch_gauss_gauss_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_gauss_1d *spatial_sampler, mui_chrono_sampler_gauss_1d *temporal_sampler);
double mui_fetch_gauss_gauss_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t,
		mui_sampler_gauss_1dx *spatial_sampler, mui_chrono_sampler_gauss_1dx *temporal_sampler);
double mui_fetch_gauss_gauss_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_gauss_1t *spatial_sampler, mui_chrono_sampler_gauss_1t *temporal_sampler);
float mui_fetch_gauss_mean_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_gauss_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler);
float mui_fetch_gauss_mean_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t,
		mui_sampler_gauss_1fx *spatial_sampler, mui_chrono_sampler_mean_1fx *temporal_sampler);
double mui_fetch_gauss_mean_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_gauss_1d *spatial_sampler, mui_chrono_sampler_mean_1d *temporal_sampler);
double mui_fetch_gauss_mean_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t,
		mui_sampler_gauss_1dx *spatial_sampler, mui_chrono_sampler_mean_1dx *temporal_sampler);
double mui_fetch_gauss_mean_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_gauss_1t *spatial_sampler, mui_chrono_sampler_mean_1t *temporal_sampler);
float mui_fetch_moving_average_exact_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_moving_average_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler);
float mui_fetch_moving_average_exact_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t,
		mui_sampler_moving_average_1fx *spatial_sampler, mui_chrono_sampler_exact_1fx *temporal_sampler);
double mui_fetch_moving_average_exact_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_moving_average_1d *spatial_sampler, mui_chrono_sampler_exact_1d *temporal_sampler);
double mui_fetch_moving_average_exact_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t,
		mui_sampler_moving_average_1dx *spatial_sampler, mui_chrono_sampler_exact_1dx *temporal_sampler);
double mui_fetch_moving_average_exact_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_moving_average_1t *spatial_sampler, mui_chrono_sampler_exact_1t *temporal_sampler);
float mui_fetch_moving_average_gauss_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_moving_average_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler);
float mui_fetch_moving_average_gauss_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t,
		mui_sampler_moving_average_1fx *spatial_sampler, mui_chrono_sampler_gauss_1fx *temporal_sampler);
double mui_fetch_moving_average_gauss_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_moving_average_1d *spatial_sampler, mui_chrono_sampler_gauss_1d *temporal_sampler);
double mui_fetch_moving_average_gauss_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t,
		mui_sampler_moving_average_1dx *spatial_sampler, mui_chrono_sampler_gauss_1dx *temporal_sampler);
double mui_fetch_moving_average_gauss_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_moving_average_1t *spatial_sampler, mui_chrono_sampler_gauss_1t *temporal_sampler);
float mui_fetch_moving_average_mean_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_moving_average_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler);
float mui_fetch_moving_average_mean_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t,
		mui_sampler_moving_average_1fx *spatial_sampler, mui_chrono_sampler_mean_1fx *temporal_sampler);
double mui_fetch_moving_average_mean_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_moving_average_1d *spatial_sampler, mui_chrono_sampler_mean_1d *temporal_sampler);
double mui_fetch_moving_average_mean_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t,
		mui_sampler_moving_average_1dx *spatial_sampler, mui_chrono_sampler_mean_1dx *temporal_sampler);
double mui_fetch_moving_average_mean_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_moving_average_1t *spatial_sampler, mui_chrono_sampler_mean_1t *temporal_sampler);
float mui_fetch_moving_average_sum_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_moving_average_1f *spatial_sampler, mui_chrono_sampler_sum_1f *temporal_sampler);
float mui_fetch_moving_average_sum_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t,
		mui_sampler_moving_average_1fx *spatial_sampler, mui_chrono_sampler_sum_1fx *temporal_sampler);
double mui_fetch_moving_average_sum_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_moving_average_1d *spatial_sampler, mui_chrono_sampler_sum_1d *temporal_sampler);
double mui_fetch_moving_average_sum_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t,
		mui_sampler_moving_average_1dx *spatial_sampler, mui_chrono_sampler_sum_1dx *temporal_sampler);
double mui_fetch_moving_average_sum_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_moving_average_1t *spatial_sampler, mui_chrono_sampler_sum_1t *temporal_sampler);
float mui_fetch_nearest_neighbor_exact_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_nearest_neighbor_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler);
float mui_fetch_nearest_neighbor_exact_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t,
		mui_sampler_nearest_neighbor_1fx *spatial_sampler, mui_chrono_sampler_exact_1fx *temporal_sampler);
double mui_fetch_nearest_neighbor_exact_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_nearest_neighbor_1d *spatial_sampler, mui_chrono_sampler_exact_1d *temporal_sampler);
double mui_fetch_nearest_neighbor_exact_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t,
		mui_sampler_nearest_neighbor_1dx *spatial_sampler, mui_chrono_sampler_exact_1dx *temporal_sampler);
double mui_fetch_nearest_neighbor_exact_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_nearest_neighbor_1t *spatial_sampler, mui_chrono_sampler_exact_1t *temporal_sampler);
float mui_fetch_nearest_neighbor_gauss_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_nearest_neighbor_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler);
float mui_fetch_nearest_neighbor_gauss_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t,
		mui_sampler_nearest_neighbor_1fx *spatial_sampler, mui_chrono_sampler_gauss_1fx *temporal_sampler);
double mui_fetch_nearest_neighbor_gauss_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_nearest_neighbor_1d *spatial_sampler, mui_chrono_sampler_gauss_1d *temporal_sampler);
double mui_fetch_nearest_neighbor_gauss_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t,
		mui_sampler_nearest_neighbor_1dx *spatial_sampler, mui_chrono_sampler_gauss_1dx *temporal_sampler);
double mui_fetch_nearest_neighbor_gauss_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_nearest_neighbor_1t *spatial_sampler, mui_chrono_sampler_gauss_1t *temporal_sampler);
float mui_fetch_nearest_neighbor_mean_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_nearest_neighbor_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler);
float mui_fetch_nearest_neighbor_mean_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t,
		mui_sampler_nearest_neighbor_1fx *spatial_sampler, mui_chrono_sampler_mean_1fx *temporal_sampler);
double mui_fetch_nearest_neighbor_mean_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_nearest_neighbor_1d *spatial_sampler, mui_chrono_sampler_mean_1d *temporal_sampler);
double mui_fetch_nearest_neighbor_mean_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t,
		mui_sampler_nearest_neighbor_1dx *spatial_sampler, mui_chrono_sampler_mean_1dx *temporal_sampler);
double mui_fetch_nearest_neighbor_mean_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_nearest_neighbor_1t *spatial_sampler, mui_chrono_sampler_mean_1t *temporal_sampler);
float mui_fetch_nearest_neighbor_sum_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_nearest_neighbor_1f *spatial_sampler, mui_chrono_sampler_sum_1f *temporal_sampler);
float mui_fetch_nearest_neighbor_sum_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t,
		mui_sampler_nearest_neighbor_1fx *spatial_sampler, mui_chrono_sampler_sum_1fx *temporal_sampler);
double mui_fetch_nearest_neighbor_sum_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_nearest_neighbor_1d *spatial_sampler, mui_chrono_sampler_sum_1d *temporal_sampler);
double mui_fetch_nearest_neighbor_sum_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t,
		mui_sampler_nearest_neighbor_1dx *spatial_sampler, mui_chrono_sampler_sum_1dx *temporal_sampler);
double mui_fetch_nearest_neighbor_sum_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_nearest_neighbor_1t *spatial_sampler, mui_chrono_sampler_sum_1t *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_exact_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_pseudo_nearest_neighbor_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_exact_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point,
		float t, mui_sampler_pseudo_nearest_neighbor_1fx *spatial_sampler,
		mui_chrono_sampler_exact_1fx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_exact_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point,
		double t, mui_sampler_pseudo_nearest_neighbor_1d *spatial_sampler,
		mui_chrono_sampler_exact_1d *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_exact_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point,
		double t, mui_sampler_pseudo_nearest_neighbor_1dx *spatial_sampler,
		mui_chrono_sampler_exact_1dx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_exact_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point,
		double t, mui_sampler_pseudo_nearest_neighbor_1t *spatial_sampler,
		mui_chrono_sampler_exact_1t *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_gauss_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_pseudo_nearest_neighbor_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_gauss_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point,
		float t, mui_sampler_pseudo_nearest_neighbor_1fx *spatial_sampler,
		mui_chrono_sampler_gauss_1fx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_gauss_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point,
		double t, mui_sampler_pseudo_nearest_neighbor_1d *spatial_sampler,
		mui_chrono_sampler_gauss_1d *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_gauss_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point,
		double t, mui_sampler_pseudo_nearest_neighbor_1dx *spatial_sampler,
		mui_chrono_sampler_gauss_1dx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_gauss_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point,
		double t, mui_sampler_pseudo_nearest_neighbor_1t *spatial_sampler,
		mui_chrono_sampler_gauss_1t *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_mean_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_pseudo_nearest_neighbor_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_mean_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point,
		float t, mui_sampler_pseudo_nearest_neighbor_1fx *spatial_sampler,
		mui_chrono_sampler_mean_1fx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_mean_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point,
		double t, mui_sampler_pseudo_nearest_neighbor_1d *spatial_sampler,
		mui_chrono_sampler_mean_1d *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_mean_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point,
		double t, mui_sampler_pseudo_nearest_neighbor_1dx *spatial_sampler,
		mui_chrono_sampler_mean_1dx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_mean_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point,
		double t, mui_sampler_pseudo_nearest_neighbor_1t *spatial_sampler,
		mui_chrono_sampler_mean_1t *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_sum_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_pseudo_nearest_neighbor_1f *spatial_sampler, mui_chrono_sampler_sum_1f *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_sum_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point,
		float t, mui_sampler_pseudo_nearest_neighbor_1fx *spatial_sampler,
		mui_chrono_sampler_sum_1fx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_sum_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_pseudo_nearest_neighbor_1d *spatial_sampler, mui_chrono_sampler_sum_1d *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_sum_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point,
		double t, mui_sampler_pseudo_nearest_neighbor_1dx *spatial_sampler,
		mui_chrono_sampler_sum_1dx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_sum_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_pseudo_nearest_neighbor_1t *spatial_sampler, mui_chrono_sampler_sum_1t *temporal_sampler);
float mui_fetch_shepard_quintic_exact_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_shepard_quintic_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler);
float mui_fetch_shepard_quintic_exact_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t,
		mui_sampler_shepard_quintic_1fx *spatial_sampler, mui_chrono_sampler_exact_1fx *temporal_sampler);
double mui_fetch_shepard_quintic_exact_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_shepard_quintic_1d *spatial_sampler, mui_chrono_sampler_exact_1d *temporal_sampler);
double mui_fetch_shepard_quintic_exact_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t,
		mui_sampler_shepard_quintic_1dx *spatial_sampler, mui_chrono_sampler_exact_1dx *temporal_sampler);
double mui_fetch_shepard_quintic_exact_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_shepard_quintic_1t *spatial_sampler, mui_chrono_sampler_exact_1t *temporal_sampler);
float mui_fetch_shepard_quintic_gauss_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_shepard_quintic_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler);
float mui_fetch_shepard_quintic_gauss_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t,
		mui_sampler_shepard_quintic_1fx *spatial_sampler, mui_chrono_sampler_gauss_1fx *temporal_sampler);
double mui_fetch_shepard_quintic_gauss_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_shepard_quintic_1d *spatial_sampler, mui_chrono_sampler_gauss_1d *temporal_sampler);
double mui_fetch_shepard_quintic_gauss_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t,
		mui_sampler_shepard_quintic_1dx *spatial_sampler, mui_chrono_sampler_gauss_1dx *temporal_sampler);
double mui_fetch_shepard_quintic_gauss_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_shepard_quintic_1t *spatial_sampler, mui_chrono_sampler_gauss_1t *temporal_sampler);
float mui_fetch_shepard_quintic_mean_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_shepard_quintic_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler);
float mui_fetch_shepard_quintic_mean_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t,
		mui_sampler_shepard_quintic_1fx *spatial_sampler, mui_chrono_sampler_mean_1fx *temporal_sampler);
double mui_fetch_shepard_quintic_mean_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_shepard_quintic_1d *spatial_sampler, mui_chrono_sampler_mean_1d *temporal_sampler);
double mui_fetch_shepard_quintic_mean_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t,
		mui_sampler_shepard_quintic_1dx *spatial_sampler, mui_chrono_sampler_mean_1dx *temporal_sampler);
double mui_fetch_shepard_quintic_mean_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_shepard_quintic_1t *spatial_sampler, mui_chrono_sampler_mean_1t *temporal_sampler);
float mui_fetch_shepard_quintic_sum_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_shepard_quintic_1f *spatial_sampler, mui_chrono_sampler_sum_1f *temporal_sampler);
float mui_fetch_shepard_quintic_sum_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t,
		mui_sampler_shepard_quintic_1fx *spatial_sampler, mui_chrono_sampler_sum_1fx *temporal_sampler);
double mui_fetch_shepard_quintic_sum_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_shepard_quintic_1d *spatial_sampler, mui_chrono_sampler_sum_1d *temporal_sampler);
double mui_fetch_shepard_quintic_sum_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t,
		mui_sampler_shepard_quintic_1dx *spatial_sampler, mui_chrono_sampler_sum_1dx *temporal_sampler);
double mui_fetch_shepard_quintic_sum_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_shepard_quintic_1t *spatial_sampler, mui_chrono_sampler_sum_1t *temporal_sampler);
float mui_fetch_sph_quintic_exact_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_sph_quintic_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler);
float mui_fetch_sph_quintic_exact_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t,
		mui_sampler_sph_quintic_1fx *spatial_sampler, mui_chrono_sampler_exact_1fx *temporal_sampler);
double mui_fetch_sph_quintic_exact_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_sph_quintic_1d *spatial_sampler, mui_chrono_sampler_exact_1d *temporal_sampler);
double mui_fetch_sph_quintic_exact_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t,
		mui_sampler_sph_quintic_1dx *spatial_sampler, mui_chrono_sampler_exact_1dx *temporal_sampler);
double mui_fetch_sph_quintic_exact_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_sph_quintic_1t *spatial_sampler, mui_chrono_sampler_exact_1t *temporal_sampler);
float mui_fetch_sph_quintic_gauss_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_sph_quintic_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler);
float mui_fetch_sph_quintic_gauss_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t,
		mui_sampler_sph_quintic_1fx *spatial_sampler, mui_chrono_sampler_gauss_1fx *temporal_sampler);
double mui_fetch_sph_quintic_gauss_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_sph_quintic_1d *spatial_sampler, mui_chrono_sampler_gauss_1d *temporal_sampler);
double mui_fetch_sph_quintic_gauss_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t,
		mui_sampler_sph_quintic_1dx *spatial_sampler, mui_chrono_sampler_gauss_1dx *temporal_sampler);
double mui_fetch_sph_quintic_gauss_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_sph_quintic_1t *spatial_sampler, mui_chrono_sampler_gauss_1t *temporal_sampler);
float mui_fetch_sph_quintic_mean_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_sph_quintic_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler);
float mui_fetch_sph_quintic_mean_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t,
		mui_sampler_sph_quintic_1fx *spatial_sampler, mui_chrono_sampler_mean_1fx *temporal_sampler);
double mui_fetch_sph_quintic_mean_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_sph_quintic_1d *spatial_sampler, mui_chrono_sampler_mean_1d *temporal_sampler);
double mui_fetch_sph_quintic_mean_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t,
		mui_sampler_sph_quintic_1dx *spatial_sampler, mui_chrono_sampler_mean_1dx *temporal_sampler);
double mui_fetch_sph_quintic_mean_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_sph_quintic_1t *spatial_sampler, mui_chrono_sampler_mean_1t *temporal_sampler);
float mui_fetch_sph_quintic_sum_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_sph_quintic_1f *spatial_sampler, mui_chrono_sampler_sum_1f *temporal_sampler);
float mui_fetch_sph_quintic_sum_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t,
		mui_sampler_sph_quintic_1fx *spatial_sampler, mui_chrono_sampler_sum_1fx *temporal_sampler);
double mui_fetch_sph_quintic_sum_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_sph_quintic_1d *spatial_sampler, mui_chrono_sampler_sum_1d *temporal_sampler);
double mui_fetch_sph_quintic_sum_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t,
		mui_sampler_sph_quintic_1dx *spatial_sampler, mui_chrono_sampler_sum_1dx *temporal_sampler);
double mui_fetch_sph_quintic_sum_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_sph_quintic_1t *spatial_sampler, mui_chrono_sampler_sum_1t *temporal_sampler);
float mui_fetch_sum_quintic_exact_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_sum_quintic_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler);
float mui_fetch_sum_quintic_exact_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t,
		mui_sampler_sum_quintic_1fx *spatial_sampler, mui_chrono_sampler_exact_1fx *temporal_sampler);
double mui_fetch_sum_quintic_exact_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_sum_quintic_1d *spatial_sampler, mui_chrono_sampler_exact_1d *temporal_sampler);
double mui_fetch_sum_quintic_exact_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t,
		mui_sampler_sum_quintic_1dx *spatial_sampler, mui_chrono_sampler_exact_1dx *temporal_sampler);
double mui_fetch_sum_quintic_exact_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_sum_quintic_1t *spatial_sampler, mui_chrono_sampler_exact_1t *temporal_sampler);
float mui_fetch_sum_quintic_gauss_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_sum_quintic_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler);
float mui_fetch_sum_quintic_gauss_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t,
		mui_sampler_sum_quintic_1fx *spatial_sampler, mui_chrono_sampler_gauss_1fx *temporal_sampler);
double mui_fetch_sum_quintic_gauss_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_sum_quintic_1d *spatial_sampler, mui_chrono_sampler_gauss_1d *temporal_sampler);
double mui_fetch_sum_quintic_gauss_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t,
		mui_sampler_sum_quintic_1dx *spatial_sampler, mui_chrono_sampler_gauss_1dx *temporal_sampler);
double mui_fetch_sum_quintic_gauss_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_sum_quintic_1t *spatial_sampler, mui_chrono_sampler_gauss_1t *temporal_sampler);
float mui_fetch_sum_quintic_mean_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_sum_quintic_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler);
float mui_fetch_sum_quintic_mean_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t,
		mui_sampler_sum_quintic_1fx *spatial_sampler, mui_chrono_sampler_mean_1fx *temporal_sampler);
double mui_fetch_sum_quintic_mean_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_sum_quintic_1d *spatial_sampler, mui_chrono_sampler_mean_1d *temporal_sampler);
double mui_fetch_sum_quintic_mean_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t,
		mui_sampler_sum_quintic_1dx *spatial_sampler, mui_chrono_sampler_mean_1dx *temporal_sampler);
double mui_fetch_sum_quintic_mean_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_sum_quintic_1t *spatial_sampler, mui_chrono_sampler_mean_1t *temporal_sampler);
float mui_fetch_sum_quintic_sum_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_sum_quintic_1f *spatial_sampler, mui_chrono_sampler_sum_1f *temporal_sampler);
float mui_fetch_sum_quintic_sum_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t,
		mui_sampler_sum_quintic_1fx *spatial_sampler, mui_chrono_sampler_sum_1fx *temporal_sampler);
double mui_fetch_sum_quintic_sum_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_sum_quintic_1d *spatial_sampler, mui_chrono_sampler_sum_1d *temporal_sampler);
double mui_fetch_sum_quintic_sum_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t,
		mui_sampler_sum_quintic_1dx *spatial_sampler, mui_chrono_sampler_sum_1dx *temporal_sampler);
double mui_fetch_sum_quintic_sum_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_sum_quintic_1t *spatial_sampler, mui_chrono_sampler_sum_1t *temporal_sampler);
#ifdef USE_RBF
float mui_fetch_rbf_exact_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_rbf_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler);
float mui_fetch_rbf_exact_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t,
		mui_sampler_rbf_1fx *spatial_sampler, mui_chrono_sampler_exact_1fx *temporal_sampler);
double mui_fetch_rbf_exact_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_rbf_1d *spatial_sampler, mui_chrono_sampler_exact_1d *temporal_sampler);
double mui_fetch_rbf_exact_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t,
		mui_sampler_rbf_1dx *spatial_sampler, mui_chrono_sampler_exact_1dx *temporal_sampler);
double mui_fetch_rbf_exact_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_rbf_1t *spatial_sampler, mui_chrono_sampler_exact_1t *temporal_sampler);
float mui_fetch_rbf_gauss_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_rbf_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler);
float mui_fetch_rbf_gauss_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t,
		mui_sampler_rbf_1fx *spatial_sampler, mui_chrono_sampler_gauss_1fx *temporal_sampler);
double mui_fetch_rbf_gauss_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_rbf_1d *spatial_sampler, mui_chrono_sampler_gauss_1d *temporal_sampler);
double mui_fetch_rbf_gauss_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t,
		mui_sampler_rbf_1dx *spatial_sampler, mui_chrono_sampler_gauss_1dx *temporal_sampler);
double mui_fetch_rbf_gauss_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_rbf_1t *spatial_sampler, mui_chrono_sampler_gauss_1t *temporal_sampler);
float mui_fetch_rbf_mean_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_rbf_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler);
float mui_fetch_rbf_mean_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t,
		mui_sampler_rbf_1fx *spatial_sampler, mui_chrono_sampler_mean_1fx *temporal_sampler);
double mui_fetch_rbf_mean_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_rbf_1d *spatial_sampler, mui_chrono_sampler_mean_1d *temporal_sampler);
double mui_fetch_rbf_mean_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t,
		mui_sampler_rbf_1dx *spatial_sampler, mui_chrono_sampler_mean_1dx *temporal_sampler);
double mui_fetch_rbf_mean_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_rbf_1t *spatial_sampler, mui_chrono_sampler_mean_1t *temporal_sampler);
float mui_fetch_rbf_sum_1f(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t,
		mui_sampler_rbf_1f *spatial_sampler, mui_chrono_sampler_sum_1f *temporal_sampler);
float mui_fetch_rbf_sum_1fx(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t,
		mui_sampler_rbf_1fx *spatial_sampler, mui_chrono_sampler_sum_1fx *temporal_sampler);
double mui_fetch_rbf_sum_1d(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t,
		mui_sampler_rbf_1d *spatial_sampler, mui_chrono_sampler_sum_1d *temporal_sampler);
double mui_fetch_rbf_sum_1dx(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t,
		mui_sampler_rbf_1dx *spatial_sampler, mui_chrono_sampler_sum_1dx *temporal_sampler);
double mui_fetch_rbf_sum_1t(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t,
		mui_sampler_rbf_1t *spatial_sampler, mui_chrono_sampler_sum_1t *temporal_sampler);
#endif

/********************************************************
 * MUI functions for 1D data fetch using two time values *
 *********************************************************/

float mui_fetch_exact_exact_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t_1, float t_2,
		mui_sampler_exact_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler);
float mui_fetch_exact_exact_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t_1,
		float t_2, mui_sampler_exact_1fx *spatial_sampler, mui_chrono_sampler_exact_1fx *temporal_sampler);
double mui_fetch_exact_exact_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t_1,
		double t_2, mui_sampler_exact_1d *spatial_sampler, mui_chrono_sampler_exact_1d *temporal_sampler);
double mui_fetch_exact_exact_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t_1,
		double t_2, mui_sampler_exact_1dx *spatial_sampler, mui_chrono_sampler_exact_1dx *temporal_sampler);
double mui_fetch_exact_exact_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t_1,
		double t_2, mui_sampler_exact_1t *spatial_sampler, mui_chrono_sampler_exact_1t *temporal_sampler);
float mui_fetch_exact_gauss_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t_1, float t_2,
		mui_sampler_exact_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler);
float mui_fetch_exact_gauss_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t_1,
		float t_2, mui_sampler_exact_1fx *spatial_sampler, mui_chrono_sampler_gauss_1fx *temporal_sampler);
double mui_fetch_exact_gauss_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t_1,
		double t_2, mui_sampler_exact_1d *spatial_sampler, mui_chrono_sampler_gauss_1d *temporal_sampler);
double mui_fetch_exact_gauss_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t_1,
		double t_2, mui_sampler_exact_1dx *spatial_sampler, mui_chrono_sampler_gauss_1dx *temporal_sampler);
double mui_fetch_exact_gauss_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t_1,
		double t_2, mui_sampler_exact_1t *spatial_sampler, mui_chrono_sampler_gauss_1t *temporal_sampler);
float mui_fetch_exact_mean_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t_1, float t_2,
		mui_sampler_exact_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler);
float mui_fetch_exact_mean_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t_1,
		float t_2, mui_sampler_exact_1fx *spatial_sampler, mui_chrono_sampler_mean_1fx *temporal_sampler);
double mui_fetch_exact_mean_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t_1,
		double t_2, mui_sampler_exact_1d *spatial_sampler, mui_chrono_sampler_mean_1d *temporal_sampler);
double mui_fetch_exact_mean_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t_1,
		double t_2, mui_sampler_exact_1dx *spatial_sampler, mui_chrono_sampler_mean_1dx *temporal_sampler);
double mui_fetch_exact_mean_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t_1,
		double t_2, mui_sampler_exact_1t *spatial_sampler, mui_chrono_sampler_mean_1t *temporal_sampler);
float mui_fetch_exact_sum_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t_1, float t_2,
		mui_sampler_exact_1f *spatial_sampler, mui_chrono_sampler_sum_1f *temporal_sampler);
float mui_fetch_exact_sum_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t_1,
		float t_2, mui_sampler_exact_1fx *spatial_sampler, mui_chrono_sampler_sum_1fx *temporal_sampler);
double mui_fetch_exact_sum_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t_1,
		double t_2, mui_sampler_exact_1d *spatial_sampler, mui_chrono_sampler_sum_1d *temporal_sampler);
double mui_fetch_exact_sum_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t_1,
		double t_2, mui_sampler_exact_1dx *spatial_sampler, mui_chrono_sampler_sum_1dx *temporal_sampler);
double mui_fetch_exact_sum_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t_1,
		double t_2, mui_sampler_exact_1t *spatial_sampler, mui_chrono_sampler_sum_1t *temporal_sampler);
float mui_fetch_gauss_exact_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t_1, float t_2,
		mui_sampler_gauss_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler);
float mui_fetch_gauss_exact_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t_1,
		float t_2, mui_sampler_gauss_1fx *spatial_sampler, mui_chrono_sampler_exact_1fx *temporal_sampler);
double mui_fetch_gauss_exact_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t_1,
		double t_2, mui_sampler_gauss_1d *spatial_sampler, mui_chrono_sampler_exact_1d *temporal_sampler);
double mui_fetch_gauss_exact_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t_1,
		double t_2, mui_sampler_gauss_1dx *spatial_sampler, mui_chrono_sampler_exact_1dx *temporal_sampler);
double mui_fetch_gauss_exact_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t_1,
		double t_2, mui_sampler_gauss_1t *spatial_sampler, mui_chrono_sampler_exact_1t *temporal_sampler);
float mui_fetch_gauss_gauss_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t_1, float t_2,
		mui_sampler_gauss_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler);
float mui_fetch_gauss_gauss_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t_1,
		float t_2, mui_sampler_gauss_1fx *spatial_sampler, mui_chrono_sampler_gauss_1fx *temporal_sampler);
double mui_fetch_gauss_gauss_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t_1,
		double t_2, mui_sampler_gauss_1d *spatial_sampler, mui_chrono_sampler_gauss_1d *temporal_sampler);
double mui_fetch_gauss_gauss_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t_1,
		double t_2, mui_sampler_gauss_1dx *spatial_sampler, mui_chrono_sampler_gauss_1dx *temporal_sampler);
double mui_fetch_gauss_gauss_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t_1,
		double t_2, mui_sampler_gauss_1t *spatial_sampler, mui_chrono_sampler_gauss_1t *temporal_sampler);
float mui_fetch_gauss_mean_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t_1, float t_2,
		mui_sampler_gauss_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler);
float mui_fetch_gauss_mean_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t_1,
		float t_2, mui_sampler_gauss_1fx *spatial_sampler, mui_chrono_sampler_mean_1fx *temporal_sampler);
double mui_fetch_gauss_mean_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t_1,
		double t_2, mui_sampler_gauss_1d *spatial_sampler, mui_chrono_sampler_mean_1d *temporal_sampler);
double mui_fetch_gauss_mean_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t_1,
		double t_2, mui_sampler_gauss_1dx *spatial_sampler, mui_chrono_sampler_mean_1dx *temporal_sampler);
double mui_fetch_gauss_mean_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t_1,
		double t_2, mui_sampler_gauss_1t *spatial_sampler, mui_chrono_sampler_mean_1t *temporal_sampler);
float mui_fetch_moving_average_exact_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t_1,
		float t_2, mui_sampler_moving_average_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler);
float mui_fetch_moving_average_exact_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point,
		float t_1, float t_2, mui_sampler_moving_average_1fx *spatial_sampler,
		mui_chrono_sampler_exact_1fx *temporal_sampler);
double mui_fetch_moving_average_exact_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t_1,
		double t_2, mui_sampler_moving_average_1d *spatial_sampler, mui_chrono_sampler_exact_1d *temporal_sampler);
double mui_fetch_moving_average_exact_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point,
		double t_1, double t_2, mui_sampler_moving_average_1dx *spatial_sampler,
		mui_chrono_sampler_exact_1dx *temporal_sampler);
double mui_fetch_moving_average_exact_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t_1,
		double t_2, mui_sampler_moving_average_1t *spatial_sampler, mui_chrono_sampler_exact_1t *temporal_sampler);
float mui_fetch_moving_average_gauss_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t_1,
		float t_2, mui_sampler_moving_average_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler);
float mui_fetch_moving_average_gauss_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point,
		float t_1, float t_2, mui_sampler_moving_average_1fx *spatial_sampler,
		mui_chrono_sampler_gauss_1fx *temporal_sampler);
double mui_fetch_moving_average_gauss_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t_1,
		double t_2, mui_sampler_moving_average_1d *spatial_sampler, mui_chrono_sampler_gauss_1d *temporal_sampler);
double mui_fetch_moving_average_gauss_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point,
		double t_1, double t_2, mui_sampler_moving_average_1dx *spatial_sampler,
		mui_chrono_sampler_gauss_1dx *temporal_sampler);
double mui_fetch_moving_average_gauss_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t_1,
		double t_2, mui_sampler_moving_average_1t *spatial_sampler, mui_chrono_sampler_gauss_1t *temporal_sampler);
float mui_fetch_moving_average_mean_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t_1,
		float t_2, mui_sampler_moving_average_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler);
float mui_fetch_moving_average_mean_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t_1,
		float t_2, mui_sampler_moving_average_1fx *spatial_sampler, mui_chrono_sampler_mean_1fx *temporal_sampler);
double mui_fetch_moving_average_mean_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t_1,
		double t_2, mui_sampler_moving_average_1d *spatial_sampler, mui_chrono_sampler_mean_1d *temporal_sampler);
double mui_fetch_moving_average_mean_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point,
		double t_1, double t_2, mui_sampler_moving_average_1dx *spatial_sampler,
		mui_chrono_sampler_mean_1dx *temporal_sampler);
double mui_fetch_moving_average_mean_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t_1,
		double t_2, mui_sampler_moving_average_1t *spatial_sampler, mui_chrono_sampler_mean_1t *temporal_sampler);
float mui_fetch_moving_average_sum_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t_1,
		float t_2, mui_sampler_moving_average_1f *spatial_sampler, mui_chrono_sampler_sum_1f *temporal_sampler);
float mui_fetch_moving_average_sum_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t_1,
		float t_2, mui_sampler_moving_average_1fx *spatial_sampler, mui_chrono_sampler_sum_1fx *temporal_sampler);
double mui_fetch_moving_average_sum_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t_1,
		double t_2, mui_sampler_moving_average_1d *spatial_sampler, mui_chrono_sampler_sum_1d *temporal_sampler);
double mui_fetch_moving_average_sum_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point,
		double t_1, double t_2, mui_sampler_moving_average_1dx *spatial_sampler,
		mui_chrono_sampler_sum_1dx *temporal_sampler);
double mui_fetch_moving_average_sum_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t_1,
		double t_2, mui_sampler_moving_average_1t *spatial_sampler, mui_chrono_sampler_sum_1t *temporal_sampler);
float mui_fetch_nearest_neighbor_exact_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t_1,
		float t_2, mui_sampler_nearest_neighbor_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler);
float mui_fetch_nearest_neighbor_exact_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point,
		float t_1, float t_2, mui_sampler_nearest_neighbor_1fx *spatial_sampler,
		mui_chrono_sampler_exact_1fx *temporal_sampler);
double mui_fetch_nearest_neighbor_exact_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point,
		double t_1, double t_2, mui_sampler_nearest_neighbor_1d *spatial_sampler,
		mui_chrono_sampler_exact_1d *temporal_sampler);
double mui_fetch_nearest_neighbor_exact_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point,
		double t_1, double t_2, mui_sampler_nearest_neighbor_1dx *spatial_sampler,
		mui_chrono_sampler_exact_1dx *temporal_sampler);
double mui_fetch_nearest_neighbor_exact_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point,
		double t_1, double t_2, mui_sampler_nearest_neighbor_1t *spatial_sampler,
		mui_chrono_sampler_exact_1t *temporal_sampler);
float mui_fetch_nearest_neighbor_gauss_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t_1,
		float t_2, mui_sampler_nearest_neighbor_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler);
float mui_fetch_nearest_neighbor_gauss_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point,
		float t_1, float t_2, mui_sampler_nearest_neighbor_1fx *spatial_sampler,
		mui_chrono_sampler_gauss_1fx *temporal_sampler);
double mui_fetch_nearest_neighbor_gauss_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point,
		double t_1, double t_2, mui_sampler_nearest_neighbor_1d *spatial_sampler,
		mui_chrono_sampler_gauss_1d *temporal_sampler);
double mui_fetch_nearest_neighbor_gauss_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point,
		double t_1, double t_2, mui_sampler_nearest_neighbor_1dx *spatial_sampler,
		mui_chrono_sampler_gauss_1dx *temporal_sampler);
double mui_fetch_nearest_neighbor_gauss_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point,
		double t_1, double t_2, mui_sampler_nearest_neighbor_1t *spatial_sampler,
		mui_chrono_sampler_gauss_1t *temporal_sampler);
float mui_fetch_nearest_neighbor_mean_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t_1,
		float t_2, mui_sampler_nearest_neighbor_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler);
float mui_fetch_nearest_neighbor_mean_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point,
		float t_1, float t_2, mui_sampler_nearest_neighbor_1fx *spatial_sampler,
		mui_chrono_sampler_mean_1fx *temporal_sampler);
double mui_fetch_nearest_neighbor_mean_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point,
		double t_1, double t_2, mui_sampler_nearest_neighbor_1d *spatial_sampler,
		mui_chrono_sampler_mean_1d *temporal_sampler);
double mui_fetch_nearest_neighbor_mean_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point,
		double t_1, double t_2, mui_sampler_nearest_neighbor_1dx *spatial_sampler,
		mui_chrono_sampler_mean_1dx *temporal_sampler);
double mui_fetch_nearest_neighbor_mean_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point,
		double t_1, double t_2, mui_sampler_nearest_neighbor_1t *spatial_sampler,
		mui_chrono_sampler_mean_1t *temporal_sampler);
float mui_fetch_nearest_neighbor_sum_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t_1,
		float t_2, mui_sampler_nearest_neighbor_1f *spatial_sampler, mui_chrono_sampler_sum_1f *temporal_sampler);
float mui_fetch_nearest_neighbor_sum_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point,
		float t_1, float t_2, mui_sampler_nearest_neighbor_1fx *spatial_sampler,
		mui_chrono_sampler_sum_1fx *temporal_sampler);
double mui_fetch_nearest_neighbor_sum_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t_1,
		double t_2, mui_sampler_nearest_neighbor_1d *spatial_sampler, mui_chrono_sampler_sum_1d *temporal_sampler);
double mui_fetch_nearest_neighbor_sum_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point,
		double t_1, double t_2, mui_sampler_nearest_neighbor_1dx *spatial_sampler,
		mui_chrono_sampler_sum_1dx *temporal_sampler);
double mui_fetch_nearest_neighbor_sum_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t_1,
		double t_2, mui_sampler_nearest_neighbor_1t *spatial_sampler, mui_chrono_sampler_sum_1t *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_exact_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point,
		float t_1, float t_2, mui_sampler_pseudo_nearest_neighbor_1f *spatial_sampler,
		mui_chrono_sampler_exact_1f *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_exact_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point,
		float t_1, float t_2, mui_sampler_pseudo_nearest_neighbor_1fx *spatial_sampler,
		mui_chrono_sampler_exact_1fx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_exact_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_1d *spatial_sampler,
		mui_chrono_sampler_exact_1d *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_exact_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_1dx *spatial_sampler,
		mui_chrono_sampler_exact_1dx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_exact_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_1t *spatial_sampler,
		mui_chrono_sampler_exact_1t *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_gauss_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point,
		float t_1, float t_2, mui_sampler_pseudo_nearest_neighbor_1f *spatial_sampler,
		mui_chrono_sampler_gauss_1f *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_gauss_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point,
		float t_1, float t_2, mui_sampler_pseudo_nearest_neighbor_1fx *spatial_sampler,
		mui_chrono_sampler_gauss_1fx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_gauss_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_1d *spatial_sampler,
		mui_chrono_sampler_gauss_1d *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_gauss_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_1dx *spatial_sampler,
		mui_chrono_sampler_gauss_1dx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_gauss_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_1t *spatial_sampler,
		mui_chrono_sampler_gauss_1t *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_mean_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point,
		float t_1, float t_2, mui_sampler_pseudo_nearest_neighbor_1f *spatial_sampler,
		mui_chrono_sampler_mean_1f *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_mean_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point,
		float t_1, float t_2, mui_sampler_pseudo_nearest_neighbor_1fx *spatial_sampler,
		mui_chrono_sampler_mean_1fx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_mean_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_1d *spatial_sampler,
		mui_chrono_sampler_mean_1d *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_mean_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_1dx *spatial_sampler,
		mui_chrono_sampler_mean_1dx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_mean_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_1t *spatial_sampler,
		mui_chrono_sampler_mean_1t *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_sum_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point,
		float t_1, float t_2, mui_sampler_pseudo_nearest_neighbor_1f *spatial_sampler,
		mui_chrono_sampler_sum_1f *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_sum_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point,
		float t_1, float t_2, mui_sampler_pseudo_nearest_neighbor_1fx *spatial_sampler,
		mui_chrono_sampler_sum_1fx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_sum_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_1d *spatial_sampler,
		mui_chrono_sampler_sum_1d *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_sum_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_1dx *spatial_sampler,
		mui_chrono_sampler_sum_1dx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_sum_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_1t *spatial_sampler,
		mui_chrono_sampler_sum_1t *temporal_sampler);
float mui_fetch_shepard_quintic_exact_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t_1,
		float t_2, mui_sampler_shepard_quintic_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler);
float mui_fetch_shepard_quintic_exact_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point,
		float t_1, float t_2, mui_sampler_shepard_quintic_1fx *spatial_sampler,
		mui_chrono_sampler_exact_1fx *temporal_sampler);
double mui_fetch_shepard_quintic_exact_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point,
		double t_1, double t_2, mui_sampler_shepard_quintic_1d *spatial_sampler,
		mui_chrono_sampler_exact_1d *temporal_sampler);
double mui_fetch_shepard_quintic_exact_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point,
		double t_1, double t_2, mui_sampler_shepard_quintic_1dx *spatial_sampler,
		mui_chrono_sampler_exact_1dx *temporal_sampler);
double mui_fetch_shepard_quintic_exact_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point,
		double t_1, double t_2, mui_sampler_shepard_quintic_1t *spatial_sampler,
		mui_chrono_sampler_exact_1t *temporal_sampler);
float mui_fetch_shepard_quintic_gauss_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t_1,
		float t_2, mui_sampler_shepard_quintic_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler);
float mui_fetch_shepard_quintic_gauss_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point,
		float t_1, float t_2, mui_sampler_shepard_quintic_1fx *spatial_sampler,
		mui_chrono_sampler_gauss_1fx *temporal_sampler);
double mui_fetch_shepard_quintic_gauss_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point,
		double t_1, double t_2, mui_sampler_shepard_quintic_1d *spatial_sampler,
		mui_chrono_sampler_gauss_1d *temporal_sampler);
double mui_fetch_shepard_quintic_gauss_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point,
		double t_1, double t_2, mui_sampler_shepard_quintic_1dx *spatial_sampler,
		mui_chrono_sampler_gauss_1dx *temporal_sampler);
double mui_fetch_shepard_quintic_gauss_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point,
		double t_1, double t_2, mui_sampler_shepard_quintic_1t *spatial_sampler,
		mui_chrono_sampler_gauss_1t *temporal_sampler);
float mui_fetch_shepard_quintic_mean_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t_1,
		float t_2, mui_sampler_shepard_quintic_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler);
float mui_fetch_shepard_quintic_mean_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point,
		float t_1, float t_2, mui_sampler_shepard_quintic_1fx *spatial_sampler,
		mui_chrono_sampler_mean_1fx *temporal_sampler);
double mui_fetch_shepard_quintic_mean_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t_1,
		double t_2, mui_sampler_shepard_quintic_1d *spatial_sampler, mui_chrono_sampler_mean_1d *temporal_sampler);
double mui_fetch_shepard_quintic_mean_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point,
		double t_1, double t_2, mui_sampler_shepard_quintic_1dx *spatial_sampler,
		mui_chrono_sampler_mean_1dx *temporal_sampler);
double mui_fetch_shepard_quintic_mean_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t_1,
		double t_2, mui_sampler_shepard_quintic_1t *spatial_sampler, mui_chrono_sampler_mean_1t *temporal_sampler);
float mui_fetch_shepard_quintic_sum_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t_1,
		float t_2, mui_sampler_shepard_quintic_1f *spatial_sampler, mui_chrono_sampler_sum_1f *temporal_sampler);
float mui_fetch_shepard_quintic_sum_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t_1,
		float t_2, mui_sampler_shepard_quintic_1fx *spatial_sampler, mui_chrono_sampler_sum_1fx *temporal_sampler);
double mui_fetch_shepard_quintic_sum_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t_1,
		double t_2, mui_sampler_shepard_quintic_1d *spatial_sampler, mui_chrono_sampler_sum_1d *temporal_sampler);
double mui_fetch_shepard_quintic_sum_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point,
		double t_1, double t_2, mui_sampler_shepard_quintic_1dx *spatial_sampler,
		mui_chrono_sampler_sum_1dx *temporal_sampler);
double mui_fetch_shepard_quintic_sum_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t_1,
		double t_2, mui_sampler_shepard_quintic_1t *spatial_sampler, mui_chrono_sampler_sum_1t *temporal_sampler);
float mui_fetch_sph_quintic_exact_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t_1,
		float t_2, mui_sampler_sph_quintic_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler);
float mui_fetch_sph_quintic_exact_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t_1,
		float t_2, mui_sampler_sph_quintic_1fx *spatial_sampler, mui_chrono_sampler_exact_1fx *temporal_sampler);
double mui_fetch_sph_quintic_exact_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t_1,
		double t_2, mui_sampler_sph_quintic_1d *spatial_sampler, mui_chrono_sampler_exact_1d *temporal_sampler);
double mui_fetch_sph_quintic_exact_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t_1,
		double t_2, mui_sampler_sph_quintic_1dx *spatial_sampler, mui_chrono_sampler_exact_1dx *temporal_sampler);
double mui_fetch_sph_quintic_exact_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t_1,
		double t_2, mui_sampler_sph_quintic_1t *spatial_sampler, mui_chrono_sampler_exact_1t *temporal_sampler);
float mui_fetch_sph_quintic_gauss_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t_1,
		float t_2, mui_sampler_sph_quintic_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler);
float mui_fetch_sph_quintic_gauss_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t_1,
		float t_2, mui_sampler_sph_quintic_1fx *spatial_sampler, mui_chrono_sampler_gauss_1fx *temporal_sampler);
double mui_fetch_sph_quintic_gauss_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t_1,
		double t_2, mui_sampler_sph_quintic_1d *spatial_sampler, mui_chrono_sampler_gauss_1d *temporal_sampler);
double mui_fetch_sph_quintic_gauss_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t_1,
		double t_2, mui_sampler_sph_quintic_1dx *spatial_sampler, mui_chrono_sampler_gauss_1dx *temporal_sampler);
double mui_fetch_sph_quintic_gauss_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t_1,
		double t_2, mui_sampler_sph_quintic_1t *spatial_sampler, mui_chrono_sampler_gauss_1t *temporal_sampler);
float mui_fetch_sph_quintic_mean_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t_1,
		float t_2, mui_sampler_sph_quintic_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler);
float mui_fetch_sph_quintic_mean_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t_1,
		float t_2, mui_sampler_sph_quintic_1fx *spatial_sampler, mui_chrono_sampler_mean_1fx *temporal_sampler);
double mui_fetch_sph_quintic_mean_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t_1,
		double t_2, mui_sampler_sph_quintic_1d *spatial_sampler, mui_chrono_sampler_mean_1d *temporal_sampler);
double mui_fetch_sph_quintic_mean_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t_1,
		double t_2, mui_sampler_sph_quintic_1dx *spatial_sampler, mui_chrono_sampler_mean_1dx *temporal_sampler);
double mui_fetch_sph_quintic_mean_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t_1,
		double t_2, mui_sampler_sph_quintic_1t *spatial_sampler, mui_chrono_sampler_mean_1t *temporal_sampler);
float mui_fetch_sph_quintic_sum_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t_1,
		float t_2, mui_sampler_sph_quintic_1f *spatial_sampler, mui_chrono_sampler_sum_1f *temporal_sampler);
float mui_fetch_sph_quintic_sum_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t_1,
		float t_2, mui_sampler_sph_quintic_1fx *spatial_sampler, mui_chrono_sampler_sum_1fx *temporal_sampler);
double mui_fetch_sph_quintic_sum_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t_1,
		double t_2, mui_sampler_sph_quintic_1d *spatial_sampler, mui_chrono_sampler_sum_1d *temporal_sampler);
double mui_fetch_sph_quintic_sum_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t_1,
		double t_2, mui_sampler_sph_quintic_1dx *spatial_sampler, mui_chrono_sampler_sum_1dx *temporal_sampler);
double mui_fetch_sph_quintic_sum_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t_1,
		double t_2, mui_sampler_sph_quintic_1t *spatial_sampler, mui_chrono_sampler_sum_1t *temporal_sampler);
float mui_fetch_sum_quintic_exact_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t_1,
		float t_2, mui_sampler_sum_quintic_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler);
float mui_fetch_sum_quintic_exact_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t_1,
		float t_2, mui_sampler_sum_quintic_1fx *spatial_sampler, mui_chrono_sampler_exact_1fx *temporal_sampler);
double mui_fetch_sum_quintic_exact_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t_1,
		double t_2, mui_sampler_sum_quintic_1d *spatial_sampler, mui_chrono_sampler_exact_1d *temporal_sampler);
double mui_fetch_sum_quintic_exact_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t_1,
		double t_2, mui_sampler_sum_quintic_1dx *spatial_sampler, mui_chrono_sampler_exact_1dx *temporal_sampler);
double mui_fetch_sum_quintic_exact_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t_1,
		double t_2, mui_sampler_sum_quintic_1t *spatial_sampler, mui_chrono_sampler_exact_1t *temporal_sampler);
float mui_fetch_sum_quintic_gauss_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t_1,
		float t_2, mui_sampler_sum_quintic_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler);
float mui_fetch_sum_quintic_gauss_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t_1,
		float t_2, mui_sampler_sum_quintic_1fx *spatial_sampler, mui_chrono_sampler_gauss_1fx *temporal_sampler);
double mui_fetch_sum_quintic_gauss_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t_1,
		double t_2, mui_sampler_sum_quintic_1d *spatial_sampler, mui_chrono_sampler_gauss_1d *temporal_sampler);
double mui_fetch_sum_quintic_gauss_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t_1,
		double t_2, mui_sampler_sum_quintic_1dx *spatial_sampler, mui_chrono_sampler_gauss_1dx *temporal_sampler);
double mui_fetch_sum_quintic_gauss_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t_1,
		double t_2, mui_sampler_sum_quintic_1t *spatial_sampler, mui_chrono_sampler_gauss_1t *temporal_sampler);
float mui_fetch_sum_quintic_mean_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t_1,
		float t_2, mui_sampler_sum_quintic_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler);
float mui_fetch_sum_quintic_mean_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t_1,
		float t_2, mui_sampler_sum_quintic_1fx *spatial_sampler, mui_chrono_sampler_mean_1fx *temporal_sampler);
double mui_fetch_sum_quintic_mean_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t_1,
		double t_2, mui_sampler_sum_quintic_1d *spatial_sampler, mui_chrono_sampler_mean_1d *temporal_sampler);
double mui_fetch_sum_quintic_mean_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t_1,
		double t_2, mui_sampler_sum_quintic_1dx *spatial_sampler, mui_chrono_sampler_mean_1dx *temporal_sampler);
double mui_fetch_sum_quintic_mean_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t_1,
		double t_2, mui_sampler_sum_quintic_1t *spatial_sampler, mui_chrono_sampler_mean_1t *temporal_sampler);
float mui_fetch_sum_quintic_sum_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t_1,
		float t_2, mui_sampler_sum_quintic_1f *spatial_sampler, mui_chrono_sampler_sum_1f *temporal_sampler);
float mui_fetch_sum_quintic_sum_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t_1,
		float t_2, mui_sampler_sum_quintic_1fx *spatial_sampler, mui_chrono_sampler_sum_1fx *temporal_sampler);
double mui_fetch_sum_quintic_sum_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t_1,
		double t_2, mui_sampler_sum_quintic_1d *spatial_sampler, mui_chrono_sampler_sum_1d *temporal_sampler);
double mui_fetch_sum_quintic_sum_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t_1,
		double t_2, mui_sampler_sum_quintic_1dx *spatial_sampler, mui_chrono_sampler_sum_1dx *temporal_sampler);
double mui_fetch_sum_quintic_sum_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t_1,
		double t_2, mui_sampler_sum_quintic_1t *spatial_sampler, mui_chrono_sampler_sum_1t *temporal_sampler);
#ifdef USE_RBF
float mui_fetch_rbf_exact_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t_1, float t_2,
		mui_sampler_rbf_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler);
float mui_fetch_rbf_exact_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t_1,
		float t_2, mui_sampler_rbf_1fx *spatial_sampler, mui_chrono_sampler_exact_1fx *temporal_sampler);
double mui_fetch_rbf_exact_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t_1,
		double t_2, mui_sampler_rbf_1d *spatial_sampler, mui_chrono_sampler_exact_1d *temporal_sampler);
double mui_fetch_rbf_exact_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t_1,
		double t_2, mui_sampler_rbf_1dx *spatial_sampler, mui_chrono_sampler_exact_1dx *temporal_sampler);
double mui_fetch_rbf_exact_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t_1,
		double t_2, mui_sampler_rbf_1t *spatial_sampler, mui_chrono_sampler_exact_1t *temporal_sampler);
float mui_fetch_rbf_gauss_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t_1, float t_2,
		mui_sampler_rbf_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler);
float mui_fetch_rbf_gauss_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t_1,
		float t_2, mui_sampler_rbf_1fx *spatial_sampler, mui_chrono_sampler_gauss_1fx *temporal_sampler);
double mui_fetch_rbf_gauss_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t_1,
		double t_2, mui_sampler_rbf_1d *spatial_sampler, mui_chrono_sampler_gauss_1d *temporal_sampler);
double mui_fetch_rbf_gauss_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t_1,
		double t_2, mui_sampler_rbf_1dx *spatial_sampler, mui_chrono_sampler_gauss_1dx *temporal_sampler);
double mui_fetch_rbf_gauss_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t_1,
		double t_2, mui_sampler_rbf_1t *spatial_sampler, mui_chrono_sampler_gauss_1t *temporal_sampler);
float mui_fetch_rbf_mean_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t_1, float t_2,
		mui_sampler_rbf_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler);
float mui_fetch_rbf_mean_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t_1, float t_2,
		mui_sampler_rbf_1fx *spatial_sampler, mui_chrono_sampler_mean_1fx *temporal_sampler);
double mui_fetch_rbf_mean_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t_1, double t_2,
		mui_sampler_rbf_1d *spatial_sampler, mui_chrono_sampler_mean_1d *temporal_sampler);
double mui_fetch_rbf_mean_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t_1,
		double t_2, mui_sampler_rbf_1dx *spatial_sampler, mui_chrono_sampler_mean_1dx *temporal_sampler);
double mui_fetch_rbf_mean_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t_1, double t_2,
		mui_sampler_rbf_1t *spatial_sampler, mui_chrono_sampler_mean_1t *temporal_sampler);
float mui_fetch_rbf_sum_1f_pair(mui_uniface_1f *uniface, const char *attr, mui_point_1f point, float t_1, float t_2,
		mui_sampler_rbf_1f *spatial_sampler, mui_chrono_sampler_sum_1f *temporal_sampler);
float mui_fetch_rbf_sum_1fx_pair(mui_uniface_1fx *uniface, const char *attr, mui_point_1fx point, float t_1, float t_2,
		mui_sampler_rbf_1fx *spatial_sampler, mui_chrono_sampler_sum_1fx *temporal_sampler);
double mui_fetch_rbf_sum_1d_pair(mui_uniface_1d *uniface, const char *attr, mui_point_1d point, double t_1, double t_2,
		mui_sampler_rbf_1d *spatial_sampler, mui_chrono_sampler_sum_1d *temporal_sampler);
double mui_fetch_rbf_sum_1dx_pair(mui_uniface_1dx *uniface, const char *attr, mui_point_1dx point, double t_1,
		double t_2, mui_sampler_rbf_1dx *spatial_sampler, mui_chrono_sampler_sum_1dx *temporal_sampler);
double mui_fetch_rbf_sum_1t_pair(mui_uniface_1t *uniface, const char *attr, mui_point_1t point, double t_1, double t_2,
		mui_sampler_rbf_1t *spatial_sampler, mui_chrono_sampler_sum_1t *temporal_sampler);
#endif

// MUI functions for 1D data point only fetch using single time value
void mui_fetch_points_exact_1f(mui_uniface_1f *uniface, const char *attr, float t,
		mui_chrono_sampler_exact_1f *temporal_sampler, mui_point_1f *ret_points, int *num_points);
void mui_fetch_points_exact_1fx(mui_uniface_1fx *uniface, const char *attr, float t,
		mui_chrono_sampler_exact_1fx *temporal_sampler, mui_point_1fx *ret_points, int *num_points);
void mui_fetch_points_exact_1d(mui_uniface_1d *uniface, const char *attr, double t,
		mui_chrono_sampler_exact_1d *temporal_sampler, mui_point_1d *ret_points, int *num_points);
void mui_fetch_points_exact_1dx(mui_uniface_1dx *uniface, const char *attr, double t,
		mui_chrono_sampler_exact_1dx *temporal_sampler, mui_point_1dx *ret_points, int *num_points);
void mui_fetch_points_exact_1t(mui_uniface_1t *uniface, const char *attr, double t,
		mui_chrono_sampler_exact_1t *temporal_sampler, mui_point_1t *ret_points, int *num_points);
void mui_fetch_points_gauss_1f(mui_uniface_1f *uniface, const char *attr, float t,
		mui_chrono_sampler_gauss_1f *temporal_sampler, mui_point_1f *ret_points, int *num_points);
void mui_fetch_points_gauss_1fx(mui_uniface_1fx *uniface, const char *attr, float t,
		mui_chrono_sampler_gauss_1fx *temporal_sampler, mui_point_1fx *ret_points, int *num_points);
void mui_fetch_points_gauss_1d(mui_uniface_1d *uniface, const char *attr, double t,
		mui_chrono_sampler_gauss_1d *temporal_sampler, mui_point_1d *ret_points, int *num_points);
void mui_fetch_points_gauss_1dx(mui_uniface_1dx *uniface, const char *attr, double t,
		mui_chrono_sampler_gauss_1dx *temporal_sampler, mui_point_1dx *ret_points, int *num_points);
void mui_fetch_points_gauss_1t(mui_uniface_1t *uniface, const char *attr, double t,
		mui_chrono_sampler_gauss_1t *temporal_sampler, mui_point_1t *ret_points, int *num_points);
void mui_fetch_points_mean_1f(mui_uniface_1f *uniface, const char *attr, float t,
		mui_chrono_sampler_mean_1f *temporal_sampler, mui_point_1f *ret_points, int *num_points);
void mui_fetch_points_mean_1fx(mui_uniface_1fx *uniface, const char *attr, float t,
		mui_chrono_sampler_mean_1fx *temporal_sampler, mui_point_1fx *ret_points, int *num_points);
void mui_fetch_points_mean_1d(mui_uniface_1d *uniface, const char *attr, double t,
		mui_chrono_sampler_mean_1d *temporal_sampler, mui_point_1d *ret_points, int *num_points);
void mui_fetch_points_mean_1dx(mui_uniface_1dx *uniface, const char *attr, double t,
		mui_chrono_sampler_mean_1dx *temporal_sampler, mui_point_1dx *ret_points, int *num_points);
void mui_fetch_points_mean_1t(mui_uniface_1t *uniface, const char *attr, double t,
		mui_chrono_sampler_mean_1t *temporal_sampler, mui_point_1t *ret_points, int *num_points);
void mui_fetch_points_sum_1f(mui_uniface_1f *uniface, const char *attr, float t,
		mui_chrono_sampler_sum_1f *temporal_sampler, mui_point_1f *ret_points, int *num_points);
void mui_fetch_points_sum_1fx(mui_uniface_1fx *uniface, const char *attr, float t,
		mui_chrono_sampler_sum_1fx *temporal_sampler, mui_point_1fx *ret_points, int *num_points);
void mui_fetch_points_sum_1d(mui_uniface_1d *uniface, const char *attr, double t,
		mui_chrono_sampler_sum_1d *temporal_sampler, mui_point_1d *ret_points, int *num_points);
void mui_fetch_points_sum_1dx(mui_uniface_1dx *uniface, const char *attr, double t,
		mui_chrono_sampler_sum_1dx *temporal_sampler, mui_point_1dx *ret_points, int *num_points);
void mui_fetch_points_sum_1t(mui_uniface_1t *uniface, const char *attr, double t,
		mui_chrono_sampler_sum_1t *temporal_sampler, mui_point_1t *ret_points, int *num_points);

// MUI functions for 1D data point only fetch using two time values
void mui_fetch_values_exact_1f_pair(mui_uniface_1f *uniface, const char *attr, float t_1, float t_2,
		mui_chrono_sampler_exact_1f *temporal_sampler, float *values, int *num_values);
void mui_fetch_values_exact_1fx_pair(mui_uniface_1fx *uniface, const char *attr, float t_1, float t_2,
		mui_chrono_sampler_exact_1fx *temporal_sampler, float *values, int *num_values);
void mui_fetch_values_exact_1d_pair(mui_uniface_1d *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_exact_1d *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_exact_1dx_pair(mui_uniface_1dx *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_exact_1dx *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_exact_1t_pair(mui_uniface_1t *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_exact_1t *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_gauss_1f_pair(mui_uniface_1f *uniface, const char *attr, float t_1, float t_2,
		mui_chrono_sampler_gauss_1f *temporal_sampler, float *values, int *num_values);
void mui_fetch_values_gauss_1fx_pair(mui_uniface_1fx *uniface, const char *attr, float t_1, float t_2,
		mui_chrono_sampler_gauss_1fx *temporal_sampler, float *values, int *num_values);
void mui_fetch_values_gauss_1d_pair(mui_uniface_1d *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_gauss_1d *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_gauss_1dx_pair(mui_uniface_1dx *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_gauss_1dx *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_gauss_1t_pair(mui_uniface_1t *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_gauss_1t *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_mean_1f_pair(mui_uniface_1f *uniface, const char *attr, float t_1, float t_2,
		mui_chrono_sampler_mean_1f *temporal_sampler, float *values, int *num_values);
void mui_fetch_values_mean_1fx_pair(mui_uniface_1fx *uniface, const char *attr, float t_1, float t_2,
		mui_chrono_sampler_mean_1fx *temporal_sampler, float *values, int *num_values);
void mui_fetch_values_mean_1d_pair(mui_uniface_1d *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_mean_1d *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_mean_1dx_pair(mui_uniface_1dx *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_mean_1dx *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_mean_1t_pair(mui_uniface_1t *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_mean_1t *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_sum_1f_pair(mui_uniface_1f *uniface, const char *attr, float t_1, float t_2,
		mui_chrono_sampler_sum_1f *temporal_sampler, float *values, int *num_values);
void mui_fetch_values_sum_1fx_pair(mui_uniface_1fx *uniface, const char *attr, float t_1, float t_2,
		mui_chrono_sampler_sum_1fx *temporal_sampler, float *values, int *num_values);
void mui_fetch_values_sum_1d_pair(mui_uniface_1d *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_sum_1d *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_sum_1dx_pair(mui_uniface_1dx *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_sum_1dx *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_sum_1t_pair(mui_uniface_1t *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_sum_1t *temporal_sampler, double *values, int *num_values);

// MUI functions for single parameter fetch
float mui_fetch_1f_param(mui_uniface_1f *uniface, const char *attr);
float mui_fetch_1fx_param(mui_uniface_1fx *uniface, const char *attr);
double mui_fetch_1d_param(mui_uniface_1d *uniface, const char *attr);
double mui_fetch_1dx_param(mui_uniface_1dx *uniface, const char *attr);
double mui_fetch_1t_param(mui_uniface_1t *uniface, const char *attr);

// MUI data receive test functions
int mui_is_ready_1f(mui_uniface_1f *uniface, const char *attr, float t);
int mui_is_ready_1fx(mui_uniface_1fx *uniface, const char *attr, float t);
int mui_is_ready_1d(mui_uniface_1d *uniface, const char *attr, double t);
int mui_is_ready_1dx(mui_uniface_1dx *uniface, const char *attr, double t);
int mui_is_ready_1t(mui_uniface_1t *uniface, const char *attr, double t);
int mui_is_ready_1f_pair(mui_uniface_1f *uniface, const char *attr, float t_1, float t_2);
int mui_is_ready_1fx_pair(mui_uniface_1fx *uniface, const char *attr, float t_1, float t_2);
int mui_is_ready_1d_pair(mui_uniface_1d *uniface, const char *attr, double t_1, double t_2);
int mui_is_ready_1dx_pair(mui_uniface_1dx *uniface, const char *attr, double t_1, double t_2);
int mui_is_ready_1t_pair(mui_uniface_1t *uniface, const char *attr, double t_1, double t_2);

// MUI Smart Send functions
void mui_announce_send_span_1f_box(mui_uniface_1f *uniface, float box_1_1, float box_1_2, float t_start,
		float t_timeout, int synchronised);
void mui_announce_send_span_1fx_box(mui_uniface_1fx *uniface, float box_1_1, float box_1_2, float t_start,
		float t_timeout, int synchronised);
void mui_announce_send_span_1d_box(mui_uniface_1d *uniface, double box_1_1, double box_1_2, double t_start,
		double t_timeout, int synchronised);
void mui_announce_send_span_1dx_box(mui_uniface_1dx *uniface, double box_1_1, double box_1_2, double t_start,
		double t_timeout, int synchronised);
void mui_announce_send_span_1t_box(mui_uniface_1t *uniface, double box_1_1, double box_1_2, double t_start,
		double t_timeout, int synchronised);
void mui_announce_send_span_1f_sphere(mui_uniface_1f *uniface, mui_point_1f centre, float radius, float t_start,
		float t_timeout, int synchronised);
void mui_announce_send_span_1fx_sphere(mui_uniface_1fx *uniface, mui_point_1fx centre, float radius, float t_start,
		float t_timeout, int synchronised);
void mui_announce_send_span_1d_sphere(mui_uniface_1d *uniface, mui_point_1d centre, double radius, double t_start,
		double t_timeout, int synchronised);
void mui_announce_send_span_1dx_sphere(mui_uniface_1dx *uniface, mui_point_1dx centre, double radius, double t_start,
		double t_timeout, int synchronised);
void mui_announce_send_span_1t_sphere(mui_uniface_1t *uniface, mui_point_1t centre, double radius, double t_start,
		double t_timeout, int synchronised);
void mui_announce_recv_span_1f_box(mui_uniface_1f *uniface, float box_1_1, float box_1_2, float t_start,
		float t_timeout, int synchronised);
void mui_announce_recv_span_1fx_box(mui_uniface_1fx *uniface, float box_1_1, float box_1_2, float t_start,
		float t_timeout, int synchronised);
void mui_announce_recv_span_1d_box(mui_uniface_1d *uniface, double box_1_1, double box_1_2, double t_start,
		double t_timeout, int synchronised);
void mui_announce_recv_span_1dx_box(mui_uniface_1dx *uniface, double box_1_1, double box_1_2, double t_start,
		double t_timeout, int synchronised);
void mui_announce_recv_span_1t_box(mui_uniface_1t *uniface, double box_1_1, double box_1_2, double t_start,
		double t_timeout, int synchronised);
void mui_announce_recv_span_1f_sphere(mui_uniface_1f *uniface, mui_point_1f centre, float radius, float t_start,
		float t_timeout, int synchronised);
void mui_announce_recv_span_1fx_sphere(mui_uniface_1fx *uniface, mui_point_1fx centre, float radius, float t_start,
		float t_timeout, int synchronised);
void mui_announce_recv_span_1d_sphere(mui_uniface_1d *uniface, mui_point_1d centre, double radius, double t_start,
		double t_timeout, int synchronised);
void mui_announce_recv_span_1dx_sphere(mui_uniface_1dx *uniface, mui_point_1dx centre, double radius, double t_start,
		double t_timeout, int synchronised);
void mui_announce_recv_span_1t_sphere(mui_uniface_1t *uniface, mui_point_1t centre, double radius, double t_start,
		double t_timeout, int synchronised);
void mui_announce_send_disable_1f(mui_uniface_1f *uniface, int synchronised);
void mui_announce_send_disable_1fx(mui_uniface_1fx *uniface, int synchronised);
void mui_announce_send_disable_1d(mui_uniface_1d *uniface, int synchronised);
void mui_announce_send_disable_1dx(mui_uniface_1dx *uniface, int synchronised);
void mui_announce_send_disable_1t(mui_uniface_1t *uniface, int synchronised);
void mui_announce_recv_disable_1f(mui_uniface_1f *uniface, int synchronised);
void mui_announce_recv_disable_1fx(mui_uniface_1fx *uniface, int synchronised);
void mui_announce_recv_disable_1d(mui_uniface_1d *uniface, int synchronised);
void mui_announce_recv_disable_1dx(mui_uniface_1dx *uniface, int synchronised);
void mui_announce_recv_disable_1t(mui_uniface_1t *uniface, int synchronised);

// MUI barrier functions
void mui_barrier_1f(mui_uniface_1f *uniface, float t);
void mui_barrier_1fx(mui_uniface_1fx *uniface, float t);
void mui_barrier_1d(mui_uniface_1d *uniface, double t);
void mui_barrier_1dx(mui_uniface_1dx *uniface, double t);
void mui_barrier_1t(mui_uniface_1t *uniface, double t);
void mui_barrier_1f_pair(mui_uniface_1f *uniface, float t_1, float t_2);
void mui_barrier_1fx_pair(mui_uniface_1fx *uniface, float t_1, float t_2);
void mui_barrier_1d_pair(mui_uniface_1d *uniface, double t_1, double t_2);
void mui_barrier_1dx_pair(mui_uniface_1dx *uniface, double t_1, double t_2);
void mui_barrier_1t_pair(mui_uniface_1t *uniface, double t_1, double t_2);

// MUI forget functions
void mui_forget_upper_1f(mui_uniface_1f *uniface, float upper, int reset_log);
void mui_forget_upper_1fx(mui_uniface_1fx *uniface, float upper, int reset_log);
void mui_forget_upper_1d(mui_uniface_1d *uniface, double upper, int reset_log);
void mui_forget_upper_1dx(mui_uniface_1dx *uniface, double upper, int reset_log);
void mui_forget_upper_1t(mui_uniface_1t *uniface, double upper, int reset_log);
void mui_forget_upper_1f_pair(mui_uniface_1f *uniface, float upper_1, float upper_2, int reset_log);
void mui_forget_upper_1fx_pair(mui_uniface_1fx *uniface, float upper_1, float upper_2, int reset_log);
void mui_forget_upper_1d_pair(mui_uniface_1d *uniface, double upper_1, double upper_2, int reset_log);
void mui_forget_upper_1dx_pair(mui_uniface_1dx *uniface, double upper_1, double upper_2, int reset_log);
void mui_forget_upper_1t_pair(mui_uniface_1t *uniface, double upper_1, double upper_2, int reset_log);
void mui_forget_lower_upper_1f(mui_uniface_1f *uniface, float lower, float upper, int reset_log);
void mui_forget_lower_upper_1fx(mui_uniface_1fx *uniface, float lower, float upper, int reset_log);
void mui_forget_lower_upper_1d(mui_uniface_1d *uniface, double lower, double upper, int reset_log);
void mui_forget_lower_upper_1dx(mui_uniface_1dx *uniface, double lower, double upper, int reset_log);
void mui_forget_lower_upper_1t(mui_uniface_1t *uniface, double lower, double upper, int reset_log);
void mui_forget_lower_upper_1f_pair(mui_uniface_1f *uniface, float lower_1, float lower_2, float upper_1, float upper_2,
		int reset_log);
void mui_forget_lower_upper_1fx_pair(mui_uniface_1fx *uniface, float lower_1, float lower_2, float upper_1,
		float upper_2, int reset_log);
void mui_forget_lower_upper_1d_pair(mui_uniface_1d *uniface, double lower_1, double lower_2, double upper_1,
		double upper_2, int reset_log);
void mui_forget_lower_upper_1dx_pair(mui_uniface_1dx *uniface, double lower_1, double lower_2, double upper_1,
		double upper_2, int reset_log);
void mui_forget_lower_upper_1t_pair(mui_uniface_1t *uniface, double lower_1, double lower_2, double upper_1,
		double upper_2, int reset_log);
void mui_set_forget_length_1f(mui_uniface_1f *uniface, float length);
void mui_set_forget_length_1fx(mui_uniface_1fx *uniface, float length);
void mui_set_forget_length_1d(mui_uniface_1d *uniface, double length);
void mui_set_forget_length_1dx(mui_uniface_1dx *uniface, double length);
void mui_set_forget_length_1t(mui_uniface_1t *uniface, double length);
const char* mui_uri_host_1f(mui_uniface_1f *uniface);
const char* mui_uri_host_1fx(mui_uniface_1fx *uniface);
const char* mui_uri_host_1d(mui_uniface_1d *uniface);
const char* mui_uri_host_1dx(mui_uniface_1dx *uniface);
const char* mui_uri_host_1t(mui_uniface_1t *uniface);
const char* mui_uri_path_1f(mui_uniface_1f *uniface);
const char* mui_uri_path_1fx(mui_uniface_1fx *uniface);
const char* mui_uri_path_1d(mui_uniface_1d *uniface);
const char* mui_uri_path_1dx(mui_uniface_1dx *uniface);
const char* mui_uri_path_1t(mui_uniface_1t *uniface);
const char* mui_uri_protocol_1f(mui_uniface_1f *uniface);
const char* mui_uri_protocol_1fx(mui_uniface_1fx *uniface);
const char* mui_uri_protocol_1d(mui_uniface_1d *uniface);
const char* mui_uri_protocol_1dx(mui_uniface_1dx *uniface);
const char* mui_uri_protocol_1t(mui_uniface_1t *uniface);

#endif /* MUI_C_WRAPPER_1D_H_ */
