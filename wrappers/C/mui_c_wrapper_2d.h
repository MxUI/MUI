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
 * @file mui_c_wrapper_2d.h
 * @author S. M. Longshaw (derived from original 3D wrapper by Y. H. Tang)
 * @date Aug 4, 2021
 * @brief Header for C wrapper to create and manage 2D MUI interfaces and
 *        associated sampler objects
 *
 *        NOTE: Any point co-ordinates are enumerated rather than assuming
 *              Cartesian form, i.e. {1, 2, 3} rather than {x, y, z}.
 */

#ifndef MUI_C_WRAPPER_2D_H_
#define MUI_C_WRAPPER_2D_H_

// C-defined simple struct versions of MUI point types
typedef struct mui_point_2f {
	float point_1;
	float point_2;
} mui_point_2f;

typedef struct mui_point_2fx {
	float point_1;
	float point_2;
} mui_point_2fx;

typedef struct mui_point_2d {
	double point_1;
	double point_2;
} mui_point_2d;

typedef struct mui_point_2dx {
	double point_1;
	double point_2;
} mui_point_2dx;

typedef struct mui_point_2t {
	double point_1;
	double point_2;
} mui_point_2t;

// C access typedefs for uniface and sampler types
typedef struct mui_uniface_2f mui_uniface_2f;
typedef struct mui_uniface_2fx mui_uniface_2fx;
typedef struct mui_uniface_2d mui_uniface_2d;
typedef struct mui_uniface_2dx mui_uniface_2dx;
typedef struct mui_uniface_2t mui_uniface_2t;

typedef struct mui_sampler_exact_2f mui_sampler_exact_2f;
typedef struct mui_sampler_exact_2fx mui_sampler_exact_2fx;
typedef struct mui_sampler_exact_2d mui_sampler_exact_2d;
typedef struct mui_sampler_exact_2dx mui_sampler_exact_2dx;
typedef struct mui_sampler_exact_2t mui_sampler_exact_2t;

typedef struct mui_sampler_gauss_2f mui_sampler_gauss_2f;
typedef struct mui_sampler_gauss_2fx mui_sampler_gauss_2fx;
typedef struct mui_sampler_gauss_2d mui_sampler_gauss_2d;
typedef struct mui_sampler_gauss_2dx mui_sampler_gauss_2dx;
typedef struct mui_sampler_gauss_2t mui_sampler_gauss_2t;

typedef struct mui_sampler_moving_average_2f mui_sampler_moving_average_2f;
typedef struct mui_sampler_moving_average_2fx mui_sampler_moving_average_2fx;
typedef struct mui_sampler_moving_average_2d mui_sampler_moving_average_2d;
typedef struct mui_sampler_moving_average_2dx mui_sampler_moving_average_2dx;
typedef struct mui_sampler_moving_average_2t mui_sampler_moving_average_2t;

typedef struct mui_sampler_nearest_neighbor_2f mui_sampler_nearest_neighbor_2f;
typedef struct mui_sampler_nearest_neighbor_2fx mui_sampler_nearest_neighbor_2fx;
typedef struct mui_sampler_nearest_neighbor_2d mui_sampler_nearest_neighbor_2d;
typedef struct mui_sampler_nearest_neighbor_2dx mui_sampler_nearest_neighbor_2dx;
typedef struct mui_sampler_nearest_neighbor_2t mui_sampler_nearest_neighbor_2t;

typedef struct mui_sampler_pseudo_n2_linear_2f mui_sampler_pseudo_n2_linear_2f;
typedef struct mui_sampler_pseudo_n2_linear_2fx mui_sampler_pseudo_n2_linear_2fx;
typedef struct mui_sampler_pseudo_n2_linear_2d mui_sampler_pseudo_n2_linear_2d;
typedef struct mui_sampler_pseudo_n2_linear_2dx mui_sampler_pseudo_n2_linear_2dx;
typedef struct mui_sampler_pseudo_n2_linear_2t mui_sampler_pseudo_n2_linear_2t;

typedef struct mui_sampler_pseudo_nearest_neighbor_2f mui_sampler_pseudo_nearest_neighbor_2f;
typedef struct mui_sampler_pseudo_nearest_neighbor_2fx mui_sampler_pseudo_nearest_neighbor_2fx;
typedef struct mui_sampler_pseudo_nearest_neighbor_2d mui_sampler_pseudo_nearest_neighbor_2d;
typedef struct mui_sampler_pseudo_nearest_neighbor_2dx mui_sampler_pseudo_nearest_neighbor_2dx;
typedef struct mui_sampler_pseudo_nearest_neighbor_2t mui_sampler_pseudo_nearest_neighbor_2t;

typedef struct mui_sampler_shepard_quintic_2f mui_sampler_shepard_quintic_2f;
typedef struct mui_sampler_shepard_quintic_2fx mui_sampler_shepard_quintic_2fx;
typedef struct mui_sampler_shepard_quintic_2d mui_sampler_shepard_quintic_2d;
typedef struct mui_sampler_shepard_quintic_2dx mui_sampler_shepard_quintic_2dx;
typedef struct mui_sampler_shepard_quintic_2t mui_sampler_shepard_quintic_2t;

typedef struct mui_sampler_sph_quintic_2f mui_sampler_sph_quintic_2f;
typedef struct mui_sampler_sph_quintic_2fx mui_sampler_sph_quintic_2fx;
typedef struct mui_sampler_sph_quintic_2d mui_sampler_sph_quintic_2d;
typedef struct mui_sampler_sph_quintic_2dx mui_sampler_sph_quintic_2dx;
typedef struct mui_sampler_sph_quintic_2t mui_sampler_sph_quintic_2t;

typedef struct mui_sampler_sum_quintic_2f mui_sampler_sum_quintic_2f;
typedef struct mui_sampler_sum_quintic_2fx mui_sampler_sum_quintic_2fx;
typedef struct mui_sampler_sum_quintic_2d mui_sampler_sum_quintic_2d;
typedef struct mui_sampler_sum_quintic_2dx mui_sampler_sum_quintic_2dx;
typedef struct mui_sampler_sum_quintic_2t mui_sampler_sum_quintic_2t;

#ifdef USE_RBF
typedef struct mui_sampler_rbf_2f mui_sampler_rbf_2f;
typedef struct mui_sampler_rbf_2fx mui_sampler_rbf_2fx;
typedef struct mui_sampler_rbf_2d mui_sampler_rbf_2d;
typedef struct mui_sampler_rbf_2dx mui_sampler_rbf_2dx;
typedef struct mui_sampler_rbf_2t mui_sampler_rbf_2t;
#endif

typedef struct mui_chrono_sampler_exact_2f mui_chrono_sampler_exact_2f;
typedef struct mui_chrono_sampler_exact_2fx mui_chrono_sampler_exact_2fx;
typedef struct mui_chrono_sampler_exact_2d mui_chrono_sampler_exact_2d;
typedef struct mui_chrono_sampler_exact_2dx mui_chrono_sampler_exact_2dx;
typedef struct mui_chrono_sampler_exact_2t mui_chrono_sampler_exact_2t;

typedef struct mui_chrono_sampler_gauss_2f mui_chrono_sampler_gauss_2f;
typedef struct mui_chrono_sampler_gauss_2fx mui_chrono_sampler_gauss_2fx;
typedef struct mui_chrono_sampler_gauss_2d mui_chrono_sampler_gauss_2d;
typedef struct mui_chrono_sampler_gauss_2dx mui_chrono_sampler_gauss_2dx;
typedef struct mui_chrono_sampler_gauss_2t mui_chrono_sampler_gauss_2t;

typedef struct mui_chrono_sampler_mean_2f mui_chrono_sampler_mean_2f;
typedef struct mui_chrono_sampler_mean_2fx mui_chrono_sampler_mean_2fx;
typedef struct mui_chrono_sampler_mean_2d mui_chrono_sampler_mean_2d;
typedef struct mui_chrono_sampler_mean_2dx mui_chrono_sampler_mean_2dx;
typedef struct mui_chrono_sampler_mean_2t mui_chrono_sampler_mean_2t;

typedef struct mui_chrono_sampler_sum_2f mui_chrono_sampler_sum_2f;
typedef struct mui_chrono_sampler_sum_2fx mui_chrono_sampler_sum_2fx;
typedef struct mui_chrono_sampler_sum_2d mui_chrono_sampler_sum_2d;
typedef struct mui_chrono_sampler_sum_2dx mui_chrono_sampler_sum_2dx;
typedef struct mui_chrono_sampler_sum_2t mui_chrono_sampler_sum_2t;

// MUI uniface creation
mui_uniface_2f* mui_create_uniface_2f(const char *URI);
mui_uniface_2fx* mui_create_uniface_2fx(const char *URI);
mui_uniface_2d* mui_create_uniface_2d(const char *URI);
mui_uniface_2dx* mui_create_uniface_2dx(const char *URI);
mui_uniface_2t* mui_create_uniface_2t(const char *URI);

// MUI multi uniface creation
mui_uniface_2f** mui_create_uniface_multi_2f( const char *domain, const char **interfaces, int interface_count );
mui_uniface_2fx** mui_create_uniface_multi_2fx( const char *domain, const char **interfaces, int interface_count );
mui_uniface_2d** mui_create_uniface_multi_2d( const char *domain, const char **interfaces, int interface_count );
mui_uniface_2dx** mui_create_uniface_multi_2dx( const char *domain, const char **interfaces, int interface_count );
mui_uniface_2t** mui_create_uniface_multi_2t( const char *domain, const char **interfaces, int interface_count );

// MUI uniface destruction
void mui_destroy_uniface_2f(mui_uniface_2f *uniface);
void mui_destroy_uniface_2fx(mui_uniface_2fx *uniface);
void mui_destroy_uniface_2d(mui_uniface_2d *uniface);
void mui_destroy_uniface_2dx(mui_uniface_2dx *uniface);
void mui_destroy_uniface_2t(mui_uniface_2t *uniface);

// MUI spatial samplers creation
mui_sampler_exact_2f* mui_create_sampler_exact_2f(float tolerance);
mui_sampler_exact_2fx* mui_create_sampler_exact_2fx(float tolerance);
mui_sampler_exact_2d* mui_create_sampler_exact_2d(double tolerance);
mui_sampler_exact_2dx* mui_create_sampler_exact_2dx(double tolerance);
mui_sampler_exact_2t* mui_create_sampler_exact_2t(double tolerance);
mui_sampler_gauss_2f* mui_create_sampler_gauss_2f(float r, float h);
mui_sampler_gauss_2fx* mui_create_sampler_gauss_2fx(float r, float h);
mui_sampler_gauss_2d* mui_create_sampler_gauss_2d(double r, double h);
mui_sampler_gauss_2dx* mui_create_sampler_gauss_2dx(double r, double h);
mui_sampler_gauss_2t* mui_create_sampler_gauss_2t(double r, double h);
mui_sampler_moving_average_2f* mui_create_sampler_moving_average_2f(float bbox_1, float bbox_2);
mui_sampler_moving_average_2fx* mui_create_sampler_moving_average_2fx(float bbox_1, float bbox_2);
mui_sampler_moving_average_2d* mui_create_sampler_moving_average_2d(double bbox_1, double bbox_2);
mui_sampler_moving_average_2dx* mui_create_sampler_moving_average_2dx(double bbox_1, double bbox_2);
mui_sampler_moving_average_2t* mui_create_sampler_moving_average_2t(double bbox_1, double bbox_2);
mui_sampler_nearest_neighbor_2f* mui_create_sampler_nearest_neighbor_2f();
mui_sampler_nearest_neighbor_2fx* mui_create_sampler_nearest_neighbor_2fx();
mui_sampler_nearest_neighbor_2d* mui_create_sampler_nearest_neighbor_2d();
mui_sampler_nearest_neighbor_2dx* mui_create_sampler_nearest_neighbor_2dx();
mui_sampler_nearest_neighbor_2t* mui_create_sampler_nearest_neighbor_2t();
mui_sampler_pseudo_n2_linear_2f* mui_create_sampler_pseudo_n2_linear_2f(float r);
mui_sampler_pseudo_n2_linear_2fx* mui_create_sampler_pseudo_n2_linear_2fx(float r);
mui_sampler_pseudo_n2_linear_2d* mui_create_sampler_pseudo_n2_linear_2d(double r);
mui_sampler_pseudo_n2_linear_2dx* mui_create_sampler_pseudo_n2_linear_2dx(double r);
mui_sampler_pseudo_n2_linear_2t* mui_create_sampler_pseudo_n2_linear_2t(double r);
mui_sampler_pseudo_nearest_neighbor_2f* mui_create_sampler_pseudo_nearest_neighbor_2f(float h);
mui_sampler_pseudo_nearest_neighbor_2fx* mui_create_sampler_pseudo_nearest_neighbor_2fx(float h);
mui_sampler_pseudo_nearest_neighbor_2d* mui_create_sampler_pseudo_nearest_neighbor_2d(double h);
mui_sampler_pseudo_nearest_neighbor_2dx* mui_create_sampler_pseudo_nearest_neighbor_2dx(double h);
mui_sampler_pseudo_nearest_neighbor_2t* mui_create_sampler_pseudo_nearest_neighbor_2t(double h);
mui_sampler_shepard_quintic_2f* mui_create_sampler_shepard_quintic_2f(float r);
mui_sampler_shepard_quintic_2fx* mui_create_sampler_shepard_quintic_2fx(float r);
mui_sampler_shepard_quintic_2d* mui_create_sampler_shepard_quintic_2d(double r);
mui_sampler_shepard_quintic_2dx* mui_create_sampler_shepard_quintic_2dx(double r);
mui_sampler_shepard_quintic_2t* mui_create_sampler_shepard_quintic_2t(double r);
mui_sampler_sph_quintic_2f* mui_create_sampler_sph_quintic_2f(float r);
mui_sampler_sph_quintic_2fx* mui_create_sampler_sph_quintic_2fx(float r);
mui_sampler_sph_quintic_2d* mui_create_sampler_sph_quintic_2d(double r);
mui_sampler_sph_quintic_2dx* mui_create_sampler_sph_quintic_2dx(double r);
mui_sampler_sph_quintic_2t* mui_create_sampler_sph_quintic_2t(double r);
mui_sampler_sum_quintic_2f* mui_create_sampler_sum_quintic_2f(float r);
mui_sampler_sum_quintic_2fx* mui_create_sampler_sum_quintic_2fx(float r);
mui_sampler_sum_quintic_2d* mui_create_sampler_sum_quintic_2d(double r);
mui_sampler_sum_quintic_2dx* mui_create_sampler_sum_quintic_2dx(double r);
mui_sampler_sum_quintic_2t* mui_create_sampler_sum_quintic_2t(double r);
#ifdef USE_RBF
mui_sampler_rbf_2f* mui_create_sampler_rbf_2f(float r, mui_point_2f *points, int points_count, int basis_func,
		int conservative, int polynomial, int smoothFunc, int readMatrix, int writeMatrix, const char *file_address,
		float cutoff, float cg_solve_tol, int cg_solve_it, int pou_size);
mui_sampler_rbf_2fx* mui_create_sampler_rbf_2fx(float r, mui_point_2fx *points, int points_count, int basis_func,
		int conservative, int polynomial, int smoothFunc, int readMatrix, int writeMatrix, const char *file_address,
		float cutoff, float cg_solve_tol, int cg_solve_it, int pou_size);
mui_sampler_rbf_2d* mui_create_sampler_rbf_2d(double r, mui_point_2d *points, int points_count, int basis_func,
		int conservative, int polynomial, int smoothFunc, int readMatrix, int writeMatrix, const char *file_address,
		double cutoff, double cg_solve_tol, int cg_solve_it, int pou_size);
mui_sampler_rbf_2dx* mui_create_sampler_rbf_2dx(double r, mui_point_2dx *points, int points_count, int basis_func,
		int conservative, int polynomial, int smoothFunc, int readMatrix, int writeMatrix, const char *file_address,
		double cutoff, double cg_solve_tol, int cg_solve_it, int pou_size);
mui_sampler_rbf_2t* mui_create_sampler_rbf_2t(double r, mui_point_2t *points, int points_count, int basis_func,
		int conservative, int polynomial, int smoothFunc, int readMatrix, int writeMatrix, const char *file_address,
		double cutoff, double cg_solve_tol, int cg_solve_it, int pou_size);
#endif

// MUI spatial samplers destruction
void mui_destroy_sampler_exact_2f(mui_sampler_exact_2f *sampler);
void mui_destroy_sampler_exact_2fx(mui_sampler_exact_2fx *sampler);
void mui_destroy_sampler_exact_2d(mui_sampler_exact_2d *sampler);
void mui_destroy_sampler_exact_2dx(mui_sampler_exact_2dx *sampler);
void mui_destroy_sampler_exact_2t(mui_sampler_exact_2t *sampler);
void mui_destroy_sampler_gauss_2f(mui_sampler_gauss_2f *sampler);
void mui_destroy_sampler_gauss_2fx(mui_sampler_gauss_2fx *sampler);
void mui_destroy_sampler_gauss_2d(mui_sampler_gauss_2d *sampler);
void mui_destroy_sampler_gauss_2dx(mui_sampler_gauss_2dx *sampler);
void mui_destroy_sampler_gauss_2t(mui_sampler_gauss_2t *sampler);
void mui_destroy_sampler_moving_average_2f(mui_sampler_moving_average_2f *sampler);
void mui_destroy_sampler_moving_average_2fx(mui_sampler_moving_average_2fx *sampler);
void mui_destroy_sampler_moving_average_2d(mui_sampler_moving_average_2d *sampler);
void mui_destroy_sampler_moving_average_2dx(mui_sampler_moving_average_2dx *sampler);
void mui_destroy_sampler_moving_average_2t(mui_sampler_moving_average_2t *sampler);
void mui_destroy_sampler_nearest_neighbor_2f(mui_sampler_nearest_neighbor_2f *sampler);
void mui_destroy_sampler_nearest_neighbor_2fx(mui_sampler_nearest_neighbor_2fx *sampler);
void mui_destroy_sampler_nearest_neighbor_2d(mui_sampler_nearest_neighbor_2d *sampler);
void mui_destroy_sampler_nearest_neighbor_2dx(mui_sampler_nearest_neighbor_2dx *sampler);
void mui_destroy_sampler_nearest_neighbor_2t(mui_sampler_nearest_neighbor_2t *sampler);
void mui_destroy_sampler_pseudo_nearest2_linear_2f(mui_sampler_pseudo_nearest_neighbor_2f *sampler);
void mui_destroy_sampler_pseudo_nearest2_linear_2fx(mui_sampler_pseudo_nearest_neighbor_2fx *sampler);
void mui_destroy_sampler_pseudo_nearest2_linear_2d(mui_sampler_pseudo_nearest_neighbor_2d *sampler);
void mui_destroy_sampler_pseudo_nearest2_linear_2dx(mui_sampler_pseudo_nearest_neighbor_2dx *sampler);
void mui_destroy_sampler_pseudo_nearest2_linear_2t(mui_sampler_pseudo_nearest_neighbor_2t *sampler);
void mui_destroy_sampler_pseudo_nearest_neighbor_2f(mui_sampler_pseudo_nearest_neighbor_2f *sampler);
void mui_destroy_sampler_pseudo_nearest_neighbor_2fx(mui_sampler_pseudo_nearest_neighbor_2fx *sampler);
void mui_destroy_sampler_pseudo_nearest_neighbor_2d(mui_sampler_pseudo_nearest_neighbor_2d *sampler);
void mui_destroy_sampler_pseudo_nearest_neighbor_2dx(mui_sampler_pseudo_nearest_neighbor_2dx *sampler);
void mui_destroy_sampler_pseudo_nearest_neighbor_2t(mui_sampler_pseudo_nearest_neighbor_2t *sampler);
void mui_destroy_sampler_shepard_quintic_2f(mui_sampler_shepard_quintic_2f *sampler);
void mui_destroy_sampler_shepard_quintic_2fx(mui_sampler_shepard_quintic_2fx *sampler);
void mui_destroy_sampler_shepard_quintic_2d(mui_sampler_shepard_quintic_2d *sampler);
void mui_destroy_sampler_shepard_quintic_2dx(mui_sampler_shepard_quintic_2dx *sampler);
void mui_destroy_sampler_shepard_quintic_2t(mui_sampler_shepard_quintic_2t *sampler);
void mui_destroy_sampler_sph_quintic_2f(mui_sampler_sph_quintic_2f *sampler);
void mui_destroy_sampler_sph_quintic_2fx(mui_sampler_sph_quintic_2fx *sampler);
void mui_destroy_sampler_sph_quintic_2d(mui_sampler_sph_quintic_2d *sampler);
void mui_destroy_sampler_sph_quintic_2dx(mui_sampler_sph_quintic_2dx *sampler);
void mui_destroy_sampler_sph_quintic_2t(mui_sampler_sph_quintic_2t *sampler);
void mui_destroy_sampler_sum_quintic_2f(mui_sampler_sum_quintic_2f *sampler);
void mui_destroy_sampler_sum_quintic_2fx(mui_sampler_sum_quintic_2fx *sampler);
void mui_destroy_sampler_sum_quintic_2d(mui_sampler_sum_quintic_2d *sampler);
void mui_destroy_sampler_sum_quintic_2dx(mui_sampler_sum_quintic_2dx *sampler);
void mui_destroy_sampler_sum_quintic_2t(mui_sampler_sum_quintic_2t *sampler);
#ifdef USE_RBF
void mui_destroy_sampler_rbf_2f(mui_sampler_rbf_2f *sampler);
void mui_destroy_sampler_rbf_2fx(mui_sampler_rbf_2fx *sampler);
void mui_destroy_sampler_rbf_2d(mui_sampler_rbf_2d *sampler);
void mui_destroy_sampler_rbf_2dx(mui_sampler_rbf_2dx *sampler);
void mui_destroy_sampler_rbf_2t(mui_sampler_rbf_2t *sampler);
#endif

// MUI temporal samplers creation
mui_chrono_sampler_exact_2f* mui_create_chrono_sampler_exact_2f(float tolerance);
mui_chrono_sampler_exact_2fx* mui_create_chrono_sampler_exact_2fx(float tolerance);
mui_chrono_sampler_exact_2d* mui_create_chrono_sampler_exact_2d(double tolerance);
mui_chrono_sampler_exact_2dx* mui_create_chrono_sampler_exact_2dx(double tolerance);
mui_chrono_sampler_exact_2t* mui_create_chrono_sampler_exact_2t(double tolerance);
mui_chrono_sampler_gauss_2f* mui_create_chrono_sampler_gauss_2f(float cutoff, float sigma);
mui_chrono_sampler_gauss_2fx* mui_create_chrono_sampler_gauss_2fx(float cutoff, float sigma);
mui_chrono_sampler_gauss_2d* mui_create_chrono_sampler_gauss_2d(double cutoff, double sigma);
mui_chrono_sampler_gauss_2dx* mui_create_chrono_sampler_gauss_2dx(double cutoff, double sigma);
mui_chrono_sampler_gauss_2t* mui_create_chrono_sampler_gauss_2t(double cutoff, double sigma);
mui_chrono_sampler_mean_2f* mui_create_chrono_sampler_mean_2f(float lower, float upper);
mui_chrono_sampler_mean_2fx* mui_create_chrono_sampler_mean_2fx(float lower, float upper);
mui_chrono_sampler_mean_2d* mui_create_chrono_sampler_mean_2d(double lower, double upper);
mui_chrono_sampler_mean_2dx* mui_create_chrono_sampler_mean_2dx(double lower, double upper);
mui_chrono_sampler_mean_2t* mui_create_chrono_sampler_mean_2t(double lower, double upper);
mui_chrono_sampler_sum_2f* mui_create_chrono_sampler_sum_2f(float lower, float upper);
mui_chrono_sampler_sum_2fx* mui_create_chrono_sampler_sum_2fx(float lower, float upper);
mui_chrono_sampler_sum_2d* mui_create_chrono_sampler_sum_2d(double lower, double upper);
mui_chrono_sampler_sum_2dx* mui_create_chrono_sampler_sum_2dx(double lower, double upper);
mui_chrono_sampler_sum_2t* mui_create_chrono_sampler_sum_2t(double lower, double upper);

// MUI temporal samplers destruction
void mui_destroy_chrono_sampler_exact_2f(mui_chrono_sampler_exact_2f *sampler);
void mui_destroy_chrono_sampler_exact_2fx(mui_chrono_sampler_exact_2fx *sampler);
void mui_destroy_chrono_sampler_exact_2d(mui_chrono_sampler_exact_2d *sampler);
void mui_destroy_chrono_sampler_exact_2dx(mui_chrono_sampler_exact_2dx *sampler);
void mui_destroy_chrono_sampler_exact_2t(mui_chrono_sampler_exact_2t *sampler);
void mui_destroy_chrono_sampler_gauss_2f(mui_chrono_sampler_gauss_2f *sampler);
void mui_destroy_chrono_sampler_gauss_2fx(mui_chrono_sampler_gauss_2fx *sampler);
void mui_destroy_chrono_sampler_gauss_2d(mui_chrono_sampler_gauss_2d *sampler);
void mui_destroy_chrono_sampler_gauss_2dx(mui_chrono_sampler_gauss_2dx *sampler);
void mui_destroy_chrono_sampler_gauss_2t(mui_chrono_sampler_gauss_2t *sampler);
void mui_destroy_chrono_sampler_mean_2f(mui_chrono_sampler_mean_2f *sampler);
void mui_destroy_chrono_sampler_mean_2fx(mui_chrono_sampler_mean_2fx *sampler);
void mui_destroy_chrono_sampler_mean_2d(mui_chrono_sampler_mean_2d *sampler);
void mui_destroy_chrono_sampler_mean_2dx(mui_chrono_sampler_mean_2dx *sampler);
void mui_destroy_chrono_sampler_mean_2t(mui_chrono_sampler_mean_2t *sampler);
void mui_destroy_chrono_sampler_sum_2f(mui_chrono_sampler_sum_2f *sampler);
void mui_destroy_chrono_sampler_sum_2fx(mui_chrono_sampler_sum_2fx *sampler);
void mui_destroy_chrono_sampler_sum_2d(mui_chrono_sampler_sum_2d *sampler);
void mui_destroy_chrono_sampler_sum_2dx(mui_chrono_sampler_sum_2dx *sampler);
void mui_destroy_chrono_sampler_sum_2t(mui_chrono_sampler_sum_2t *sampler);

// MUI functions for data push
void mui_push_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float value);
void mui_push_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float value);
void mui_push_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double value);
void mui_push_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double value);
void mui_push_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double value);
void mui_push_2f_param(mui_uniface_2f *uniface, const char *attr, float value);
void mui_push_2fx_param(mui_uniface_2fx *uniface, const char *attr, float value);
void mui_push_2d_param(mui_uniface_2d *uniface, const char *attr, double value);
void mui_push_2dx_param(mui_uniface_2dx *uniface, const char *attr, double value);
void mui_push_2t_param(mui_uniface_2t *uniface, const char *attr, double value);

// MUI functions for data commit
void mui_commit_2f(mui_uniface_2f *uniface, float t);
void mui_commit_2fx(mui_uniface_2fx *uniface, float t);
void mui_commit_2d(mui_uniface_2d *uniface, double t);
void mui_commit_2dx(mui_uniface_2dx *uniface, double t);
void mui_commit_2t(mui_uniface_2t *uniface, double t);
void mui_commit_2f_pair(mui_uniface_2f *uniface, float t_1, float t_2);
void mui_commit_2fx_pair(mui_uniface_2fx *uniface, float t_1, float t_2);
void mui_commit_2d_pair(mui_uniface_2d *uniface, double t_1, double t_2);
void mui_commit_2dx_pair(mui_uniface_2dx *uniface, double t_1, double t_2);
void mui_commit_2t_pair(mui_uniface_2t *uniface, double t_1, double t_2);

// MUI functions for data forecast
void mui_forecast_2f(mui_uniface_2f *uniface, float t);
void mui_forecast_2fx(mui_uniface_2fx *uniface, float t);
void mui_forecast_2d(mui_uniface_2d *uniface, double t);
void mui_forecast_2dx(mui_uniface_2dx *uniface, double t);
void mui_forecast_2t(mui_uniface_2t *uniface, double t);
void mui_forecast_2f_pair(mui_uniface_2f *uniface, float t_1, float t_2);
void mui_forecast_2fx_pair(mui_uniface_2fx *uniface, float t_1, float t_2);
void mui_forecast_2d_pair(mui_uniface_2d *uniface, double t_1, double t_2);
void mui_forecast_2dx_pair(mui_uniface_2dx *uniface, double t_1, double t_2);
void mui_forecast_2t_pair(mui_uniface_2t *uniface, double t_1, double t_2);

// MUI functions for data fetch
float mui_fetch_exact_exact_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_exact_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler);
float mui_fetch_exact_exact_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t,
		mui_sampler_exact_2fx *spatial_sampler, mui_chrono_sampler_exact_2fx *temporal_sampler);
double mui_fetch_exact_exact_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_exact_2d *spatial_sampler, mui_chrono_sampler_exact_2d *temporal_sampler);
double mui_fetch_exact_exact_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t,
		mui_sampler_exact_2dx *spatial_sampler, mui_chrono_sampler_exact_2dx *temporal_sampler);
double mui_fetch_exact_exact_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_exact_2t *spatial_sampler, mui_chrono_sampler_exact_2t *temporal_sampler);
float mui_fetch_exact_gauss_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_exact_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler);
float mui_fetch_exact_gauss_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t,
		mui_sampler_exact_2fx *spatial_sampler, mui_chrono_sampler_gauss_2fx *temporal_sampler);
double mui_fetch_exact_gauss_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_exact_2d *spatial_sampler, mui_chrono_sampler_gauss_2d *temporal_sampler);
double mui_fetch_exact_gauss_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t,
		mui_sampler_exact_2dx *spatial_sampler, mui_chrono_sampler_gauss_2dx *temporal_sampler);
double mui_fetch_exact_gauss_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_exact_2t *spatial_sampler, mui_chrono_sampler_gauss_2t *temporal_sampler);
float mui_fetch_exact_mean_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_exact_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler);
float mui_fetch_exact_mean_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t,
		mui_sampler_exact_2fx *spatial_sampler, mui_chrono_sampler_mean_2fx *temporal_sampler);
double mui_fetch_exact_mean_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_exact_2d *spatial_sampler, mui_chrono_sampler_mean_2d *temporal_sampler);
double mui_fetch_exact_mean_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t,
		mui_sampler_exact_2dx *spatial_sampler, mui_chrono_sampler_mean_2dx *temporal_sampler);
double mui_fetch_exact_mean_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_exact_2t *spatial_sampler, mui_chrono_sampler_mean_2t *temporal_sampler);
float mui_fetch_exact_sum_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_exact_2f *spatial_sampler, mui_chrono_sampler_sum_2f *temporal_sampler);
float mui_fetch_exact_sum_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t,
		mui_sampler_exact_2fx *spatial_sampler, mui_chrono_sampler_sum_2fx *temporal_sampler);
double mui_fetch_exact_sum_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_exact_2d *spatial_sampler, mui_chrono_sampler_sum_2d *temporal_sampler);
double mui_fetch_exact_sum_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t,
		mui_sampler_exact_2dx *spatial_sampler, mui_chrono_sampler_sum_2dx *temporal_sampler);
double mui_fetch_exact_sum_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_exact_2t *spatial_sampler, mui_chrono_sampler_sum_2t *temporal_sampler);
float mui_fetch_gauss_exact_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_gauss_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler);
float mui_fetch_gauss_exact_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t,
		mui_sampler_gauss_2fx *spatial_sampler, mui_chrono_sampler_exact_2fx *temporal_sampler);
double mui_fetch_gauss_exact_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_gauss_2d *spatial_sampler, mui_chrono_sampler_exact_2d *temporal_sampler);
double mui_fetch_gauss_exact_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t,
		mui_sampler_gauss_2dx *spatial_sampler, mui_chrono_sampler_exact_2dx *temporal_sampler);
double mui_fetch_gauss_exact_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_gauss_2t *spatial_sampler, mui_chrono_sampler_exact_2t *temporal_sampler);
float mui_fetch_gauss_gauss_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_gauss_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler);
float mui_fetch_gauss_gauss_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t,
		mui_sampler_gauss_2fx *spatial_sampler, mui_chrono_sampler_gauss_2fx *temporal_sampler);
double mui_fetch_gauss_gauss_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_gauss_2d *spatial_sampler, mui_chrono_sampler_gauss_2d *temporal_sampler);
double mui_fetch_gauss_gauss_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t,
		mui_sampler_gauss_2dx *spatial_sampler, mui_chrono_sampler_gauss_2dx *temporal_sampler);
double mui_fetch_gauss_gauss_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_gauss_2t *spatial_sampler, mui_chrono_sampler_gauss_2t *temporal_sampler);
float mui_fetch_gauss_mean_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_gauss_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler);
float mui_fetch_gauss_mean_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t,
		mui_sampler_gauss_2fx *spatial_sampler, mui_chrono_sampler_mean_2fx *temporal_sampler);
double mui_fetch_gauss_mean_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_gauss_2d *spatial_sampler, mui_chrono_sampler_mean_2d *temporal_sampler);
double mui_fetch_gauss_mean_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t,
		mui_sampler_gauss_2dx *spatial_sampler, mui_chrono_sampler_mean_2dx *temporal_sampler);
double mui_fetch_gauss_mean_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_gauss_2t *spatial_sampler, mui_chrono_sampler_mean_2t *temporal_sampler);
float mui_fetch_moving_average_exact_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_moving_average_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler);
float mui_fetch_moving_average_exact_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t,
		mui_sampler_moving_average_2fx *spatial_sampler, mui_chrono_sampler_exact_2fx *temporal_sampler);
double mui_fetch_moving_average_exact_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_moving_average_2d *spatial_sampler, mui_chrono_sampler_exact_2d *temporal_sampler);
double mui_fetch_moving_average_exact_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t,
		mui_sampler_moving_average_2dx *spatial_sampler, mui_chrono_sampler_exact_2dx *temporal_sampler);
double mui_fetch_moving_average_exact_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_moving_average_2t *spatial_sampler, mui_chrono_sampler_exact_2t *temporal_sampler);
float mui_fetch_moving_average_gauss_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_moving_average_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler);
float mui_fetch_moving_average_gauss_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t,
		mui_sampler_moving_average_2fx *spatial_sampler, mui_chrono_sampler_gauss_2fx *temporal_sampler);
double mui_fetch_moving_average_gauss_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_moving_average_2d *spatial_sampler, mui_chrono_sampler_gauss_2d *temporal_sampler);
double mui_fetch_moving_average_gauss_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t,
		mui_sampler_moving_average_2dx *spatial_sampler, mui_chrono_sampler_gauss_2dx *temporal_sampler);
double mui_fetch_moving_average_gauss_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_moving_average_2t *spatial_sampler, mui_chrono_sampler_gauss_2t *temporal_sampler);
float mui_fetch_moving_average_mean_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_moving_average_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler);
float mui_fetch_moving_average_mean_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t,
		mui_sampler_moving_average_2fx *spatial_sampler, mui_chrono_sampler_mean_2fx *temporal_sampler);
double mui_fetch_moving_average_mean_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_moving_average_2d *spatial_sampler, mui_chrono_sampler_mean_2d *temporal_sampler);
double mui_fetch_moving_average_mean_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t,
		mui_sampler_moving_average_2dx *spatial_sampler, mui_chrono_sampler_mean_2dx *temporal_sampler);
double mui_fetch_moving_average_mean_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_moving_average_2t *spatial_sampler, mui_chrono_sampler_mean_2t *temporal_sampler);
float mui_fetch_moving_average_sum_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_moving_average_2f *spatial_sampler, mui_chrono_sampler_sum_2f *temporal_sampler);
float mui_fetch_moving_average_sum_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t,
		mui_sampler_moving_average_2fx *spatial_sampler, mui_chrono_sampler_sum_2fx *temporal_sampler);
double mui_fetch_moving_average_sum_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_moving_average_2d *spatial_sampler, mui_chrono_sampler_sum_2d *temporal_sampler);
double mui_fetch_moving_average_sum_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t,
		mui_sampler_moving_average_2dx *spatial_sampler, mui_chrono_sampler_sum_2dx *temporal_sampler);
double mui_fetch_moving_average_sum_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_moving_average_2t *spatial_sampler, mui_chrono_sampler_sum_2t *temporal_sampler);
float mui_fetch_nearest_neighbor_exact_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_nearest_neighbor_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler);
float mui_fetch_nearest_neighbor_exact_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t,
		mui_sampler_nearest_neighbor_2fx *spatial_sampler, mui_chrono_sampler_exact_2fx *temporal_sampler);
double mui_fetch_nearest_neighbor_exact_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_nearest_neighbor_2d *spatial_sampler, mui_chrono_sampler_exact_2d *temporal_sampler);
double mui_fetch_nearest_neighbor_exact_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t,
		mui_sampler_nearest_neighbor_2dx *spatial_sampler, mui_chrono_sampler_exact_2dx *temporal_sampler);
double mui_fetch_nearest_neighbor_exact_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_nearest_neighbor_2t *spatial_sampler, mui_chrono_sampler_exact_2t *temporal_sampler);
float mui_fetch_nearest_neighbor_gauss_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_nearest_neighbor_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler);
float mui_fetch_nearest_neighbor_gauss_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t,
		mui_sampler_nearest_neighbor_2fx *spatial_sampler, mui_chrono_sampler_gauss_2fx *temporal_sampler);
double mui_fetch_nearest_neighbor_gauss_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_nearest_neighbor_2d *spatial_sampler, mui_chrono_sampler_gauss_2d *temporal_sampler);
double mui_fetch_nearest_neighbor_gauss_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t,
		mui_sampler_nearest_neighbor_2dx *spatial_sampler, mui_chrono_sampler_gauss_2dx *temporal_sampler);
double mui_fetch_nearest_neighbor_gauss_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_nearest_neighbor_2t *spatial_sampler, mui_chrono_sampler_gauss_2t *temporal_sampler);
float mui_fetch_nearest_neighbor_mean_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_nearest_neighbor_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler);
float mui_fetch_nearest_neighbor_mean_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t,
		mui_sampler_nearest_neighbor_2fx *spatial_sampler, mui_chrono_sampler_mean_2fx *temporal_sampler);
double mui_fetch_nearest_neighbor_mean_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_nearest_neighbor_2d *spatial_sampler, mui_chrono_sampler_mean_2d *temporal_sampler);
double mui_fetch_nearest_neighbor_mean_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t,
		mui_sampler_nearest_neighbor_2dx *spatial_sampler, mui_chrono_sampler_mean_2dx *temporal_sampler);
double mui_fetch_nearest_neighbor_mean_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_nearest_neighbor_2t *spatial_sampler, mui_chrono_sampler_mean_2t *temporal_sampler);
float mui_fetch_nearest_neighbor_sum_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_nearest_neighbor_2f *spatial_sampler, mui_chrono_sampler_sum_2f *temporal_sampler);
float mui_fetch_nearest_neighbor_sum_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t,
		mui_sampler_nearest_neighbor_2fx *spatial_sampler, mui_chrono_sampler_sum_2fx *temporal_sampler);
double mui_fetch_nearest_neighbor_sum_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_nearest_neighbor_2d *spatial_sampler, mui_chrono_sampler_sum_2d *temporal_sampler);
double mui_fetch_nearest_neighbor_sum_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t,
		mui_sampler_nearest_neighbor_2dx *spatial_sampler, mui_chrono_sampler_sum_2dx *temporal_sampler);
double mui_fetch_nearest_neighbor_sum_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_nearest_neighbor_2t *spatial_sampler, mui_chrono_sampler_sum_2t *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_exact_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_pseudo_nearest_neighbor_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_exact_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point,
		float t, mui_sampler_pseudo_nearest_neighbor_2fx *spatial_sampler,
		mui_chrono_sampler_exact_2fx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_exact_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point,
		double t, mui_sampler_pseudo_nearest_neighbor_2d *spatial_sampler,
		mui_chrono_sampler_exact_2d *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_exact_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point,
		double t, mui_sampler_pseudo_nearest_neighbor_2dx *spatial_sampler,
		mui_chrono_sampler_exact_2dx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_exact_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point,
		double t, mui_sampler_pseudo_nearest_neighbor_2t *spatial_sampler,
		mui_chrono_sampler_exact_2t *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_gauss_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_pseudo_nearest_neighbor_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_gauss_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point,
		float t, mui_sampler_pseudo_nearest_neighbor_2fx *spatial_sampler,
		mui_chrono_sampler_gauss_2fx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_gauss_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point,
		double t, mui_sampler_pseudo_nearest_neighbor_2d *spatial_sampler,
		mui_chrono_sampler_gauss_2d *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_gauss_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point,
		double t, mui_sampler_pseudo_nearest_neighbor_2dx *spatial_sampler,
		mui_chrono_sampler_gauss_2dx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_gauss_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point,
		double t, mui_sampler_pseudo_nearest_neighbor_2t *spatial_sampler,
		mui_chrono_sampler_gauss_2t *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_mean_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_pseudo_nearest_neighbor_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_mean_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point,
		float t, mui_sampler_pseudo_nearest_neighbor_2fx *spatial_sampler,
		mui_chrono_sampler_mean_2fx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_mean_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point,
		double t, mui_sampler_pseudo_nearest_neighbor_2d *spatial_sampler,
		mui_chrono_sampler_mean_2d *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_mean_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point,
		double t, mui_sampler_pseudo_nearest_neighbor_2dx *spatial_sampler,
		mui_chrono_sampler_mean_2dx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_mean_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point,
		double t, mui_sampler_pseudo_nearest_neighbor_2t *spatial_sampler,
		mui_chrono_sampler_mean_2t *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_sum_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_pseudo_nearest_neighbor_2f *spatial_sampler, mui_chrono_sampler_sum_2f *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_sum_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point,
		float t, mui_sampler_pseudo_nearest_neighbor_2fx *spatial_sampler,
		mui_chrono_sampler_sum_2fx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_sum_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_pseudo_nearest_neighbor_2d *spatial_sampler, mui_chrono_sampler_sum_2d *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_sum_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point,
		double t, mui_sampler_pseudo_nearest_neighbor_2dx *spatial_sampler,
		mui_chrono_sampler_sum_2dx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_sum_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_pseudo_nearest_neighbor_2t *spatial_sampler, mui_chrono_sampler_sum_2t *temporal_sampler);
float mui_fetch_shepard_quintic_exact_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_shepard_quintic_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler);
float mui_fetch_shepard_quintic_exact_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t,
		mui_sampler_shepard_quintic_2fx *spatial_sampler, mui_chrono_sampler_exact_2fx *temporal_sampler);
double mui_fetch_shepard_quintic_exact_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_shepard_quintic_2d *spatial_sampler, mui_chrono_sampler_exact_2d *temporal_sampler);
double mui_fetch_shepard_quintic_exact_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t,
		mui_sampler_shepard_quintic_2dx *spatial_sampler, mui_chrono_sampler_exact_2dx *temporal_sampler);
double mui_fetch_shepard_quintic_exact_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_shepard_quintic_2t *spatial_sampler, mui_chrono_sampler_exact_2t *temporal_sampler);
float mui_fetch_shepard_quintic_gauss_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_shepard_quintic_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler);
float mui_fetch_shepard_quintic_gauss_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t,
		mui_sampler_shepard_quintic_2fx *spatial_sampler, mui_chrono_sampler_gauss_2fx *temporal_sampler);
double mui_fetch_shepard_quintic_gauss_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_shepard_quintic_2d *spatial_sampler, mui_chrono_sampler_gauss_2d *temporal_sampler);
double mui_fetch_shepard_quintic_gauss_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t,
		mui_sampler_shepard_quintic_2dx *spatial_sampler, mui_chrono_sampler_gauss_2dx *temporal_sampler);
double mui_fetch_shepard_quintic_gauss_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_shepard_quintic_2t *spatial_sampler, mui_chrono_sampler_gauss_2t *temporal_sampler);
float mui_fetch_shepard_quintic_mean_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_shepard_quintic_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler);
float mui_fetch_shepard_quintic_mean_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t,
		mui_sampler_shepard_quintic_2fx *spatial_sampler, mui_chrono_sampler_mean_2fx *temporal_sampler);
double mui_fetch_shepard_quintic_mean_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_shepard_quintic_2d *spatial_sampler, mui_chrono_sampler_mean_2d *temporal_sampler);
double mui_fetch_shepard_quintic_mean_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t,
		mui_sampler_shepard_quintic_2dx *spatial_sampler, mui_chrono_sampler_mean_2dx *temporal_sampler);
double mui_fetch_shepard_quintic_mean_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_shepard_quintic_2t *spatial_sampler, mui_chrono_sampler_mean_2t *temporal_sampler);
float mui_fetch_shepard_quintic_sum_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_shepard_quintic_2f *spatial_sampler, mui_chrono_sampler_sum_2f *temporal_sampler);
float mui_fetch_shepard_quintic_sum_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t,
		mui_sampler_shepard_quintic_2fx *spatial_sampler, mui_chrono_sampler_sum_2fx *temporal_sampler);
double mui_fetch_shepard_quintic_sum_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_shepard_quintic_2d *spatial_sampler, mui_chrono_sampler_sum_2d *temporal_sampler);
double mui_fetch_shepard_quintic_sum_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t,
		mui_sampler_shepard_quintic_2dx *spatial_sampler, mui_chrono_sampler_sum_2dx *temporal_sampler);
double mui_fetch_shepard_quintic_sum_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_shepard_quintic_2t *spatial_sampler, mui_chrono_sampler_sum_2t *temporal_sampler);
float mui_fetch_sph_quintic_exact_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_sph_quintic_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler);
float mui_fetch_sph_quintic_exact_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t,
		mui_sampler_sph_quintic_2fx *spatial_sampler, mui_chrono_sampler_exact_2fx *temporal_sampler);
double mui_fetch_sph_quintic_exact_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_sph_quintic_2d *spatial_sampler, mui_chrono_sampler_exact_2d *temporal_sampler);
double mui_fetch_sph_quintic_exact_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t,
		mui_sampler_sph_quintic_2dx *spatial_sampler, mui_chrono_sampler_exact_2dx *temporal_sampler);
double mui_fetch_sph_quintic_exact_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_sph_quintic_2t *spatial_sampler, mui_chrono_sampler_exact_2t *temporal_sampler);
float mui_fetch_sph_quintic_gauss_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_sph_quintic_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler);
float mui_fetch_sph_quintic_gauss_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t,
		mui_sampler_sph_quintic_2fx *spatial_sampler, mui_chrono_sampler_gauss_2fx *temporal_sampler);
double mui_fetch_sph_quintic_gauss_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_sph_quintic_2d *spatial_sampler, mui_chrono_sampler_gauss_2d *temporal_sampler);
double mui_fetch_sph_quintic_gauss_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t,
		mui_sampler_sph_quintic_2dx *spatial_sampler, mui_chrono_sampler_gauss_2dx *temporal_sampler);
double mui_fetch_sph_quintic_gauss_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_sph_quintic_2t *spatial_sampler, mui_chrono_sampler_gauss_2t *temporal_sampler);
float mui_fetch_sph_quintic_mean_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_sph_quintic_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler);
float mui_fetch_sph_quintic_mean_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t,
		mui_sampler_sph_quintic_2fx *spatial_sampler, mui_chrono_sampler_mean_2fx *temporal_sampler);
double mui_fetch_sph_quintic_mean_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_sph_quintic_2d *spatial_sampler, mui_chrono_sampler_mean_2d *temporal_sampler);
double mui_fetch_sph_quintic_mean_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t,
		mui_sampler_sph_quintic_2dx *spatial_sampler, mui_chrono_sampler_mean_2dx *temporal_sampler);
double mui_fetch_sph_quintic_mean_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_sph_quintic_2t *spatial_sampler, mui_chrono_sampler_mean_2t *temporal_sampler);
float mui_fetch_sph_quintic_sum_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_sph_quintic_2f *spatial_sampler, mui_chrono_sampler_sum_2f *temporal_sampler);
float mui_fetch_sph_quintic_sum_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t,
		mui_sampler_sph_quintic_2fx *spatial_sampler, mui_chrono_sampler_sum_2fx *temporal_sampler);
double mui_fetch_sph_quintic_sum_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_sph_quintic_2d *spatial_sampler, mui_chrono_sampler_sum_2d *temporal_sampler);
double mui_fetch_sph_quintic_sum_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t,
		mui_sampler_sph_quintic_2dx *spatial_sampler, mui_chrono_sampler_sum_2dx *temporal_sampler);
double mui_fetch_sph_quintic_sum_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_sph_quintic_2t *spatial_sampler, mui_chrono_sampler_sum_2t *temporal_sampler);
float mui_fetch_sum_quintic_exact_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_sum_quintic_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler);
float mui_fetch_sum_quintic_exact_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t,
		mui_sampler_sum_quintic_2fx *spatial_sampler, mui_chrono_sampler_exact_2fx *temporal_sampler);
double mui_fetch_sum_quintic_exact_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_sum_quintic_2d *spatial_sampler, mui_chrono_sampler_exact_2d *temporal_sampler);
double mui_fetch_sum_quintic_exact_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t,
		mui_sampler_sum_quintic_2dx *spatial_sampler, mui_chrono_sampler_exact_2dx *temporal_sampler);
double mui_fetch_sum_quintic_exact_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_sum_quintic_2t *spatial_sampler, mui_chrono_sampler_exact_2t *temporal_sampler);
float mui_fetch_sum_quintic_gauss_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_sum_quintic_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler);
float mui_fetch_sum_quintic_gauss_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t,
		mui_sampler_sum_quintic_2fx *spatial_sampler, mui_chrono_sampler_gauss_2fx *temporal_sampler);
double mui_fetch_sum_quintic_gauss_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_sum_quintic_2d *spatial_sampler, mui_chrono_sampler_gauss_2d *temporal_sampler);
double mui_fetch_sum_quintic_gauss_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t,
		mui_sampler_sum_quintic_2dx *spatial_sampler, mui_chrono_sampler_gauss_2dx *temporal_sampler);
double mui_fetch_sum_quintic_gauss_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_sum_quintic_2t *spatial_sampler, mui_chrono_sampler_gauss_2t *temporal_sampler);
float mui_fetch_sum_quintic_mean_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_sum_quintic_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler);
float mui_fetch_sum_quintic_mean_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t,
		mui_sampler_sum_quintic_2fx *spatial_sampler, mui_chrono_sampler_mean_2fx *temporal_sampler);
double mui_fetch_sum_quintic_mean_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_sum_quintic_2d *spatial_sampler, mui_chrono_sampler_mean_2d *temporal_sampler);
double mui_fetch_sum_quintic_mean_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t,
		mui_sampler_sum_quintic_2dx *spatial_sampler, mui_chrono_sampler_mean_2dx *temporal_sampler);
double mui_fetch_sum_quintic_mean_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_sum_quintic_2t *spatial_sampler, mui_chrono_sampler_mean_2t *temporal_sampler);
float mui_fetch_sum_quintic_sum_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_sum_quintic_2f *spatial_sampler, mui_chrono_sampler_sum_2f *temporal_sampler);
float mui_fetch_sum_quintic_sum_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t,
		mui_sampler_sum_quintic_2fx *spatial_sampler, mui_chrono_sampler_sum_2fx *temporal_sampler);
double mui_fetch_sum_quintic_sum_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_sum_quintic_2d *spatial_sampler, mui_chrono_sampler_sum_2d *temporal_sampler);
double mui_fetch_sum_quintic_sum_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t,
		mui_sampler_sum_quintic_2dx *spatial_sampler, mui_chrono_sampler_sum_2dx *temporal_sampler);
double mui_fetch_sum_quintic_sum_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_sum_quintic_2t *spatial_sampler, mui_chrono_sampler_sum_2t *temporal_sampler);
#ifdef USE_RBF
float mui_fetch_rbf_exact_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_rbf_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler);
float mui_fetch_rbf_exact_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t,
		mui_sampler_rbf_2fx *spatial_sampler, mui_chrono_sampler_exact_2fx *temporal_sampler);
double mui_fetch_rbf_exact_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_rbf_2d *spatial_sampler, mui_chrono_sampler_exact_2d *temporal_sampler);
double mui_fetch_rbf_exact_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t,
		mui_sampler_rbf_2dx *spatial_sampler, mui_chrono_sampler_exact_2dx *temporal_sampler);
double mui_fetch_rbf_exact_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_rbf_2t *spatial_sampler, mui_chrono_sampler_exact_2t *temporal_sampler);
float mui_fetch_rbf_gauss_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_rbf_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler);
float mui_fetch_rbf_gauss_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t,
		mui_sampler_rbf_2fx *spatial_sampler, mui_chrono_sampler_gauss_2fx *temporal_sampler);
double mui_fetch_rbf_gauss_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_rbf_2d *spatial_sampler, mui_chrono_sampler_gauss_2d *temporal_sampler);
double mui_fetch_rbf_gauss_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t,
		mui_sampler_rbf_2dx *spatial_sampler, mui_chrono_sampler_gauss_2dx *temporal_sampler);
double mui_fetch_rbf_gauss_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_rbf_2t *spatial_sampler, mui_chrono_sampler_gauss_2t *temporal_sampler);
float mui_fetch_rbf_mean_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_rbf_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler);
float mui_fetch_rbf_mean_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t,
		mui_sampler_rbf_2fx *spatial_sampler, mui_chrono_sampler_mean_2fx *temporal_sampler);
double mui_fetch_rbf_mean_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_rbf_2d *spatial_sampler, mui_chrono_sampler_mean_2d *temporal_sampler);
double mui_fetch_rbf_mean_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t,
		mui_sampler_rbf_2dx *spatial_sampler, mui_chrono_sampler_mean_2dx *temporal_sampler);
double mui_fetch_rbf_mean_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_rbf_2t *spatial_sampler, mui_chrono_sampler_mean_2t *temporal_sampler);
float mui_fetch_rbf_sum_2f(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t,
		mui_sampler_rbf_2f *spatial_sampler, mui_chrono_sampler_sum_2f *temporal_sampler);
float mui_fetch_rbf_sum_2fx(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t,
		mui_sampler_rbf_2fx *spatial_sampler, mui_chrono_sampler_sum_2fx *temporal_sampler);
double mui_fetch_rbf_sum_2d(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t,
		mui_sampler_rbf_2d *spatial_sampler, mui_chrono_sampler_sum_2d *temporal_sampler);
double mui_fetch_rbf_sum_2dx(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t,
		mui_sampler_rbf_2dx *spatial_sampler, mui_chrono_sampler_sum_2dx *temporal_sampler);
double mui_fetch_rbf_sum_2t(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t,
		mui_sampler_rbf_2t *spatial_sampler, mui_chrono_sampler_sum_2t *temporal_sampler);
#endif

/********************************************************
 * MUI functions for 2D data fetch using two time values *
 *********************************************************/

float mui_fetch_exact_exact_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t_1, float t_2,
		mui_sampler_exact_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler);
float mui_fetch_exact_exact_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t_1,
		float t_2, mui_sampler_exact_2fx *spatial_sampler, mui_chrono_sampler_exact_2fx *temporal_sampler);
double mui_fetch_exact_exact_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t_1,
		double t_2, mui_sampler_exact_2d *spatial_sampler, mui_chrono_sampler_exact_2d *temporal_sampler);
double mui_fetch_exact_exact_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t_1,
		double t_2, mui_sampler_exact_2dx *spatial_sampler, mui_chrono_sampler_exact_2dx *temporal_sampler);
double mui_fetch_exact_exact_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t_1,
		double t_2, mui_sampler_exact_2t *spatial_sampler, mui_chrono_sampler_exact_2t *temporal_sampler);
float mui_fetch_exact_gauss_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t_1, float t_2,
		mui_sampler_exact_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler);
float mui_fetch_exact_gauss_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t_1,
		float t_2, mui_sampler_exact_2fx *spatial_sampler, mui_chrono_sampler_gauss_2fx *temporal_sampler);
double mui_fetch_exact_gauss_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t_1,
		double t_2, mui_sampler_exact_2d *spatial_sampler, mui_chrono_sampler_gauss_2d *temporal_sampler);
double mui_fetch_exact_gauss_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t_1,
		double t_2, mui_sampler_exact_2dx *spatial_sampler, mui_chrono_sampler_gauss_2dx *temporal_sampler);
double mui_fetch_exact_gauss_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t_1,
		double t_2, mui_sampler_exact_2t *spatial_sampler, mui_chrono_sampler_gauss_2t *temporal_sampler);
float mui_fetch_exact_mean_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t_1, float t_2,
		mui_sampler_exact_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler);
float mui_fetch_exact_mean_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t_1,
		float t_2, mui_sampler_exact_2fx *spatial_sampler, mui_chrono_sampler_mean_2fx *temporal_sampler);
double mui_fetch_exact_mean_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t_1,
		double t_2, mui_sampler_exact_2d *spatial_sampler, mui_chrono_sampler_mean_2d *temporal_sampler);
double mui_fetch_exact_mean_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t_1,
		double t_2, mui_sampler_exact_2dx *spatial_sampler, mui_chrono_sampler_mean_2dx *temporal_sampler);
double mui_fetch_exact_mean_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t_1,
		double t_2, mui_sampler_exact_2t *spatial_sampler, mui_chrono_sampler_mean_2t *temporal_sampler);
float mui_fetch_exact_sum_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t_1, float t_2,
		mui_sampler_exact_2f *spatial_sampler, mui_chrono_sampler_sum_2f *temporal_sampler);
float mui_fetch_exact_sum_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t_1,
		float t_2, mui_sampler_exact_2fx *spatial_sampler, mui_chrono_sampler_sum_2fx *temporal_sampler);
double mui_fetch_exact_sum_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t_1,
		double t_2, mui_sampler_exact_2d *spatial_sampler, mui_chrono_sampler_sum_2d *temporal_sampler);
double mui_fetch_exact_sum_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t_1,
		double t_2, mui_sampler_exact_2dx *spatial_sampler, mui_chrono_sampler_sum_2dx *temporal_sampler);
double mui_fetch_exact_sum_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t_1,
		double t_2, mui_sampler_exact_2t *spatial_sampler, mui_chrono_sampler_sum_2t *temporal_sampler);
float mui_fetch_gauss_exact_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t_1, float t_2,
		mui_sampler_gauss_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler);
float mui_fetch_gauss_exact_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t_1,
		float t_2, mui_sampler_gauss_2fx *spatial_sampler, mui_chrono_sampler_exact_2fx *temporal_sampler);
double mui_fetch_gauss_exact_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t_1,
		double t_2, mui_sampler_gauss_2d *spatial_sampler, mui_chrono_sampler_exact_2d *temporal_sampler);
double mui_fetch_gauss_exact_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t_1,
		double t_2, mui_sampler_gauss_2dx *spatial_sampler, mui_chrono_sampler_exact_2dx *temporal_sampler);
double mui_fetch_gauss_exact_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t_1,
		double t_2, mui_sampler_gauss_2t *spatial_sampler, mui_chrono_sampler_exact_2t *temporal_sampler);
float mui_fetch_gauss_gauss_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t_1, float t_2,
		mui_sampler_gauss_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler);
float mui_fetch_gauss_gauss_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t_1,
		float t_2, mui_sampler_gauss_2fx *spatial_sampler, mui_chrono_sampler_gauss_2fx *temporal_sampler);
double mui_fetch_gauss_gauss_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t_1,
		double t_2, mui_sampler_gauss_2d *spatial_sampler, mui_chrono_sampler_gauss_2d *temporal_sampler);
double mui_fetch_gauss_gauss_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t_1,
		double t_2, mui_sampler_gauss_2dx *spatial_sampler, mui_chrono_sampler_gauss_2dx *temporal_sampler);
double mui_fetch_gauss_gauss_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t_1,
		double t_2, mui_sampler_gauss_2t *spatial_sampler, mui_chrono_sampler_gauss_2t *temporal_sampler);
float mui_fetch_gauss_mean_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t_1, float t_2,
		mui_sampler_gauss_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler);
float mui_fetch_gauss_mean_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t_1,
		float t_2, mui_sampler_gauss_2fx *spatial_sampler, mui_chrono_sampler_mean_2fx *temporal_sampler);
double mui_fetch_gauss_mean_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t_1,
		double t_2, mui_sampler_gauss_2d *spatial_sampler, mui_chrono_sampler_mean_2d *temporal_sampler);
double mui_fetch_gauss_mean_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t_1,
		double t_2, mui_sampler_gauss_2dx *spatial_sampler, mui_chrono_sampler_mean_2dx *temporal_sampler);
double mui_fetch_gauss_mean_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t_1,
		double t_2, mui_sampler_gauss_2t *spatial_sampler, mui_chrono_sampler_mean_2t *temporal_sampler);
float mui_fetch_moving_average_exact_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t_1,
		float t_2, mui_sampler_moving_average_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler);
float mui_fetch_moving_average_exact_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point,
		float t_1, float t_2, mui_sampler_moving_average_2fx *spatial_sampler,
		mui_chrono_sampler_exact_2fx *temporal_sampler);
double mui_fetch_moving_average_exact_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t_1,
		double t_2, mui_sampler_moving_average_2d *spatial_sampler, mui_chrono_sampler_exact_2d *temporal_sampler);
double mui_fetch_moving_average_exact_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point,
		double t_1, double t_2, mui_sampler_moving_average_2dx *spatial_sampler,
		mui_chrono_sampler_exact_2dx *temporal_sampler);
double mui_fetch_moving_average_exact_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t_1,
		double t_2, mui_sampler_moving_average_2t *spatial_sampler, mui_chrono_sampler_exact_2t *temporal_sampler);
float mui_fetch_moving_average_gauss_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t_1,
		float t_2, mui_sampler_moving_average_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler);
float mui_fetch_moving_average_gauss_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point,
		float t_1, float t_2, mui_sampler_moving_average_2fx *spatial_sampler,
		mui_chrono_sampler_gauss_2fx *temporal_sampler);
double mui_fetch_moving_average_gauss_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t_1,
		double t_2, mui_sampler_moving_average_2d *spatial_sampler, mui_chrono_sampler_gauss_2d *temporal_sampler);
double mui_fetch_moving_average_gauss_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point,
		double t_1, double t_2, mui_sampler_moving_average_2dx *spatial_sampler,
		mui_chrono_sampler_gauss_2dx *temporal_sampler);
double mui_fetch_moving_average_gauss_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t_1,
		double t_2, mui_sampler_moving_average_2t *spatial_sampler, mui_chrono_sampler_gauss_2t *temporal_sampler);
float mui_fetch_moving_average_mean_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t_1,
		float t_2, mui_sampler_moving_average_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler);
float mui_fetch_moving_average_mean_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t_1,
		float t_2, mui_sampler_moving_average_2fx *spatial_sampler, mui_chrono_sampler_mean_2fx *temporal_sampler);
double mui_fetch_moving_average_mean_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t_1,
		double t_2, mui_sampler_moving_average_2d *spatial_sampler, mui_chrono_sampler_mean_2d *temporal_sampler);
double mui_fetch_moving_average_mean_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point,
		double t_1, double t_2, mui_sampler_moving_average_2dx *spatial_sampler,
		mui_chrono_sampler_mean_2dx *temporal_sampler);
double mui_fetch_moving_average_mean_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t_1,
		double t_2, mui_sampler_moving_average_2t *spatial_sampler, mui_chrono_sampler_mean_2t *temporal_sampler);
float mui_fetch_moving_average_sum_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t_1,
		float t_2, mui_sampler_moving_average_2f *spatial_sampler, mui_chrono_sampler_sum_2f *temporal_sampler);
float mui_fetch_moving_average_sum_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t_1,
		float t_2, mui_sampler_moving_average_2fx *spatial_sampler, mui_chrono_sampler_sum_2fx *temporal_sampler);
double mui_fetch_moving_average_sum_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t_1,
		double t_2, mui_sampler_moving_average_2d *spatial_sampler, mui_chrono_sampler_sum_2d *temporal_sampler);
double mui_fetch_moving_average_sum_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point,
		double t_1, double t_2, mui_sampler_moving_average_2dx *spatial_sampler,
		mui_chrono_sampler_sum_2dx *temporal_sampler);
double mui_fetch_moving_average_sum_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t_1,
		double t_2, mui_sampler_moving_average_2t *spatial_sampler, mui_chrono_sampler_sum_2t *temporal_sampler);
float mui_fetch_nearest_neighbor_exact_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t_1,
		float t_2, mui_sampler_nearest_neighbor_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler);
float mui_fetch_nearest_neighbor_exact_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point,
		float t_1, float t_2, mui_sampler_nearest_neighbor_2fx *spatial_sampler,
		mui_chrono_sampler_exact_2fx *temporal_sampler);
double mui_fetch_nearest_neighbor_exact_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point,
		double t_1, double t_2, mui_sampler_nearest_neighbor_2d *spatial_sampler,
		mui_chrono_sampler_exact_2d *temporal_sampler);
double mui_fetch_nearest_neighbor_exact_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point,
		double t_1, double t_2, mui_sampler_nearest_neighbor_2dx *spatial_sampler,
		mui_chrono_sampler_exact_2dx *temporal_sampler);
double mui_fetch_nearest_neighbor_exact_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point,
		double t_1, double t_2, mui_sampler_nearest_neighbor_2t *spatial_sampler,
		mui_chrono_sampler_exact_2t *temporal_sampler);
float mui_fetch_nearest_neighbor_gauss_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t_1,
		float t_2, mui_sampler_nearest_neighbor_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler);
float mui_fetch_nearest_neighbor_gauss_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point,
		float t_1, float t_2, mui_sampler_nearest_neighbor_2fx *spatial_sampler,
		mui_chrono_sampler_gauss_2fx *temporal_sampler);
double mui_fetch_nearest_neighbor_gauss_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point,
		double t_1, double t_2, mui_sampler_nearest_neighbor_2d *spatial_sampler,
		mui_chrono_sampler_gauss_2d *temporal_sampler);
double mui_fetch_nearest_neighbor_gauss_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point,
		double t_1, double t_2, mui_sampler_nearest_neighbor_2dx *spatial_sampler,
		mui_chrono_sampler_gauss_2dx *temporal_sampler);
double mui_fetch_nearest_neighbor_gauss_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point,
		double t_1, double t_2, mui_sampler_nearest_neighbor_2t *spatial_sampler,
		mui_chrono_sampler_gauss_2t *temporal_sampler);
float mui_fetch_nearest_neighbor_mean_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t_1,
		float t_2, mui_sampler_nearest_neighbor_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler);
float mui_fetch_nearest_neighbor_mean_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point,
		float t_1, float t_2, mui_sampler_nearest_neighbor_2fx *spatial_sampler,
		mui_chrono_sampler_mean_2fx *temporal_sampler);
double mui_fetch_nearest_neighbor_mean_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point,
		double t_1, double t_2, mui_sampler_nearest_neighbor_2d *spatial_sampler,
		mui_chrono_sampler_mean_2d *temporal_sampler);
double mui_fetch_nearest_neighbor_mean_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point,
		double t_1, double t_2, mui_sampler_nearest_neighbor_2dx *spatial_sampler,
		mui_chrono_sampler_mean_2dx *temporal_sampler);
double mui_fetch_nearest_neighbor_mean_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point,
		double t_1, double t_2, mui_sampler_nearest_neighbor_2t *spatial_sampler,
		mui_chrono_sampler_mean_2t *temporal_sampler);
float mui_fetch_nearest_neighbor_sum_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t_1,
		float t_2, mui_sampler_nearest_neighbor_2f *spatial_sampler, mui_chrono_sampler_sum_2f *temporal_sampler);
float mui_fetch_nearest_neighbor_sum_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point,
		float t_1, float t_2, mui_sampler_nearest_neighbor_2fx *spatial_sampler,
		mui_chrono_sampler_sum_2fx *temporal_sampler);
double mui_fetch_nearest_neighbor_sum_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t_1,
		double t_2, mui_sampler_nearest_neighbor_2d *spatial_sampler, mui_chrono_sampler_sum_2d *temporal_sampler);
double mui_fetch_nearest_neighbor_sum_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point,
		double t_1, double t_2, mui_sampler_nearest_neighbor_2dx *spatial_sampler,
		mui_chrono_sampler_sum_2dx *temporal_sampler);
double mui_fetch_nearest_neighbor_sum_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t_1,
		double t_2, mui_sampler_nearest_neighbor_2t *spatial_sampler, mui_chrono_sampler_sum_2t *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_exact_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point,
		float t_1, float t_2, mui_sampler_pseudo_nearest_neighbor_2f *spatial_sampler,
		mui_chrono_sampler_exact_2f *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_exact_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point,
		float t_1, float t_2, mui_sampler_pseudo_nearest_neighbor_2fx *spatial_sampler,
		mui_chrono_sampler_exact_2fx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_exact_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_2d *spatial_sampler,
		mui_chrono_sampler_exact_2d *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_exact_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_2dx *spatial_sampler,
		mui_chrono_sampler_exact_2dx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_exact_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_2t *spatial_sampler,
		mui_chrono_sampler_exact_2t *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_gauss_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point,
		float t_1, float t_2, mui_sampler_pseudo_nearest_neighbor_2f *spatial_sampler,
		mui_chrono_sampler_gauss_2f *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_gauss_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point,
		float t_1, float t_2, mui_sampler_pseudo_nearest_neighbor_2fx *spatial_sampler,
		mui_chrono_sampler_gauss_2fx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_gauss_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_2d *spatial_sampler,
		mui_chrono_sampler_gauss_2d *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_gauss_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_2dx *spatial_sampler,
		mui_chrono_sampler_gauss_2dx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_gauss_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_2t *spatial_sampler,
		mui_chrono_sampler_gauss_2t *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_mean_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point,
		float t_1, float t_2, mui_sampler_pseudo_nearest_neighbor_2f *spatial_sampler,
		mui_chrono_sampler_mean_2f *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_mean_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point,
		float t_1, float t_2, mui_sampler_pseudo_nearest_neighbor_2fx *spatial_sampler,
		mui_chrono_sampler_mean_2fx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_mean_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_2d *spatial_sampler,
		mui_chrono_sampler_mean_2d *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_mean_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_2dx *spatial_sampler,
		mui_chrono_sampler_mean_2dx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_mean_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_2t *spatial_sampler,
		mui_chrono_sampler_mean_2t *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_sum_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point,
		float t_1, float t_2, mui_sampler_pseudo_nearest_neighbor_2f *spatial_sampler,
		mui_chrono_sampler_sum_2f *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_sum_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point,
		float t_1, float t_2, mui_sampler_pseudo_nearest_neighbor_2fx *spatial_sampler,
		mui_chrono_sampler_sum_2fx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_sum_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_2d *spatial_sampler,
		mui_chrono_sampler_sum_2d *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_sum_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_2dx *spatial_sampler,
		mui_chrono_sampler_sum_2dx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_sum_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_2t *spatial_sampler,
		mui_chrono_sampler_sum_2t *temporal_sampler);
float mui_fetch_shepard_quintic_exact_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t_1,
		float t_2, mui_sampler_shepard_quintic_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler);
float mui_fetch_shepard_quintic_exact_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point,
		float t_1, float t_2, mui_sampler_shepard_quintic_2fx *spatial_sampler,
		mui_chrono_sampler_exact_2fx *temporal_sampler);
double mui_fetch_shepard_quintic_exact_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point,
		double t_1, double t_2, mui_sampler_shepard_quintic_2d *spatial_sampler,
		mui_chrono_sampler_exact_2d *temporal_sampler);
double mui_fetch_shepard_quintic_exact_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point,
		double t_1, double t_2, mui_sampler_shepard_quintic_2dx *spatial_sampler,
		mui_chrono_sampler_exact_2dx *temporal_sampler);
double mui_fetch_shepard_quintic_exact_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point,
		double t_1, double t_2, mui_sampler_shepard_quintic_2t *spatial_sampler,
		mui_chrono_sampler_exact_2t *temporal_sampler);
float mui_fetch_shepard_quintic_gauss_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t_1,
		float t_2, mui_sampler_shepard_quintic_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler);
float mui_fetch_shepard_quintic_gauss_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point,
		float t_1, float t_2, mui_sampler_shepard_quintic_2fx *spatial_sampler,
		mui_chrono_sampler_gauss_2fx *temporal_sampler);
double mui_fetch_shepard_quintic_gauss_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point,
		double t_1, double t_2, mui_sampler_shepard_quintic_2d *spatial_sampler,
		mui_chrono_sampler_gauss_2d *temporal_sampler);
double mui_fetch_shepard_quintic_gauss_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point,
		double t_1, double t_2, mui_sampler_shepard_quintic_2dx *spatial_sampler,
		mui_chrono_sampler_gauss_2dx *temporal_sampler);
double mui_fetch_shepard_quintic_gauss_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point,
		double t_1, double t_2, mui_sampler_shepard_quintic_2t *spatial_sampler,
		mui_chrono_sampler_gauss_2t *temporal_sampler);
float mui_fetch_shepard_quintic_mean_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t_1,
		float t_2, mui_sampler_shepard_quintic_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler);
float mui_fetch_shepard_quintic_mean_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point,
		float t_1, float t_2, mui_sampler_shepard_quintic_2fx *spatial_sampler,
		mui_chrono_sampler_mean_2fx *temporal_sampler);
double mui_fetch_shepard_quintic_mean_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t_1,
		double t_2, mui_sampler_shepard_quintic_2d *spatial_sampler, mui_chrono_sampler_mean_2d *temporal_sampler);
double mui_fetch_shepard_quintic_mean_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point,
		double t_1, double t_2, mui_sampler_shepard_quintic_2dx *spatial_sampler,
		mui_chrono_sampler_mean_2dx *temporal_sampler);
double mui_fetch_shepard_quintic_mean_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t_1,
		double t_2, mui_sampler_shepard_quintic_2t *spatial_sampler, mui_chrono_sampler_mean_2t *temporal_sampler);
float mui_fetch_shepard_quintic_sum_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t_1,
		float t_2, mui_sampler_shepard_quintic_2f *spatial_sampler, mui_chrono_sampler_sum_2f *temporal_sampler);
float mui_fetch_shepard_quintic_sum_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t_1,
		float t_2, mui_sampler_shepard_quintic_2fx *spatial_sampler, mui_chrono_sampler_sum_2fx *temporal_sampler);
double mui_fetch_shepard_quintic_sum_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t_1,
		double t_2, mui_sampler_shepard_quintic_2d *spatial_sampler, mui_chrono_sampler_sum_2d *temporal_sampler);
double mui_fetch_shepard_quintic_sum_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point,
		double t_1, double t_2, mui_sampler_shepard_quintic_2dx *spatial_sampler,
		mui_chrono_sampler_sum_2dx *temporal_sampler);
double mui_fetch_shepard_quintic_sum_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t_1,
		double t_2, mui_sampler_shepard_quintic_2t *spatial_sampler, mui_chrono_sampler_sum_2t *temporal_sampler);
float mui_fetch_sph_quintic_exact_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t_1,
		float t_2, mui_sampler_sph_quintic_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler);
float mui_fetch_sph_quintic_exact_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t_1,
		float t_2, mui_sampler_sph_quintic_2fx *spatial_sampler, mui_chrono_sampler_exact_2fx *temporal_sampler);
double mui_fetch_sph_quintic_exact_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t_1,
		double t_2, mui_sampler_sph_quintic_2d *spatial_sampler, mui_chrono_sampler_exact_2d *temporal_sampler);
double mui_fetch_sph_quintic_exact_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t_1,
		double t_2, mui_sampler_sph_quintic_2dx *spatial_sampler, mui_chrono_sampler_exact_2dx *temporal_sampler);
double mui_fetch_sph_quintic_exact_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t_1,
		double t_2, mui_sampler_sph_quintic_2t *spatial_sampler, mui_chrono_sampler_exact_2t *temporal_sampler);
float mui_fetch_sph_quintic_gauss_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t_1,
		float t_2, mui_sampler_sph_quintic_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler);
float mui_fetch_sph_quintic_gauss_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t_1,
		float t_2, mui_sampler_sph_quintic_2fx *spatial_sampler, mui_chrono_sampler_gauss_2fx *temporal_sampler);
double mui_fetch_sph_quintic_gauss_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t_1,
		double t_2, mui_sampler_sph_quintic_2d *spatial_sampler, mui_chrono_sampler_gauss_2d *temporal_sampler);
double mui_fetch_sph_quintic_gauss_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t_1,
		double t_2, mui_sampler_sph_quintic_2dx *spatial_sampler, mui_chrono_sampler_gauss_2dx *temporal_sampler);
double mui_fetch_sph_quintic_gauss_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t_1,
		double t_2, mui_sampler_sph_quintic_2t *spatial_sampler, mui_chrono_sampler_gauss_2t *temporal_sampler);
float mui_fetch_sph_quintic_mean_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t_1,
		float t_2, mui_sampler_sph_quintic_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler);
float mui_fetch_sph_quintic_mean_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t_1,
		float t_2, mui_sampler_sph_quintic_2fx *spatial_sampler, mui_chrono_sampler_mean_2fx *temporal_sampler);
double mui_fetch_sph_quintic_mean_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t_1,
		double t_2, mui_sampler_sph_quintic_2d *spatial_sampler, mui_chrono_sampler_mean_2d *temporal_sampler);
double mui_fetch_sph_quintic_mean_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t_1,
		double t_2, mui_sampler_sph_quintic_2dx *spatial_sampler, mui_chrono_sampler_mean_2dx *temporal_sampler);
double mui_fetch_sph_quintic_mean_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t_1,
		double t_2, mui_sampler_sph_quintic_2t *spatial_sampler, mui_chrono_sampler_mean_2t *temporal_sampler);
float mui_fetch_sph_quintic_sum_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t_1,
		float t_2, mui_sampler_sph_quintic_2f *spatial_sampler, mui_chrono_sampler_sum_2f *temporal_sampler);
float mui_fetch_sph_quintic_sum_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t_1,
		float t_2, mui_sampler_sph_quintic_2fx *spatial_sampler, mui_chrono_sampler_sum_2fx *temporal_sampler);
double mui_fetch_sph_quintic_sum_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t_1,
		double t_2, mui_sampler_sph_quintic_2d *spatial_sampler, mui_chrono_sampler_sum_2d *temporal_sampler);
double mui_fetch_sph_quintic_sum_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t_1,
		double t_2, mui_sampler_sph_quintic_2dx *spatial_sampler, mui_chrono_sampler_sum_2dx *temporal_sampler);
double mui_fetch_sph_quintic_sum_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t_1,
		double t_2, mui_sampler_sph_quintic_2t *spatial_sampler, mui_chrono_sampler_sum_2t *temporal_sampler);
float mui_fetch_sum_quintic_exact_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t_1,
		float t_2, mui_sampler_sum_quintic_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler);
float mui_fetch_sum_quintic_exact_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t_1,
		float t_2, mui_sampler_sum_quintic_2fx *spatial_sampler, mui_chrono_sampler_exact_2fx *temporal_sampler);
double mui_fetch_sum_quintic_exact_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t_1,
		double t_2, mui_sampler_sum_quintic_2d *spatial_sampler, mui_chrono_sampler_exact_2d *temporal_sampler);
double mui_fetch_sum_quintic_exact_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t_1,
		double t_2, mui_sampler_sum_quintic_2dx *spatial_sampler, mui_chrono_sampler_exact_2dx *temporal_sampler);
double mui_fetch_sum_quintic_exact_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t_1,
		double t_2, mui_sampler_sum_quintic_2t *spatial_sampler, mui_chrono_sampler_exact_2t *temporal_sampler);
float mui_fetch_sum_quintic_gauss_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t_1,
		float t_2, mui_sampler_sum_quintic_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler);
float mui_fetch_sum_quintic_gauss_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t_1,
		float t_2, mui_sampler_sum_quintic_2fx *spatial_sampler, mui_chrono_sampler_gauss_2fx *temporal_sampler);
double mui_fetch_sum_quintic_gauss_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t_1,
		double t_2, mui_sampler_sum_quintic_2d *spatial_sampler, mui_chrono_sampler_gauss_2d *temporal_sampler);
double mui_fetch_sum_quintic_gauss_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t_1,
		double t_2, mui_sampler_sum_quintic_2dx *spatial_sampler, mui_chrono_sampler_gauss_2dx *temporal_sampler);
double mui_fetch_sum_quintic_gauss_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t_1,
		double t_2, mui_sampler_sum_quintic_2t *spatial_sampler, mui_chrono_sampler_gauss_2t *temporal_sampler);
float mui_fetch_sum_quintic_mean_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t_1,
		float t_2, mui_sampler_sum_quintic_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler);
float mui_fetch_sum_quintic_mean_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t_1,
		float t_2, mui_sampler_sum_quintic_2fx *spatial_sampler, mui_chrono_sampler_mean_2fx *temporal_sampler);
double mui_fetch_sum_quintic_mean_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t_1,
		double t_2, mui_sampler_sum_quintic_2d *spatial_sampler, mui_chrono_sampler_mean_2d *temporal_sampler);
double mui_fetch_sum_quintic_mean_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t_1,
		double t_2, mui_sampler_sum_quintic_2dx *spatial_sampler, mui_chrono_sampler_mean_2dx *temporal_sampler);
double mui_fetch_sum_quintic_mean_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t_1,
		double t_2, mui_sampler_sum_quintic_2t *spatial_sampler, mui_chrono_sampler_mean_2t *temporal_sampler);
float mui_fetch_sum_quintic_sum_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t_1,
		float t_2, mui_sampler_sum_quintic_2f *spatial_sampler, mui_chrono_sampler_sum_2f *temporal_sampler);
float mui_fetch_sum_quintic_sum_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t_1,
		float t_2, mui_sampler_sum_quintic_2fx *spatial_sampler, mui_chrono_sampler_sum_2fx *temporal_sampler);
double mui_fetch_sum_quintic_sum_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t_1,
		double t_2, mui_sampler_sum_quintic_2d *spatial_sampler, mui_chrono_sampler_sum_2d *temporal_sampler);
double mui_fetch_sum_quintic_sum_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t_1,
		double t_2, mui_sampler_sum_quintic_2dx *spatial_sampler, mui_chrono_sampler_sum_2dx *temporal_sampler);
double mui_fetch_sum_quintic_sum_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t_1,
		double t_2, mui_sampler_sum_quintic_2t *spatial_sampler, mui_chrono_sampler_sum_2t *temporal_sampler);
#ifdef USE_RBF
float mui_fetch_rbf_exact_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t_1, float t_2,
		mui_sampler_rbf_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler);
float mui_fetch_rbf_exact_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t_1,
		float t_2, mui_sampler_rbf_2fx *spatial_sampler, mui_chrono_sampler_exact_2fx *temporal_sampler);
double mui_fetch_rbf_exact_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t_1,
		double t_2, mui_sampler_rbf_2d *spatial_sampler, mui_chrono_sampler_exact_2d *temporal_sampler);
double mui_fetch_rbf_exact_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t_1,
		double t_2, mui_sampler_rbf_2dx *spatial_sampler, mui_chrono_sampler_exact_2dx *temporal_sampler);
double mui_fetch_rbf_exact_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t_1,
		double t_2, mui_sampler_rbf_2t *spatial_sampler, mui_chrono_sampler_exact_2t *temporal_sampler);
float mui_fetch_rbf_gauss_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t_1, float t_2,
		mui_sampler_rbf_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler);
float mui_fetch_rbf_gauss_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t_1,
		float t_2, mui_sampler_rbf_2fx *spatial_sampler, mui_chrono_sampler_gauss_2fx *temporal_sampler);
double mui_fetch_rbf_gauss_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t_1,
		double t_2, mui_sampler_rbf_2d *spatial_sampler, mui_chrono_sampler_gauss_2d *temporal_sampler);
double mui_fetch_rbf_gauss_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t_1,
		double t_2, mui_sampler_rbf_2dx *spatial_sampler, mui_chrono_sampler_gauss_2dx *temporal_sampler);
double mui_fetch_rbf_gauss_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t_1,
		double t_2, mui_sampler_rbf_2t *spatial_sampler, mui_chrono_sampler_gauss_2t *temporal_sampler);
float mui_fetch_rbf_mean_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t_1, float t_2,
		mui_sampler_rbf_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler);
float mui_fetch_rbf_mean_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t_1, float t_2,
		mui_sampler_rbf_2fx *spatial_sampler, mui_chrono_sampler_mean_2fx *temporal_sampler);
double mui_fetch_rbf_mean_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t_1, double t_2,
		mui_sampler_rbf_2d *spatial_sampler, mui_chrono_sampler_mean_2d *temporal_sampler);
double mui_fetch_rbf_mean_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t_1,
		double t_2, mui_sampler_rbf_2dx *spatial_sampler, mui_chrono_sampler_mean_2dx *temporal_sampler);
double mui_fetch_rbf_mean_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t_1, double t_2,
		mui_sampler_rbf_2t *spatial_sampler, mui_chrono_sampler_mean_2t *temporal_sampler);
float mui_fetch_rbf_sum_2f_pair(mui_uniface_2f *uniface, const char *attr, mui_point_2f point, float t_1, float t_2,
		mui_sampler_rbf_2f *spatial_sampler, mui_chrono_sampler_sum_2f *temporal_sampler);
float mui_fetch_rbf_sum_2fx_pair(mui_uniface_2fx *uniface, const char *attr, mui_point_2fx point, float t_1, float t_2,
		mui_sampler_rbf_2fx *spatial_sampler, mui_chrono_sampler_sum_2fx *temporal_sampler);
double mui_fetch_rbf_sum_2d_pair(mui_uniface_2d *uniface, const char *attr, mui_point_2d point, double t_1, double t_2,
		mui_sampler_rbf_2d *spatial_sampler, mui_chrono_sampler_sum_2d *temporal_sampler);
double mui_fetch_rbf_sum_2dx_pair(mui_uniface_2dx *uniface, const char *attr, mui_point_2dx point, double t_1,
		double t_2, mui_sampler_rbf_2dx *spatial_sampler, mui_chrono_sampler_sum_2dx *temporal_sampler);
double mui_fetch_rbf_sum_2t_pair(mui_uniface_2t *uniface, const char *attr, mui_point_2t point, double t_1, double t_2,
		mui_sampler_rbf_2t *spatial_sampler, mui_chrono_sampler_sum_2t *temporal_sampler);
#endif

// MUI functions for 2D data point only fetch using single time value
void mui_fetch_points_exact_2f(mui_uniface_2f *uniface, const char *attr, float t,
		mui_chrono_sampler_exact_2f *temporal_sampler, mui_point_2f *ret_points, int *num_points);
void mui_fetch_points_exact_2fx(mui_uniface_2fx *uniface, const char *attr, float t,
		mui_chrono_sampler_exact_2fx *temporal_sampler, mui_point_2fx *ret_points, int *num_points);
void mui_fetch_points_exact_2d(mui_uniface_2d *uniface, const char *attr, double t,
		mui_chrono_sampler_exact_2d *temporal_sampler, mui_point_2d *ret_points, int *num_points);
void mui_fetch_points_exact_2dx(mui_uniface_2dx *uniface, const char *attr, double t,
		mui_chrono_sampler_exact_2dx *temporal_sampler, mui_point_2dx *ret_points, int *num_points);
void mui_fetch_points_exact_2t(mui_uniface_2t *uniface, const char *attr, double t,
		mui_chrono_sampler_exact_2t *temporal_sampler, mui_point_2t *ret_points, int *num_points);
void mui_fetch_points_gauss_2f(mui_uniface_2f *uniface, const char *attr, float t,
		mui_chrono_sampler_gauss_2f *temporal_sampler, mui_point_2f *ret_points, int *num_points);
void mui_fetch_points_gauss_2fx(mui_uniface_2fx *uniface, const char *attr, float t,
		mui_chrono_sampler_gauss_2fx *temporal_sampler, mui_point_2fx *ret_points, int *num_points);
void mui_fetch_points_gauss_2d(mui_uniface_2d *uniface, const char *attr, double t,
		mui_chrono_sampler_gauss_2d *temporal_sampler, mui_point_2d *ret_points, int *num_points);
void mui_fetch_points_gauss_2dx(mui_uniface_2dx *uniface, const char *attr, double t,
		mui_chrono_sampler_gauss_2dx *temporal_sampler, mui_point_2dx *ret_points, int *num_points);
void mui_fetch_points_gauss_2t(mui_uniface_2t *uniface, const char *attr, double t,
		mui_chrono_sampler_gauss_2t *temporal_sampler, mui_point_2t *ret_points, int *num_points);
void mui_fetch_points_mean_2f(mui_uniface_2f *uniface, const char *attr, float t,
		mui_chrono_sampler_mean_2f *temporal_sampler, mui_point_2f *ret_points, int *num_points);
void mui_fetch_points_mean_2fx(mui_uniface_2fx *uniface, const char *attr, float t,
		mui_chrono_sampler_mean_2fx *temporal_sampler, mui_point_2fx *ret_points, int *num_points);
void mui_fetch_points_mean_2d(mui_uniface_2d *uniface, const char *attr, double t,
		mui_chrono_sampler_mean_2d *temporal_sampler, mui_point_2d *ret_points, int *num_points);
void mui_fetch_points_mean_2dx(mui_uniface_2dx *uniface, const char *attr, double t,
		mui_chrono_sampler_mean_2dx *temporal_sampler, mui_point_2dx *ret_points, int *num_points);
void mui_fetch_points_mean_2t(mui_uniface_2t *uniface, const char *attr, double t,
		mui_chrono_sampler_mean_2t *temporal_sampler, mui_point_2t *ret_points, int *num_points);
void mui_fetch_points_sum_2f(mui_uniface_2f *uniface, const char *attr, float t,
		mui_chrono_sampler_sum_2f *temporal_sampler, mui_point_2f *ret_points, int *num_points);
void mui_fetch_points_sum_2fx(mui_uniface_2fx *uniface, const char *attr, float t,
		mui_chrono_sampler_sum_2fx *temporal_sampler, mui_point_2fx *ret_points, int *num_points);
void mui_fetch_points_sum_2d(mui_uniface_2d *uniface, const char *attr, double t,
		mui_chrono_sampler_sum_2d *temporal_sampler, mui_point_2d *ret_points, int *num_points);
void mui_fetch_points_sum_2dx(mui_uniface_2dx *uniface, const char *attr, double t,
		mui_chrono_sampler_sum_2dx *temporal_sampler, mui_point_2dx *ret_points, int *num_points);
void mui_fetch_points_sum_2t(mui_uniface_2t *uniface, const char *attr, double t,
		mui_chrono_sampler_sum_2t *temporal_sampler, mui_point_2t *ret_points, int *num_points);

// MUI functions for 2D data point only fetch using two time values
void mui_fetch_values_exact_2f_pair(mui_uniface_2f *uniface, const char *attr, float t_1, float t_2,
		mui_chrono_sampler_exact_2f *temporal_sampler, float *values, int *num_values);
void mui_fetch_values_exact_2fx_pair(mui_uniface_2fx *uniface, const char *attr, float t_1, float t_2,
		mui_chrono_sampler_exact_2fx *temporal_sampler, float *values, int *num_values);
void mui_fetch_values_exact_2d_pair(mui_uniface_2d *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_exact_2d *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_exact_2dx_pair(mui_uniface_2dx *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_exact_2dx *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_exact_2t_pair(mui_uniface_2t *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_exact_2t *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_gauss_2f_pair(mui_uniface_2f *uniface, const char *attr, float t_1, float t_2,
		mui_chrono_sampler_gauss_2f *temporal_sampler, float *values, int *num_values);
void mui_fetch_values_gauss_2fx_pair(mui_uniface_2fx *uniface, const char *attr, float t_1, float t_2,
		mui_chrono_sampler_gauss_2fx *temporal_sampler, float *values, int *num_values);
void mui_fetch_values_gauss_2d_pair(mui_uniface_2d *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_gauss_2d *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_gauss_2dx_pair(mui_uniface_2dx *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_gauss_2dx *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_gauss_2t_pair(mui_uniface_2t *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_gauss_2t *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_mean_2f_pair(mui_uniface_2f *uniface, const char *attr, float t_1, float t_2,
		mui_chrono_sampler_mean_2f *temporal_sampler, float *values, int *num_values);
void mui_fetch_values_mean_2fx_pair(mui_uniface_2fx *uniface, const char *attr, float t_1, float t_2,
		mui_chrono_sampler_mean_2fx *temporal_sampler, float *values, int *num_values);
void mui_fetch_values_mean_2d_pair(mui_uniface_2d *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_mean_2d *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_mean_2dx_pair(mui_uniface_2dx *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_mean_2dx *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_mean_2t_pair(mui_uniface_2t *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_mean_2t *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_sum_2f_pair(mui_uniface_2f *uniface, const char *attr, float t_1, float t_2,
		mui_chrono_sampler_sum_2f *temporal_sampler, float *values, int *num_values);
void mui_fetch_values_sum_2fx_pair(mui_uniface_2fx *uniface, const char *attr, float t_1, float t_2,
		mui_chrono_sampler_sum_2fx *temporal_sampler, float *values, int *num_values);
void mui_fetch_values_sum_2d_pair(mui_uniface_2d *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_sum_2d *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_sum_2dx_pair(mui_uniface_2dx *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_sum_2dx *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_sum_2t_pair(mui_uniface_2t *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_sum_2t *temporal_sampler, double *values, int *num_values);

// MUI functions for single parameter fetch
float mui_fetch_2f_param(mui_uniface_2f *uniface, const char *attr);
float mui_fetch_2fx_param(mui_uniface_2fx *uniface, const char *attr);
double mui_fetch_2d_param(mui_uniface_2d *uniface, const char *attr);
double mui_fetch_2dx_param(mui_uniface_2dx *uniface, const char *attr);
double mui_fetch_2t_param(mui_uniface_2t *uniface, const char *attr);

// MUI data receive test functions
int mui_is_ready_2f(mui_uniface_2f *uniface, const char *attr, float t);
int mui_is_ready_2fx(mui_uniface_2fx *uniface, const char *attr, float t);
int mui_is_ready_2d(mui_uniface_2d *uniface, const char *attr, double t);
int mui_is_ready_2dx(mui_uniface_2dx *uniface, const char *attr, double t);
int mui_is_ready_2t(mui_uniface_2t *uniface, const char *attr, double t);
int mui_is_ready_2f_pair(mui_uniface_2f *uniface, const char *attr, float t_1, float t_2);
int mui_is_ready_2fx_pair(mui_uniface_2fx *uniface, const char *attr, float t_1, float t_2);
int mui_is_ready_2d_pair(mui_uniface_2d *uniface, const char *attr, double t_1, double t_2);
int mui_is_ready_2dx_pair(mui_uniface_2dx *uniface, const char *attr, double t_1, double t_2);
int mui_is_ready_2t_pair(mui_uniface_2t *uniface, const char *attr, double t_1, double t_2);

// MUI Smart Send functions
void mui_announce_send_span_2f_box(mui_uniface_2f *uniface, float box_1_1, float box_1_2, float box_2_1, float box_2_2,
		float t_start, float t_timeout, int synchronised);
void mui_announce_send_span_2fx_box(mui_uniface_2fx *uniface, float box_1_1, float box_1_2, float box_2_1,
		float box_2_2, float t_start, float t_timeout, int synchronised);
void mui_announce_send_span_2d_box(mui_uniface_2d *uniface, double box_1_1, double box_1_2, double box_2_1,
		double box_2_2, double t_start, double t_timeout, int synchronised);
void mui_announce_send_span_2dx_box(mui_uniface_2dx *uniface, double box_1_1, double box_1_2, double box_2_1,
		double box_2_2, double t_start, double t_timeout, int synchronised);
void mui_announce_send_span_2t_box(mui_uniface_2t *uniface, double box_1_1, double box_1_2, double box_2_1,
		double box_2_2, double t_start, double t_timeout, int synchronised);
void mui_announce_send_span_2f_sphere(mui_uniface_2f *uniface, mui_point_2f centre, float radius, float t_start,
		float t_timeout, int synchronised);
void mui_announce_send_span_2fx_sphere(mui_uniface_2fx *uniface, mui_point_2fx centre, float radius, float t_start,
		float t_timeout, int synchronised);
void mui_announce_send_span_2d_sphere(mui_uniface_2d *uniface, mui_point_2d centre, double radius, double t_start,
		double t_timeout, int synchronised);
void mui_announce_send_span_2dx_sphere(mui_uniface_2dx *uniface, mui_point_2dx centre, double radius, double t_start,
		double t_timeout, int synchronised);
void mui_announce_send_span_2t_sphere(mui_uniface_2t *uniface, mui_point_2t centre, double radius, double t_start,
		double t_timeout, int synchronised);
void mui_announce_recv_span_2f_box(mui_uniface_2f *uniface, float box_1_1, float box_1_2, float box_2_1, float box_2_2,
		float t_start, float t_timeout, int synchronised);
void mui_announce_recv_span_2fx_box(mui_uniface_2fx *uniface, float box_1_1, float box_1_2, float box_2_1,
		float box_2_2, float t_start, float t_timeout, int synchronised);
void mui_announce_recv_span_2d_box(mui_uniface_2d *uniface, double box_1_1, double box_1_2, double box_2_1,
		double box_2_2, double t_start, double t_timeout, int synchronised);
void mui_announce_recv_span_2dx_box(mui_uniface_2dx *uniface, double box_1_1, double box_1_2, double box_2_1,
		double box_2_2, double t_start, double t_timeout, int synchronised);
void mui_announce_recv_span_2t_box(mui_uniface_2t *uniface, double box_1_1, double box_1_2, double box_2_1,
		double box_2_2, double t_start, double t_timeout, int synchronised);
void mui_announce_recv_span_2f_sphere(mui_uniface_2f *uniface, mui_point_2f centre, float radius, float t_start,
		float t_timeout, int synchronised);
void mui_announce_recv_span_2fx_sphere(mui_uniface_2fx *uniface, mui_point_2fx centre, float radius, float t_start,
		float t_timeout, int synchronised);
void mui_announce_recv_span_2d_sphere(mui_uniface_2d *uniface, mui_point_2d centre, double radius, double t_start,
		double t_timeout, int synchronised);
void mui_announce_recv_span_2dx_sphere(mui_uniface_2dx *uniface, mui_point_2dx centre, double radius, double t_start,
		double t_timeout, int synchronised);
void mui_announce_recv_span_2t_sphere(mui_uniface_2t *uniface, mui_point_2t centre, double radius, double t_start,
		double t_timeout, int synchronised);
void mui_announce_send_disable_2f(mui_uniface_2f *uniface, int synchronised);
void mui_announce_send_disable_2fx(mui_uniface_2fx *uniface, int synchronised);
void mui_announce_send_disable_2d(mui_uniface_2d *uniface, int synchronised);
void mui_announce_send_disable_2dx(mui_uniface_2dx *uniface, int synchronised);
void mui_announce_send_disable_2t(mui_uniface_2t *uniface, int synchronised);
void mui_announce_recv_disable_2f(mui_uniface_2f *uniface, int synchronised);
void mui_announce_recv_disable_2fx(mui_uniface_2fx *uniface, int synchronised);
void mui_announce_recv_disable_2d(mui_uniface_2d *uniface, int synchronised);
void mui_announce_recv_disable_2dx(mui_uniface_2dx *uniface, int synchronised);
void mui_announce_recv_disable_2t(mui_uniface_2t *uniface, int synchronised);

// MUI barrier functions
void mui_barrier_2f(mui_uniface_2f *uniface, float t);
void mui_barrier_2fx(mui_uniface_2fx *uniface, float t);
void mui_barrier_2d(mui_uniface_2d *uniface, double t);
void mui_barrier_2dx(mui_uniface_2dx *uniface, double t);
void mui_barrier_2t(mui_uniface_2t *uniface, double t);
void mui_barrier_2f_pair(mui_uniface_2f *uniface, float t_1, float t_2);
void mui_barrier_2fx_pair(mui_uniface_2fx *uniface, float t_1, float t_2);
void mui_barrier_2d_pair(mui_uniface_2d *uniface, double t_1, double t_2);
void mui_barrier_2dx_pair(mui_uniface_2dx *uniface, double t_1, double t_2);
void mui_barrier_2t_pair(mui_uniface_2t *uniface, double t_1, double t_2);

// MUI forget functions
void mui_forget_upper_2f(mui_uniface_2f *uniface, float upper, int reset_log);
void mui_forget_upper_2fx(mui_uniface_2fx *uniface, float upper, int reset_log);
void mui_forget_upper_2d(mui_uniface_2d *uniface, double upper, int reset_log);
void mui_forget_upper_2dx(mui_uniface_2dx *uniface, double upper, int reset_log);
void mui_forget_upper_2t(mui_uniface_2t *uniface, double upper, int reset_log);
void mui_forget_upper_2f_pair(mui_uniface_2f *uniface, float upper_1, float upper_2, int reset_log);
void mui_forget_upper_2fx_pair(mui_uniface_2fx *uniface, float upper_1, float upper_2, int reset_log);
void mui_forget_upper_2d_pair(mui_uniface_2d *uniface, double upper_1, double upper_2, int reset_log);
void mui_forget_upper_2dx_pair(mui_uniface_2dx *uniface, double upper_1, double upper_2, int reset_log);
void mui_forget_upper_2t_pair(mui_uniface_2t *uniface, double upper_1, double upper_2, int reset_log);
void mui_forget_lower_upper_2f(mui_uniface_2f *uniface, float lower, float upper, int reset_log);
void mui_forget_lower_upper_2fx(mui_uniface_2fx *uniface, float lower, float upper, int reset_log);
void mui_forget_lower_upper_2d(mui_uniface_2d *uniface, double lower, double upper, int reset_log);
void mui_forget_lower_upper_2dx(mui_uniface_2dx *uniface, double lower, double upper, int reset_log);
void mui_forget_lower_upper_2t(mui_uniface_2t *uniface, double lower, double upper, int reset_log);
void mui_forget_lower_upper_2f_pair(mui_uniface_2f *uniface, float lower_1, float lower_2, float upper_1, float upper_2,
		int reset_log);
void mui_forget_lower_upper_2fx_pair(mui_uniface_2fx *uniface, float lower_1, float lower_2, float upper_1,
		float upper_2, int reset_log);
void mui_forget_lower_upper_2d_pair(mui_uniface_2d *uniface, double lower_1, double lower_2, double upper_1,
		double upper_2, int reset_log);
void mui_forget_lower_upper_2dx_pair(mui_uniface_2dx *uniface, double lower_1, double lower_2, double upper_1,
		double upper_2, int reset_log);
void mui_forget_lower_upper_2t_pair(mui_uniface_2t *uniface, double lower_1, double lower_2, double upper_1,
		double upper_2, int reset_log);
void mui_set_forget_length_2f(mui_uniface_2f *uniface, float length);
void mui_set_forget_length_2fx(mui_uniface_2fx *uniface, float length);
void mui_set_forget_length_2d(mui_uniface_2d *uniface, double length);
void mui_set_forget_length_2dx(mui_uniface_2dx *uniface, double length);
void mui_set_forget_length_2t(mui_uniface_2t *uniface, double length);
const char* mui_uri_host_2f(mui_uniface_2f *uniface);
const char* mui_uri_host_2fx(mui_uniface_2fx *uniface);
const char* mui_uri_host_2d(mui_uniface_2d *uniface);
const char* mui_uri_host_2dx(mui_uniface_2dx *uniface);
const char* mui_uri_host_2t(mui_uniface_2t *uniface);
const char* mui_uri_path_2f(mui_uniface_2f *uniface);
const char* mui_uri_path_2fx(mui_uniface_2fx *uniface);
const char* mui_uri_path_2d(mui_uniface_2d *uniface);
const char* mui_uri_path_2dx(mui_uniface_2dx *uniface);
const char* mui_uri_path_2t(mui_uniface_2t *uniface);
const char* mui_uri_protocol_2f(mui_uniface_2f *uniface);
const char* mui_uri_protocol_2fx(mui_uniface_2fx *uniface);
const char* mui_uri_protocol_2d(mui_uniface_2d *uniface);
const char* mui_uri_protocol_2dx(mui_uniface_2dx *uniface);
const char* mui_uri_protocol_2t(mui_uniface_2t *uniface);

#endif /* MUI_C_WRAPPER_2D_H_ */
