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
 * @file mui_c_wrapper_3d.h
 * @author S. M. Longshaw (derived from original 3D wrapper by Y. H. Tang)
 * @date Aug 4, 2021
 * @brief Header for C wrapper to create and manage 3D MUI interfaces and
 *        associated sampler objects
 *
 *        NOTE: Any point co-ordinates are enumerated rather than assuming
 *              Cartesian form, i.e. {1, 2, 3} rather than {x, y, z}.
 */

#ifndef MUI_C_WRAPPER_3D_H_
#define MUI_C_WRAPPER_3D_H_

// C-defined simple struct versions of MUI point types
typedef struct mui_point_3f {
	float point_1;
	float point_2;
	float point_3;
} mui_point_3f;

typedef struct mui_point_3fx {
	float point_1;
	float point_2;
	float point_3;
} mui_point_3fx;

typedef struct mui_point_3d {
	double point_1;
	double point_2;
	double point_3;
} mui_point_3d;

typedef struct mui_point_3dx {
	double point_1;
	double point_2;
	double point_3;
} mui_point_3dx;

typedef struct mui_point_3t {
	double point_1;
	double point_2;
	double point_3;
} mui_point_3t;

// C access typedefs for uniface and sampler types
typedef struct mui_uniface_3f mui_uniface_3f;
typedef struct mui_uniface_3fx mui_uniface_3fx;
typedef struct mui_uniface_3d mui_uniface_3d;
typedef struct mui_uniface_3dx mui_uniface_3dx;
typedef struct mui_uniface_3t mui_uniface_3t;

typedef struct mui_sampler_exact_3f mui_sampler_exact_3f;
typedef struct mui_sampler_exact_3fx mui_sampler_exact_3fx;
typedef struct mui_sampler_exact_3d mui_sampler_exact_3d;
typedef struct mui_sampler_exact_3dx mui_sampler_exact_3dx;
typedef struct mui_sampler_exact_3t mui_sampler_exact_3t;

typedef struct mui_sampler_gauss_3f mui_sampler_gauss_3f;
typedef struct mui_sampler_gauss_3fx mui_sampler_gauss_3fx;
typedef struct mui_sampler_gauss_3d mui_sampler_gauss_3d;
typedef struct mui_sampler_gauss_3dx mui_sampler_gauss_3dx;
typedef struct mui_sampler_gauss_3t mui_sampler_gauss_3t;

typedef struct mui_sampler_moving_average_3f mui_sampler_moving_average_3f;
typedef struct mui_sampler_moving_average_3fx mui_sampler_moving_average_3fx;
typedef struct mui_sampler_moving_average_3d mui_sampler_moving_average_3d;
typedef struct mui_sampler_moving_average_3dx mui_sampler_moving_average_3dx;
typedef struct mui_sampler_moving_average_3t mui_sampler_moving_average_3t;

typedef struct mui_sampler_nearest_neighbor_3f mui_sampler_nearest_neighbor_3f;
typedef struct mui_sampler_nearest_neighbor_3fx mui_sampler_nearest_neighbor_3fx;
typedef struct mui_sampler_nearest_neighbor_3d mui_sampler_nearest_neighbor_3d;
typedef struct mui_sampler_nearest_neighbor_3dx mui_sampler_nearest_neighbor_3dx;
typedef struct mui_sampler_nearest_neighbor_3t mui_sampler_nearest_neighbor_3t;

typedef struct mui_sampler_pseudo_n2_linear_3f mui_sampler_pseudo_n2_linear_3f;
typedef struct mui_sampler_pseudo_n2_linear_3fx mui_sampler_pseudo_n2_linear_3fx;
typedef struct mui_sampler_pseudo_n2_linear_3d mui_sampler_pseudo_n2_linear_3d;
typedef struct mui_sampler_pseudo_n2_linear_3dx mui_sampler_pseudo_n2_linear_3dx;
typedef struct mui_sampler_pseudo_n2_linear_3t mui_sampler_pseudo_n2_linear_3t;

typedef struct mui_sampler_pseudo_nearest_neighbor_3f mui_sampler_pseudo_nearest_neighbor_3f;
typedef struct mui_sampler_pseudo_nearest_neighbor_3fx mui_sampler_pseudo_nearest_neighbor_3fx;
typedef struct mui_sampler_pseudo_nearest_neighbor_3d mui_sampler_pseudo_nearest_neighbor_3d;
typedef struct mui_sampler_pseudo_nearest_neighbor_3dx mui_sampler_pseudo_nearest_neighbor_3dx;
typedef struct mui_sampler_pseudo_nearest_neighbor_3t mui_sampler_pseudo_nearest_neighbor_3t;

typedef struct mui_sampler_shepard_quintic_3f mui_sampler_shepard_quintic_3f;
typedef struct mui_sampler_shepard_quintic_3fx mui_sampler_shepard_quintic_3fx;
typedef struct mui_sampler_shepard_quintic_3d mui_sampler_shepard_quintic_3d;
typedef struct mui_sampler_shepard_quintic_3dx mui_sampler_shepard_quintic_3dx;
typedef struct mui_sampler_shepard_quintic_3t mui_sampler_shepard_quintic_3t;

typedef struct mui_sampler_sph_quintic_3f mui_sampler_sph_quintic_3f;
typedef struct mui_sampler_sph_quintic_3fx mui_sampler_sph_quintic_3fx;
typedef struct mui_sampler_sph_quintic_3d mui_sampler_sph_quintic_3d;
typedef struct mui_sampler_sph_quintic_3dx mui_sampler_sph_quintic_3dx;
typedef struct mui_sampler_sph_quintic_3t mui_sampler_sph_quintic_3t;

typedef struct mui_sampler_sum_quintic_3f mui_sampler_sum_quintic_3f;
typedef struct mui_sampler_sum_quintic_3fx mui_sampler_sum_quintic_3fx;
typedef struct mui_sampler_sum_quintic_3d mui_sampler_sum_quintic_3d;
typedef struct mui_sampler_sum_quintic_3dx mui_sampler_sum_quintic_3dx;
typedef struct mui_sampler_sum_quintic_3t mui_sampler_sum_quintic_3t;

#ifdef USE_RBF
typedef struct mui_sampler_rbf_3f mui_sampler_rbf_3f;
typedef struct mui_sampler_rbf_3fx mui_sampler_rbf_3fx;
typedef struct mui_sampler_rbf_3d mui_sampler_rbf_3d;
typedef struct mui_sampler_rbf_3dx mui_sampler_rbf_3dx;
typedef struct mui_sampler_rbf_3t mui_sampler_rbf_3t;
#endif

typedef struct mui_chrono_sampler_exact_3f mui_chrono_sampler_exact_3f;
typedef struct mui_chrono_sampler_exact_3fx mui_chrono_sampler_exact_3fx;
typedef struct mui_chrono_sampler_exact_3d mui_chrono_sampler_exact_3d;
typedef struct mui_chrono_sampler_exact_3dx mui_chrono_sampler_exact_3dx;
typedef struct mui_chrono_sampler_exact_3t mui_chrono_sampler_exact_3t;

typedef struct mui_chrono_sampler_gauss_3f mui_chrono_sampler_gauss_3f;
typedef struct mui_chrono_sampler_gauss_3fx mui_chrono_sampler_gauss_3fx;
typedef struct mui_chrono_sampler_gauss_3d mui_chrono_sampler_gauss_3d;
typedef struct mui_chrono_sampler_gauss_3dx mui_chrono_sampler_gauss_3dx;
typedef struct mui_chrono_sampler_gauss_3t mui_chrono_sampler_gauss_3t;

typedef struct mui_chrono_sampler_mean_3f mui_chrono_sampler_mean_3f;
typedef struct mui_chrono_sampler_mean_3fx mui_chrono_sampler_mean_3fx;
typedef struct mui_chrono_sampler_mean_3d mui_chrono_sampler_mean_3d;
typedef struct mui_chrono_sampler_mean_3dx mui_chrono_sampler_mean_3dx;
typedef struct mui_chrono_sampler_mean_3t mui_chrono_sampler_mean_3t;

typedef struct mui_chrono_sampler_sum_3f mui_chrono_sampler_sum_3f;
typedef struct mui_chrono_sampler_sum_3fx mui_chrono_sampler_sum_3fx;
typedef struct mui_chrono_sampler_sum_3d mui_chrono_sampler_sum_3d;
typedef struct mui_chrono_sampler_sum_3dx mui_chrono_sampler_sum_3dx;
typedef struct mui_chrono_sampler_sum_3t mui_chrono_sampler_sum_3t;

// MUI uniface creation
mui_uniface_3f* mui_create_uniface_3f(const char *URI);
mui_uniface_3fx* mui_create_uniface_3fx(const char *URI);
mui_uniface_3d* mui_create_uniface_3d(const char *URI);
mui_uniface_3dx* mui_create_uniface_3dx(const char *URI);
mui_uniface_3t* mui_create_uniface_3t(const char *URI);

// MUI multi uniface creation
mui_uniface_3f** mui_create_uniface_multi_3f( const char *domain, const char **interfaces, int interface_count );
mui_uniface_3fx** mui_create_uniface_multi_3fx( const char *domain, const char **interfaces, int interface_count );
mui_uniface_3d** mui_create_uniface_multi_3d( const char *domain, const char **interfaces, int interface_count );
mui_uniface_3dx** mui_create_uniface_multi_3dx( const char *domain, const char **interfaces, int interface_count );
mui_uniface_3t** mui_create_uniface_multi_3t( const char *domain, const char **interfaces, int interface_count );

// MUI uniface destruction
void mui_destroy_uniface_3f(mui_uniface_3f *uniface);
void mui_destroy_uniface_3fx(mui_uniface_3fx *uniface);
void mui_destroy_uniface_3d(mui_uniface_3d *uniface);
void mui_destroy_uniface_3dx(mui_uniface_3dx *uniface);
void mui_destroy_uniface_3t(mui_uniface_3t *uniface);

// MUI spatial samplers creation
mui_sampler_exact_3f* mui_create_sampler_exact_3f(float tolerance);
mui_sampler_exact_3fx* mui_create_sampler_exact_3fx(float tolerance);
mui_sampler_exact_3d* mui_create_sampler_exact_3d(double tolerance);
mui_sampler_exact_3dx* mui_create_sampler_exact_3dx(double tolerance);
mui_sampler_exact_3t* mui_create_sampler_exact_3t(double tolerance);
mui_sampler_gauss_3f* mui_create_sampler_gauss_3f(float r, float h);
mui_sampler_gauss_3fx* mui_create_sampler_gauss_3fx(float r, float h);
mui_sampler_gauss_3d* mui_create_sampler_gauss_3d(double r, double h);
mui_sampler_gauss_3dx* mui_create_sampler_gauss_3dx(double r, double h);
mui_sampler_gauss_3t* mui_create_sampler_gauss_3t(double r, double h);
mui_sampler_moving_average_3f* mui_create_sampler_moving_average_3f(float bbox_1, float bbox_2, float bbox_3);
mui_sampler_moving_average_3fx* mui_create_sampler_moving_average_3fx(float bbox_1, float bbox_2, float bbox_3);
mui_sampler_moving_average_3d* mui_create_sampler_moving_average_3d(double bbox_1, double bbox_2, double bbox_3);
mui_sampler_moving_average_3dx* mui_create_sampler_moving_average_3dx(double bbox_1, double bbox_2, double bbox_3);
mui_sampler_moving_average_3t* mui_create_sampler_moving_average_3t(double bbox_1, double bbox_2, double bbox_3);
mui_sampler_nearest_neighbor_3f* mui_create_sampler_nearest_neighbor_3f();
mui_sampler_nearest_neighbor_3fx* mui_create_sampler_nearest_neighbor_3fx();
mui_sampler_nearest_neighbor_3d* mui_create_sampler_nearest_neighbor_3d();
mui_sampler_nearest_neighbor_3dx* mui_create_sampler_nearest_neighbor_3dx();
mui_sampler_nearest_neighbor_3t* mui_create_sampler_nearest_neighbor_3t();
mui_sampler_pseudo_n2_linear_3f* mui_create_sampler_pseudo_n2_linear_3f(float r);
mui_sampler_pseudo_n2_linear_3fx* mui_create_sampler_pseudo_n2_linear_3fx(float r);
mui_sampler_pseudo_n2_linear_3d* mui_create_sampler_pseudo_n2_linear_3d(double r);
mui_sampler_pseudo_n2_linear_3dx* mui_create_sampler_pseudo_n2_linear_3dx(double r);
mui_sampler_pseudo_n2_linear_3t* mui_create_sampler_pseudo_n2_linear_3t(double r);
mui_sampler_pseudo_nearest_neighbor_3f* mui_create_sampler_pseudo_nearest_neighbor_3f(float h);
mui_sampler_pseudo_nearest_neighbor_3fx* mui_create_sampler_pseudo_nearest_neighbor_3fx(float h);
mui_sampler_pseudo_nearest_neighbor_3d* mui_create_sampler_pseudo_nearest_neighbor_3d(double h);
mui_sampler_pseudo_nearest_neighbor_3dx* mui_create_sampler_pseudo_nearest_neighbor_3dx(double h);
mui_sampler_pseudo_nearest_neighbor_3t* mui_create_sampler_pseudo_nearest_neighbor_3t(double h);
mui_sampler_shepard_quintic_3f* mui_create_sampler_shepard_quintic_3f(float r);
mui_sampler_shepard_quintic_3fx* mui_create_sampler_shepard_quintic_3fx(float r);
mui_sampler_shepard_quintic_3d* mui_create_sampler_shepard_quintic_3d(double r);
mui_sampler_shepard_quintic_3dx* mui_create_sampler_shepard_quintic_3dx(double r);
mui_sampler_shepard_quintic_3t* mui_create_sampler_shepard_quintic_3t(double r);
mui_sampler_sph_quintic_3f* mui_create_sampler_sph_quintic_3f(float r);
mui_sampler_sph_quintic_3fx* mui_create_sampler_sph_quintic_3fx(float r);
mui_sampler_sph_quintic_3d* mui_create_sampler_sph_quintic_3d(double r);
mui_sampler_sph_quintic_3dx* mui_create_sampler_sph_quintic_3dx(double r);
mui_sampler_sph_quintic_3t* mui_create_sampler_sph_quintic_3t(double r);
mui_sampler_sum_quintic_3f* mui_create_sampler_sum_quintic_3f(float r);
mui_sampler_sum_quintic_3fx* mui_create_sampler_sum_quintic_3fx(float r);
mui_sampler_sum_quintic_3d* mui_create_sampler_sum_quintic_3d(double r);
mui_sampler_sum_quintic_3dx* mui_create_sampler_sum_quintic_3dx(double r);
mui_sampler_sum_quintic_3t* mui_create_sampler_sum_quintic_3t(double r);
#ifdef USE_RBF
mui_sampler_rbf_3f* mui_create_sampler_rbf_3f(float r, mui_point_3f *points, int points_count, int basis_func,
		int conservative, int polynomial, int smoothFunc, int readMatrix, int writeMatrix, const char *file_address,
		float cutoff, float cg_solve_tol, int cg_solve_it, int pou_size);
mui_sampler_rbf_3fx* mui_create_sampler_rbf_3fx(float r, mui_point_3fx *points, int points_count, int basis_func,
		int conservative, int polynomial, int smoothFunc, int readMatrix, int writeMatrix, const char *file_address,
		float cutoff, float cg_solve_tol, int cg_solve_it, int pou_size);
mui_sampler_rbf_3d* mui_create_sampler_rbf_3d(double r, mui_point_3d *points, int points_count, int basis_func,
		int conservative, int polynomial, int smoothFunc, int readMatrix, int writeMatrix, const char *file_address,
		double cutoff, double cg_solve_tol, int cg_solve_it, int pou_size);
mui_sampler_rbf_3dx* mui_create_sampler_rbf_3dx(double r, mui_point_3dx *points, int points_count, int basis_func,
		int conservative, int polynomial, int smoothFunc, int readMatrix, int writeMatrix, const char *file_address,
		double cutoff, double cg_solve_tol, int cg_solve_it, int pou_size);
mui_sampler_rbf_3t* mui_create_sampler_rbf_3t(double r, mui_point_3t *points, int points_count, int basis_func,
		int conservative, int polynomial, int smoothFunc, int readMatrix, int writeMatrix, const char *file_address,
		double cutoff, double cg_solve_tol, int cg_solve_it, int pou_size);
#endif

// MUI spatial samplers destruction
void mui_destroy_sampler_exact_3f(mui_sampler_exact_3f *sampler);
void mui_destroy_sampler_exact_3fx(mui_sampler_exact_3fx *sampler);
void mui_destroy_sampler_exact_3d(mui_sampler_exact_3d *sampler);
void mui_destroy_sampler_exact_3dx(mui_sampler_exact_3dx *sampler);
void mui_destroy_sampler_exact_3t(mui_sampler_exact_3t *sampler);
void mui_destroy_sampler_gauss_3f(mui_sampler_gauss_3f *sampler);
void mui_destroy_sampler_gauss_3fx(mui_sampler_gauss_3fx *sampler);
void mui_destroy_sampler_gauss_3d(mui_sampler_gauss_3d *sampler);
void mui_destroy_sampler_gauss_3dx(mui_sampler_gauss_3dx *sampler);
void mui_destroy_sampler_gauss_3t(mui_sampler_gauss_3t *sampler);
void mui_destroy_sampler_moving_average_3f(mui_sampler_moving_average_3f *sampler);
void mui_destroy_sampler_moving_average_3fx(mui_sampler_moving_average_3fx *sampler);
void mui_destroy_sampler_moving_average_3d(mui_sampler_moving_average_3d *sampler);
void mui_destroy_sampler_moving_average_3dx(mui_sampler_moving_average_3dx *sampler);
void mui_destroy_sampler_moving_average_3t(mui_sampler_moving_average_3t *sampler);
void mui_destroy_sampler_nearest_neighbor_3f(mui_sampler_nearest_neighbor_3f *sampler);
void mui_destroy_sampler_nearest_neighbor_3fx(mui_sampler_nearest_neighbor_3fx *sampler);
void mui_destroy_sampler_nearest_neighbor_3d(mui_sampler_nearest_neighbor_3d *sampler);
void mui_destroy_sampler_nearest_neighbor_3dx(mui_sampler_nearest_neighbor_3dx *sampler);
void mui_destroy_sampler_nearest_neighbor_3t(mui_sampler_nearest_neighbor_3t *sampler);
void mui_destroy_sampler_pseudo_nearest2_linear_3f(mui_sampler_pseudo_nearest_neighbor_3f *sampler);
void mui_destroy_sampler_pseudo_nearest2_linear_3fx(mui_sampler_pseudo_nearest_neighbor_3fx *sampler);
void mui_destroy_sampler_pseudo_nearest2_linear_3d(mui_sampler_pseudo_nearest_neighbor_3d *sampler);
void mui_destroy_sampler_pseudo_nearest2_linear_3dx(mui_sampler_pseudo_nearest_neighbor_3dx *sampler);
void mui_destroy_sampler_pseudo_nearest2_linear_3t(mui_sampler_pseudo_nearest_neighbor_3t *sampler);
void mui_destroy_sampler_pseudo_nearest_neighbor_3f(mui_sampler_pseudo_nearest_neighbor_3f *sampler);
void mui_destroy_sampler_pseudo_nearest_neighbor_3fx(mui_sampler_pseudo_nearest_neighbor_3fx *sampler);
void mui_destroy_sampler_pseudo_nearest_neighbor_3d(mui_sampler_pseudo_nearest_neighbor_3d *sampler);
void mui_destroy_sampler_pseudo_nearest_neighbor_3dx(mui_sampler_pseudo_nearest_neighbor_3dx *sampler);
void mui_destroy_sampler_pseudo_nearest_neighbor_3t(mui_sampler_pseudo_nearest_neighbor_3t *sampler);
void mui_destroy_sampler_shepard_quintic_3f(mui_sampler_shepard_quintic_3f *sampler);
void mui_destroy_sampler_shepard_quintic_3fx(mui_sampler_shepard_quintic_3fx *sampler);
void mui_destroy_sampler_shepard_quintic_3d(mui_sampler_shepard_quintic_3d *sampler);
void mui_destroy_sampler_shepard_quintic_3dx(mui_sampler_shepard_quintic_3dx *sampler);
void mui_destroy_sampler_shepard_quintic_3t(mui_sampler_shepard_quintic_3t *sampler);
void mui_destroy_sampler_sph_quintic_3f(mui_sampler_sph_quintic_3f *sampler);
void mui_destroy_sampler_sph_quintic_3fx(mui_sampler_sph_quintic_3fx *sampler);
void mui_destroy_sampler_sph_quintic_3d(mui_sampler_sph_quintic_3d *sampler);
void mui_destroy_sampler_sph_quintic_3dx(mui_sampler_sph_quintic_3dx *sampler);
void mui_destroy_sampler_sph_quintic_3t(mui_sampler_sph_quintic_3t *sampler);
void mui_destroy_sampler_sum_quintic_3f(mui_sampler_sum_quintic_3f *sampler);
void mui_destroy_sampler_sum_quintic_3fx(mui_sampler_sum_quintic_3fx *sampler);
void mui_destroy_sampler_sum_quintic_3d(mui_sampler_sum_quintic_3d *sampler);
void mui_destroy_sampler_sum_quintic_3dx(mui_sampler_sum_quintic_3dx *sampler);
void mui_destroy_sampler_sum_quintic_3t(mui_sampler_sum_quintic_3t *sampler);
#ifdef USE_RBF
void mui_destroy_sampler_rbf_3f(mui_sampler_rbf_3f *sampler);
void mui_destroy_sampler_rbf_3fx(mui_sampler_rbf_3fx *sampler);
void mui_destroy_sampler_rbf_3d(mui_sampler_rbf_3d *sampler);
void mui_destroy_sampler_rbf_3dx(mui_sampler_rbf_3dx *sampler);
void mui_destroy_sampler_rbf_3t(mui_sampler_rbf_3t *sampler);
#endif

// MUI temporal samplers creation
mui_chrono_sampler_exact_3f* mui_create_chrono_sampler_exact_3f(float tolerance);
mui_chrono_sampler_exact_3fx* mui_create_chrono_sampler_exact_3fx(float tolerance);
mui_chrono_sampler_exact_3d* mui_create_chrono_sampler_exact_3d(double tolerance);
mui_chrono_sampler_exact_3dx* mui_create_chrono_sampler_exact_3dx(double tolerance);
mui_chrono_sampler_exact_3t* mui_create_chrono_sampler_exact_3t(double tolerance);
mui_chrono_sampler_gauss_3f* mui_create_chrono_sampler_gauss_3f(float cutoff, float sigma);
mui_chrono_sampler_gauss_3fx* mui_create_chrono_sampler_gauss_3fx(float cutoff, float sigma);
mui_chrono_sampler_gauss_3d* mui_create_chrono_sampler_gauss_3d(double cutoff, double sigma);
mui_chrono_sampler_gauss_3dx* mui_create_chrono_sampler_gauss_3dx(double cutoff, double sigma);
mui_chrono_sampler_gauss_3t* mui_create_chrono_sampler_gauss_3t(double cutoff, double sigma);
mui_chrono_sampler_mean_3f* mui_create_chrono_sampler_mean_3f(float lower, float upper);
mui_chrono_sampler_mean_3fx* mui_create_chrono_sampler_mean_3fx(float lower, float upper);
mui_chrono_sampler_mean_3d* mui_create_chrono_sampler_mean_3d(double lower, double upper);
mui_chrono_sampler_mean_3dx* mui_create_chrono_sampler_mean_3dx(double lower, double upper);
mui_chrono_sampler_mean_3t* mui_create_chrono_sampler_mean_3t(double lower, double upper);
mui_chrono_sampler_sum_3f* mui_create_chrono_sampler_sum_3f(float lower, float upper);
mui_chrono_sampler_sum_3fx* mui_create_chrono_sampler_sum_3fx(float lower, float upper);
mui_chrono_sampler_sum_3d* mui_create_chrono_sampler_sum_3d(double lower, double upper);
mui_chrono_sampler_sum_3dx* mui_create_chrono_sampler_sum_3dx(double lower, double upper);
mui_chrono_sampler_sum_3t* mui_create_chrono_sampler_sum_3t(double lower, double upper);

// MUI temporal samplers destruction
void mui_destroy_chrono_sampler_exact_3f(mui_chrono_sampler_exact_3f *sampler);
void mui_destroy_chrono_sampler_exact_3fx(mui_chrono_sampler_exact_3fx *sampler);
void mui_destroy_chrono_sampler_exact_3d(mui_chrono_sampler_exact_3d *sampler);
void mui_destroy_chrono_sampler_exact_3dx(mui_chrono_sampler_exact_3dx *sampler);
void mui_destroy_chrono_sampler_exact_3t(mui_chrono_sampler_exact_3t *sampler);
void mui_destroy_chrono_sampler_gauss_3f(mui_chrono_sampler_gauss_3f *sampler);
void mui_destroy_chrono_sampler_gauss_3fx(mui_chrono_sampler_gauss_3fx *sampler);
void mui_destroy_chrono_sampler_gauss_3d(mui_chrono_sampler_gauss_3d *sampler);
void mui_destroy_chrono_sampler_gauss_3dx(mui_chrono_sampler_gauss_3dx *sampler);
void mui_destroy_chrono_sampler_gauss_3t(mui_chrono_sampler_gauss_3t *sampler);
void mui_destroy_chrono_sampler_mean_3f(mui_chrono_sampler_mean_3f *sampler);
void mui_destroy_chrono_sampler_mean_3fx(mui_chrono_sampler_mean_3fx *sampler);
void mui_destroy_chrono_sampler_mean_3d(mui_chrono_sampler_mean_3d *sampler);
void mui_destroy_chrono_sampler_mean_3dx(mui_chrono_sampler_mean_3dx *sampler);
void mui_destroy_chrono_sampler_mean_3t(mui_chrono_sampler_mean_3t *sampler);
void mui_destroy_chrono_sampler_sum_3f(mui_chrono_sampler_sum_3f *sampler);
void mui_destroy_chrono_sampler_sum_3fx(mui_chrono_sampler_sum_3fx *sampler);
void mui_destroy_chrono_sampler_sum_3d(mui_chrono_sampler_sum_3d *sampler);
void mui_destroy_chrono_sampler_sum_3dx(mui_chrono_sampler_sum_3dx *sampler);
void mui_destroy_chrono_sampler_sum_3t(mui_chrono_sampler_sum_3t *sampler);

// MUI functions for data push
void mui_push_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float value);
void mui_push_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float value);
void mui_push_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double value);
void mui_push_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double value);
void mui_push_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double value);
void mui_push_3f_param(mui_uniface_3f *uniface, const char *attr, float value);
void mui_push_3fx_param(mui_uniface_3fx *uniface, const char *attr, float value);
void mui_push_3d_param(mui_uniface_3d *uniface, const char *attr, double value);
void mui_push_3dx_param(mui_uniface_3dx *uniface, const char *attr, double value);
void mui_push_3t_param(mui_uniface_3t *uniface, const char *attr, double value);

// MUI functions for data commit
void mui_commit_3f(mui_uniface_3f *uniface, float t);
void mui_commit_3fx(mui_uniface_3fx *uniface, float t);
void mui_commit_3d(mui_uniface_3d *uniface, double t);
void mui_commit_3dx(mui_uniface_3dx *uniface, double t);
void mui_commit_3t(mui_uniface_3t *uniface, double t);
void mui_commit_3f_pair(mui_uniface_3f *uniface, float t_1, float t_2);
void mui_commit_3fx_pair(mui_uniface_3fx *uniface, float t_1, float t_2);
void mui_commit_3d_pair(mui_uniface_3d *uniface, double t_1, double t_2);
void mui_commit_3dx_pair(mui_uniface_3dx *uniface, double t_1, double t_2);
void mui_commit_3t_pair(mui_uniface_3t *uniface, double t_1, double t_2);

// MUI functions for data forecast
void mui_forecast_3f(mui_uniface_3f *uniface, float t);
void mui_forecast_3fx(mui_uniface_3fx *uniface, float t);
void mui_forecast_3d(mui_uniface_3d *uniface, double t);
void mui_forecast_3dx(mui_uniface_3dx *uniface, double t);
void mui_forecast_3t(mui_uniface_3t *uniface, double t);
void mui_forecast_3f_pair(mui_uniface_3f *uniface, float t_1, float t_2);
void mui_forecast_3fx_pair(mui_uniface_3fx *uniface, float t_1, float t_2);
void mui_forecast_3d_pair(mui_uniface_3d *uniface, double t_1, double t_2);
void mui_forecast_3dx_pair(mui_uniface_3dx *uniface, double t_1, double t_2);
void mui_forecast_3t_pair(mui_uniface_3t *uniface, double t_1, double t_2);

// MUI functions for data fetch
float mui_fetch_exact_exact_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_exact_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler);
float mui_fetch_exact_exact_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t,
		mui_sampler_exact_3fx *spatial_sampler, mui_chrono_sampler_exact_3fx *temporal_sampler);
double mui_fetch_exact_exact_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_exact_3d *spatial_sampler, mui_chrono_sampler_exact_3d *temporal_sampler);
double mui_fetch_exact_exact_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t,
		mui_sampler_exact_3dx *spatial_sampler, mui_chrono_sampler_exact_3dx *temporal_sampler);
double mui_fetch_exact_exact_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_exact_3t *spatial_sampler, mui_chrono_sampler_exact_3t *temporal_sampler);
float mui_fetch_exact_gauss_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_exact_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler);
float mui_fetch_exact_gauss_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t,
		mui_sampler_exact_3fx *spatial_sampler, mui_chrono_sampler_gauss_3fx *temporal_sampler);
double mui_fetch_exact_gauss_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_exact_3d *spatial_sampler, mui_chrono_sampler_gauss_3d *temporal_sampler);
double mui_fetch_exact_gauss_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t,
		mui_sampler_exact_3dx *spatial_sampler, mui_chrono_sampler_gauss_3dx *temporal_sampler);
double mui_fetch_exact_gauss_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_exact_3t *spatial_sampler, mui_chrono_sampler_gauss_3t *temporal_sampler);
float mui_fetch_exact_mean_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_exact_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler);
float mui_fetch_exact_mean_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t,
		mui_sampler_exact_3fx *spatial_sampler, mui_chrono_sampler_mean_3fx *temporal_sampler);
double mui_fetch_exact_mean_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_exact_3d *spatial_sampler, mui_chrono_sampler_mean_3d *temporal_sampler);
double mui_fetch_exact_mean_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t,
		mui_sampler_exact_3dx *spatial_sampler, mui_chrono_sampler_mean_3dx *temporal_sampler);
double mui_fetch_exact_mean_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_exact_3t *spatial_sampler, mui_chrono_sampler_mean_3t *temporal_sampler);
float mui_fetch_exact_sum_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_exact_3f *spatial_sampler, mui_chrono_sampler_sum_3f *temporal_sampler);
float mui_fetch_exact_sum_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t,
		mui_sampler_exact_3fx *spatial_sampler, mui_chrono_sampler_sum_3fx *temporal_sampler);
double mui_fetch_exact_sum_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_exact_3d *spatial_sampler, mui_chrono_sampler_sum_3d *temporal_sampler);
double mui_fetch_exact_sum_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t,
		mui_sampler_exact_3dx *spatial_sampler, mui_chrono_sampler_sum_3dx *temporal_sampler);
double mui_fetch_exact_sum_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_exact_3t *spatial_sampler, mui_chrono_sampler_sum_3t *temporal_sampler);
float mui_fetch_gauss_exact_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_gauss_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler);
float mui_fetch_gauss_exact_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t,
		mui_sampler_gauss_3fx *spatial_sampler, mui_chrono_sampler_exact_3fx *temporal_sampler);
double mui_fetch_gauss_exact_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_gauss_3d *spatial_sampler, mui_chrono_sampler_exact_3d *temporal_sampler);
double mui_fetch_gauss_exact_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t,
		mui_sampler_gauss_3dx *spatial_sampler, mui_chrono_sampler_exact_3dx *temporal_sampler);
double mui_fetch_gauss_exact_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_gauss_3t *spatial_sampler, mui_chrono_sampler_exact_3t *temporal_sampler);
float mui_fetch_gauss_gauss_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_gauss_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler);
float mui_fetch_gauss_gauss_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t,
		mui_sampler_gauss_3fx *spatial_sampler, mui_chrono_sampler_gauss_3fx *temporal_sampler);
double mui_fetch_gauss_gauss_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_gauss_3d *spatial_sampler, mui_chrono_sampler_gauss_3d *temporal_sampler);
double mui_fetch_gauss_gauss_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t,
		mui_sampler_gauss_3dx *spatial_sampler, mui_chrono_sampler_gauss_3dx *temporal_sampler);
double mui_fetch_gauss_gauss_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_gauss_3t *spatial_sampler, mui_chrono_sampler_gauss_3t *temporal_sampler);
float mui_fetch_gauss_mean_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_gauss_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler);
float mui_fetch_gauss_mean_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t,
		mui_sampler_gauss_3fx *spatial_sampler, mui_chrono_sampler_mean_3fx *temporal_sampler);
double mui_fetch_gauss_mean_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_gauss_3d *spatial_sampler, mui_chrono_sampler_mean_3d *temporal_sampler);
double mui_fetch_gauss_mean_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t,
		mui_sampler_gauss_3dx *spatial_sampler, mui_chrono_sampler_mean_3dx *temporal_sampler);
double mui_fetch_gauss_mean_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_gauss_3t *spatial_sampler, mui_chrono_sampler_mean_3t *temporal_sampler);
float mui_fetch_moving_average_exact_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_moving_average_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler);
float mui_fetch_moving_average_exact_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t,
		mui_sampler_moving_average_3fx *spatial_sampler, mui_chrono_sampler_exact_3fx *temporal_sampler);
double mui_fetch_moving_average_exact_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_moving_average_3d *spatial_sampler, mui_chrono_sampler_exact_3d *temporal_sampler);
double mui_fetch_moving_average_exact_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t,
		mui_sampler_moving_average_3dx *spatial_sampler, mui_chrono_sampler_exact_3dx *temporal_sampler);
double mui_fetch_moving_average_exact_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_moving_average_3t *spatial_sampler, mui_chrono_sampler_exact_3t *temporal_sampler);
float mui_fetch_moving_average_gauss_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_moving_average_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler);
float mui_fetch_moving_average_gauss_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t,
		mui_sampler_moving_average_3fx *spatial_sampler, mui_chrono_sampler_gauss_3fx *temporal_sampler);
double mui_fetch_moving_average_gauss_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_moving_average_3d *spatial_sampler, mui_chrono_sampler_gauss_3d *temporal_sampler);
double mui_fetch_moving_average_gauss_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t,
		mui_sampler_moving_average_3dx *spatial_sampler, mui_chrono_sampler_gauss_3dx *temporal_sampler);
double mui_fetch_moving_average_gauss_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_moving_average_3t *spatial_sampler, mui_chrono_sampler_gauss_3t *temporal_sampler);
float mui_fetch_moving_average_mean_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_moving_average_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler);
float mui_fetch_moving_average_mean_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t,
		mui_sampler_moving_average_3fx *spatial_sampler, mui_chrono_sampler_mean_3fx *temporal_sampler);
double mui_fetch_moving_average_mean_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_moving_average_3d *spatial_sampler, mui_chrono_sampler_mean_3d *temporal_sampler);
double mui_fetch_moving_average_mean_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t,
		mui_sampler_moving_average_3dx *spatial_sampler, mui_chrono_sampler_mean_3dx *temporal_sampler);
double mui_fetch_moving_average_mean_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_moving_average_3t *spatial_sampler, mui_chrono_sampler_mean_3t *temporal_sampler);
float mui_fetch_moving_average_sum_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_moving_average_3f *spatial_sampler, mui_chrono_sampler_sum_3f *temporal_sampler);
float mui_fetch_moving_average_sum_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t,
		mui_sampler_moving_average_3fx *spatial_sampler, mui_chrono_sampler_sum_3fx *temporal_sampler);
double mui_fetch_moving_average_sum_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_moving_average_3d *spatial_sampler, mui_chrono_sampler_sum_3d *temporal_sampler);
double mui_fetch_moving_average_sum_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t,
		mui_sampler_moving_average_3dx *spatial_sampler, mui_chrono_sampler_sum_3dx *temporal_sampler);
double mui_fetch_moving_average_sum_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_moving_average_3t *spatial_sampler, mui_chrono_sampler_sum_3t *temporal_sampler);
float mui_fetch_nearest_neighbor_exact_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_nearest_neighbor_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler);
float mui_fetch_nearest_neighbor_exact_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t,
		mui_sampler_nearest_neighbor_3fx *spatial_sampler, mui_chrono_sampler_exact_3fx *temporal_sampler);
double mui_fetch_nearest_neighbor_exact_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_nearest_neighbor_3d *spatial_sampler, mui_chrono_sampler_exact_3d *temporal_sampler);
double mui_fetch_nearest_neighbor_exact_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t,
		mui_sampler_nearest_neighbor_3dx *spatial_sampler, mui_chrono_sampler_exact_3dx *temporal_sampler);
double mui_fetch_nearest_neighbor_exact_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_nearest_neighbor_3t *spatial_sampler, mui_chrono_sampler_exact_3t *temporal_sampler);
float mui_fetch_nearest_neighbor_gauss_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_nearest_neighbor_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler);
float mui_fetch_nearest_neighbor_gauss_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t,
		mui_sampler_nearest_neighbor_3fx *spatial_sampler, mui_chrono_sampler_gauss_3fx *temporal_sampler);
double mui_fetch_nearest_neighbor_gauss_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_nearest_neighbor_3d *spatial_sampler, mui_chrono_sampler_gauss_3d *temporal_sampler);
double mui_fetch_nearest_neighbor_gauss_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t,
		mui_sampler_nearest_neighbor_3dx *spatial_sampler, mui_chrono_sampler_gauss_3dx *temporal_sampler);
double mui_fetch_nearest_neighbor_gauss_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_nearest_neighbor_3t *spatial_sampler, mui_chrono_sampler_gauss_3t *temporal_sampler);
float mui_fetch_nearest_neighbor_mean_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_nearest_neighbor_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler);
float mui_fetch_nearest_neighbor_mean_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t,
		mui_sampler_nearest_neighbor_3fx *spatial_sampler, mui_chrono_sampler_mean_3fx *temporal_sampler);
double mui_fetch_nearest_neighbor_mean_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_nearest_neighbor_3d *spatial_sampler, mui_chrono_sampler_mean_3d *temporal_sampler);
double mui_fetch_nearest_neighbor_mean_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t,
		mui_sampler_nearest_neighbor_3dx *spatial_sampler, mui_chrono_sampler_mean_3dx *temporal_sampler);
double mui_fetch_nearest_neighbor_mean_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_nearest_neighbor_3t *spatial_sampler, mui_chrono_sampler_mean_3t *temporal_sampler);
float mui_fetch_nearest_neighbor_sum_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_nearest_neighbor_3f *spatial_sampler, mui_chrono_sampler_sum_3f *temporal_sampler);
float mui_fetch_nearest_neighbor_sum_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t,
		mui_sampler_nearest_neighbor_3fx *spatial_sampler, mui_chrono_sampler_sum_3fx *temporal_sampler);
double mui_fetch_nearest_neighbor_sum_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_nearest_neighbor_3d *spatial_sampler, mui_chrono_sampler_sum_3d *temporal_sampler);
double mui_fetch_nearest_neighbor_sum_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t,
		mui_sampler_nearest_neighbor_3dx *spatial_sampler, mui_chrono_sampler_sum_3dx *temporal_sampler);
double mui_fetch_nearest_neighbor_sum_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_nearest_neighbor_3t *spatial_sampler, mui_chrono_sampler_sum_3t *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_exact_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_pseudo_nearest_neighbor_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_exact_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point,
		float t, mui_sampler_pseudo_nearest_neighbor_3fx *spatial_sampler,
		mui_chrono_sampler_exact_3fx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_exact_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point,
		double t, mui_sampler_pseudo_nearest_neighbor_3d *spatial_sampler,
		mui_chrono_sampler_exact_3d *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_exact_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point,
		double t, mui_sampler_pseudo_nearest_neighbor_3dx *spatial_sampler,
		mui_chrono_sampler_exact_3dx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_exact_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point,
		double t, mui_sampler_pseudo_nearest_neighbor_3t *spatial_sampler,
		mui_chrono_sampler_exact_3t *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_gauss_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_pseudo_nearest_neighbor_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_gauss_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point,
		float t, mui_sampler_pseudo_nearest_neighbor_3fx *spatial_sampler,
		mui_chrono_sampler_gauss_3fx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_gauss_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point,
		double t, mui_sampler_pseudo_nearest_neighbor_3d *spatial_sampler,
		mui_chrono_sampler_gauss_3d *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_gauss_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point,
		double t, mui_sampler_pseudo_nearest_neighbor_3dx *spatial_sampler,
		mui_chrono_sampler_gauss_3dx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_gauss_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point,
		double t, mui_sampler_pseudo_nearest_neighbor_3t *spatial_sampler,
		mui_chrono_sampler_gauss_3t *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_mean_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_pseudo_nearest_neighbor_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_mean_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point,
		float t, mui_sampler_pseudo_nearest_neighbor_3fx *spatial_sampler,
		mui_chrono_sampler_mean_3fx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_mean_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point,
		double t, mui_sampler_pseudo_nearest_neighbor_3d *spatial_sampler,
		mui_chrono_sampler_mean_3d *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_mean_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point,
		double t, mui_sampler_pseudo_nearest_neighbor_3dx *spatial_sampler,
		mui_chrono_sampler_mean_3dx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_mean_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point,
		double t, mui_sampler_pseudo_nearest_neighbor_3t *spatial_sampler,
		mui_chrono_sampler_mean_3t *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_sum_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_pseudo_nearest_neighbor_3f *spatial_sampler, mui_chrono_sampler_sum_3f *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_sum_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point,
		float t, mui_sampler_pseudo_nearest_neighbor_3fx *spatial_sampler,
		mui_chrono_sampler_sum_3fx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_sum_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_pseudo_nearest_neighbor_3d *spatial_sampler, mui_chrono_sampler_sum_3d *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_sum_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point,
		double t, mui_sampler_pseudo_nearest_neighbor_3dx *spatial_sampler,
		mui_chrono_sampler_sum_3dx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_sum_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_pseudo_nearest_neighbor_3t *spatial_sampler, mui_chrono_sampler_sum_3t *temporal_sampler);
float mui_fetch_shepard_quintic_exact_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_shepard_quintic_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler);
float mui_fetch_shepard_quintic_exact_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t,
		mui_sampler_shepard_quintic_3fx *spatial_sampler, mui_chrono_sampler_exact_3fx *temporal_sampler);
double mui_fetch_shepard_quintic_exact_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_shepard_quintic_3d *spatial_sampler, mui_chrono_sampler_exact_3d *temporal_sampler);
double mui_fetch_shepard_quintic_exact_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t,
		mui_sampler_shepard_quintic_3dx *spatial_sampler, mui_chrono_sampler_exact_3dx *temporal_sampler);
double mui_fetch_shepard_quintic_exact_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_shepard_quintic_3t *spatial_sampler, mui_chrono_sampler_exact_3t *temporal_sampler);
float mui_fetch_shepard_quintic_gauss_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_shepard_quintic_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler);
float mui_fetch_shepard_quintic_gauss_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t,
		mui_sampler_shepard_quintic_3fx *spatial_sampler, mui_chrono_sampler_gauss_3fx *temporal_sampler);
double mui_fetch_shepard_quintic_gauss_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_shepard_quintic_3d *spatial_sampler, mui_chrono_sampler_gauss_3d *temporal_sampler);
double mui_fetch_shepard_quintic_gauss_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t,
		mui_sampler_shepard_quintic_3dx *spatial_sampler, mui_chrono_sampler_gauss_3dx *temporal_sampler);
double mui_fetch_shepard_quintic_gauss_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_shepard_quintic_3t *spatial_sampler, mui_chrono_sampler_gauss_3t *temporal_sampler);
float mui_fetch_shepard_quintic_mean_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_shepard_quintic_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler);
float mui_fetch_shepard_quintic_mean_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t,
		mui_sampler_shepard_quintic_3fx *spatial_sampler, mui_chrono_sampler_mean_3fx *temporal_sampler);
double mui_fetch_shepard_quintic_mean_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_shepard_quintic_3d *spatial_sampler, mui_chrono_sampler_mean_3d *temporal_sampler);
double mui_fetch_shepard_quintic_mean_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t,
		mui_sampler_shepard_quintic_3dx *spatial_sampler, mui_chrono_sampler_mean_3dx *temporal_sampler);
double mui_fetch_shepard_quintic_mean_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_shepard_quintic_3t *spatial_sampler, mui_chrono_sampler_mean_3t *temporal_sampler);
float mui_fetch_shepard_quintic_sum_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_shepard_quintic_3f *spatial_sampler, mui_chrono_sampler_sum_3f *temporal_sampler);
float mui_fetch_shepard_quintic_sum_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t,
		mui_sampler_shepard_quintic_3fx *spatial_sampler, mui_chrono_sampler_sum_3fx *temporal_sampler);
double mui_fetch_shepard_quintic_sum_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_shepard_quintic_3d *spatial_sampler, mui_chrono_sampler_sum_3d *temporal_sampler);
double mui_fetch_shepard_quintic_sum_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t,
		mui_sampler_shepard_quintic_3dx *spatial_sampler, mui_chrono_sampler_sum_3dx *temporal_sampler);
double mui_fetch_shepard_quintic_sum_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_shepard_quintic_3t *spatial_sampler, mui_chrono_sampler_sum_3t *temporal_sampler);
float mui_fetch_sph_quintic_exact_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_sph_quintic_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler);
float mui_fetch_sph_quintic_exact_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t,
		mui_sampler_sph_quintic_3fx *spatial_sampler, mui_chrono_sampler_exact_3fx *temporal_sampler);
double mui_fetch_sph_quintic_exact_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_sph_quintic_3d *spatial_sampler, mui_chrono_sampler_exact_3d *temporal_sampler);
double mui_fetch_sph_quintic_exact_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t,
		mui_sampler_sph_quintic_3dx *spatial_sampler, mui_chrono_sampler_exact_3dx *temporal_sampler);
double mui_fetch_sph_quintic_exact_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_sph_quintic_3t *spatial_sampler, mui_chrono_sampler_exact_3t *temporal_sampler);
float mui_fetch_sph_quintic_gauss_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_sph_quintic_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler);
float mui_fetch_sph_quintic_gauss_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t,
		mui_sampler_sph_quintic_3fx *spatial_sampler, mui_chrono_sampler_gauss_3fx *temporal_sampler);
double mui_fetch_sph_quintic_gauss_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_sph_quintic_3d *spatial_sampler, mui_chrono_sampler_gauss_3d *temporal_sampler);
double mui_fetch_sph_quintic_gauss_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t,
		mui_sampler_sph_quintic_3dx *spatial_sampler, mui_chrono_sampler_gauss_3dx *temporal_sampler);
double mui_fetch_sph_quintic_gauss_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_sph_quintic_3t *spatial_sampler, mui_chrono_sampler_gauss_3t *temporal_sampler);
float mui_fetch_sph_quintic_mean_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_sph_quintic_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler);
float mui_fetch_sph_quintic_mean_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t,
		mui_sampler_sph_quintic_3fx *spatial_sampler, mui_chrono_sampler_mean_3fx *temporal_sampler);
double mui_fetch_sph_quintic_mean_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_sph_quintic_3d *spatial_sampler, mui_chrono_sampler_mean_3d *temporal_sampler);
double mui_fetch_sph_quintic_mean_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t,
		mui_sampler_sph_quintic_3dx *spatial_sampler, mui_chrono_sampler_mean_3dx *temporal_sampler);
double mui_fetch_sph_quintic_mean_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_sph_quintic_3t *spatial_sampler, mui_chrono_sampler_mean_3t *temporal_sampler);
float mui_fetch_sph_quintic_sum_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_sph_quintic_3f *spatial_sampler, mui_chrono_sampler_sum_3f *temporal_sampler);
float mui_fetch_sph_quintic_sum_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t,
		mui_sampler_sph_quintic_3fx *spatial_sampler, mui_chrono_sampler_sum_3fx *temporal_sampler);
double mui_fetch_sph_quintic_sum_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_sph_quintic_3d *spatial_sampler, mui_chrono_sampler_sum_3d *temporal_sampler);
double mui_fetch_sph_quintic_sum_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t,
		mui_sampler_sph_quintic_3dx *spatial_sampler, mui_chrono_sampler_sum_3dx *temporal_sampler);
double mui_fetch_sph_quintic_sum_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_sph_quintic_3t *spatial_sampler, mui_chrono_sampler_sum_3t *temporal_sampler);
float mui_fetch_sum_quintic_exact_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_sum_quintic_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler);
float mui_fetch_sum_quintic_exact_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t,
		mui_sampler_sum_quintic_3fx *spatial_sampler, mui_chrono_sampler_exact_3fx *temporal_sampler);
double mui_fetch_sum_quintic_exact_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_sum_quintic_3d *spatial_sampler, mui_chrono_sampler_exact_3d *temporal_sampler);
double mui_fetch_sum_quintic_exact_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t,
		mui_sampler_sum_quintic_3dx *spatial_sampler, mui_chrono_sampler_exact_3dx *temporal_sampler);
double mui_fetch_sum_quintic_exact_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_sum_quintic_3t *spatial_sampler, mui_chrono_sampler_exact_3t *temporal_sampler);
float mui_fetch_sum_quintic_gauss_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_sum_quintic_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler);
float mui_fetch_sum_quintic_gauss_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t,
		mui_sampler_sum_quintic_3fx *spatial_sampler, mui_chrono_sampler_gauss_3fx *temporal_sampler);
double mui_fetch_sum_quintic_gauss_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_sum_quintic_3d *spatial_sampler, mui_chrono_sampler_gauss_3d *temporal_sampler);
double mui_fetch_sum_quintic_gauss_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t,
		mui_sampler_sum_quintic_3dx *spatial_sampler, mui_chrono_sampler_gauss_3dx *temporal_sampler);
double mui_fetch_sum_quintic_gauss_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_sum_quintic_3t *spatial_sampler, mui_chrono_sampler_gauss_3t *temporal_sampler);
float mui_fetch_sum_quintic_mean_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_sum_quintic_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler);
float mui_fetch_sum_quintic_mean_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t,
		mui_sampler_sum_quintic_3fx *spatial_sampler, mui_chrono_sampler_mean_3fx *temporal_sampler);
double mui_fetch_sum_quintic_mean_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_sum_quintic_3d *spatial_sampler, mui_chrono_sampler_mean_3d *temporal_sampler);
double mui_fetch_sum_quintic_mean_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t,
		mui_sampler_sum_quintic_3dx *spatial_sampler, mui_chrono_sampler_mean_3dx *temporal_sampler);
double mui_fetch_sum_quintic_mean_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_sum_quintic_3t *spatial_sampler, mui_chrono_sampler_mean_3t *temporal_sampler);
float mui_fetch_sum_quintic_sum_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_sum_quintic_3f *spatial_sampler, mui_chrono_sampler_sum_3f *temporal_sampler);
float mui_fetch_sum_quintic_sum_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t,
		mui_sampler_sum_quintic_3fx *spatial_sampler, mui_chrono_sampler_sum_3fx *temporal_sampler);
double mui_fetch_sum_quintic_sum_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_sum_quintic_3d *spatial_sampler, mui_chrono_sampler_sum_3d *temporal_sampler);
double mui_fetch_sum_quintic_sum_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t,
		mui_sampler_sum_quintic_3dx *spatial_sampler, mui_chrono_sampler_sum_3dx *temporal_sampler);
double mui_fetch_sum_quintic_sum_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_sum_quintic_3t *spatial_sampler, mui_chrono_sampler_sum_3t *temporal_sampler);
#ifdef USE_RBF
float mui_fetch_rbf_exact_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_rbf_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler);
float mui_fetch_rbf_exact_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t,
		mui_sampler_rbf_3fx *spatial_sampler, mui_chrono_sampler_exact_3fx *temporal_sampler);
double mui_fetch_rbf_exact_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_rbf_3d *spatial_sampler, mui_chrono_sampler_exact_3d *temporal_sampler);
double mui_fetch_rbf_exact_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t,
		mui_sampler_rbf_3dx *spatial_sampler, mui_chrono_sampler_exact_3dx *temporal_sampler);
double mui_fetch_rbf_exact_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_rbf_3t *spatial_sampler, mui_chrono_sampler_exact_3t *temporal_sampler);
float mui_fetch_rbf_gauss_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_rbf_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler);
float mui_fetch_rbf_gauss_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t,
		mui_sampler_rbf_3fx *spatial_sampler, mui_chrono_sampler_gauss_3fx *temporal_sampler);
double mui_fetch_rbf_gauss_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_rbf_3d *spatial_sampler, mui_chrono_sampler_gauss_3d *temporal_sampler);
double mui_fetch_rbf_gauss_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t,
		mui_sampler_rbf_3dx *spatial_sampler, mui_chrono_sampler_gauss_3dx *temporal_sampler);
double mui_fetch_rbf_gauss_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_rbf_3t *spatial_sampler, mui_chrono_sampler_gauss_3t *temporal_sampler);
float mui_fetch_rbf_mean_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_rbf_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler);
float mui_fetch_rbf_mean_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t,
		mui_sampler_rbf_3fx *spatial_sampler, mui_chrono_sampler_mean_3fx *temporal_sampler);
double mui_fetch_rbf_mean_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_rbf_3d *spatial_sampler, mui_chrono_sampler_mean_3d *temporal_sampler);
double mui_fetch_rbf_mean_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t,
		mui_sampler_rbf_3dx *spatial_sampler, mui_chrono_sampler_mean_3dx *temporal_sampler);
double mui_fetch_rbf_mean_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_rbf_3t *spatial_sampler, mui_chrono_sampler_mean_3t *temporal_sampler);
float mui_fetch_rbf_sum_3f(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t,
		mui_sampler_rbf_3f *spatial_sampler, mui_chrono_sampler_sum_3f *temporal_sampler);
float mui_fetch_rbf_sum_3fx(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t,
		mui_sampler_rbf_3fx *spatial_sampler, mui_chrono_sampler_sum_3fx *temporal_sampler);
double mui_fetch_rbf_sum_3d(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t,
		mui_sampler_rbf_3d *spatial_sampler, mui_chrono_sampler_sum_3d *temporal_sampler);
double mui_fetch_rbf_sum_3dx(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t,
		mui_sampler_rbf_3dx *spatial_sampler, mui_chrono_sampler_sum_3dx *temporal_sampler);
double mui_fetch_rbf_sum_3t(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t,
		mui_sampler_rbf_3t *spatial_sampler, mui_chrono_sampler_sum_3t *temporal_sampler);
#endif

/********************************************************
 * MUI functions for 3D data fetch using two time values *
 *********************************************************/

float mui_fetch_exact_exact_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t_1, float t_2,
		mui_sampler_exact_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler);
float mui_fetch_exact_exact_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t_1,
		float t_2, mui_sampler_exact_3fx *spatial_sampler, mui_chrono_sampler_exact_3fx *temporal_sampler);
double mui_fetch_exact_exact_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t_1,
		double t_2, mui_sampler_exact_3d *spatial_sampler, mui_chrono_sampler_exact_3d *temporal_sampler);
double mui_fetch_exact_exact_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t_1,
		double t_2, mui_sampler_exact_3dx *spatial_sampler, mui_chrono_sampler_exact_3dx *temporal_sampler);
double mui_fetch_exact_exact_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t_1,
		double t_2, mui_sampler_exact_3t *spatial_sampler, mui_chrono_sampler_exact_3t *temporal_sampler);
float mui_fetch_exact_gauss_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t_1, float t_2,
		mui_sampler_exact_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler);
float mui_fetch_exact_gauss_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t_1,
		float t_2, mui_sampler_exact_3fx *spatial_sampler, mui_chrono_sampler_gauss_3fx *temporal_sampler);
double mui_fetch_exact_gauss_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t_1,
		double t_2, mui_sampler_exact_3d *spatial_sampler, mui_chrono_sampler_gauss_3d *temporal_sampler);
double mui_fetch_exact_gauss_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t_1,
		double t_2, mui_sampler_exact_3dx *spatial_sampler, mui_chrono_sampler_gauss_3dx *temporal_sampler);
double mui_fetch_exact_gauss_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t_1,
		double t_2, mui_sampler_exact_3t *spatial_sampler, mui_chrono_sampler_gauss_3t *temporal_sampler);
float mui_fetch_exact_mean_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t_1, float t_2,
		mui_sampler_exact_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler);
float mui_fetch_exact_mean_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t_1,
		float t_2, mui_sampler_exact_3fx *spatial_sampler, mui_chrono_sampler_mean_3fx *temporal_sampler);
double mui_fetch_exact_mean_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t_1,
		double t_2, mui_sampler_exact_3d *spatial_sampler, mui_chrono_sampler_mean_3d *temporal_sampler);
double mui_fetch_exact_mean_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t_1,
		double t_2, mui_sampler_exact_3dx *spatial_sampler, mui_chrono_sampler_mean_3dx *temporal_sampler);
double mui_fetch_exact_mean_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t_1,
		double t_2, mui_sampler_exact_3t *spatial_sampler, mui_chrono_sampler_mean_3t *temporal_sampler);
float mui_fetch_exact_sum_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t_1, float t_2,
		mui_sampler_exact_3f *spatial_sampler, mui_chrono_sampler_sum_3f *temporal_sampler);
float mui_fetch_exact_sum_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t_1,
		float t_2, mui_sampler_exact_3fx *spatial_sampler, mui_chrono_sampler_sum_3fx *temporal_sampler);
double mui_fetch_exact_sum_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t_1,
		double t_2, mui_sampler_exact_3d *spatial_sampler, mui_chrono_sampler_sum_3d *temporal_sampler);
double mui_fetch_exact_sum_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t_1,
		double t_2, mui_sampler_exact_3dx *spatial_sampler, mui_chrono_sampler_sum_3dx *temporal_sampler);
double mui_fetch_exact_sum_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t_1,
		double t_2, mui_sampler_exact_3t *spatial_sampler, mui_chrono_sampler_sum_3t *temporal_sampler);
float mui_fetch_gauss_exact_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t_1, float t_2,
		mui_sampler_gauss_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler);
float mui_fetch_gauss_exact_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t_1,
		float t_2, mui_sampler_gauss_3fx *spatial_sampler, mui_chrono_sampler_exact_3fx *temporal_sampler);
double mui_fetch_gauss_exact_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t_1,
		double t_2, mui_sampler_gauss_3d *spatial_sampler, mui_chrono_sampler_exact_3d *temporal_sampler);
double mui_fetch_gauss_exact_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t_1,
		double t_2, mui_sampler_gauss_3dx *spatial_sampler, mui_chrono_sampler_exact_3dx *temporal_sampler);
double mui_fetch_gauss_exact_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t_1,
		double t_2, mui_sampler_gauss_3t *spatial_sampler, mui_chrono_sampler_exact_3t *temporal_sampler);
float mui_fetch_gauss_gauss_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t_1, float t_2,
		mui_sampler_gauss_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler);
float mui_fetch_gauss_gauss_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t_1,
		float t_2, mui_sampler_gauss_3fx *spatial_sampler, mui_chrono_sampler_gauss_3fx *temporal_sampler);
double mui_fetch_gauss_gauss_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t_1,
		double t_2, mui_sampler_gauss_3d *spatial_sampler, mui_chrono_sampler_gauss_3d *temporal_sampler);
double mui_fetch_gauss_gauss_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t_1,
		double t_2, mui_sampler_gauss_3dx *spatial_sampler, mui_chrono_sampler_gauss_3dx *temporal_sampler);
double mui_fetch_gauss_gauss_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t_1,
		double t_2, mui_sampler_gauss_3t *spatial_sampler, mui_chrono_sampler_gauss_3t *temporal_sampler);
float mui_fetch_gauss_mean_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t_1, float t_2,
		mui_sampler_gauss_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler);
float mui_fetch_gauss_mean_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t_1,
		float t_2, mui_sampler_gauss_3fx *spatial_sampler, mui_chrono_sampler_mean_3fx *temporal_sampler);
double mui_fetch_gauss_mean_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t_1,
		double t_2, mui_sampler_gauss_3d *spatial_sampler, mui_chrono_sampler_mean_3d *temporal_sampler);
double mui_fetch_gauss_mean_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t_1,
		double t_2, mui_sampler_gauss_3dx *spatial_sampler, mui_chrono_sampler_mean_3dx *temporal_sampler);
double mui_fetch_gauss_mean_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t_1,
		double t_2, mui_sampler_gauss_3t *spatial_sampler, mui_chrono_sampler_mean_3t *temporal_sampler);
float mui_fetch_moving_average_exact_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t_1,
		float t_2, mui_sampler_moving_average_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler);
float mui_fetch_moving_average_exact_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point,
		float t_1, float t_2, mui_sampler_moving_average_3fx *spatial_sampler,
		mui_chrono_sampler_exact_3fx *temporal_sampler);
double mui_fetch_moving_average_exact_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t_1,
		double t_2, mui_sampler_moving_average_3d *spatial_sampler, mui_chrono_sampler_exact_3d *temporal_sampler);
double mui_fetch_moving_average_exact_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point,
		double t_1, double t_2, mui_sampler_moving_average_3dx *spatial_sampler,
		mui_chrono_sampler_exact_3dx *temporal_sampler);
double mui_fetch_moving_average_exact_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t_1,
		double t_2, mui_sampler_moving_average_3t *spatial_sampler, mui_chrono_sampler_exact_3t *temporal_sampler);
float mui_fetch_moving_average_gauss_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t_1,
		float t_2, mui_sampler_moving_average_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler);
float mui_fetch_moving_average_gauss_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point,
		float t_1, float t_2, mui_sampler_moving_average_3fx *spatial_sampler,
		mui_chrono_sampler_gauss_3fx *temporal_sampler);
double mui_fetch_moving_average_gauss_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t_1,
		double t_2, mui_sampler_moving_average_3d *spatial_sampler, mui_chrono_sampler_gauss_3d *temporal_sampler);
double mui_fetch_moving_average_gauss_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point,
		double t_1, double t_2, mui_sampler_moving_average_3dx *spatial_sampler,
		mui_chrono_sampler_gauss_3dx *temporal_sampler);
double mui_fetch_moving_average_gauss_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t_1,
		double t_2, mui_sampler_moving_average_3t *spatial_sampler, mui_chrono_sampler_gauss_3t *temporal_sampler);
float mui_fetch_moving_average_mean_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t_1,
		float t_2, mui_sampler_moving_average_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler);
float mui_fetch_moving_average_mean_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t_1,
		float t_2, mui_sampler_moving_average_3fx *spatial_sampler, mui_chrono_sampler_mean_3fx *temporal_sampler);
double mui_fetch_moving_average_mean_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t_1,
		double t_2, mui_sampler_moving_average_3d *spatial_sampler, mui_chrono_sampler_mean_3d *temporal_sampler);
double mui_fetch_moving_average_mean_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point,
		double t_1, double t_2, mui_sampler_moving_average_3dx *spatial_sampler,
		mui_chrono_sampler_mean_3dx *temporal_sampler);
double mui_fetch_moving_average_mean_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t_1,
		double t_2, mui_sampler_moving_average_3t *spatial_sampler, mui_chrono_sampler_mean_3t *temporal_sampler);
float mui_fetch_moving_average_sum_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t_1,
		float t_2, mui_sampler_moving_average_3f *spatial_sampler, mui_chrono_sampler_sum_3f *temporal_sampler);
float mui_fetch_moving_average_sum_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t_1,
		float t_2, mui_sampler_moving_average_3fx *spatial_sampler, mui_chrono_sampler_sum_3fx *temporal_sampler);
double mui_fetch_moving_average_sum_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t_1,
		double t_2, mui_sampler_moving_average_3d *spatial_sampler, mui_chrono_sampler_sum_3d *temporal_sampler);
double mui_fetch_moving_average_sum_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point,
		double t_1, double t_2, mui_sampler_moving_average_3dx *spatial_sampler,
		mui_chrono_sampler_sum_3dx *temporal_sampler);
double mui_fetch_moving_average_sum_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t_1,
		double t_2, mui_sampler_moving_average_3t *spatial_sampler, mui_chrono_sampler_sum_3t *temporal_sampler);
float mui_fetch_nearest_neighbor_exact_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t_1,
		float t_2, mui_sampler_nearest_neighbor_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler);
float mui_fetch_nearest_neighbor_exact_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point,
		float t_1, float t_2, mui_sampler_nearest_neighbor_3fx *spatial_sampler,
		mui_chrono_sampler_exact_3fx *temporal_sampler);
double mui_fetch_nearest_neighbor_exact_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point,
		double t_1, double t_2, mui_sampler_nearest_neighbor_3d *spatial_sampler,
		mui_chrono_sampler_exact_3d *temporal_sampler);
double mui_fetch_nearest_neighbor_exact_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point,
		double t_1, double t_2, mui_sampler_nearest_neighbor_3dx *spatial_sampler,
		mui_chrono_sampler_exact_3dx *temporal_sampler);
double mui_fetch_nearest_neighbor_exact_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point,
		double t_1, double t_2, mui_sampler_nearest_neighbor_3t *spatial_sampler,
		mui_chrono_sampler_exact_3t *temporal_sampler);
float mui_fetch_nearest_neighbor_gauss_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t_1,
		float t_2, mui_sampler_nearest_neighbor_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler);
float mui_fetch_nearest_neighbor_gauss_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point,
		float t_1, float t_2, mui_sampler_nearest_neighbor_3fx *spatial_sampler,
		mui_chrono_sampler_gauss_3fx *temporal_sampler);
double mui_fetch_nearest_neighbor_gauss_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point,
		double t_1, double t_2, mui_sampler_nearest_neighbor_3d *spatial_sampler,
		mui_chrono_sampler_gauss_3d *temporal_sampler);
double mui_fetch_nearest_neighbor_gauss_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point,
		double t_1, double t_2, mui_sampler_nearest_neighbor_3dx *spatial_sampler,
		mui_chrono_sampler_gauss_3dx *temporal_sampler);
double mui_fetch_nearest_neighbor_gauss_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point,
		double t_1, double t_2, mui_sampler_nearest_neighbor_3t *spatial_sampler,
		mui_chrono_sampler_gauss_3t *temporal_sampler);
float mui_fetch_nearest_neighbor_mean_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t_1,
		float t_2, mui_sampler_nearest_neighbor_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler);
float mui_fetch_nearest_neighbor_mean_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point,
		float t_1, float t_2, mui_sampler_nearest_neighbor_3fx *spatial_sampler,
		mui_chrono_sampler_mean_3fx *temporal_sampler);
double mui_fetch_nearest_neighbor_mean_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point,
		double t_1, double t_2, mui_sampler_nearest_neighbor_3d *spatial_sampler,
		mui_chrono_sampler_mean_3d *temporal_sampler);
double mui_fetch_nearest_neighbor_mean_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point,
		double t_1, double t_2, mui_sampler_nearest_neighbor_3dx *spatial_sampler,
		mui_chrono_sampler_mean_3dx *temporal_sampler);
double mui_fetch_nearest_neighbor_mean_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point,
		double t_1, double t_2, mui_sampler_nearest_neighbor_3t *spatial_sampler,
		mui_chrono_sampler_mean_3t *temporal_sampler);
float mui_fetch_nearest_neighbor_sum_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t_1,
		float t_2, mui_sampler_nearest_neighbor_3f *spatial_sampler, mui_chrono_sampler_sum_3f *temporal_sampler);
float mui_fetch_nearest_neighbor_sum_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point,
		float t_1, float t_2, mui_sampler_nearest_neighbor_3fx *spatial_sampler,
		mui_chrono_sampler_sum_3fx *temporal_sampler);
double mui_fetch_nearest_neighbor_sum_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t_1,
		double t_2, mui_sampler_nearest_neighbor_3d *spatial_sampler, mui_chrono_sampler_sum_3d *temporal_sampler);
double mui_fetch_nearest_neighbor_sum_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point,
		double t_1, double t_2, mui_sampler_nearest_neighbor_3dx *spatial_sampler,
		mui_chrono_sampler_sum_3dx *temporal_sampler);
double mui_fetch_nearest_neighbor_sum_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t_1,
		double t_2, mui_sampler_nearest_neighbor_3t *spatial_sampler, mui_chrono_sampler_sum_3t *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_exact_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point,
		float t_1, float t_2, mui_sampler_pseudo_nearest_neighbor_3f *spatial_sampler,
		mui_chrono_sampler_exact_3f *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_exact_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point,
		float t_1, float t_2, mui_sampler_pseudo_nearest_neighbor_3fx *spatial_sampler,
		mui_chrono_sampler_exact_3fx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_exact_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_3d *spatial_sampler,
		mui_chrono_sampler_exact_3d *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_exact_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_3dx *spatial_sampler,
		mui_chrono_sampler_exact_3dx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_exact_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_3t *spatial_sampler,
		mui_chrono_sampler_exact_3t *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_gauss_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point,
		float t_1, float t_2, mui_sampler_pseudo_nearest_neighbor_3f *spatial_sampler,
		mui_chrono_sampler_gauss_3f *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_gauss_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point,
		float t_1, float t_2, mui_sampler_pseudo_nearest_neighbor_3fx *spatial_sampler,
		mui_chrono_sampler_gauss_3fx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_gauss_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_3d *spatial_sampler,
		mui_chrono_sampler_gauss_3d *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_gauss_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_3dx *spatial_sampler,
		mui_chrono_sampler_gauss_3dx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_gauss_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_3t *spatial_sampler,
		mui_chrono_sampler_gauss_3t *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_mean_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point,
		float t_1, float t_2, mui_sampler_pseudo_nearest_neighbor_3f *spatial_sampler,
		mui_chrono_sampler_mean_3f *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_mean_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point,
		float t_1, float t_2, mui_sampler_pseudo_nearest_neighbor_3fx *spatial_sampler,
		mui_chrono_sampler_mean_3fx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_mean_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_3d *spatial_sampler,
		mui_chrono_sampler_mean_3d *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_mean_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_3dx *spatial_sampler,
		mui_chrono_sampler_mean_3dx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_mean_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_3t *spatial_sampler,
		mui_chrono_sampler_mean_3t *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_sum_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point,
		float t_1, float t_2, mui_sampler_pseudo_nearest_neighbor_3f *spatial_sampler,
		mui_chrono_sampler_sum_3f *temporal_sampler);
float mui_fetch_pseudo_nearest_neighbor_sum_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point,
		float t_1, float t_2, mui_sampler_pseudo_nearest_neighbor_3fx *spatial_sampler,
		mui_chrono_sampler_sum_3fx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_sum_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_3d *spatial_sampler,
		mui_chrono_sampler_sum_3d *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_sum_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_3dx *spatial_sampler,
		mui_chrono_sampler_sum_3dx *temporal_sampler);
double mui_fetch_pseudo_nearest_neighbor_sum_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point,
		double t_1, double t_2, mui_sampler_pseudo_nearest_neighbor_3t *spatial_sampler,
		mui_chrono_sampler_sum_3t *temporal_sampler);
float mui_fetch_shepard_quintic_exact_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t_1,
		float t_2, mui_sampler_shepard_quintic_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler);
float mui_fetch_shepard_quintic_exact_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point,
		float t_1, float t_2, mui_sampler_shepard_quintic_3fx *spatial_sampler,
		mui_chrono_sampler_exact_3fx *temporal_sampler);
double mui_fetch_shepard_quintic_exact_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point,
		double t_1, double t_2, mui_sampler_shepard_quintic_3d *spatial_sampler,
		mui_chrono_sampler_exact_3d *temporal_sampler);
double mui_fetch_shepard_quintic_exact_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point,
		double t_1, double t_2, mui_sampler_shepard_quintic_3dx *spatial_sampler,
		mui_chrono_sampler_exact_3dx *temporal_sampler);
double mui_fetch_shepard_quintic_exact_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point,
		double t_1, double t_2, mui_sampler_shepard_quintic_3t *spatial_sampler,
		mui_chrono_sampler_exact_3t *temporal_sampler);
float mui_fetch_shepard_quintic_gauss_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t_1,
		float t_2, mui_sampler_shepard_quintic_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler);
float mui_fetch_shepard_quintic_gauss_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point,
		float t_1, float t_2, mui_sampler_shepard_quintic_3fx *spatial_sampler,
		mui_chrono_sampler_gauss_3fx *temporal_sampler);
double mui_fetch_shepard_quintic_gauss_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point,
		double t_1, double t_2, mui_sampler_shepard_quintic_3d *spatial_sampler,
		mui_chrono_sampler_gauss_3d *temporal_sampler);
double mui_fetch_shepard_quintic_gauss_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point,
		double t_1, double t_2, mui_sampler_shepard_quintic_3dx *spatial_sampler,
		mui_chrono_sampler_gauss_3dx *temporal_sampler);
double mui_fetch_shepard_quintic_gauss_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point,
		double t_1, double t_2, mui_sampler_shepard_quintic_3t *spatial_sampler,
		mui_chrono_sampler_gauss_3t *temporal_sampler);
float mui_fetch_shepard_quintic_mean_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t_1,
		float t_2, mui_sampler_shepard_quintic_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler);
float mui_fetch_shepard_quintic_mean_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point,
		float t_1, float t_2, mui_sampler_shepard_quintic_3fx *spatial_sampler,
		mui_chrono_sampler_mean_3fx *temporal_sampler);
double mui_fetch_shepard_quintic_mean_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t_1,
		double t_2, mui_sampler_shepard_quintic_3d *spatial_sampler, mui_chrono_sampler_mean_3d *temporal_sampler);
double mui_fetch_shepard_quintic_mean_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point,
		double t_1, double t_2, mui_sampler_shepard_quintic_3dx *spatial_sampler,
		mui_chrono_sampler_mean_3dx *temporal_sampler);
double mui_fetch_shepard_quintic_mean_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t_1,
		double t_2, mui_sampler_shepard_quintic_3t *spatial_sampler, mui_chrono_sampler_mean_3t *temporal_sampler);
float mui_fetch_shepard_quintic_sum_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t_1,
		float t_2, mui_sampler_shepard_quintic_3f *spatial_sampler, mui_chrono_sampler_sum_3f *temporal_sampler);
float mui_fetch_shepard_quintic_sum_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t_1,
		float t_2, mui_sampler_shepard_quintic_3fx *spatial_sampler, mui_chrono_sampler_sum_3fx *temporal_sampler);
double mui_fetch_shepard_quintic_sum_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t_1,
		double t_2, mui_sampler_shepard_quintic_3d *spatial_sampler, mui_chrono_sampler_sum_3d *temporal_sampler);
double mui_fetch_shepard_quintic_sum_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point,
		double t_1, double t_2, mui_sampler_shepard_quintic_3dx *spatial_sampler,
		mui_chrono_sampler_sum_3dx *temporal_sampler);
double mui_fetch_shepard_quintic_sum_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t_1,
		double t_2, mui_sampler_shepard_quintic_3t *spatial_sampler, mui_chrono_sampler_sum_3t *temporal_sampler);
float mui_fetch_sph_quintic_exact_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t_1,
		float t_2, mui_sampler_sph_quintic_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler);
float mui_fetch_sph_quintic_exact_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t_1,
		float t_2, mui_sampler_sph_quintic_3fx *spatial_sampler, mui_chrono_sampler_exact_3fx *temporal_sampler);
double mui_fetch_sph_quintic_exact_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t_1,
		double t_2, mui_sampler_sph_quintic_3d *spatial_sampler, mui_chrono_sampler_exact_3d *temporal_sampler);
double mui_fetch_sph_quintic_exact_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t_1,
		double t_2, mui_sampler_sph_quintic_3dx *spatial_sampler, mui_chrono_sampler_exact_3dx *temporal_sampler);
double mui_fetch_sph_quintic_exact_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t_1,
		double t_2, mui_sampler_sph_quintic_3t *spatial_sampler, mui_chrono_sampler_exact_3t *temporal_sampler);
float mui_fetch_sph_quintic_gauss_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t_1,
		float t_2, mui_sampler_sph_quintic_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler);
float mui_fetch_sph_quintic_gauss_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t_1,
		float t_2, mui_sampler_sph_quintic_3fx *spatial_sampler, mui_chrono_sampler_gauss_3fx *temporal_sampler);
double mui_fetch_sph_quintic_gauss_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t_1,
		double t_2, mui_sampler_sph_quintic_3d *spatial_sampler, mui_chrono_sampler_gauss_3d *temporal_sampler);
double mui_fetch_sph_quintic_gauss_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t_1,
		double t_2, mui_sampler_sph_quintic_3dx *spatial_sampler, mui_chrono_sampler_gauss_3dx *temporal_sampler);
double mui_fetch_sph_quintic_gauss_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t_1,
		double t_2, mui_sampler_sph_quintic_3t *spatial_sampler, mui_chrono_sampler_gauss_3t *temporal_sampler);
float mui_fetch_sph_quintic_mean_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t_1,
		float t_2, mui_sampler_sph_quintic_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler);
float mui_fetch_sph_quintic_mean_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t_1,
		float t_2, mui_sampler_sph_quintic_3fx *spatial_sampler, mui_chrono_sampler_mean_3fx *temporal_sampler);
double mui_fetch_sph_quintic_mean_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t_1,
		double t_2, mui_sampler_sph_quintic_3d *spatial_sampler, mui_chrono_sampler_mean_3d *temporal_sampler);
double mui_fetch_sph_quintic_mean_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t_1,
		double t_2, mui_sampler_sph_quintic_3dx *spatial_sampler, mui_chrono_sampler_mean_3dx *temporal_sampler);
double mui_fetch_sph_quintic_mean_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t_1,
		double t_2, mui_sampler_sph_quintic_3t *spatial_sampler, mui_chrono_sampler_mean_3t *temporal_sampler);
float mui_fetch_sph_quintic_sum_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t_1,
		float t_2, mui_sampler_sph_quintic_3f *spatial_sampler, mui_chrono_sampler_sum_3f *temporal_sampler);
float mui_fetch_sph_quintic_sum_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t_1,
		float t_2, mui_sampler_sph_quintic_3fx *spatial_sampler, mui_chrono_sampler_sum_3fx *temporal_sampler);
double mui_fetch_sph_quintic_sum_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t_1,
		double t_2, mui_sampler_sph_quintic_3d *spatial_sampler, mui_chrono_sampler_sum_3d *temporal_sampler);
double mui_fetch_sph_quintic_sum_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t_1,
		double t_2, mui_sampler_sph_quintic_3dx *spatial_sampler, mui_chrono_sampler_sum_3dx *temporal_sampler);
double mui_fetch_sph_quintic_sum_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t_1,
		double t_2, mui_sampler_sph_quintic_3t *spatial_sampler, mui_chrono_sampler_sum_3t *temporal_sampler);
float mui_fetch_sum_quintic_exact_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t_1,
		float t_2, mui_sampler_sum_quintic_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler);
float mui_fetch_sum_quintic_exact_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t_1,
		float t_2, mui_sampler_sum_quintic_3fx *spatial_sampler, mui_chrono_sampler_exact_3fx *temporal_sampler);
double mui_fetch_sum_quintic_exact_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t_1,
		double t_2, mui_sampler_sum_quintic_3d *spatial_sampler, mui_chrono_sampler_exact_3d *temporal_sampler);
double mui_fetch_sum_quintic_exact_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t_1,
		double t_2, mui_sampler_sum_quintic_3dx *spatial_sampler, mui_chrono_sampler_exact_3dx *temporal_sampler);
double mui_fetch_sum_quintic_exact_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t_1,
		double t_2, mui_sampler_sum_quintic_3t *spatial_sampler, mui_chrono_sampler_exact_3t *temporal_sampler);
float mui_fetch_sum_quintic_gauss_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t_1,
		float t_2, mui_sampler_sum_quintic_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler);
float mui_fetch_sum_quintic_gauss_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t_1,
		float t_2, mui_sampler_sum_quintic_3fx *spatial_sampler, mui_chrono_sampler_gauss_3fx *temporal_sampler);
double mui_fetch_sum_quintic_gauss_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t_1,
		double t_2, mui_sampler_sum_quintic_3d *spatial_sampler, mui_chrono_sampler_gauss_3d *temporal_sampler);
double mui_fetch_sum_quintic_gauss_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t_1,
		double t_2, mui_sampler_sum_quintic_3dx *spatial_sampler, mui_chrono_sampler_gauss_3dx *temporal_sampler);
double mui_fetch_sum_quintic_gauss_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t_1,
		double t_2, mui_sampler_sum_quintic_3t *spatial_sampler, mui_chrono_sampler_gauss_3t *temporal_sampler);
float mui_fetch_sum_quintic_mean_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t_1,
		float t_2, mui_sampler_sum_quintic_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler);
float mui_fetch_sum_quintic_mean_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t_1,
		float t_2, mui_sampler_sum_quintic_3fx *spatial_sampler, mui_chrono_sampler_mean_3fx *temporal_sampler);
double mui_fetch_sum_quintic_mean_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t_1,
		double t_2, mui_sampler_sum_quintic_3d *spatial_sampler, mui_chrono_sampler_mean_3d *temporal_sampler);
double mui_fetch_sum_quintic_mean_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t_1,
		double t_2, mui_sampler_sum_quintic_3dx *spatial_sampler, mui_chrono_sampler_mean_3dx *temporal_sampler);
double mui_fetch_sum_quintic_mean_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t_1,
		double t_2, mui_sampler_sum_quintic_3t *spatial_sampler, mui_chrono_sampler_mean_3t *temporal_sampler);
float mui_fetch_sum_quintic_sum_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t_1,
		float t_2, mui_sampler_sum_quintic_3f *spatial_sampler, mui_chrono_sampler_sum_3f *temporal_sampler);
float mui_fetch_sum_quintic_sum_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t_1,
		float t_2, mui_sampler_sum_quintic_3fx *spatial_sampler, mui_chrono_sampler_sum_3fx *temporal_sampler);
double mui_fetch_sum_quintic_sum_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t_1,
		double t_2, mui_sampler_sum_quintic_3d *spatial_sampler, mui_chrono_sampler_sum_3d *temporal_sampler);
double mui_fetch_sum_quintic_sum_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t_1,
		double t_2, mui_sampler_sum_quintic_3dx *spatial_sampler, mui_chrono_sampler_sum_3dx *temporal_sampler);
double mui_fetch_sum_quintic_sum_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t_1,
		double t_2, mui_sampler_sum_quintic_3t *spatial_sampler, mui_chrono_sampler_sum_3t *temporal_sampler);
#ifdef USE_RBF
float mui_fetch_rbf_exact_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t_1, float t_2,
		mui_sampler_rbf_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler);
float mui_fetch_rbf_exact_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t_1,
		float t_2, mui_sampler_rbf_3fx *spatial_sampler, mui_chrono_sampler_exact_3fx *temporal_sampler);
double mui_fetch_rbf_exact_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t_1,
		double t_2, mui_sampler_rbf_3d *spatial_sampler, mui_chrono_sampler_exact_3d *temporal_sampler);
double mui_fetch_rbf_exact_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t_1,
		double t_2, mui_sampler_rbf_3dx *spatial_sampler, mui_chrono_sampler_exact_3dx *temporal_sampler);
double mui_fetch_rbf_exact_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t_1,
		double t_2, mui_sampler_rbf_3t *spatial_sampler, mui_chrono_sampler_exact_3t *temporal_sampler);
float mui_fetch_rbf_gauss_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t_1, float t_2,
		mui_sampler_rbf_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler);
float mui_fetch_rbf_gauss_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t_1,
		float t_2, mui_sampler_rbf_3fx *spatial_sampler, mui_chrono_sampler_gauss_3fx *temporal_sampler);
double mui_fetch_rbf_gauss_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t_1,
		double t_2, mui_sampler_rbf_3d *spatial_sampler, mui_chrono_sampler_gauss_3d *temporal_sampler);
double mui_fetch_rbf_gauss_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t_1,
		double t_2, mui_sampler_rbf_3dx *spatial_sampler, mui_chrono_sampler_gauss_3dx *temporal_sampler);
double mui_fetch_rbf_gauss_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t_1,
		double t_2, mui_sampler_rbf_3t *spatial_sampler, mui_chrono_sampler_gauss_3t *temporal_sampler);
float mui_fetch_rbf_mean_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t_1, float t_2,
		mui_sampler_rbf_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler);
float mui_fetch_rbf_mean_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t_1, float t_2,
		mui_sampler_rbf_3fx *spatial_sampler, mui_chrono_sampler_mean_3fx *temporal_sampler);
double mui_fetch_rbf_mean_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t_1, double t_2,
		mui_sampler_rbf_3d *spatial_sampler, mui_chrono_sampler_mean_3d *temporal_sampler);
double mui_fetch_rbf_mean_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t_1,
		double t_2, mui_sampler_rbf_3dx *spatial_sampler, mui_chrono_sampler_mean_3dx *temporal_sampler);
double mui_fetch_rbf_mean_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t_1, double t_2,
		mui_sampler_rbf_3t *spatial_sampler, mui_chrono_sampler_mean_3t *temporal_sampler);
float mui_fetch_rbf_sum_3f_pair(mui_uniface_3f *uniface, const char *attr, mui_point_3f point, float t_1, float t_2,
		mui_sampler_rbf_3f *spatial_sampler, mui_chrono_sampler_sum_3f *temporal_sampler);
float mui_fetch_rbf_sum_3fx_pair(mui_uniface_3fx *uniface, const char *attr, mui_point_3fx point, float t_1, float t_2,
		mui_sampler_rbf_3fx *spatial_sampler, mui_chrono_sampler_sum_3fx *temporal_sampler);
double mui_fetch_rbf_sum_3d_pair(mui_uniface_3d *uniface, const char *attr, mui_point_3d point, double t_1, double t_2,
		mui_sampler_rbf_3d *spatial_sampler, mui_chrono_sampler_sum_3d *temporal_sampler);
double mui_fetch_rbf_sum_3dx_pair(mui_uniface_3dx *uniface, const char *attr, mui_point_3dx point, double t_1,
		double t_2, mui_sampler_rbf_3dx *spatial_sampler, mui_chrono_sampler_sum_3dx *temporal_sampler);
double mui_fetch_rbf_sum_3t_pair(mui_uniface_3t *uniface, const char *attr, mui_point_3t point, double t_1, double t_2,
		mui_sampler_rbf_3t *spatial_sampler, mui_chrono_sampler_sum_3t *temporal_sampler);
#endif

// MUI functions for 3D data point only fetch using single time value
void mui_fetch_points_exact_3f(mui_uniface_3f *uniface, const char *attr, float t,
		mui_chrono_sampler_exact_3f *temporal_sampler, mui_point_3f *ret_points, int *num_points);
void mui_fetch_points_exact_3fx(mui_uniface_3fx *uniface, const char *attr, float t,
		mui_chrono_sampler_exact_3fx *temporal_sampler, mui_point_3fx *ret_points, int *num_points);
void mui_fetch_points_exact_3d(mui_uniface_3d *uniface, const char *attr, double t,
		mui_chrono_sampler_exact_3d *temporal_sampler, mui_point_3d *ret_points, int *num_points);
void mui_fetch_points_exact_3dx(mui_uniface_3dx *uniface, const char *attr, double t,
		mui_chrono_sampler_exact_3dx *temporal_sampler, mui_point_3dx *ret_points, int *num_points);
void mui_fetch_points_exact_3t(mui_uniface_3t *uniface, const char *attr, double t,
		mui_chrono_sampler_exact_3t *temporal_sampler, mui_point_3t *ret_points, int *num_points);
void mui_fetch_points_gauss_3f(mui_uniface_3f *uniface, const char *attr, float t,
		mui_chrono_sampler_gauss_3f *temporal_sampler, mui_point_3f *ret_points, int *num_points);
void mui_fetch_points_gauss_3fx(mui_uniface_3fx *uniface, const char *attr, float t,
		mui_chrono_sampler_gauss_3fx *temporal_sampler, mui_point_3fx *ret_points, int *num_points);
void mui_fetch_points_gauss_3d(mui_uniface_3d *uniface, const char *attr, double t,
		mui_chrono_sampler_gauss_3d *temporal_sampler, mui_point_3d *ret_points, int *num_points);
void mui_fetch_points_gauss_3dx(mui_uniface_3dx *uniface, const char *attr, double t,
		mui_chrono_sampler_gauss_3dx *temporal_sampler, mui_point_3dx *ret_points, int *num_points);
void mui_fetch_points_gauss_3t(mui_uniface_3t *uniface, const char *attr, double t,
		mui_chrono_sampler_gauss_3t *temporal_sampler, mui_point_3t *ret_points, int *num_points);
void mui_fetch_points_mean_3f(mui_uniface_3f *uniface, const char *attr, float t,
		mui_chrono_sampler_mean_3f *temporal_sampler, mui_point_3f *ret_points, int *num_points);
void mui_fetch_points_mean_3fx(mui_uniface_3fx *uniface, const char *attr, float t,
		mui_chrono_sampler_mean_3fx *temporal_sampler, mui_point_3fx *ret_points, int *num_points);
void mui_fetch_points_mean_3d(mui_uniface_3d *uniface, const char *attr, double t,
		mui_chrono_sampler_mean_3d *temporal_sampler, mui_point_3d *ret_points, int *num_points);
void mui_fetch_points_mean_3dx(mui_uniface_3dx *uniface, const char *attr, double t,
		mui_chrono_sampler_mean_3dx *temporal_sampler, mui_point_3dx *ret_points, int *num_points);
void mui_fetch_points_mean_3t(mui_uniface_3t *uniface, const char *attr, double t,
		mui_chrono_sampler_mean_3t *temporal_sampler, mui_point_3t *ret_points, int *num_points);
void mui_fetch_points_sum_3f(mui_uniface_3f *uniface, const char *attr, float t,
		mui_chrono_sampler_sum_3f *temporal_sampler, mui_point_3f *ret_points, int *num_points);
void mui_fetch_points_sum_3fx(mui_uniface_3fx *uniface, const char *attr, float t,
		mui_chrono_sampler_sum_3fx *temporal_sampler, mui_point_3fx *ret_points, int *num_points);
void mui_fetch_points_sum_3d(mui_uniface_3d *uniface, const char *attr, double t,
		mui_chrono_sampler_sum_3d *temporal_sampler, mui_point_3d *ret_points, int *num_points);
void mui_fetch_points_sum_3dx(mui_uniface_3dx *uniface, const char *attr, double t,
		mui_chrono_sampler_sum_3dx *temporal_sampler, mui_point_3dx *ret_points, int *num_points);
void mui_fetch_points_sum_3t(mui_uniface_3t *uniface, const char *attr, double t,
		mui_chrono_sampler_sum_3t *temporal_sampler, mui_point_3t *ret_points, int *num_points);

// MUI functions for 3D data point only fetch using two time values
void mui_fetch_values_exact_3f_pair(mui_uniface_3f *uniface, const char *attr, float t_1, float t_2,
		mui_chrono_sampler_exact_3f *temporal_sampler, float *values, int *num_values);
void mui_fetch_values_exact_3fx_pair(mui_uniface_3fx *uniface, const char *attr, float t_1, float t_2,
		mui_chrono_sampler_exact_3fx *temporal_sampler, float *values, int *num_values);
void mui_fetch_values_exact_3d_pair(mui_uniface_3d *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_exact_3d *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_exact_3dx_pair(mui_uniface_3dx *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_exact_3dx *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_exact_3t_pair(mui_uniface_3t *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_exact_3t *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_gauss_3f_pair(mui_uniface_3f *uniface, const char *attr, float t_1, float t_2,
		mui_chrono_sampler_gauss_3f *temporal_sampler, float *values, int *num_values);
void mui_fetch_values_gauss_3fx_pair(mui_uniface_3fx *uniface, const char *attr, float t_1, float t_2,
		mui_chrono_sampler_gauss_3fx *temporal_sampler, float *values, int *num_values);
void mui_fetch_values_gauss_3d_pair(mui_uniface_3d *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_gauss_3d *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_gauss_3dx_pair(mui_uniface_3dx *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_gauss_3dx *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_gauss_3t_pair(mui_uniface_3t *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_gauss_3t *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_mean_3f_pair(mui_uniface_3f *uniface, const char *attr, float t_1, float t_2,
		mui_chrono_sampler_mean_3f *temporal_sampler, float *values, int *num_values);
void mui_fetch_values_mean_3fx_pair(mui_uniface_3fx *uniface, const char *attr, float t_1, float t_2,
		mui_chrono_sampler_mean_3fx *temporal_sampler, float *values, int *num_values);
void mui_fetch_values_mean_3d_pair(mui_uniface_3d *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_mean_3d *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_mean_3dx_pair(mui_uniface_3dx *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_mean_3dx *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_mean_3t_pair(mui_uniface_3t *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_mean_3t *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_sum_3f_pair(mui_uniface_3f *uniface, const char *attr, float t_1, float t_2,
		mui_chrono_sampler_sum_3f *temporal_sampler, float *values, int *num_values);
void mui_fetch_values_sum_3fx_pair(mui_uniface_3fx *uniface, const char *attr, float t_1, float t_2,
		mui_chrono_sampler_sum_3fx *temporal_sampler, float *values, int *num_values);
void mui_fetch_values_sum_3d_pair(mui_uniface_3d *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_sum_3d *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_sum_3dx_pair(mui_uniface_3dx *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_sum_3dx *temporal_sampler, double *values, int *num_values);
void mui_fetch_values_sum_3t_pair(mui_uniface_3t *uniface, const char *attr, double t_1, double t_2,
		mui_chrono_sampler_sum_3t *temporal_sampler, double *values, int *num_values);

// MUI functions for single parameter fetch
float mui_fetch_3f_param(mui_uniface_3f *uniface, const char *attr);
float mui_fetch_3fx_param(mui_uniface_3fx *uniface, const char *attr);
double mui_fetch_3d_param(mui_uniface_3d *uniface, const char *attr);
double mui_fetch_3dx_param(mui_uniface_3dx *uniface, const char *attr);
double mui_fetch_3t_param(mui_uniface_3t *uniface, const char *attr);

// MUI data receive test functions
int mui_is_ready_3f(mui_uniface_3f *uniface, const char *attr, float t);
int mui_is_ready_3fx(mui_uniface_3fx *uniface, const char *attr, float t);
int mui_is_ready_3d(mui_uniface_3d *uniface, const char *attr, double t);
int mui_is_ready_3dx(mui_uniface_3dx *uniface, const char *attr, double t);
int mui_is_ready_3t(mui_uniface_3t *uniface, const char *attr, double t);
int mui_is_ready_3f_pair(mui_uniface_3f *uniface, const char *attr, float t_1, float t_2);
int mui_is_ready_3fx_pair(mui_uniface_3fx *uniface, const char *attr, float t_1, float t_2);
int mui_is_ready_3d_pair(mui_uniface_3d *uniface, const char *attr, double t_1, double t_2);
int mui_is_ready_3dx_pair(mui_uniface_3dx *uniface, const char *attr, double t_1, double t_2);
int mui_is_ready_3t_pair(mui_uniface_3t *uniface, const char *attr, double t_1, double t_2);

// MUI Smart Send functions
void mui_announce_send_span_3f_box(mui_uniface_3f *uniface, float box_1_1, float box_1_2, float box_1_3, float box_2_1,
		float box_2_2, float box_2_3, float t_start, float t_timeout, int synchronised);
void mui_announce_send_span_3fx_box(mui_uniface_3fx *uniface, float box_1_1, float box_1_2, float box_1_3,
		float box_2_1, float box_2_2, float box_2_3, float t_start, float t_timeout, int synchronised);
void mui_announce_send_span_3d_box(mui_uniface_3d *uniface, double box_1_1, double box_1_2, double box_1_3,
		double box_2_1, double box_2_2, double box_2_3, double t_start, double t_timeout, int synchronised);
void mui_announce_send_span_3dx_box(mui_uniface_3dx *uniface, double box_1_1, double box_1_2, double box_1_3,
		double box_2_1, double box_2_2, double box_2_3, double t_start, double t_timeout, int synchronised);
void mui_announce_send_span_3t_box(mui_uniface_3t *uniface, double box_1_1, double box_1_2, double box_1_3,
		double box_2_1, double box_2_2, double box_2_3, double t_start, double t_timeout, int synchronised);
void mui_announce_send_span_3f_sphere(mui_uniface_3f *uniface, mui_point_3f centre, float radius, float t_start,
		float t_timeout, int synchronised);
void mui_announce_send_span_3fx_sphere(mui_uniface_3fx *uniface, mui_point_3fx centre, float radius, float t_start,
		float t_timeout, int synchronised);
void mui_announce_send_span_3d_sphere(mui_uniface_3d *uniface, mui_point_3d centre, double radius, double t_start,
		double t_timeout, int synchronised);
void mui_announce_send_span_3dx_sphere(mui_uniface_3dx *uniface, mui_point_3dx centre, double radius, double t_start,
		double t_timeout, int synchronised);
void mui_announce_send_span_3t_sphere(mui_uniface_3t *uniface, mui_point_3t centre, double radius, double t_start,
		double t_timeout, int synchronised);
void mui_announce_recv_span_3f_box(mui_uniface_3f *uniface, float box_1_1, float box_1_2, float box_1_3, float box_2_1,
		float box_2_2, float box_2_3, float t_start, float t_timeout, int synchronised);
void mui_announce_recv_span_3fx_box(mui_uniface_3fx *uniface, float box_1_1, float box_1_2, float box_1_3,
		float box_2_1, float box_2_2, float box_2_3, float t_start, float t_timeout, int synchronised);
void mui_announce_recv_span_3d_box(mui_uniface_3d *uniface, double box_1_1, double box_1_2, double box_1_3,
		double box_2_1, double box_2_2, double box_2_3, double t_start, double t_timeout, int synchronised);
void mui_announce_recv_span_3dx_box(mui_uniface_3dx *uniface, double box_1_1, double box_1_2, double box_1_3,
		double box_2_1, double box_2_2, double box_2_3, double t_start, double t_timeout, int synchronised);
void mui_announce_recv_span_3t_box(mui_uniface_3t *uniface, double box_1_1, double box_1_2, double box_1_3,
		double box_2_1, double box_2_2, double box_2_3, double t_start, double t_timeout, int synchronised);
void mui_announce_recv_span_3f_sphere(mui_uniface_3f *uniface, mui_point_3f centre, float radius, float t_start,
		float t_timeout, int synchronised);
void mui_announce_recv_span_3fx_sphere(mui_uniface_3fx *uniface, mui_point_3fx centre, float radius, float t_start,
		float t_timeout, int synchronised);
void mui_announce_recv_span_3d_sphere(mui_uniface_3d *uniface, mui_point_3d centre, double radius, double t_start,
		double t_timeout, int synchronised);
void mui_announce_recv_span_3dx_sphere(mui_uniface_3dx *uniface, mui_point_3dx centre, double radius, double t_start,
		double t_timeout, int synchronised);
void mui_announce_recv_span_3t_sphere(mui_uniface_3t *uniface, mui_point_3t centre, double radius, double t_start,
		double t_timeout, int synchronised);
void mui_announce_send_disable_3f(mui_uniface_3f *uniface, int synchronised);
void mui_announce_send_disable_3fx(mui_uniface_3fx *uniface, int synchronised);
void mui_announce_send_disable_3d(mui_uniface_3d *uniface, int synchronised);
void mui_announce_send_disable_3dx(mui_uniface_3dx *uniface, int synchronised);
void mui_announce_send_disable_3t(mui_uniface_3t *uniface, int synchronised);
void mui_announce_recv_disable_3f(mui_uniface_3f *uniface, int synchronised);
void mui_announce_recv_disable_3fx(mui_uniface_3fx *uniface, int synchronised);
void mui_announce_recv_disable_3d(mui_uniface_3d *uniface, int synchronised);
void mui_announce_recv_disable_3dx(mui_uniface_3dx *uniface, int synchronised);
void mui_announce_recv_disable_3t(mui_uniface_3t *uniface, int synchronised);

// MUI barrier functions
void mui_barrier_3f(mui_uniface_3f *uniface, float t);
void mui_barrier_3fx(mui_uniface_3fx *uniface, float t);
void mui_barrier_3d(mui_uniface_3d *uniface, double t);
void mui_barrier_3dx(mui_uniface_3dx *uniface, double t);
void mui_barrier_3t(mui_uniface_3t *uniface, double t);
void mui_barrier_3f_pair(mui_uniface_3f *uniface, float t_1, float t_2);
void mui_barrier_3fx_pair(mui_uniface_3fx *uniface, float t_1, float t_2);
void mui_barrier_3d_pair(mui_uniface_3d *uniface, double t_1, double t_2);
void mui_barrier_3dx_pair(mui_uniface_3dx *uniface, double t_1, double t_2);
void mui_barrier_3t_pair(mui_uniface_3t *uniface, double t_1, double t_2);

// MUI forget functions
void mui_forget_upper_3f(mui_uniface_3f *uniface, float upper, int reset_log);
void mui_forget_upper_3fx(mui_uniface_3fx *uniface, float upper, int reset_log);
void mui_forget_upper_3d(mui_uniface_3d *uniface, double upper, int reset_log);
void mui_forget_upper_3dx(mui_uniface_3dx *uniface, double upper, int reset_log);
void mui_forget_upper_3t(mui_uniface_3t *uniface, double upper, int reset_log);
void mui_forget_upper_3f_pair(mui_uniface_3f *uniface, float upper_1, float upper_2, int reset_log);
void mui_forget_upper_3fx_pair(mui_uniface_3fx *uniface, float upper_1, float upper_2, int reset_log);
void mui_forget_upper_3d_pair(mui_uniface_3d *uniface, double upper_1, double upper_2, int reset_log);
void mui_forget_upper_3dx_pair(mui_uniface_3dx *uniface, double upper_1, double upper_2, int reset_log);
void mui_forget_upper_3t_pair(mui_uniface_3t *uniface, double upper_1, double upper_2, int reset_log);
void mui_forget_lower_upper_3f(mui_uniface_3f *uniface, float lower, float upper, int reset_log);
void mui_forget_lower_upper_3fx(mui_uniface_3fx *uniface, float lower, float upper, int reset_log);
void mui_forget_lower_upper_3d(mui_uniface_3d *uniface, double lower, double upper, int reset_log);
void mui_forget_lower_upper_3dx(mui_uniface_3dx *uniface, double lower, double upper, int reset_log);
void mui_forget_lower_upper_3t(mui_uniface_3t *uniface, double lower, double upper, int reset_log);
void mui_forget_lower_upper_3f_pair(mui_uniface_3f *uniface, float lower_1, float lower_2, float upper_1, float upper_2,
		int reset_log);
void mui_forget_lower_upper_3fx_pair(mui_uniface_3fx *uniface, float lower_1, float lower_2, float upper_1,
		float upper_2, int reset_log);
void mui_forget_lower_upper_3d_pair(mui_uniface_3d *uniface, double lower_1, double lower_2, double upper_1,
		double upper_2, int reset_log);
void mui_forget_lower_upper_3dx_pair(mui_uniface_3dx *uniface, double lower_1, double lower_2, double upper_1,
		double upper_2, int reset_log);
void mui_forget_lower_upper_3t_pair(mui_uniface_3t *uniface, double lower_1, double lower_2, double upper_1,
		double upper_2, int reset_log);
void mui_set_forget_length_3f(mui_uniface_3f *uniface, float length);
void mui_set_forget_length_3fx(mui_uniface_3fx *uniface, float length);
void mui_set_forget_length_3d(mui_uniface_3d *uniface, double length);
void mui_set_forget_length_3dx(mui_uniface_3dx *uniface, double length);
void mui_set_forget_length_3t(mui_uniface_3t *uniface, double length);
const char* mui_uri_host_3f(mui_uniface_3f *uniface);
const char* mui_uri_host_3fx(mui_uniface_3fx *uniface);
const char* mui_uri_host_3d(mui_uniface_3d *uniface);
const char* mui_uri_host_3dx(mui_uniface_3dx *uniface);
const char* mui_uri_host_3t(mui_uniface_3t *uniface);
const char* mui_uri_path_3f(mui_uniface_3f *uniface);
const char* mui_uri_path_3fx(mui_uniface_3fx *uniface);
const char* mui_uri_path_3d(mui_uniface_3d *uniface);
const char* mui_uri_path_3dx(mui_uniface_3dx *uniface);
const char* mui_uri_path_3t(mui_uniface_3t *uniface);
const char* mui_uri_protocol_3f(mui_uniface_3f *uniface);
const char* mui_uri_protocol_3fx(mui_uniface_3fx *uniface);
const char* mui_uri_protocol_3d(mui_uniface_3d *uniface);
const char* mui_uri_protocol_3dx(mui_uniface_3dx *uniface);
const char* mui_uri_protocol_3t(mui_uniface_3t *uniface);

#endif /* MUI_C_WRAPPER_3D_H_ */
