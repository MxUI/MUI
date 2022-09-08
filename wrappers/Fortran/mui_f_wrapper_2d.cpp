/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library                       *
*                                                                            *
* Copyright (C) 2021 Y. H. Tang, S. Kudo, X. Bian, Z. Li, G. E. Karniadakis, *
*                    S. M. Longshaw, W. Liu                                  *
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
 * @file mui_f_wrapper_2d.cpp
 * @author S. M. Longshaw (derived from original 3D wrapper by S. Kudo)
 * @date Dec 07, 2021
 * @brief C interface for Fortran wrapper to create and manage 2D MUI interfaces
 *        and associated sampler objects
 *
 *        NOTE: Any point co-ordinates are enumerated rather than assuming
 *              Cartesian form, i.e. {1, 2, 3} rather than {x, y, z}.
 */

// Main MUI header include (contains any other needed includes)
#include "../../mui.h"
// Include config header local to C wrapper
#include "config_f_wrapper.h"

extern "C" {

// MUI Interface typedefs for specialism creation
typedef mui::uniface2f mui_uniface_2f;
typedef mui::uniface2fx mui_uniface_2fx;
typedef mui::uniface2d mui_uniface_2d;
typedef mui::uniface2dx mui_uniface_2dx;

// MUI Interface typedefs for template creation (recommended)
typedef mui::uniface<mui::mui_f_wrapper_2D> mui_uniface_2t;

// Exact spatial sampler typedefs for specialism creation
typedef mui::sampler_exact2f<float> mui_sampler_exact_2f;
typedef mui::sampler_exact2fx<float> mui_sampler_exact_2fx;
typedef mui::sampler_exact2d<double> mui_sampler_exact_2d;
typedef mui::sampler_exact2dx<double> mui_sampler_exact_2dx;

// Exact spatial sampler typedef for template creation (recommended)
typedef mui::sampler_exact<mui::mui_f_wrapper_2D> mui_sampler_exact_2t;

// Gaussian spatial sampler typedefs for specialism creation
typedef mui::sampler_gauss2f<float> mui_sampler_gauss_2f;
typedef mui::sampler_gauss2fx<float> mui_sampler_gauss_2fx;
typedef mui::sampler_gauss2d<double> mui_sampler_gauss_2d;
typedef mui::sampler_gauss2dx<double> mui_sampler_gauss_2dx;

// Gaussian spatial sampler typedef for template creation (recommended)
typedef mui::sampler_gauss<mui::mui_f_wrapper_2D> mui_sampler_gauss_2t;

// Moving average spatial sampler typedefs for specialism creation
typedef mui::sampler_moving_average2f<float> mui_sampler_moving_average_2f;
typedef mui::sampler_moving_average2fx<float> mui_sampler_moving_average_2fx;
typedef mui::sampler_moving_average2d<double> mui_sampler_moving_average_2d;
typedef mui::sampler_moving_average2dx<double> mui_sampler_moving_average_2dx;

// Moving average spatial sampler typedef for template creation (recommended)
typedef mui::sampler_moving_average<mui::mui_f_wrapper_2D> mui_sampler_moving_average_2t;

// Nearest neighbor spatial sampler typedefs for specialism creation
typedef mui::sampler_nearest_neighbor2f<float> mui_sampler_nearest_neighbor_2f;
typedef mui::sampler_nearest_neighbor2fx<float> mui_sampler_nearest_neighbor_2fx;
typedef mui::sampler_nearest_neighbor2d<double> mui_sampler_nearest_neighbor_2d;
typedef mui::sampler_nearest_neighbor2dx<double> mui_sampler_nearest_neighbor_2dx;

// Moving average spatial sampler typedef for template creation (recommended)
typedef mui::sampler_nearest_neighbor<mui::mui_f_wrapper_2D> mui_sampler_nearest_neighbor_2t;

// Pseudo-linear n^2 spatial sampler typedefs for specialism creation
typedef mui::sampler_pseudo_n2_linear2f<float> mui_sampler_pseudo_n2_linear_2f;
typedef mui::sampler_pseudo_n2_linear2fx<float> mui_sampler_pseudo_n2_linear_2fx;
typedef mui::sampler_pseudo_n2_linear2d<double> mui_sampler_pseudo_n2_linear_2d;
typedef mui::sampler_pseudo_n2_linear2dx<double> mui_sampler_pseudo_n2_linear_2dx;

// Pseudo-linear n^2 spatial sampler typedef for template creation (recommended)
typedef mui::sampler_pseudo_n2_linear<mui::mui_f_wrapper_2D> mui_sampler_pseudo_n2_linear_2t;

// Pseudo-nearest neighbor spatial sampler typedefs for specialism creation
typedef mui::sampler_pseudo_nearest_neighbor2f<float> mui_sampler_pseudo_nearest_neighbor_2f;
typedef mui::sampler_pseudo_nearest_neighbor2fx<float> mui_sampler_pseudo_nearest_neighbor_2fx;
typedef mui::sampler_pseudo_nearest_neighbor2d<double> mui_sampler_pseudo_nearest_neighbor_2d;
typedef mui::sampler_pseudo_nearest_neighbor2dx<double> mui_sampler_pseudo_nearest_neighbor_2dx;

// Pseudo-nearest neighbor spatial sampler typedef for template creation (recommended)
typedef mui::sampler_pseudo_nearest_neighbor<mui::mui_f_wrapper_2D> mui_sampler_pseudo_nearest_neighbor_2t;

// Shepard interpolation with quintic kernel spatial sampler typedefs for specialism creation
typedef mui::sampler_shepard_quintic2f<float> mui_sampler_shepard_quintic_2f;
typedef mui::sampler_shepard_quintic2fx<float> mui_sampler_shepard_quintic_2fx;
typedef mui::sampler_shepard_quintic2d<double> mui_sampler_shepard_quintic_2d;
typedef mui::sampler_shepard_quintic2dx<double> mui_sampler_shepard_quintic_2dx;

// Shepard interpolation with quintic kernel spatial sampler typedef for template creation (recommended)
typedef mui::sampler_shepard_quintic<mui::mui_f_wrapper_2D> mui_sampler_shepard_quintic_2t;

// Spatial sampler that provides a value at a point using a Smoothed Particle Hydrodynamics (SPH) derived
// interpolation method with a quintic spline kernel spatial sampler typedefs for specialism creation
typedef mui::sampler_sph_quintic2f<float> mui_sampler_sph_quintic_2f;
typedef mui::sampler_sph_quintic2fx<float> mui_sampler_sph_quintic_2fx;
typedef mui::sampler_sph_quintic2d<double> mui_sampler_sph_quintic_2d;
typedef mui::sampler_sph_quintic2dx<double> mui_sampler_sph_quintic_2dx;

// Spatial sampler that provides a value at a point using a Smoothed Particle Hydrodynamics (SPH) derived
// interpolation method with a quintic spline kernel spatial sampler typedefs for template creation (recommended)
typedef mui::sampler_sph_quintic<mui::mui_f_wrapper_2D> mui_sampler_sph_quintic_2t;

// Summation with quintic kernel spatial sampler typedefs for specialism creation
typedef mui::sampler_sum_quintic2f<float> mui_sampler_sum_quintic_2f;
typedef mui::sampler_sum_quintic2fx<float> mui_sampler_sum_quintic_2fx;
typedef mui::sampler_sum_quintic2d<double> mui_sampler_sum_quintic_2d;
typedef mui::sampler_sum_quintic2dx<double> mui_sampler_sum_quintic_2dx;

// Summation with quintic kernel spatial sampler typedef for template creation (recommended)
typedef mui::sampler_sum_quintic<mui::mui_f_wrapper_2D> mui_sampler_sum_quintic_2t;

#ifdef USE_RBF
// Radial Basis Function (RBF) spatial sampler typedefs for specialism creation
typedef mui::sampler_rbf2f<float> mui_sampler_rbf_2f;
typedef mui::sampler_rbf2fx<float> mui_sampler_rbf_2fx;
typedef mui::sampler_rbf2d<double> mui_sampler_rbf_2d;
typedef mui::sampler_rbf2dx<double> mui_sampler_rbf_2dx;

// Radial Basis Function (RBF) spatial sampler typedef for template creation (recommended)
typedef mui::sampler_rbf<mui::mui_f_wrapper_2D> mui_sampler_rbf_2t;
#endif

// Exact temporal sampler typedefs for specialism creation
typedef mui::chrono_sampler_exact2f mui_chrono_sampler_exact_2f;
typedef mui::chrono_sampler_exact2fx mui_chrono_sampler_exact_2fx;
typedef mui::chrono_sampler_exact2d mui_chrono_sampler_exact_2d;
typedef mui::chrono_sampler_exact2dx mui_chrono_sampler_exact_2dx;

// Exact temporal sampler typedefs for template creation (recommended)
typedef mui::chrono_sampler_exact<mui::mui_f_wrapper_2D> mui_chrono_sampler_exact_2t;

// Gaussian temporal sampler typedefs for specialism creation
typedef mui::chrono_sampler_gauss2f mui_chrono_sampler_gauss_2f;
typedef mui::chrono_sampler_gauss2fx mui_chrono_sampler_gauss_2fx;
typedef mui::chrono_sampler_gauss2d mui_chrono_sampler_gauss_2d;
typedef mui::chrono_sampler_gauss2dx mui_chrono_sampler_gauss_2dx;

// Gaussian temporal sampler typedefs for template creation (recommended)
typedef mui::chrono_sampler_gauss<mui::mui_f_wrapper_2D> mui_chrono_sampler_gauss_2t;

// Mean average temporal sampler typedefs for specialism creation
typedef mui::chrono_sampler_mean2f mui_chrono_sampler_mean_2f;
typedef mui::chrono_sampler_mean2fx mui_chrono_sampler_mean_2fx;
typedef mui::chrono_sampler_mean2d mui_chrono_sampler_mean_2d;
typedef mui::chrono_sampler_mean2dx mui_chrono_sampler_mean_2dx;

// Mean average temporal sampler typedefs for template creation (recommended)
typedef mui::chrono_sampler_mean<mui::mui_f_wrapper_2D> mui_chrono_sampler_mean_2t;

// Summation temporal sampler typedefs for specialism creation
typedef mui::chrono_sampler_sum2f mui_chrono_sampler_sum_2f;
typedef mui::chrono_sampler_sum2fx mui_chrono_sampler_sum_2fx;
typedef mui::chrono_sampler_sum2d mui_chrono_sampler_sum_2d;
typedef mui::chrono_sampler_sum2dx mui_chrono_sampler_sum_2dx;

// Summation temporal sampler typedefs for template creation (recommended)
typedef mui::chrono_sampler_sum<mui::mui_f_wrapper_2D> mui_chrono_sampler_sum_2t;

// MUI set of specialism interface on multi-domain coupling
mui_uniface_2f** mui_uniface_multi_2f;
mui_uniface_2fx** mui_uniface_multi_2fx;
mui_uniface_2d** mui_uniface_multi_2d;
mui_uniface_2dx** mui_uniface_multi_2dx;

// MUI set of template interface on multi-domain coupling
mui_uniface_2t** mui_uniface_multi_2t;

/****************************************
 * Create MUI interfaces                 *
 ****************************************/

// 2d interface with float=single and int=int32
void mui_create_uniface_2f_f(mui_uniface_2f **ret, const char *URI) {
    *ret = new mui_uniface_2f(URI);
}

// 2d interface with float=single and int=int64
void mui_create_uniface_2fx_f(mui_uniface_2fx **ret, const char *URI) {
    *ret = new mui_uniface_2fx(URI);
}

// 2d interface with float=double and int=int32
void mui_create_uniface_2d_f(mui_uniface_2d **ret, const char *URI) {
    *ret = new mui_uniface_2d(URI);
}

// 2d interface with float=double and int=int64
void mui_create_uniface_2dx_f(mui_uniface_2dx **ret, const char *URI) {
    *ret = new mui_uniface_2dx(URI);
}

// 2d interface using config from config_f_wrapper.h
void mui_create_uniface_2t_f(mui_uniface_2t **ret, const char *URI) {
    *ret = new mui_uniface_2t(URI);
}

// Set of 2D interfaces with float=single and int=int32
void mui_create_uniface_multi_2f_f(const char *domain, const char **interfaces, int interface_count ) {

    std::vector<std::string> interface_names;
    size_t interface_element_size = sizeof(interfaces)/sizeof(interfaces[0]);

    if (interface_element_size == 1) {
        char interfaces_jointed[strlen(interfaces[0])];
        strcpy(interfaces_jointed, interfaces[0]);
        char delim[] = " ";

        char * pchar = strtok (interfaces_jointed,delim);
        int interface_names_count = 0;
        for(size_t i=0; i<interface_count; i++) {
            if(pchar == NULL) {
                std::cerr << "MUI Error [mui_f_wrapper_2d.cpp]: Error MUI Fortran C binding: NULL interfaces at the " << i << "th interface." << std::endl;
                std::abort();
            }
            interface_names.push_back(std::string(pchar));
            pchar = strtok (NULL, delim);
            interface_names_count = interface_names_count + 1;
        }

    } else if(interface_element_size > 1) {
        for(size_t i=0; i<interface_count; i++)
            interface_names.push_back(std::string(interfaces[i]));
    } else {
        std::cerr << "MUI Error [mui_f_wrapper_2d.cpp]: Error MUI Fortran C binding: The size of interfaces array must larger or equals to one." << std::endl;
        std::abort();
    }

    auto created_unifaces = mui::create_uniface<mui::config_2f>(domain, interface_names);

    mui_uniface_multi_2f = new mui_uniface_2f*[created_unifaces.size()];

      for(size_t i=0; i<created_unifaces.size(); i++) {
        mui_uniface_multi_2f[i] = created_unifaces[i].release();
    }
}

// Set of 2D interfaces with float=single and int=int64
void mui_create_uniface_multi_2fx_f(const char *domain, const char **interfaces, int interface_count ) {

    std::vector<std::string> interface_names;
    size_t interface_element_size = sizeof(interfaces)/sizeof(interfaces[0]);

    if (interface_element_size == 1) {
        char interfaces_jointed[strlen(interfaces[0])];
        strcpy(interfaces_jointed, interfaces[0]);
        char delim[] = " ";

        char * pchar = strtok (interfaces_jointed,delim);
        int interface_names_count = 0;
        for(size_t i=0; i<interface_count; i++) {
            if(pchar == NULL) {
                std::cerr << "MUI Error [mui_f_wrapper_2d.cpp]: Error MUI Fortran C binding: NULL interfaces at the " << i << "th interface." << std::endl;
                std::abort();
            }
            interface_names.push_back(std::string(pchar));
            pchar = strtok (NULL, delim);
            interface_names_count = interface_names_count + 1;
        }

    } else if(interface_element_size > 1) {
        for(size_t i=0; i<interface_count; i++)
            interface_names.push_back(std::string(interfaces[i]));
    } else {
        std::cerr << "MUI Error [mui_f_wrapper_2d.cpp]: Error MUI Fortran C binding: The size of interfaces array must larger or equals to one." << std::endl;
        std::abort();
    }

    auto created_unifaces = mui::create_uniface<mui::config_2fx>(domain, interface_names);

    mui_uniface_multi_2fx = new mui_uniface_2fx*[created_unifaces.size()];

      for(size_t i=0; i<created_unifaces.size(); i++) {
        mui_uniface_multi_2fx[i] = created_unifaces[i].release();
    }
}

// Set of 2D interfaces with float=double and int=int32
void mui_create_uniface_multi_2d_f(const char *domain, const char **interfaces, int interface_count ) {

    std::vector<std::string> interface_names;
    size_t interface_element_size = sizeof(interfaces)/sizeof(interfaces[0]);

    if (interface_element_size == 1) {
        char interfaces_jointed[strlen(interfaces[0])];
        strcpy(interfaces_jointed, interfaces[0]);
        char delim[] = " ";

        char * pchar = strtok (interfaces_jointed,delim);
        int interface_names_count = 0;
        for(size_t i=0; i<interface_count; i++) {
            if(pchar == NULL) {
                std::cerr << "MUI Error [mui_f_wrapper_2d.cpp]: Error MUI Fortran C binding: NULL interfaces at the " << i << "th interface." << std::endl;
                std::abort();
            }
            interface_names.push_back(std::string(pchar));
            pchar = strtok (NULL, delim);
            interface_names_count = interface_names_count + 1;
        }

    } else if(interface_element_size > 1) {
        for(size_t i=0; i<interface_count; i++)
            interface_names.push_back(std::string(interfaces[i]));
    } else {
        std::cerr << "MUI Error [mui_f_wrapper_2d.cpp]: Error MUI Fortran C binding: The size of interfaces array must larger or equals to one." << std::endl;
        std::abort();
    }

    auto created_unifaces = mui::create_uniface<mui::config_2d>(domain, interface_names);

    mui_uniface_multi_2d = new mui_uniface_2d*[created_unifaces.size()];

      for(size_t i=0; i<created_unifaces.size(); i++) {
        mui_uniface_multi_2d[i] = created_unifaces[i].release();
    }
}

// Set of 2D interfaces with float=double and int=int64
void mui_create_uniface_multi_2dx_f(const char *domain, const char **interfaces, int interface_count ) {

    std::vector<std::string> interface_names;
    size_t interface_element_size = sizeof(interfaces)/sizeof(interfaces[0]);

    if (interface_element_size == 1) {
        char interfaces_jointed[strlen(interfaces[0])];
        strcpy(interfaces_jointed, interfaces[0]);
        char delim[] = " ";

        char * pchar = strtok (interfaces_jointed,delim);
        int interface_names_count = 0;
        for(size_t i=0; i<interface_count; i++) {
            if(pchar == NULL) {
                std::cerr << "MUI Error [mui_f_wrapper_2d.cpp]: Error MUI Fortran C binding: NULL interfaces at the " << i << "th interface." << std::endl;
                std::abort();
            }
            interface_names.push_back(std::string(pchar));
            pchar = strtok (NULL, delim);
            interface_names_count = interface_names_count + 1;
        }

    } else if(interface_element_size > 1) {
        for(size_t i=0; i<interface_count; i++)
            interface_names.push_back(std::string(interfaces[i]));
    } else {
        std::cerr << "MUI Error [mui_f_wrapper_2d.cpp]: Error MUI Fortran C binding: The size of interfaces array must larger or equals to one." << std::endl;
        std::abort();
    }

    auto created_unifaces = mui::create_uniface<mui::config_2dx>(domain, interface_names);

    mui_uniface_multi_2dx = new mui_uniface_2dx*[created_unifaces.size()];

      for(size_t i=0; i<created_unifaces.size(); i++) {
        mui_uniface_multi_2dx[i] = created_unifaces[i].release();
    }
}

// Set of 2D interfaces using config from config_c_wrapper.h
void mui_create_uniface_multi_2t_f(const char *domain, const char **interfaces, int interface_count ) {

    std::vector<std::string> interface_names;
    size_t interface_element_size = sizeof(interfaces)/sizeof(interfaces[0]);

    if (interface_element_size == 1) {
        char interfaces_jointed[strlen(interfaces[0])];
        strcpy(interfaces_jointed, interfaces[0]);
        char delim[] = " ";

        char * pchar = strtok (interfaces_jointed,delim);
        int interface_names_count = 0;
        for(size_t i=0; i<interface_count; i++) {
            if(pchar == NULL) {
                std::cerr << "MUI Error [mui_f_wrapper_2d.cpp]: Error MUI Fortran C binding: NULL interfaces at the " << i << "th interface." << std::endl;
                std::abort();
            }
            interface_names.push_back(std::string(pchar));
            pchar = strtok (NULL, delim);
            interface_names_count = interface_names_count + 1;
        }

    } else if(interface_element_size > 1) {
        for(size_t i=0; i<interface_count; i++)
            interface_names.push_back(std::string(interfaces[i]));
    } else {
        std::cerr << "MUI Error [mui_f_wrapper_2d.cpp]: Error MUI Fortran C binding: The size of interfaces array must larger or equals to one." << std::endl;
        std::abort();
    }

    auto created_unifaces = mui::create_uniface<mui::mui_f_wrapper_2D>(domain, interface_names);

    mui_uniface_multi_2t = new mui_uniface_2t*[created_unifaces.size()];

      for(size_t i=0; i<created_unifaces.size(); i++) {
        mui_uniface_multi_2t[i] = created_unifaces[i].release();
    }
}

// Access to MUI set of 2D interfaces with float=single and int=int32
mui_uniface_2f* get_mui_uniface_multi_2f_f(int interface_count) {
    return mui_uniface_multi_2f[interface_count-1];
}

// Access to MUI set of 2D interfaces with float=single and int=int64
mui_uniface_2fx* get_mui_uniface_multi_2fx_f(int interface_count) {
    return mui_uniface_multi_2fx[interface_count-1];
}

// Access to MUI set of 2D interfaces with float=double and int=int32
mui_uniface_2d* get_mui_uniface_multi_2d_f(int interface_count) {
    return mui_uniface_multi_2d[interface_count-1];
}

// Access to MUI set of 2D interfaces with float=double and int=int64
mui_uniface_2dx* get_mui_uniface_multi_2dx_f(int interface_count) {
    return mui_uniface_multi_2dx[interface_count-1];
}

// Access to MUI set of 2D interfaces using config from config_f_wrapper.h
mui_uniface_2t* get_mui_uniface_multi_2t_f(int interface_count) {
    return mui_uniface_multi_2t[interface_count-1];
}

/****************************************
 * Destroy MUI interface                 *
 ****************************************/

void mui_destroy_uniface_2f_f(mui_uniface_2f *uniface) {
    delete uniface;
}

void mui_destroy_uniface_2fx_f(mui_uniface_2fx *uniface) {
    delete uniface;
}

void mui_destroy_uniface_2d_f(mui_uniface_2d *uniface) {
    delete uniface;
}

void mui_destroy_uniface_2dx_f(mui_uniface_2dx *uniface) {
    delete uniface;
}

void mui_destroy_uniface_2t_f(mui_uniface_2t *uniface) {
    delete uniface;
}

/******************************************
 * Create 2d spatial samplers              *
 ******************************************/

// Exact sampler
void mui_create_sampler_exact_2f_f(mui_sampler_exact_2f **ret, float *tolerance) {
    *ret = new mui_sampler_exact_2f(*tolerance);
}

void mui_create_sampler_exact_2fx_f(mui_sampler_exact_2fx **ret, float *tolerance) {
    *ret = new mui_sampler_exact_2fx(*tolerance);
}

void mui_create_sampler_exact_2d_f(mui_sampler_exact_2d **ret, double *tolerance) {
    *ret = new mui_sampler_exact_2d(*tolerance);
}

void mui_create_sampler_exact_2dx_f(mui_sampler_exact_2dx** ret, double *tolerance) {
    *ret = new mui_sampler_exact_2dx(*tolerance);
}

void mui_create_sampler_exact_2t_f(mui_sampler_exact_2t** ret, double* tolerance) {
    *ret = new mui_sampler_exact_2t(static_cast<mui::mui_f_wrapper_2D::REAL>(*tolerance));
}

// Gauss sampler
void mui_create_sampler_gauss_2f_f(mui_sampler_gauss_2f **ret, float* r, float* h) {
    *ret = new mui_sampler_gauss_2f(*r, *h);
}

void mui_create_sampler_gauss_2fx_f(mui_sampler_gauss_2fx **ret, float* r, float* h) {
    *ret = new mui_sampler_gauss_2fx(*r, *h);
}

void mui_create_sampler_gauss_2d_f(mui_sampler_gauss_2d **ret, double* r, double* h) {
    *ret = new mui_sampler_gauss_2d(*r, *h);
}

void mui_create_sampler_gauss_2dx_f(mui_sampler_gauss_2dx** ret, double* r, double* h) {
    *ret = new mui_sampler_gauss_2dx(*r, *h);
}

void mui_create_sampler_gauss_2t_f(mui_sampler_gauss_2t** ret, double* r, double* h) {
    *ret = new mui_sampler_gauss_2t(static_cast<mui::mui_f_wrapper_2D::REAL>(*r),
            static_cast<mui::mui_f_wrapper_2D::REAL>(*h));
}

// Moving average sampler
void mui_create_sampler_moving_average_2f_f(mui_sampler_moving_average_2f **ret, float* bbox_1, float* bbox_2) {
    mui::point2f bbox(*bbox_1, *bbox_2);
    *ret = new mui_sampler_moving_average_2f(bbox);
}

void mui_create_sampler_moving_average_2fx_f(mui_sampler_moving_average_2fx **ret, float* bbox_1, float* bbox_2) {
    mui::point2fx bbox(*bbox_1, *bbox_2);
    *ret = new mui_sampler_moving_average_2fx(bbox);
}

void mui_create_sampler_moving_average_2d_f(mui_sampler_moving_average_2d **ret, double* bbox_1, double* bbox_2) {
    mui::point2d bbox(*bbox_1, *bbox_2);
    *ret = new mui_sampler_moving_average_2d(bbox);
}

void mui_create_sampler_moving_average_2dx_f(mui_sampler_moving_average_2dx** ret, double* bbox_1, double* bbox_2) {
    mui::point2dx bbox(*bbox_1, *bbox_2);
    *ret = new mui_sampler_moving_average_2dx(bbox);
}

void mui_create_sampler_moving_average_2t_f(mui_sampler_moving_average_2t** ret, double* bbox_1, double* bbox_2) {
    mui::mui_f_wrapper_2D::point_type bbox(*bbox_1, *bbox_2);
    *ret = new mui_sampler_moving_average_2t(bbox);
}

// Nearest neighbour sampler
void mui_create_sampler_nearest_neighbor_2f_f(mui_sampler_nearest_neighbor_2f **ret) {
    *ret = new mui_sampler_nearest_neighbor_2f();
}

void mui_create_sampler_nearest_neighbor_2fx_f(mui_sampler_nearest_neighbor_2fx **ret) {
    *ret = new mui_sampler_nearest_neighbor_2fx();
}

void mui_create_sampler_nearest_neighbor_2d_f(mui_sampler_nearest_neighbor_2d **ret) {
    *ret = new mui_sampler_nearest_neighbor_2d();
}

void mui_create_sampler_nearest_neighbor_2dx_f(mui_sampler_nearest_neighbor_2dx** ret) {
    *ret = new mui_sampler_nearest_neighbor_2dx();
}

void mui_create_sampler_nearest_neighbor_2t_f(mui_sampler_nearest_neighbor_2t** ret) {
    *ret = new mui_sampler_nearest_neighbor_2t();
}

// Pseudo-linear n^2 interpolation sampler
void mui_create_sampler_pseudo_n2_linear_2f_f(mui_sampler_pseudo_n2_linear_2f **ret, float* r) {
    *ret = new mui_sampler_pseudo_n2_linear_2f(*r);
}

void mui_create_sampler_pseudo_n2_linear_2fx_f(mui_sampler_pseudo_n2_linear_2fx **ret, float* r) {
    *ret = new mui_sampler_pseudo_n2_linear_2fx(*r);
}

void mui_create_sampler_pseudo_n2_linear_2d_f(mui_sampler_pseudo_n2_linear_2d **ret, double* r) {
    *ret = new mui_sampler_pseudo_n2_linear_2d(*r);
}

void mui_create_sampler_pseudo_n2_linear_2dx_f(mui_sampler_pseudo_n2_linear_2dx** ret, double* r) {
    *ret = new mui_sampler_pseudo_n2_linear_2dx(*r);
}

void mui_create_sampler_pseudo_n2_linear_2t_f(mui_sampler_pseudo_n2_linear_2t** ret, double* r) {
    *ret = new mui_sampler_pseudo_n2_linear_2t(static_cast<mui::mui_f_wrapper_2D::REAL>(*r));
}

// Pseudo-linear nearest neighbour interpolation sampler
void mui_create_sampler_pseudo_nearest_neighbor_2f_f(mui_sampler_pseudo_nearest_neighbor_2f **ret, float* h) {
    *ret = new mui_sampler_pseudo_nearest_neighbor_2f(*h);
}

void mui_create_sampler_pseudo_nearest_neighbor_2fx_f(mui_sampler_pseudo_nearest_neighbor_2fx **ret, float* h) {
    *ret = new mui_sampler_pseudo_nearest_neighbor_2fx(*h);
}

void mui_create_sampler_pseudo_nearest_neighbor_2d_f(mui_sampler_pseudo_nearest_neighbor_2d **ret, double* h) {
    *ret = new mui_sampler_pseudo_nearest_neighbor_2d(*h);
}

void mui_create_sampler_pseudo_nearest_neighbor_2dx_f(mui_sampler_pseudo_nearest_neighbor_2dx** ret, double* h) {
    *ret = new mui_sampler_pseudo_nearest_neighbor_2dx(*h);
}

void mui_create_sampler_pseudo_nearest_neighbor_2t_f(mui_sampler_pseudo_nearest_neighbor_2t** ret, double* h) {
    *ret = new mui_sampler_pseudo_nearest_neighbor_2t(static_cast<mui::mui_f_wrapper_2D::REAL>(*h));
}

// Shepard interpolation with a quintic kernel sampler
void mui_create_sampler_shepard_quintic_2f_f(mui_sampler_shepard_quintic_2f **ret, float* r) {
    *ret = new mui_sampler_shepard_quintic_2f(*r);
}

void mui_create_sampler_shepard_quintic_2fx_f(mui_sampler_shepard_quintic_2fx **ret, float* r) {
    *ret = new mui_sampler_shepard_quintic_2fx(*r);
}

void mui_create_sampler_shepard_quintic_2d_f(mui_sampler_shepard_quintic_2d **ret, double* r) {
    *ret = new mui_sampler_shepard_quintic_2d(*r);
}

void mui_create_sampler_shepard_quintic_2dx_f(mui_sampler_shepard_quintic_2dx** ret, double* r) {
    *ret = new mui_sampler_shepard_quintic_2dx(*r);
}

void mui_create_sampler_shepard_quintic_2t_f(mui_sampler_shepard_quintic_2t** ret, double* r) {
    *ret = new mui_sampler_shepard_quintic_2t(static_cast<mui::mui_f_wrapper_2D::REAL>(*r));
}

// SPH derived interpolation method with a quintic spline kernel sampler
void mui_create_sampler_sph_quintic_2f_f(mui_sampler_sph_quintic_2f **ret, float* r) {
    *ret = new mui_sampler_sph_quintic_2f(*r);
}

void mui_create_sampler_sph_quintic_2fx_f(mui_sampler_sph_quintic_2fx **ret, float* r) {
    *ret = new mui_sampler_sph_quintic_2fx(*r);
}

void mui_create_sampler_sph_quintic_2d_f(mui_sampler_sph_quintic_2d **ret, double* r) {
    *ret = new mui_sampler_sph_quintic_2d(*r);
}

void mui_create_sampler_sph_quintic_2dx_f(mui_sampler_sph_quintic_2dx** ret, double* r) {
    *ret = new mui_sampler_sph_quintic_2dx(*r);
}

void mui_create_sampler_sph_quintic_2t_f(mui_sampler_sph_quintic_2t** ret, double* r) {
    *ret = new mui_sampler_sph_quintic_2t(static_cast<mui::mui_f_wrapper_2D::REAL>(*r));
}

// Summation with a quintic kernel sampler
void mui_create_sampler_sum_quintic_2f_f(mui_sampler_sum_quintic_2f **ret, float* r) {
    *ret = new mui_sampler_sum_quintic_2f(*r);
}

void mui_create_sampler_sum_quintic_2fx_f(mui_sampler_sum_quintic_2fx **ret, float* r) {
    *ret = new mui_sampler_sum_quintic_2fx(*r);
}

void mui_create_sampler_sum_quintic_2d_f(mui_sampler_sum_quintic_2d **ret, double* r) {
    *ret = new mui_sampler_sum_quintic_2d(*r);
}

void mui_create_sampler_sum_quintic_2dx_f(mui_sampler_sum_quintic_2dx** ret, double* r) {
    *ret = new mui_sampler_sum_quintic_2dx(*r);
}

void mui_create_sampler_sum_quintic_2t_f(mui_sampler_sum_quintic_2t** ret, double* r) {
    *ret = new mui_sampler_sum_quintic_2t(static_cast<mui::mui_f_wrapper_2D::REAL>(*r));
}

#ifdef USE_RBF
// Radial Basis Function sampler
void mui_create_sampler_rbf_2f_f(mui_sampler_rbf_2f **ret, float* r, float* points_1, float* points_2, int* points_count, int* basis_func,
        int* conservative, int* polynomial, int* smoothFunc, int* readMatrix, int* writeMatrix, const char* file_address,
    float* cutoff, float* cg_solve_tol, int* cg_solve_it, int* pou_size) {
    std::vector<mui::point2f> pts(*points_count);
    for (size_t i = 0; i < *points_count; i++) {
        pts[i][0] = points_1[i];
        pts[i][1] = points_2[i];
    }

    *ret = new mui_sampler_rbf_2f(*r, pts, *basis_func, static_cast<bool>(*conservative), static_cast<bool>(*polynomial),
      static_cast<bool>(*smoothFunc), static_cast<bool>(*readMatrix), static_cast<bool>(*writeMatrix), std::string(file_address),
      *cutoff, *cg_solve_tol, *cg_solve_it, *pou_size);
}

void mui_create_sampler_rbf_2fx_f(mui_sampler_rbf_2fx **ret, float* r, float* points_1, float* points_2, int* points_count, int* basis_func,
        int* conservative, int* polynomial, int* smoothFunc, int* readMatrix, int* writeMatrix, const char* file_address,
    float* cutoff, float* cg_solve_tol, int* cg_solve_it, int* pou_size) {
    std::vector<mui::point2fx> pts(*points_count);
    for (size_t i = 0; i < *points_count; i++) {
        pts[i][0] = points_1[i];
        pts[i][1] = points_2[i];
    }

    *ret = new mui_sampler_rbf_2fx(*r, pts, *basis_func, static_cast<bool>(*conservative), static_cast<bool>(*polynomial),
      static_cast<bool>(*smoothFunc), static_cast<bool>(*readMatrix), static_cast<bool>(*writeMatrix), std::string(file_address),
      *cutoff, *cg_solve_tol, *cg_solve_it, *pou_size);
}

void mui_create_sampler_rbf_2d_f(mui_sampler_rbf_2d **ret, double* r, double* points_1, double* points_2, int* points_count, int* basis_func,
        int* conservative, int* polynomial, int* smoothFunc, int* readMatrix, int* writeMatrix, const char* file_address,
        double* cutoff, double* cg_solve_tol, int* cg_solve_it, int* pou_size) {
    std::vector<mui::point2d> pts(*points_count);
    for (size_t i = 0; i < *points_count; i++) {
        pts[i][0] = points_1[i];
        pts[i][1] = points_2[i];
    }

    *ret = new mui_sampler_rbf_2d(*r, pts, *basis_func, static_cast<bool>(*conservative), static_cast<bool>(*polynomial),
      static_cast<bool>(*smoothFunc), static_cast<bool>(*readMatrix), static_cast<bool>(*writeMatrix), std::string(file_address),
      *cutoff, *cg_solve_tol, *cg_solve_it, *pou_size);
}

void mui_create_sampler_rbf_2dx_f(mui_sampler_rbf_2dx** ret, double* r, double* points_1, double* points_2, int* points_count, int* basis_func,
        int* conservative, int* polynomial, int* smoothFunc, int* readMatrix, int* writeMatrix, const char* file_address,
        double* cutoff, double* cg_solve_tol, int* cg_solve_it, int* pou_size) {
    std::vector<mui::point2dx> pts(*points_count);
    for (size_t i = 0; i < *points_count; i++) {
        pts[i][0] = points_1[i];
        pts[i][1] = points_2[i];
    }

    *ret = new mui_sampler_rbf_2dx(*r, pts, *basis_func, static_cast<bool>(*conservative), static_cast<bool>(*polynomial),
      static_cast<bool>(*smoothFunc), static_cast<bool>(*readMatrix), static_cast<bool>(*writeMatrix), std::string(file_address),
      *cutoff, *cg_solve_tol, *cg_solve_it, *pou_size);
}

void mui_create_sampler_rbf_2t_f(mui_sampler_rbf_2t** ret, double* r, double* points_1, double* points_2, int* points_count, int* basis_func,
        int* conservative, int* polynomial, int* smoothFunc, int* readMatrix, int* writeMatrix, const char* file_address,
        double* cutoff, double* cg_solve_tol, int* cg_solve_it, int* pou_size) {
    std::vector<mui::mui_f_wrapper_2D::point_type> pts(*points_count);
    for (size_t i = 0; i < *points_count; i++) {
        pts[i][0] = static_cast<mui::mui_f_wrapper_2D::REAL>(points_1[i]);
        pts[i][1] = static_cast<mui::mui_f_wrapper_2D::REAL>(points_2[i]);
    }

    *ret = new mui_sampler_rbf_2t(static_cast<mui::mui_f_wrapper_2D::REAL>(*r), pts, *basis_func,
      static_cast<bool>(*conservative), static_cast<bool>(*polynomial), static_cast<bool>(*smoothFunc),
      static_cast<bool>(*readMatrix), static_cast<bool>(*writeMatrix), std::string(file_address),
      static_cast<mui::mui_f_wrapper_2D::REAL>(*cutoff), static_cast<mui::mui_f_wrapper_2D::REAL>(*cg_solve_tol),
      static_cast<mui::mui_f_wrapper_2D::INT>(*cg_solve_it), static_cast<mui::mui_f_wrapper_2D::INT>(*pou_size));
}
#endif

/*******************************************
 * Destroy 2d spatial samplers              *
 *******************************************/

// Exact sampler
void mui_destroy_sampler_exact_2f_f(mui_sampler_exact_2f* sampler) {
    delete sampler;
}

void mui_destroy_sampler_exact_2fx_f(mui_sampler_exact_2fx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_exact_2d_f(mui_sampler_exact_2d* sampler) {
    delete sampler;
}

void mui_destroy_sampler_exact_2dx_f(mui_sampler_exact_2dx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_exact_2t_f(mui_sampler_exact_2t* sampler) {
    delete sampler;
}

// Gauss sampler
void mui_destroy_sampler_gauss_2f_f(mui_sampler_gauss_2f* sampler) {
    delete sampler;
}

void mui_destroy_sampler_gauss_2fx_f(mui_sampler_gauss_2fx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_gauss_2d_f(mui_sampler_gauss_2d* sampler) {
    delete sampler;
}

void mui_destroy_sampler_gauss_2dx_f(mui_sampler_gauss_2dx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_gauss_2t_f(mui_sampler_gauss_2t* sampler) {
    delete sampler;
}

// Moving average sampler
void mui_destroy_sampler_moving_average_2f_f(mui_sampler_moving_average_2f* sampler) {
    delete sampler;
}

void mui_destroy_sampler_moving_average_2fx_f(mui_sampler_moving_average_2fx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_moving_average_2d_f(mui_sampler_moving_average_2d* sampler) {
    delete sampler;
}

void mui_destroy_sampler_moving_average_2dx_f(mui_sampler_moving_average_2dx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_moving_average_2t_f(mui_sampler_moving_average_2t* sampler) {
    delete sampler;
}

// Nearest neighbour sampler
void mui_destroy_sampler_nearest_neighbor_2f_f(mui_sampler_nearest_neighbor_2f* sampler) {
    delete sampler;
}

void mui_destroy_sampler_nearest_neighbor_2fx_f(mui_sampler_nearest_neighbor_2fx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_nearest_neighbor_2d_f(mui_sampler_nearest_neighbor_2d* sampler) {
    delete sampler;
}

void mui_destroy_sampler_nearest_neighbor_2dx_f(mui_sampler_nearest_neighbor_2dx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_nearest_neighbor_2t_f(mui_sampler_nearest_neighbor_2t* sampler) {
    delete sampler;
}

// Pseudo-linear n^2 interpolation sampler
void mui_destroy_sampler_pseudo_nearest2_linear_2f_f(mui_sampler_pseudo_nearest_neighbor_2f* sampler) {
    delete sampler;
}

void mui_destroy_sampler_pseudo_nearest2_linear_2fx_f(mui_sampler_pseudo_nearest_neighbor_2fx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_pseudo_nearest2_linear_2d_f(mui_sampler_pseudo_nearest_neighbor_2d* sampler) {
    delete sampler;
}

void mui_destroy_sampler_pseudo_nearest2_linear_2dx_f(mui_sampler_pseudo_nearest_neighbor_2dx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_pseudo_nearest2_linear_2t_f(mui_sampler_pseudo_nearest_neighbor_2t* sampler) {
    delete sampler;
}

// Pseudo-linear nearest neighbour interpolation sampler
void mui_destroy_sampler_pseudo_nearest_neighbor_2f_f(mui_sampler_pseudo_nearest_neighbor_2f* sampler) {
    delete sampler;
}

void mui_destroy_sampler_pseudo_nearest_neighbor_2fx_f(mui_sampler_pseudo_nearest_neighbor_2fx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_pseudo_nearest_neighbor_2d_f(mui_sampler_pseudo_nearest_neighbor_2d* sampler) {
    delete sampler;
}

void mui_destroy_sampler_pseudo_nearest_neighbor_2dx_f(mui_sampler_pseudo_nearest_neighbor_2dx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_pseudo_nearest_neighbor_2t_f(mui_sampler_pseudo_nearest_neighbor_2t* sampler) {
    delete sampler;
}

// Shepard interpolation with a quintic kernel sampler
void mui_destroy_sampler_shepard_quintic_2f_f(mui_sampler_shepard_quintic_2f* sampler) {
    delete sampler;
}

void mui_destroy_sampler_shepard_quintic_2fx_f(mui_sampler_shepard_quintic_2fx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_shepard_quintic_2d_f(mui_sampler_shepard_quintic_2d* sampler) {
    delete sampler;
}

void mui_destroy_sampler_shepard_quintic_2dx_f(mui_sampler_shepard_quintic_2dx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_shepard_quintic_2t_f(mui_sampler_shepard_quintic_2t* sampler) {
    delete sampler;
}

// SPH derived interpolation method with a quintic spline kernel sampler
void mui_destroy_sampler_sph_quintic_2f_f(mui_sampler_sph_quintic_2f* sampler) {
    delete sampler;
}

void mui_destroy_sampler_sph_quintic_2fx_f(mui_sampler_sph_quintic_2fx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_sph_quintic_2d_f(mui_sampler_sph_quintic_2d* sampler) {
    delete sampler;
}

void mui_destroy_sampler_sph_quintic_2dx_f(mui_sampler_sph_quintic_2dx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_sph_quintic_2t_f(mui_sampler_sph_quintic_2t* sampler) {
    delete sampler;
}

// Summation with a quintic kernel sampler
void mui_destroy_sampler_sum_quintic_2f_f(mui_sampler_sum_quintic_2f* sampler) {
    delete sampler;
}

void mui_destroy_sampler_sum_quintic_2fx_f(mui_sampler_sum_quintic_2fx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_sum_quintic_2d_f(mui_sampler_sum_quintic_2d* sampler) {
    delete sampler;
}

void mui_destroy_sampler_sum_quintic_2dx_f(mui_sampler_sum_quintic_2dx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_sum_quintic_2t_f(mui_sampler_sum_quintic_2t* sampler) {
    delete sampler;
}

#ifdef USE_RBF
void mui_destroy_sampler_rbf_2f_f(mui_sampler_rbf_2f* sampler) {
    delete sampler;
}

void mui_destroy_sampler_rbf_2fx_f(mui_sampler_rbf_2fx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_rbf_2d_f(mui_sampler_rbf_2d* sampler) {
    delete sampler;
}

void mui_destroy_sampler_rbf_2dx_f(mui_sampler_rbf_2dx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_rbf_2t_f(mui_sampler_rbf_2t* sampler) {
    delete sampler;
}
#endif

/*******************************************
 * Create temporal samplers                 *
 *******************************************/

// Exact temporal sampler
void mui_create_chrono_sampler_exact_2f_f(mui_chrono_sampler_exact_2f **ret, float* tolerance) {
    *ret = new mui_chrono_sampler_exact_2f(*tolerance);
}

void mui_create_chrono_sampler_exact_2fx_f(mui_chrono_sampler_exact_2fx **ret, float* tolerance) {
    *ret = new mui_chrono_sampler_exact_2fx(*tolerance);
}

void mui_create_chrono_sampler_exact_2d_f(mui_chrono_sampler_exact_2d **ret, double* tolerance) {
    *ret = new mui_chrono_sampler_exact_2d(*tolerance);
}

void mui_create_chrono_sampler_exact_2dx_f(mui_chrono_sampler_exact_2dx** ret, double* tolerance) {
    *ret = new mui_chrono_sampler_exact_2dx(*tolerance);
}

void mui_create_chrono_sampler_exact_2t_f(mui_chrono_sampler_exact_2t** ret, double* tolerance) {
    *ret = new mui_chrono_sampler_exact_2t(static_cast<mui::mui_f_wrapper_2D::REAL>(*tolerance));
}

// Gauss temporal sampler
void mui_create_chrono_sampler_gauss_2f_f(mui_chrono_sampler_gauss_2f **ret, float* cutoff, float* sigma) {
    *ret = new mui_chrono_sampler_gauss_2f(*cutoff, *sigma);
}

void mui_create_chrono_sampler_gauss_2fx_f(mui_chrono_sampler_gauss_2fx **ret, float* cutoff, float* sigma) {
    *ret = new mui_chrono_sampler_gauss_2fx(*cutoff, *sigma);
}

void mui_create_chrono_sampler_gauss_2d_f(mui_chrono_sampler_gauss_2d **ret, double* cutoff, double* sigma) {
    *ret = new mui_chrono_sampler_gauss_2d(*cutoff, *sigma);
}

void mui_create_chrono_sampler_gauss_2dx_f(mui_chrono_sampler_gauss_2dx** ret, double* cutoff, double* sigma) {
    *ret = new mui_chrono_sampler_gauss_2dx(*cutoff, *sigma);
}

void mui_create_chrono_sampler_gauss_2t_f(mui_chrono_sampler_gauss_2t** ret, double* cutoff, double* sigma) {
    *ret = new mui_chrono_sampler_gauss_2t(static_cast<mui::mui_f_wrapper_2D::time_type>(*cutoff),
            static_cast<mui::mui_f_wrapper_2D::REAL>(*sigma));
}

// Mean temporal sampler
void mui_create_chrono_sampler_mean_2f_f(mui_chrono_sampler_mean_2f **ret, float* lower, float* upper) {
    *ret = new mui_chrono_sampler_mean_2f(*lower, *upper);
}

void mui_create_chrono_sampler_mean_2fx_f(mui_chrono_sampler_mean_2fx **ret, float* lower, float* upper) {
    *ret = new mui_chrono_sampler_mean_2fx(*lower, *upper);
}

void mui_create_chrono_sampler_mean_2d_f(mui_chrono_sampler_mean_2d **ret, double* lower, double* upper) {
    *ret = new mui_chrono_sampler_mean_2d(*lower, *upper);
}

void mui_create_chrono_sampler_mean_2dx_f(mui_chrono_sampler_mean_2dx** ret, double* lower, double* upper) {
    *ret = new mui_chrono_sampler_mean_2dx(*lower, *upper);
}

void mui_create_chrono_sampler_mean_2t_f(mui_chrono_sampler_mean_2t** ret, double* lower, double* upper) {
    *ret = new mui_chrono_sampler_mean_2t(static_cast<mui::mui_f_wrapper_2D::time_type>(*lower),
            static_cast<mui::mui_f_wrapper_2D::time_type>(*upper));
}

// Sum temporal sampler
void mui_create_chrono_sampler_sum_2f_f(mui_chrono_sampler_sum_2f **ret, float* lower, float* upper) {
    *ret = new mui_chrono_sampler_sum_2f(*lower, *upper);
}

void mui_create_chrono_sampler_sum_2fx_f(mui_chrono_sampler_sum_2fx **ret, float* lower, float* upper) {
    *ret = new mui_chrono_sampler_sum_2fx(*lower, *upper);
}

void mui_create_chrono_sampler_sum_2d_f(mui_chrono_sampler_sum_2d **ret, double* lower, double* upper) {
    *ret = new mui_chrono_sampler_sum_2d(*lower, *upper);
}

void mui_create_chrono_sampler_sum_2dx_f(mui_chrono_sampler_sum_2dx** ret, double* lower, double* upper) {
    *ret = new mui_chrono_sampler_sum_2dx(*lower, *upper);
}

void mui_create_chrono_sampler_sum_2t_f(mui_chrono_sampler_sum_2t** ret, double* lower, double* upper) {
    *ret = new mui_chrono_sampler_sum_2t(static_cast<mui::mui_f_wrapper_2D::time_type>(*lower),
            static_cast<mui::mui_f_wrapper_2D::time_type>(*upper));
}

/*******************************************
 * Destroy temporal samplers                *
 *******************************************/

// Exact temporal sampler
void mui_destroy_chrono_sampler_exact_2f_f(mui_chrono_sampler_exact_2f* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_exact_2fx_f(mui_chrono_sampler_exact_2fx* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_exact_2d_f(mui_chrono_sampler_exact_2d* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_exact_2dx_f(mui_chrono_sampler_exact_2dx* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_exact_2t_f(mui_chrono_sampler_exact_2t* sampler) {
    delete sampler;
}

// Gauss temporal sampler
void mui_destroy_chrono_sampler_gauss_2f_f(mui_chrono_sampler_gauss_2f* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_gauss_2fx_f(mui_chrono_sampler_gauss_2fx* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_gauss_2d_f(mui_chrono_sampler_gauss_2d* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_gauss_2dx_f(mui_chrono_sampler_gauss_2dx* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_gauss_2t_f(mui_chrono_sampler_gauss_2t* sampler) {
    delete sampler;
}

// Mean temporal sampler
void mui_destroy_chrono_sampler_mean_2f_f(mui_chrono_sampler_mean_2f* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_mean_2fx_f(mui_chrono_sampler_mean_2fx* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_mean_2d_f(mui_chrono_sampler_mean_2d* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_mean_2dx_f(mui_chrono_sampler_mean_2dx* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_mean_2t_f(mui_chrono_sampler_mean_2t* sampler) {
    delete sampler;
}

// Sum temporal sampler
void mui_destroy_chrono_sampler_sum_2f_f(mui_chrono_sampler_sum_2f* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_sum_2fx_f(mui_chrono_sampler_sum_2fx* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_sum_2d_f(mui_chrono_sampler_sum_2d* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_sum_2dx_f(mui_chrono_sampler_sum_2dx* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_sum_2t_f(mui_chrono_sampler_sum_2t* sampler) {
    delete sampler;
}

/******************************************
 * MUI functions for data push             *
 ******************************************/

// Standard push functions
void mui_push_2f_f(mui_uniface_2f *uniface, const char* attr, float* point_1, float* point_2, float* value) {
    uniface->push(std::string(attr), mui::point2f(*point_1,*point_2), *value);
}

void mui_push_2fx_f(mui_uniface_2fx *uniface, const char* attr, float* point_1, float* point_2, float* value) {
    uniface->push(std::string(attr), mui::point2fx(*point_1,*point_2), *value);
}

void mui_push_2d_f(mui_uniface_2d *uniface, const char* attr, double* point_1, double* point_2, double* value) {
    uniface->push(std::string(attr), mui::point2d(*point_1,*point_2), *value);
}

void mui_push_2dx_f(mui_uniface_2dx *uniface, const char* attr, double* point_1, double* point_2, double* value) {
    uniface->push(std::string(attr), mui::point2dx(*point_1,*point_2), *value);
}

void mui_push_2t_f(mui_uniface_2t *uniface, const char* attr, double* point_1, double* point_2, double* value) {
    mui::mui_f_wrapper_2D::point_type push_point(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    uniface->push(std::string(attr), push_point, static_cast<mui::mui_f_wrapper_2D::REAL>(*value));
}

// Single parameter push functions
void mui_push_2f_param_f(mui_uniface_2f *uniface, const char* attr, float* value) {
    uniface->push(std::string(attr), *value);
}

void mui_push_2fx_param_f(mui_uniface_2fx *uniface, const char* attr, float* value) {
    uniface->push(std::string(attr), *value);
}

void mui_push_2d_param_f(mui_uniface_2d *uniface, const char* attr, double* value) {
    uniface->push(std::string(attr), *value);
}

void mui_push_2dx_param_f(mui_uniface_2dx *uniface, const char* attr, double* value) {
    uniface->push(std::string(attr), *value);
}

void mui_push_2t_param_f(mui_uniface_2t *uniface, const char* attr, double* value) {
    uniface->push(std::string(attr), static_cast<mui::mui_f_wrapper_2D::REAL>(*value));
}

/******************************************
 * MUI functions for data commit           *
 ******************************************/

// Commit using one time value
void mui_commit_2f_f(mui_uniface_2f *uniface, float* t) {
    uniface->commit(*t);
}

void mui_commit_2fx_f(mui_uniface_2fx *uniface, float* t) {
    uniface->commit(*t);
}

void mui_commit_2d_f(mui_uniface_2d *uniface, double* t) {
    uniface->commit(*t);
}

void mui_commit_2dx_f(mui_uniface_2dx *uniface, double* t) {
    uniface->commit(*t);
}

void mui_commit_2t_f(mui_uniface_2t *uniface, double* t) {
    uniface->commit(static_cast<mui::mui_f_wrapper_2D::time_type>(*t));
}

// Commit using two time values
void mui_commit_2f_pair_f(mui_uniface_2f *uniface, float* t_1, float* t_2) {
    uniface->commit(*t_1, *t_2);
}

void mui_commit_2fx_pair_f(mui_uniface_2fx *uniface, float* t_1, float* t_2) {
    uniface->commit(*t_1, *t_2);
}

void mui_commit_2d_pair_f(mui_uniface_2d *uniface, double* t_1, double* t_2) {
    uniface->commit(*t_1, *t_2);
}

void mui_commit_2dx_pair_f(mui_uniface_2dx *uniface, double* t_1, double* t_2) {
    uniface->commit(*t_1, *t_2);
}

void mui_commit_2t_pair_f(mui_uniface_2t *uniface, double* t_1, double* t_2) {
    uniface->commit(static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1),
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2));
}

/******************************************
 * MUI functions for data forecast         *
 ******************************************/

// Forecast using one time value
void mui_forecast_2f_f(mui_uniface_2f *uniface, float* t) {
    uniface->forecast(*t);
}

void mui_forecast_2fx_f(mui_uniface_2fx *uniface, float* t) {
    uniface->forecast(*t);
}

void mui_forecast_2d_f(mui_uniface_2d *uniface, double* t) {
    uniface->forecast(*t);
}

void mui_forecast_2dx_f(mui_uniface_2dx *uniface, double* t) {
    uniface->forecast(*t);
}

void mui_forecast_2t_f(mui_uniface_2t *uniface, double* t) {
    uniface->forecast(static_cast<mui::mui_f_wrapper_2D::time_type>(*t));
}

// Forecast using two time values
void mui_forecast_2f_pair_f(mui_uniface_2f *uniface, float* t_1, float* t_2) {
    uniface->forecast(*t_1, *t_2);
}

void mui_forecast_2fx_pair_f(mui_uniface_2fx *uniface, float* t_1, float* t_2) {
    uniface->forecast(*t_1, *t_2);
}

void mui_forecast_2d_pair_f(mui_uniface_2d *uniface, double* t_1, double* t_2) {
    uniface->forecast(*t_1, *t_2);
}

void mui_forecast_2dx_pair_f(mui_uniface_2dx *uniface, double* t_1, double* t_2) {
    uniface->forecast(*t_1, *t_2);
}

void mui_forecast_2t_pair_f(mui_uniface_2t *uniface, double* t_1, double* t_2) {
    uniface->forecast(static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1),
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2));
}

/*******************************************************
 * MUI functions for 2d data fetch using one time value *
 ********************************************************/

// Spatial sampler: exact; temporal sampler: exact
void mui_fetch_exact_exact_2f_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t,
        mui_sampler_exact_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_exact_2fx_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2, float* t,
        mui_sampler_exact_2fx *spatial_sampler, mui_chrono_sampler_exact_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_exact_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_exact_2d *spatial_sampler, mui_chrono_sampler_exact_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_exact_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_exact_2dx *spatial_sampler, mui_chrono_sampler_exact_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_exact_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_exact_2t *spatial_sampler, mui_chrono_sampler_exact_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: exact; temporal sampler: gauss
void mui_fetch_exact_gauss_2f_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t,
        mui_sampler_exact_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_gauss_2fx_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2, float* t,
        mui_sampler_exact_2fx *spatial_sampler, mui_chrono_sampler_gauss_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_gauss_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_exact_2d *spatial_sampler, mui_chrono_sampler_gauss_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_gauss_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_exact_2dx *spatial_sampler, mui_chrono_sampler_gauss_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_gauss_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_exact_2t *spatial_sampler, mui_chrono_sampler_gauss_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: exact; temporal sampler: mean
void mui_fetch_exact_mean_2f_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t,
        mui_sampler_exact_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_mean_2fx_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2, float* t,
        mui_sampler_exact_2fx *spatial_sampler, mui_chrono_sampler_mean_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_mean_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_exact_2d *spatial_sampler, mui_chrono_sampler_mean_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_mean_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_exact_2dx *spatial_sampler, mui_chrono_sampler_mean_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_mean_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_exact_2t *spatial_sampler, mui_chrono_sampler_mean_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: exact; temporal sampler: sum
void mui_fetch_exact_sum_2f_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t,
        mui_sampler_exact_2f *spatial_sampler, mui_chrono_sampler_sum_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_sum_2fx_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2, float* t,
        mui_sampler_exact_2fx *spatial_sampler, mui_chrono_sampler_sum_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_sum_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_exact_2d *spatial_sampler, mui_chrono_sampler_sum_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_sum_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_exact_2dx *spatial_sampler, mui_chrono_sampler_sum_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_sum_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_exact_2t *spatial_sampler, mui_chrono_sampler_sum_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: gauss; temporal sampler: exact
void mui_fetch_gauss_exact_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_gauss_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_exact_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_gauss_2fx *spatial_sampler, mui_chrono_sampler_exact_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_exact_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_gauss_2d *spatial_sampler, mui_chrono_sampler_exact_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_exact_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_gauss_2dx *spatial_sampler, mui_chrono_sampler_exact_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_exact_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_gauss_2t *spatial_sampler, mui_chrono_sampler_exact_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: gauss; temporal sampler: gauss
void mui_fetch_gauss_gauss_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_gauss_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_gauss_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_gauss_2fx *spatial_sampler, mui_chrono_sampler_gauss_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_gauss_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_gauss_2d *spatial_sampler, mui_chrono_sampler_gauss_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_gauss_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_gauss_2dx *spatial_sampler, mui_chrono_sampler_gauss_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_gauss_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_gauss_2t *spatial_sampler, mui_chrono_sampler_gauss_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: gauss; temporal sampler: mean
void mui_fetch_gauss_mean_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_gauss_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_mean_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_gauss_2fx *spatial_sampler, mui_chrono_sampler_mean_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_mean_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_gauss_2d *spatial_sampler, mui_chrono_sampler_mean_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_mean_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_gauss_2dx *spatial_sampler, mui_chrono_sampler_mean_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_mean_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_gauss_2t *spatial_sampler, mui_chrono_sampler_mean_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: gauss; temporal sampler: sum
void mui_fetch_gauss_sum_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_gauss_2f *spatial_sampler, mui_chrono_sampler_sum_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_sum_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_gauss_2fx *spatial_sampler, mui_chrono_sampler_sum_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_sum_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_gauss_2d *spatial_sampler, mui_chrono_sampler_sum_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_sum_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_gauss_2dx *spatial_sampler, mui_chrono_sampler_sum_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_sum_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_gauss_2t *spatial_sampler, mui_chrono_sampler_sum_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: moving average; temporal sampler: exact
void mui_fetch_moving_average_exact_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_moving_average_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_exact_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_moving_average_2fx *spatial_sampler, mui_chrono_sampler_exact_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_exact_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_moving_average_2d *spatial_sampler, mui_chrono_sampler_exact_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_exact_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_moving_average_2dx *spatial_sampler, mui_chrono_sampler_exact_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_exact_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_moving_average_2t *spatial_sampler, mui_chrono_sampler_exact_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: moving average; temporal sampler: gauss
void mui_fetch_moving_average_gauss_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_moving_average_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_gauss_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_moving_average_2fx *spatial_sampler, mui_chrono_sampler_gauss_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_gauss_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_moving_average_2d *spatial_sampler, mui_chrono_sampler_gauss_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_gauss_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_moving_average_2dx *spatial_sampler, mui_chrono_sampler_gauss_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_gauss_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_moving_average_2t *spatial_sampler, mui_chrono_sampler_gauss_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: moving average; temporal sampler: mean
void mui_fetch_moving_average_mean_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_moving_average_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_mean_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_moving_average_2fx *spatial_sampler, mui_chrono_sampler_mean_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_mean_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_moving_average_2d *spatial_sampler, mui_chrono_sampler_mean_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_mean_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_moving_average_2dx *spatial_sampler, mui_chrono_sampler_mean_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_mean_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_moving_average_2t *spatial_sampler, mui_chrono_sampler_mean_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: moving average; temporal sampler: sum
void mui_fetch_moving_average_sum_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_moving_average_2f *spatial_sampler, mui_chrono_sampler_sum_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_sum_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_moving_average_2fx *spatial_sampler, mui_chrono_sampler_sum_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_sum_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_moving_average_2d *spatial_sampler, mui_chrono_sampler_sum_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_sum_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_moving_average_2dx *spatial_sampler, mui_chrono_sampler_sum_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_sum_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_moving_average_2t *spatial_sampler, mui_chrono_sampler_sum_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: nearest neighbor; temporal sampler: exact
void mui_fetch_nearest_neighbor_exact_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_nearest_neighbor_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_exact_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_nearest_neighbor_2fx *spatial_sampler, mui_chrono_sampler_exact_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_exact_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_nearest_neighbor_2d *spatial_sampler, mui_chrono_sampler_exact_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_exact_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_nearest_neighbor_2dx *spatial_sampler, mui_chrono_sampler_exact_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_exact_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_nearest_neighbor_2t *spatial_sampler, mui_chrono_sampler_exact_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: nearest neighbor; temporal sampler: gauss
void mui_fetch_nearest_neighbor_gauss_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_nearest_neighbor_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_gauss_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_nearest_neighbor_2fx *spatial_sampler, mui_chrono_sampler_gauss_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_gauss_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_nearest_neighbor_2d *spatial_sampler, mui_chrono_sampler_gauss_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_gauss_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_nearest_neighbor_2dx *spatial_sampler, mui_chrono_sampler_gauss_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_gauss_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_nearest_neighbor_2t *spatial_sampler, mui_chrono_sampler_gauss_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: nearest neighbor; temporal sampler: mean
void mui_fetch_nearest_neighbor_mean_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_nearest_neighbor_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_mean_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_nearest_neighbor_2fx *spatial_sampler, mui_chrono_sampler_mean_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_mean_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_nearest_neighbor_2d *spatial_sampler, mui_chrono_sampler_mean_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_mean_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_nearest_neighbor_2dx *spatial_sampler, mui_chrono_sampler_mean_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_mean_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_nearest_neighbor_2t *spatial_sampler, mui_chrono_sampler_mean_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: nearest neighbor; temporal sampler: sum
void mui_fetch_nearest_neighbor_sum_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_nearest_neighbor_2f *spatial_sampler, mui_chrono_sampler_sum_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_sum_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_nearest_neighbor_2fx *spatial_sampler, mui_chrono_sampler_sum_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_sum_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_nearest_neighbor_2d *spatial_sampler, mui_chrono_sampler_sum_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_sum_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_nearest_neighbor_2dx *spatial_sampler, mui_chrono_sampler_sum_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_sum_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_nearest_neighbor_2t *spatial_sampler, mui_chrono_sampler_sum_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: pseudo nearest neighbor; temporal sampler: exact
void mui_fetch_pseudo_nearest_neighbor_exact_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_pseudo_nearest_neighbor_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_exact_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2,
        float* t, mui_sampler_pseudo_nearest_neighbor_2fx *spatial_sampler,
        mui_chrono_sampler_exact_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_exact_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2,
        double* t, mui_sampler_pseudo_nearest_neighbor_2d *spatial_sampler,
        mui_chrono_sampler_exact_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_exact_2dx(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2,
        double* t, mui_sampler_pseudo_nearest_neighbor_2dx *spatial_sampler,
        mui_chrono_sampler_exact_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_exact_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2,
        double* t, mui_sampler_pseudo_nearest_neighbor_2t *spatial_sampler,
        mui_chrono_sampler_exact_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: nearest neighbor; temporal sampler: gauss
void mui_fetch_pseudo_nearest_neighbor_gauss_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_pseudo_nearest_neighbor_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_gauss_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2,
        float* t, mui_sampler_pseudo_nearest_neighbor_2fx *spatial_sampler,
        mui_chrono_sampler_gauss_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_gauss_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2,
        double* t, mui_sampler_pseudo_nearest_neighbor_2d *spatial_sampler,
        mui_chrono_sampler_gauss_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_gauss_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2,
        double* t, mui_sampler_pseudo_nearest_neighbor_2dx *spatial_sampler,
        mui_chrono_sampler_gauss_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_gauss_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2,
        double* t, mui_sampler_pseudo_nearest_neighbor_2t *spatial_sampler,
        mui_chrono_sampler_gauss_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: nearest neighbor; temporal sampler: mean
void mui_fetch_pseudo_nearest_neighbor_mean_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_pseudo_nearest_neighbor_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_mean_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2,
        float* t, mui_sampler_pseudo_nearest_neighbor_2fx *spatial_sampler,
        mui_chrono_sampler_mean_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_mean_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2,
        double* t, mui_sampler_pseudo_nearest_neighbor_2d *spatial_sampler,
        mui_chrono_sampler_mean_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_mean_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2,
        double* t, mui_sampler_pseudo_nearest_neighbor_2dx *spatial_sampler,
        mui_chrono_sampler_mean_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_mean_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2,
        double* t, mui_sampler_pseudo_nearest_neighbor_2t *spatial_sampler,
        mui_chrono_sampler_mean_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: nearest neighbor; temporal sampler: sum
void mui_fetch_pseudo_nearest_neighbor_sum_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_pseudo_nearest_neighbor_2f *spatial_sampler, mui_chrono_sampler_sum_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_sum_2fx(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2,
        float* t, mui_sampler_pseudo_nearest_neighbor_2fx *spatial_sampler,
        mui_chrono_sampler_sum_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_sum_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_pseudo_nearest_neighbor_2d *spatial_sampler, mui_chrono_sampler_sum_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_sum_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2,
        double* t, mui_sampler_pseudo_nearest_neighbor_2dx *spatial_sampler,
        mui_chrono_sampler_sum_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_sum_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_pseudo_nearest_neighbor_2t *spatial_sampler, mui_chrono_sampler_sum_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: shepard quintic; temporal sampler: exact
void mui_fetch_shepard_quintic_exact_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_shepard_quintic_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_exact_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_shepard_quintic_2fx *spatial_sampler, mui_chrono_sampler_exact_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_exact_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_shepard_quintic_2d *spatial_sampler, mui_chrono_sampler_exact_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_exact_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_shepard_quintic_2dx *spatial_sampler, mui_chrono_sampler_exact_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_exact_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_shepard_quintic_2t *spatial_sampler, mui_chrono_sampler_exact_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: shepard quintic; temporal sampler: gauss
void mui_fetch_shepard_quintic_gauss_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_shepard_quintic_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_gauss_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_shepard_quintic_2fx *spatial_sampler, mui_chrono_sampler_gauss_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_gauss_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_shepard_quintic_2d *spatial_sampler, mui_chrono_sampler_gauss_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_gauss_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_shepard_quintic_2dx *spatial_sampler, mui_chrono_sampler_gauss_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_gauss_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_shepard_quintic_2t *spatial_sampler, mui_chrono_sampler_gauss_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: shepard quintic; temporal sampler: mean
void mui_fetch_shepard_quintic_mean_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_shepard_quintic_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_mean_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_shepard_quintic_2fx *spatial_sampler, mui_chrono_sampler_mean_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_mean_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_shepard_quintic_2d *spatial_sampler, mui_chrono_sampler_mean_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_mean_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_shepard_quintic_2dx *spatial_sampler, mui_chrono_sampler_mean_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_mean_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_shepard_quintic_2t *spatial_sampler, mui_chrono_sampler_mean_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: shepard quintic; temporal sampler: sum
void mui_fetch_shepard_quintic_sum_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_shepard_quintic_2f *spatial_sampler, mui_chrono_sampler_sum_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_sum_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_shepard_quintic_2fx *spatial_sampler, mui_chrono_sampler_sum_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_sum_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_shepard_quintic_2d *spatial_sampler, mui_chrono_sampler_sum_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_sum_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_shepard_quintic_2dx *spatial_sampler, mui_chrono_sampler_sum_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_sum_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_shepard_quintic_2t *spatial_sampler, mui_chrono_sampler_sum_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: sph-derived quintic; temporal sampler: exact
void mui_fetch_sph_quintic_exact_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_sph_quintic_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_exact_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_sph_quintic_2fx *spatial_sampler, mui_chrono_sampler_exact_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_exact_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_sph_quintic_2d *spatial_sampler, mui_chrono_sampler_exact_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_exact_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_sph_quintic_2dx *spatial_sampler, mui_chrono_sampler_exact_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_exact_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_sph_quintic_2t *spatial_sampler, mui_chrono_sampler_exact_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: sph-derived quintic; temporal sampler: gauss
void mui_fetch_sph_quintic_gauss_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_sph_quintic_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_gauss_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_sph_quintic_2fx *spatial_sampler, mui_chrono_sampler_gauss_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_gauss_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_sph_quintic_2d *spatial_sampler, mui_chrono_sampler_gauss_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_gauss_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_sph_quintic_2dx *spatial_sampler, mui_chrono_sampler_gauss_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_gauss_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_sph_quintic_2t *spatial_sampler, mui_chrono_sampler_gauss_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: sph-derived quintic; temporal sampler: mean
void mui_fetch_sph_quintic_mean_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_sph_quintic_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_mean_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_sph_quintic_2fx *spatial_sampler, mui_chrono_sampler_mean_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_mean_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_sph_quintic_2d *spatial_sampler, mui_chrono_sampler_mean_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_mean_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_sph_quintic_2dx *spatial_sampler, mui_chrono_sampler_mean_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_mean_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_sph_quintic_2t *spatial_sampler, mui_chrono_sampler_mean_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: sph-derived quintic; temporal sampler: sum
void mui_fetch_sph_quintic_sum_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_sph_quintic_2f *spatial_sampler, mui_chrono_sampler_sum_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_sum_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_sph_quintic_2fx *spatial_sampler, mui_chrono_sampler_sum_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_sum_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_sph_quintic_2d *spatial_sampler, mui_chrono_sampler_sum_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_sum_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_sph_quintic_2dx *spatial_sampler, mui_chrono_sampler_sum_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_sum_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_sph_quintic_2t *spatial_sampler, mui_chrono_sampler_sum_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: summation quintic; temporal sampler: exact
void mui_fetch_sum_quintic_exact_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_sum_quintic_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_exact_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_sum_quintic_2fx *spatial_sampler, mui_chrono_sampler_exact_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_exact_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_sum_quintic_2d *spatial_sampler, mui_chrono_sampler_exact_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_exact_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_sum_quintic_2dx *spatial_sampler, mui_chrono_sampler_exact_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_exact_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_sum_quintic_2t *spatial_sampler, mui_chrono_sampler_exact_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: summation quintic; temporal sampler: gauss
void mui_fetch_sum_quintic_gauss_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_sum_quintic_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_gauss_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_sum_quintic_2fx *spatial_sampler, mui_chrono_sampler_gauss_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_gauss_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_sum_quintic_2d *spatial_sampler, mui_chrono_sampler_gauss_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_gauss_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_sum_quintic_2dx *spatial_sampler, mui_chrono_sampler_gauss_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_gauss_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_sum_quintic_2t *spatial_sampler, mui_chrono_sampler_gauss_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: summation quintic; temporal sampler: mean
void mui_fetch_sum_quintic_mean_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_sum_quintic_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_mean_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_sum_quintic_2fx *spatial_sampler, mui_chrono_sampler_mean_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_mean_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_sum_quintic_2d *spatial_sampler, mui_chrono_sampler_mean_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_mean_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_sum_quintic_2dx *spatial_sampler, mui_chrono_sampler_mean_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_mean_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_sum_quintic_2t *spatial_sampler, mui_chrono_sampler_mean_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: summation quintic; temporal sampler: sum
void mui_fetch_sum_quintic_sum_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_sum_quintic_2f *spatial_sampler, mui_chrono_sampler_sum_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_sum_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_sum_quintic_2fx *spatial_sampler, mui_chrono_sampler_sum_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_sum_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_sum_quintic_2d *spatial_sampler, mui_chrono_sampler_sum_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_sum_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_sum_quintic_2dx *spatial_sampler, mui_chrono_sampler_sum_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_sum_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_sum_quintic_2t *spatial_sampler, mui_chrono_sampler_sum_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

#ifdef USE_RBF
// Spatial sampler: radial basis function; temporal sampler: exact
void mui_fetch_rbf_exact_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_rbf_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_exact_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_rbf_2fx *spatial_sampler, mui_chrono_sampler_exact_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_exact_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_rbf_2d *spatial_sampler, mui_chrono_sampler_exact_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_exact_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_rbf_2dx *spatial_sampler, mui_chrono_sampler_exact_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_exact_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_rbf_2t *spatial_sampler, mui_chrono_sampler_exact_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: radial basis function; temporal sampler: gauss
void mui_fetch_rbf_gauss_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_rbf_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_gauss_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_rbf_2fx *spatial_sampler, mui_chrono_sampler_gauss_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_gauss_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_rbf_2d *spatial_sampler, mui_chrono_sampler_gauss_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_gauss_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_rbf_2dx *spatial_sampler, mui_chrono_sampler_gauss_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_gauss_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_rbf_2t *spatial_sampler, mui_chrono_sampler_gauss_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: radial basis function; temporal sampler: mean
void mui_fetch_rbf_mean_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_rbf_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_mean_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_rbf_2fx *spatial_sampler, mui_chrono_sampler_mean_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_mean_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_rbf_2d *spatial_sampler, mui_chrono_sampler_mean_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_mean_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_rbf_2dx *spatial_sampler, mui_chrono_sampler_mean_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_mean_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_rbf_2t *spatial_sampler, mui_chrono_sampler_mean_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: radial basis function; temporal sampler: sum
void mui_fetch_rbf_sum_2f_f(mui_uniface_2f *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_rbf_2f *spatial_sampler, mui_chrono_sampler_sum_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_sum_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *point_1, float *point_2, float* t,
        mui_sampler_rbf_2fx *spatial_sampler, mui_chrono_sampler_sum_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_sum_2d_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_rbf_2d *spatial_sampler, mui_chrono_sampler_sum_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_sum_2dx_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_rbf_2dx *spatial_sampler, mui_chrono_sampler_sum_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_sum_2t_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t,
        mui_sampler_rbf_2t *spatial_sampler, mui_chrono_sampler_sum_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}
#endif

/********************************************************
 * MUI functions for 2d data fetch using two time values *
 *********************************************************/

// Spatial sampler: exact; temporal sampler: exact
void mui_fetch_exact_exact_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t_1, float* t_2,
        mui_sampler_exact_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_exact_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_exact_2fx *spatial_sampler, mui_chrono_sampler_exact_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_exact_exact_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_exact_2d *spatial_sampler, mui_chrono_sampler_exact_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_exact_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_exact_2dx *spatial_sampler, mui_chrono_sampler_exact_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_exact_exact_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_exact_2t *spatial_sampler, mui_chrono_sampler_exact_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: exact; temporal sampler: gauss
void mui_fetch_exact_gauss_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t_1, float* t_2,
        mui_sampler_exact_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_gauss_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_exact_2fx *spatial_sampler, mui_chrono_sampler_gauss_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_exact_gauss_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_exact_2d *spatial_sampler, mui_chrono_sampler_gauss_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_gauss_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_exact_2dx *spatial_sampler, mui_chrono_sampler_gauss_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_exact_gauss_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_exact_2t *spatial_sampler, mui_chrono_sampler_gauss_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: exact; temporal sampler: mean
void mui_fetch_exact_mean_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t_1, float* t_2,
        mui_sampler_exact_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_mean_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_exact_2fx *spatial_sampler, mui_chrono_sampler_mean_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_exact_mean_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_exact_2d *spatial_sampler, mui_chrono_sampler_mean_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_mean_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_exact_2dx *spatial_sampler, mui_chrono_sampler_mean_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_exact_mean_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_exact_2t *spatial_sampler, mui_chrono_sampler_mean_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: exact; temporal sampler: sum
void mui_fetch_exact_sum_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t_1, float* t_2,
        mui_sampler_exact_2f *spatial_sampler, mui_chrono_sampler_sum_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_sum_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_exact_2fx *spatial_sampler, mui_chrono_sampler_sum_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_exact_sum_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_exact_2d *spatial_sampler, mui_chrono_sampler_sum_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_sum_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_exact_2dx *spatial_sampler, mui_chrono_sampler_sum_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_exact_sum_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_exact_2t *spatial_sampler, mui_chrono_sampler_sum_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: gauss; temporal sampler: exact
void mui_fetch_gauss_exact_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t_1, float* t_2,
        mui_sampler_gauss_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_exact_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_gauss_2fx *spatial_sampler, mui_chrono_sampler_exact_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_gauss_exact_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_gauss_2d *spatial_sampler, mui_chrono_sampler_exact_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_exact_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_gauss_2dx *spatial_sampler, mui_chrono_sampler_exact_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_gauss_exact_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_gauss_2t *spatial_sampler, mui_chrono_sampler_exact_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: gauss; temporal sampler: gauss
void mui_fetch_gauss_gauss_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t_1, float* t_2,
        mui_sampler_gauss_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_gauss_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_gauss_2fx *spatial_sampler, mui_chrono_sampler_gauss_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_gauss_gauss_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_gauss_2d *spatial_sampler, mui_chrono_sampler_gauss_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_gauss_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_gauss_2dx *spatial_sampler, mui_chrono_sampler_gauss_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_gauss_gauss_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_gauss_2t *spatial_sampler, mui_chrono_sampler_gauss_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: gauss; temporal sampler: mean
void mui_fetch_gauss_mean_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t_1, float* t_2,
        mui_sampler_gauss_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_mean_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_gauss_2fx *spatial_sampler, mui_chrono_sampler_mean_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_gauss_mean_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_gauss_2d *spatial_sampler, mui_chrono_sampler_mean_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_mean_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_gauss_2dx *spatial_sampler, mui_chrono_sampler_mean_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_gauss_mean_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_gauss_2t *spatial_sampler, mui_chrono_sampler_mean_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: moving average; temporal sampler: exact
void mui_fetch_moving_average_exact_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_moving_average_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_exact_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2,
        float* t_1, float* t_2, mui_sampler_moving_average_2fx *spatial_sampler,
        mui_chrono_sampler_exact_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_moving_average_exact_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_moving_average_2d *spatial_sampler, mui_chrono_sampler_exact_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_exact_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_moving_average_2dx *spatial_sampler,
        mui_chrono_sampler_exact_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_moving_average_exact_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_moving_average_2t *spatial_sampler, mui_chrono_sampler_exact_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: moving average; temporal sampler: gauss
void mui_fetch_moving_average_gauss_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_moving_average_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_gauss_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2,
        float* t_1, float* t_2, mui_sampler_moving_average_2fx *spatial_sampler,
        mui_chrono_sampler_gauss_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_moving_average_gauss_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_moving_average_2d *spatial_sampler, mui_chrono_sampler_gauss_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_gauss_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_moving_average_2dx *spatial_sampler,
        mui_chrono_sampler_gauss_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_moving_average_gauss_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_moving_average_2t *spatial_sampler, mui_chrono_sampler_gauss_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: moving average; temporal sampler: mean
void mui_fetch_moving_average_mean_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_moving_average_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_mean_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_moving_average_2fx *spatial_sampler, mui_chrono_sampler_mean_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_moving_average_mean_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_moving_average_2d *spatial_sampler, mui_chrono_sampler_mean_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_mean_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_moving_average_2dx *spatial_sampler,
        mui_chrono_sampler_mean_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_moving_average_mean_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_moving_average_2t *spatial_sampler, mui_chrono_sampler_mean_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: moving average; temporal sampler: sum
void mui_fetch_moving_average_sum_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_moving_average_2f *spatial_sampler, mui_chrono_sampler_sum_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_sum_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_moving_average_2fx *spatial_sampler, mui_chrono_sampler_sum_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_moving_average_sum_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_moving_average_2d *spatial_sampler, mui_chrono_sampler_sum_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_sum_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_moving_average_2dx *spatial_sampler,
        mui_chrono_sampler_sum_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_moving_average_sum_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_moving_average_2t *spatial_sampler, mui_chrono_sampler_sum_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: nearest neighbor; temporal sampler: exact
void mui_fetch_nearest_neighbor_exact_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_nearest_neighbor_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_exact_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2,
        float* t_1, float* t_2, mui_sampler_nearest_neighbor_2fx *spatial_sampler,
        mui_chrono_sampler_exact_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_nearest_neighbor_exact_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_nearest_neighbor_2d *spatial_sampler,
        mui_chrono_sampler_exact_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_exact_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_nearest_neighbor_2dx *spatial_sampler,
        mui_chrono_sampler_exact_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_nearest_neighbor_exact_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_nearest_neighbor_2t *spatial_sampler,
        mui_chrono_sampler_exact_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: nearest neighbor; temporal sampler: gauss
void mui_fetch_nearest_neighbor_gauss_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_nearest_neighbor_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_gauss_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2,
        float* t_1, float* t_2, mui_sampler_nearest_neighbor_2fx *spatial_sampler,
        mui_chrono_sampler_gauss_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_nearest_neighbor_gauss_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_nearest_neighbor_2d *spatial_sampler,
        mui_chrono_sampler_gauss_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_gauss_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_nearest_neighbor_2dx *spatial_sampler,
        mui_chrono_sampler_gauss_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_nearest_neighbor_gauss_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_nearest_neighbor_2t *spatial_sampler,
        mui_chrono_sampler_gauss_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: nearest neighbor; temporal sampler: mean
void mui_fetch_nearest_neighbor_mean_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_nearest_neighbor_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_mean_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2,
        float* t_1, float* t_2, mui_sampler_nearest_neighbor_2fx *spatial_sampler,
        mui_chrono_sampler_mean_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_nearest_neighbor_mean_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_nearest_neighbor_2d *spatial_sampler,
        mui_chrono_sampler_mean_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_mean_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_nearest_neighbor_2dx *spatial_sampler,
        mui_chrono_sampler_mean_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_nearest_neighbor_mean_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_nearest_neighbor_2t *spatial_sampler,
        mui_chrono_sampler_mean_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: nearest neighbor; temporal sampler: sum
void mui_fetch_nearest_neighbor_sum_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_nearest_neighbor_2f *spatial_sampler, mui_chrono_sampler_sum_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_sum_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2,
        float* t_1, float* t_2, mui_sampler_nearest_neighbor_2fx *spatial_sampler,
        mui_chrono_sampler_sum_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_nearest_neighbor_sum_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_nearest_neighbor_2d *spatial_sampler, mui_chrono_sampler_sum_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_sum_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_nearest_neighbor_2dx *spatial_sampler,
        mui_chrono_sampler_sum_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_nearest_neighbor_sum_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_nearest_neighbor_2t *spatial_sampler, mui_chrono_sampler_sum_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: nearest neighbor; temporal sampler: exact
void mui_fetch_pseudo_nearest_neighbor_exact_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2,
        float* t_1, float* t_2, mui_sampler_pseudo_nearest_neighbor_2f *spatial_sampler,
        mui_chrono_sampler_exact_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_exact_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2,
        float* t_1, float* t_2, mui_sampler_pseudo_nearest_neighbor_2fx *spatial_sampler,
        mui_chrono_sampler_exact_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_exact_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_2d *spatial_sampler,
        mui_chrono_sampler_exact_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_exact_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_2dx *spatial_sampler,
        mui_chrono_sampler_exact_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_exact_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_2t *spatial_sampler,
        mui_chrono_sampler_exact_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: nearest neighbor; temporal sampler: gauss
void mui_fetch_pseudo_nearest_neighbor_gauss_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2,
        float* t_1, float* t_2, mui_sampler_pseudo_nearest_neighbor_2f *spatial_sampler,
        mui_chrono_sampler_gauss_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_gauss_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2,
        float* t_1, float* t_2, mui_sampler_pseudo_nearest_neighbor_2fx *spatial_sampler,
        mui_chrono_sampler_gauss_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_gauss_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_2d *spatial_sampler,
        mui_chrono_sampler_gauss_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_gauss_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_2dx *spatial_sampler,
        mui_chrono_sampler_gauss_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_gauss_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_2t *spatial_sampler,
        mui_chrono_sampler_gauss_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: nearest neighbor; temporal sampler: mean
void mui_fetch_pseudo_nearest_neighbor_mean_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2,
        float* t_1, float* t_2, mui_sampler_pseudo_nearest_neighbor_2f *spatial_sampler,
        mui_chrono_sampler_mean_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_mean_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2,
        float* t_1, float* t_2, mui_sampler_pseudo_nearest_neighbor_2fx *spatial_sampler,
        mui_chrono_sampler_mean_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_mean_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_2d *spatial_sampler,
        mui_chrono_sampler_mean_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_mean_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_2dx *spatial_sampler,
        mui_chrono_sampler_mean_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_mean_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_2t *spatial_sampler,
        mui_chrono_sampler_mean_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: nearest neighbor; temporal sampler: sum
void mui_fetch_pseudo_nearest_neighbor_sum_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2,
        float* t_1, float* t_2, mui_sampler_pseudo_nearest_neighbor_2f *spatial_sampler,
        mui_chrono_sampler_sum_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_sum_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2,
        float* t_1, float* t_2, mui_sampler_pseudo_nearest_neighbor_2fx *spatial_sampler,
        mui_chrono_sampler_sum_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_sum_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_2d *spatial_sampler,
        mui_chrono_sampler_sum_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_sum_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_2dx *spatial_sampler,
        mui_chrono_sampler_sum_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_sum_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_2t *spatial_sampler,
        mui_chrono_sampler_sum_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: shepard quintic; temporal sampler: exact
void mui_fetch_shepard_quintic_exact_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_shepard_quintic_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_exact_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2,
        float* t_1, float* t_2, mui_sampler_shepard_quintic_2fx *spatial_sampler,
        mui_chrono_sampler_exact_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_shepard_quintic_exact_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_shepard_quintic_2d *spatial_sampler,
        mui_chrono_sampler_exact_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_exact_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_shepard_quintic_2dx *spatial_sampler,
        mui_chrono_sampler_exact_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_shepard_quintic_exact_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_shepard_quintic_2t *spatial_sampler,
        mui_chrono_sampler_exact_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: shepard quintic; temporal sampler: gauss
void mui_fetch_shepard_quintic_gauss_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_shepard_quintic_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_gauss_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2,
        float* t_1, float* t_2, mui_sampler_shepard_quintic_2fx *spatial_sampler,
        mui_chrono_sampler_gauss_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_shepard_quintic_gauss_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_shepard_quintic_2d *spatial_sampler,
        mui_chrono_sampler_gauss_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_gauss_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_shepard_quintic_2dx *spatial_sampler,
        mui_chrono_sampler_gauss_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_shepard_quintic_gauss_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_shepard_quintic_2t *spatial_sampler,
        mui_chrono_sampler_gauss_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: shepard quintic; temporal sampler: mean
void mui_fetch_shepard_quintic_mean_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_shepard_quintic_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_mean_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2,
        float* t_1, float* t_2, mui_sampler_shepard_quintic_2fx *spatial_sampler,
        mui_chrono_sampler_mean_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_shepard_quintic_mean_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_shepard_quintic_2d *spatial_sampler, mui_chrono_sampler_mean_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_mean_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_shepard_quintic_2dx *spatial_sampler,
        mui_chrono_sampler_mean_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_shepard_quintic_mean_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_shepard_quintic_2t *spatial_sampler, mui_chrono_sampler_mean_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: shepard quintic; temporal sampler: sum
void mui_fetch_shepard_quintic_sum_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_shepard_quintic_2f *spatial_sampler, mui_chrono_sampler_sum_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_sum_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_shepard_quintic_2fx *spatial_sampler, mui_chrono_sampler_sum_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_shepard_quintic_sum_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_shepard_quintic_2d *spatial_sampler, mui_chrono_sampler_sum_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_sum_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2,
        double* t_1, double* t_2, mui_sampler_shepard_quintic_2dx *spatial_sampler,
        mui_chrono_sampler_sum_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_shepard_quintic_sum_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_shepard_quintic_2t *spatial_sampler, mui_chrono_sampler_sum_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: sph-derived quintic; temporal sampler: exact
void mui_fetch_sph_quintic_exact_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_sph_quintic_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_exact_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_sph_quintic_2fx *spatial_sampler, mui_chrono_sampler_exact_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sph_quintic_exact_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_sph_quintic_2d *spatial_sampler, mui_chrono_sampler_exact_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_exact_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_sph_quintic_2dx *spatial_sampler, mui_chrono_sampler_exact_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sph_quintic_exact_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_sph_quintic_2t *spatial_sampler, mui_chrono_sampler_exact_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: sph-derived quintic; temporal sampler: gauss
void mui_fetch_sph_quintic_gauss_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_sph_quintic_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_gauss_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_sph_quintic_2fx *spatial_sampler, mui_chrono_sampler_gauss_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sph_quintic_gauss_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_sph_quintic_2d *spatial_sampler, mui_chrono_sampler_gauss_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_gauss_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_sph_quintic_2dx *spatial_sampler, mui_chrono_sampler_gauss_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sph_quintic_gauss_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_sph_quintic_2t *spatial_sampler, mui_chrono_sampler_gauss_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: sph-derived quintic; temporal sampler: mean
void mui_fetch_sph_quintic_mean_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_sph_quintic_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_mean_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_sph_quintic_2fx *spatial_sampler, mui_chrono_sampler_mean_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sph_quintic_mean_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_sph_quintic_2d *spatial_sampler, mui_chrono_sampler_mean_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_mean_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_sph_quintic_2dx *spatial_sampler, mui_chrono_sampler_mean_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sph_quintic_mean_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_sph_quintic_2t *spatial_sampler, mui_chrono_sampler_mean_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: sph-derived quintic; temporal sampler: sum
void mui_fetch_sph_quintic_sum_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_sph_quintic_2f *spatial_sampler, mui_chrono_sampler_sum_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_sum_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_sph_quintic_2fx *spatial_sampler, mui_chrono_sampler_sum_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sph_quintic_sum_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_sph_quintic_2d *spatial_sampler, mui_chrono_sampler_sum_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_sum_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_sph_quintic_2dx *spatial_sampler, mui_chrono_sampler_sum_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sph_quintic_sum_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_sph_quintic_2t *spatial_sampler, mui_chrono_sampler_sum_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: summation quintic; temporal sampler: exact
void mui_fetch_sum_quintic_exact_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_sum_quintic_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_exact_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_sum_quintic_2fx *spatial_sampler, mui_chrono_sampler_exact_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sum_quintic_exact_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_sum_quintic_2d *spatial_sampler, mui_chrono_sampler_exact_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_exact_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_sum_quintic_2dx *spatial_sampler, mui_chrono_sampler_exact_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sum_quintic_exact_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_sum_quintic_2t *spatial_sampler, mui_chrono_sampler_exact_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: summation quintic; temporal sampler: gauss
void mui_fetch_sum_quintic_gauss_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_sum_quintic_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_gauss_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_sum_quintic_2fx *spatial_sampler, mui_chrono_sampler_gauss_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sum_quintic_gauss_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_sum_quintic_2d *spatial_sampler, mui_chrono_sampler_gauss_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_gauss_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_sum_quintic_2dx *spatial_sampler, mui_chrono_sampler_gauss_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sum_quintic_gauss_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_sum_quintic_2t *spatial_sampler, mui_chrono_sampler_gauss_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: summation quintic; temporal sampler: mean
void mui_fetch_sum_quintic_mean_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_sum_quintic_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_mean_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_sum_quintic_2fx *spatial_sampler, mui_chrono_sampler_mean_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sum_quintic_mean_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_sum_quintic_2d *spatial_sampler, mui_chrono_sampler_mean_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_mean_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_sum_quintic_2dx *spatial_sampler, mui_chrono_sampler_mean_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sum_quintic_mean_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_sum_quintic_2t *spatial_sampler, mui_chrono_sampler_mean_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: summation quintic; temporal sampler: sum
void mui_fetch_sum_quintic_sum_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_sum_quintic_2f *spatial_sampler, mui_chrono_sampler_sum_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_sum_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_sum_quintic_2fx *spatial_sampler, mui_chrono_sampler_sum_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sum_quintic_sum_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_sum_quintic_2d *spatial_sampler, mui_chrono_sampler_sum_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_sum_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_sum_quintic_2dx *spatial_sampler, mui_chrono_sampler_sum_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sum_quintic_sum_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_sum_quintic_2t *spatial_sampler, mui_chrono_sampler_sum_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

#ifdef USE_RBF
// Spatial sampler: radial basis function; temporal sampler: exact
void mui_fetch_rbf_exact_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t_1, float* t_2,
        mui_sampler_rbf_2f *spatial_sampler, mui_chrono_sampler_exact_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_exact_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_rbf_2fx *spatial_sampler, mui_chrono_sampler_exact_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_rbf_exact_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_rbf_2d *spatial_sampler, mui_chrono_sampler_exact_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_exact_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_rbf_2dx *spatial_sampler, mui_chrono_sampler_exact_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_rbf_exact_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_rbf_2t *spatial_sampler, mui_chrono_sampler_exact_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: radial basis function; temporal sampler: gauss
void mui_fetch_rbf_gauss_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t_1, float* t_2,
        mui_sampler_rbf_2f *spatial_sampler, mui_chrono_sampler_gauss_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_gauss_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2, float* t_1,
        float* t_2, mui_sampler_rbf_2fx *spatial_sampler, mui_chrono_sampler_gauss_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_rbf_gauss_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_rbf_2d *spatial_sampler, mui_chrono_sampler_gauss_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_gauss_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_rbf_2dx *spatial_sampler, mui_chrono_sampler_gauss_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_rbf_gauss_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_rbf_2t *spatial_sampler, mui_chrono_sampler_gauss_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: radial basis function; temporal sampler: mean
void mui_fetch_rbf_mean_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t_1, float* t_2,
        mui_sampler_rbf_2f *spatial_sampler, mui_chrono_sampler_mean_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_mean_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2, float* t_1, float* t_2,
        mui_sampler_rbf_2fx *spatial_sampler, mui_chrono_sampler_mean_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_rbf_mean_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t_1, double* t_2,
        mui_sampler_rbf_2d *spatial_sampler, mui_chrono_sampler_mean_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_mean_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_rbf_2dx *spatial_sampler, mui_chrono_sampler_mean_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_rbf_mean_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t_1, double* t_2,
        mui_sampler_rbf_2t *spatial_sampler, mui_chrono_sampler_mean_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: radial basis function; temporal sampler: sum
void mui_fetch_rbf_sum_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float* point_1, float* point_2, float* t_1, float* t_2,
        mui_sampler_rbf_2f *spatial_sampler, mui_chrono_sampler_sum_2f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2f(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_sum_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float* point_1, float* point_2, float* t_1, float* t_2,
        mui_sampler_rbf_2fx *spatial_sampler, mui_chrono_sampler_sum_2fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2fx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_rbf_sum_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double* point_1, double* point_2, double* t_1, double* t_2,
        mui_sampler_rbf_2d *spatial_sampler, mui_chrono_sampler_sum_2d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2d(*point_1,*point_2), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_sum_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double* point_1, double* point_2, double* t_1,
        double* t_2, mui_sampler_rbf_2dx *spatial_sampler, mui_chrono_sampler_sum_2dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point2dx(*point_1,*point_2), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_rbf_sum_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double* point_1, double* point_2, double* t_1, double* t_2,
        mui_sampler_rbf_2t *spatial_sampler, mui_chrono_sampler_sum_2t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_2D::point_type point_fetch(static_cast<mui::mui_f_wrapper_2D::REAL>(*point_1,*point_2));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}
#endif

/*******************************************************************
 * MUI functions for 2d data point only fetch using one time value  *
 ********************************************************************/

// Temporal sampler: exact
void mui_fetch_points_exact_2f_f(mui_uniface_2f *uniface, const char *attr, float *t,
        mui_chrono_sampler_exact_2f *temporal_sampler, float** ret_points_1, float** ret_points_2, int *num_points) {
    std::vector<mui::point2f> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_exact_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *t,
        mui_chrono_sampler_exact_2fx *temporal_sampler, float** ret_points_1, float** ret_points_2, int *num_points) {
    std::vector<mui::point2fx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_exact_2d_f(mui_uniface_2d *uniface, const char *attr, double *t,
        mui_chrono_sampler_exact_2d *temporal_sampler, double** ret_points_1, double** ret_points_2, int *num_points) {
    std::vector<mui::point2d> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_exact_2dx_f(mui_uniface_2dx *uniface, const char *attr, double *t,
        mui_chrono_sampler_exact_2dx *temporal_sampler, double** ret_points_1, double** ret_points_2, int *num_points) {
    std::vector<mui::point2dx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_exact_2t_f(mui_uniface_2t *uniface, const char *attr, double *t,
        mui_chrono_sampler_exact_2t *temporal_sampler, double** ret_points_1, double** ret_points_2, int *num_points) {
    std::vector<mui::mui_f_wrapper_2D::point_type> ret_pts = uniface->fetch_points<float>(std::string(attr), *t,
            *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = static_cast<double>(ret_pts[i][0]);
        *ret_points_2[i] = static_cast<double>(ret_pts[i][1]);
    }
    *num_points = static_cast<int>(ret_pts.size());
}

// Temporal sampler: gauss
void mui_fetch_points_gauss_2f_f(mui_uniface_2f *uniface, const char *attr, float *t,
        mui_chrono_sampler_gauss_2f *temporal_sampler, float** ret_points_1, float** ret_points_2, int *num_points) {
    std::vector<mui::point2f> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_gauss_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *t,
        mui_chrono_sampler_gauss_2fx *temporal_sampler, float** ret_points_1, float** ret_points_2, int *num_points) {
    std::vector<mui::point2fx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_gauss_2d_f(mui_uniface_2d *uniface, const char *attr, double *t,
        mui_chrono_sampler_gauss_2d *temporal_sampler, double** ret_points_1, double** ret_points_2, int *num_points) {
    std::vector<mui::point2d> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_gauss_2dx_f(mui_uniface_2dx *uniface, const char *attr, double *t,
        mui_chrono_sampler_gauss_2dx *temporal_sampler, double** ret_points_1, double** ret_points_2, int *num_points) {
    std::vector<mui::point2dx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_gauss_2t_f(mui_uniface_2t *uniface, const char *attr, double *t,
        mui_chrono_sampler_gauss_2t *temporal_sampler, double** ret_points_1, double** ret_points_2, int *num_points) {
    std::vector<mui::mui_f_wrapper_2D::point_type> ret_pts = uniface->fetch_points<float>(std::string(attr), *t,
            *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = static_cast<double>(ret_pts[i][0]);
        *ret_points_2[i] = static_cast<double>(ret_pts[i][1]);
    }
    *num_points = static_cast<int>(ret_pts.size());
}

// Temporal sampler: mean
void mui_fetch_points_mean_2f_f(mui_uniface_2f *uniface, const char *attr, float *t,
        mui_chrono_sampler_mean_2f *temporal_sampler, float** ret_points_1, float** ret_points_2, int *num_points) {
    std::vector<mui::point2f> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_mean_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *t,
        mui_chrono_sampler_mean_2fx *temporal_sampler, float** ret_points_1, float** ret_points_2, int *num_points) {
    std::vector<mui::point2fx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_mean_2d_f(mui_uniface_2d *uniface, const char *attr, double *t,
        mui_chrono_sampler_mean_2d *temporal_sampler, double** ret_points_1, double** ret_points_2, int *num_points) {
    std::vector<mui::point2d> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_mean_2dx_f(mui_uniface_2dx *uniface, const char *attr, double *t,
        mui_chrono_sampler_mean_2dx *temporal_sampler, double** ret_points_1, double** ret_points_2, int *num_points) {
    std::vector<mui::point2dx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_mean_2t_f(mui_uniface_2t *uniface, const char *attr, double *t,
        mui_chrono_sampler_mean_2t *temporal_sampler, double** ret_points_1, double** ret_points_2, int *num_points) {
    std::vector<mui::mui_f_wrapper_2D::point_type> ret_pts = uniface->fetch_points<float>(std::string(attr), *t,
            *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = static_cast<double>(ret_pts[i][0]);
        *ret_points_2[i] = static_cast<double>(ret_pts[i][1]);
    }
    *num_points = static_cast<int>(ret_pts.size());
}

// Temporal sampler: sum
void mui_fetch_points_sum_2f_f(mui_uniface_2f *uniface, const char *attr, float *t,
        mui_chrono_sampler_sum_2f *temporal_sampler, float** ret_points_1, float** ret_points_2, int *num_points) {
    std::vector<mui::point2f> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_sum_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *t,
        mui_chrono_sampler_sum_2fx *temporal_sampler, float** ret_points_1, float** ret_points_2, int *num_points) {
    std::vector<mui::point2fx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_sum_2d_f(mui_uniface_2d *uniface, const char *attr, double *t,
        mui_chrono_sampler_sum_2d *temporal_sampler, double** ret_points_1, double** ret_points_2, int *num_points) {
    std::vector<mui::point2d> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_sum_2dx_f(mui_uniface_2dx *uniface, const char *attr, double *t,
        mui_chrono_sampler_sum_2dx *temporal_sampler, double** ret_points_1, double** ret_points_2, int *num_points) {
    std::vector<mui::point2dx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_sum_2t_f(mui_uniface_2t *uniface, const char *attr, double *t,
        mui_chrono_sampler_sum_2t *temporal_sampler, double** ret_points_1, double** ret_points_2, int *num_points) {
    std::vector<mui::mui_f_wrapper_2D::point_type> ret_pts = uniface->fetch_points<float>(std::string(attr), *t,
            *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = static_cast<double>(ret_pts[i][0]);
        *ret_points_2[i] = static_cast<double>(ret_pts[i][1]);
    }
    *num_points = static_cast<int>(ret_pts.size());
}

/*******************************************************************
 * MUI functions for 2d data point only fetch using two time values *
 ********************************************************************/

// Temporal sampler: exact
void mui_fetch_points_exact_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_exact_2f *temporal_sampler, float** ret_points_1, float** ret_points_2, int *num_points) {
    std::vector<mui::point2f> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_exact_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_exact_2fx *temporal_sampler, float** ret_points_1, float** ret_points_2, int *num_points) {
    std::vector<mui::point2fx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_exact_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_exact_2d *temporal_sampler, double** ret_points_1, double** ret_points_2, int *num_points) {
    std::vector<mui::point2d> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_exact_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_exact_2dx *temporal_sampler, double** ret_points_1, double** ret_points_2, int *num_points) {
    std::vector<mui::point2dx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_exact_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_exact_2t *temporal_sampler, double** ret_points_1, double** ret_points_2, int *num_points) {
    std::vector<mui::mui_f_wrapper_2D::point_type> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2,
            *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = static_cast<double>(ret_pts[i][0]);
        *ret_points_2[i] = static_cast<double>(ret_pts[i][1]);
    }
    *num_points = static_cast<int>(ret_pts.size());
}

// Temporal sampler: gauss
void mui_fetch_points_gauss_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_gauss_2f *temporal_sampler, float** ret_points_1, float** ret_points_2, int *num_points) {
    std::vector<mui::point2f> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_gauss_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_gauss_2fx *temporal_sampler, float** ret_points_1, float** ret_points_2, int *num_points) {
    std::vector<mui::point2fx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_gauss_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_gauss_2d *temporal_sampler, double** ret_points_1, double** ret_points_2, int *num_points) {
    std::vector<mui::point2d> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_gauss_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_gauss_2dx *temporal_sampler, double** ret_points_1, double** ret_points_2, int *num_points) {
    std::vector<mui::point2dx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_gauss_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_gauss_2t *temporal_sampler, double** ret_points_1, double** ret_points_2, int *num_points) {
    std::vector<mui::mui_f_wrapper_2D::point_type> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2,
            *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = static_cast<double>(ret_pts[i][0]);
        *ret_points_2[i] = static_cast<double>(ret_pts[i][1]);
    }
    *num_points = static_cast<int>(ret_pts.size());
}

// Temporal sampler: mean
void mui_fetch_points_mean_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_mean_2f *temporal_sampler, float** ret_points_1, float** ret_points_2, int *num_points) {
    std::vector<mui::point2f> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_mean_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_mean_2fx *temporal_sampler, float** ret_points_1, float** ret_points_2, int *num_points) {
    std::vector<mui::point2fx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_mean_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_mean_2d *temporal_sampler, double** ret_points_1, double** ret_points_2, int *num_points) {
    std::vector<mui::point2d> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_mean_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_mean_2dx *temporal_sampler, double** ret_points_1, double** ret_points_2, int *num_points) {
    std::vector<mui::point2dx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_mean_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_mean_2t *temporal_sampler, double** ret_points_1, double** ret_points_2, int *num_points) {
    std::vector<mui::mui_f_wrapper_2D::point_type> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2,
            *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = static_cast<double>(ret_pts[i][0]);
        *ret_points_2[i] = static_cast<double>(ret_pts[i][1]);
    }
    *num_points = static_cast<int>(ret_pts.size());
}

// Temporal sampler: sum
void mui_fetch_points_sum_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_sum_2f *temporal_sampler, float** ret_points_1, float** ret_points_2, int *num_points) {
    std::vector<mui::point2f> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_sum_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_sum_2fx *temporal_sampler, float** ret_points_1, float** ret_points_2, int *num_points) {
    std::vector<mui::point2fx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_sum_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_sum_2d *temporal_sampler, double** ret_points_1, double** ret_points_2, int *num_points) {
    std::vector<mui::point2d> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_sum_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_sum_2dx *temporal_sampler, double** ret_points_1, double** ret_points_2, int *num_points) {
    std::vector<mui::point2dx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_sum_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_sum_2t *temporal_sampler, double** ret_points_1, double** ret_points_2, int *num_points) {
    std::vector<mui::mui_f_wrapper_2D::point_type> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2,
            *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = static_cast<double>(ret_pts[i][0]);
        *ret_points_2[i] = static_cast<double>(ret_pts[i][1]);
    }
    *num_points = static_cast<int>(ret_pts.size());
}

/*******************************************************************
 * MUI functions for 2d data values only fetch using one time value  *
 ********************************************************************/

// Temporal sampler: exact
void mui_fetch_values_exact_2f_f(mui_uniface_2f *uniface, const char *attr, float *t,
        mui_chrono_sampler_exact_2f *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_exact_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *t,
        mui_chrono_sampler_exact_2fx *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_exact_2d_f(mui_uniface_2d *uniface, const char *attr, double *t,
        mui_chrono_sampler_exact_2d *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_exact_2dx_f(mui_uniface_2dx *uniface, const char *attr, double *t,
        mui_chrono_sampler_exact_2dx *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_exact_2t_f(mui_uniface_2t *uniface, const char *attr, double *t,
        mui_chrono_sampler_exact_2t *temporal_sampler, double **values, int *num_values) {
    std::vector<mui::mui_f_wrapper_2D::REAL> ret_vals = uniface->fetch_values<mui::mui_f_wrapper_2D::REAL>(
            std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = static_cast<double>(ret_vals[i]);
    }
    *num_values = static_cast<int>(ret_vals.size());
}

// Temporal sampler: gauss
void mui_fetch_values_gauss_2f_f(mui_uniface_2f *uniface, const char *attr, float *t,
        mui_chrono_sampler_gauss_2f *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_gauss_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *t,
        mui_chrono_sampler_gauss_2fx *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_gauss_2d_f(mui_uniface_2d *uniface, const char *attr, double *t,
        mui_chrono_sampler_gauss_2d *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_gauss_2dx_f(mui_uniface_2dx *uniface, const char *attr, double *t,
        mui_chrono_sampler_gauss_2dx *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_gauss_2t_f(mui_uniface_2t *uniface, const char *attr, double *t,
        mui_chrono_sampler_gauss_2t *temporal_sampler, double **values, int *num_values) {
    std::vector<mui::mui_f_wrapper_2D::REAL> ret_vals = uniface->fetch_values<mui::mui_f_wrapper_2D::REAL>(
            std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = static_cast<double>(ret_vals[i]);
    }
    *num_values = static_cast<int>(ret_vals.size());
}

// Temporal sampler: mean
void mui_fetch_values_mean_2f_f(mui_uniface_2f *uniface, const char *attr, float *t,
        mui_chrono_sampler_mean_2f *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_mean_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *t,
        mui_chrono_sampler_mean_2fx *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_mean_2d_f(mui_uniface_2d *uniface, const char *attr, double *t,
        mui_chrono_sampler_mean_2d *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_mean_2dx_f(mui_uniface_2dx *uniface, const char *attr, double *t,
        mui_chrono_sampler_mean_2dx *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_mean_2t_f(mui_uniface_2t *uniface, const char *attr, double *t,
        mui_chrono_sampler_mean_2t *temporal_sampler, double **values, int *num_values) {
    std::vector<mui::mui_f_wrapper_2D::REAL> ret_vals = uniface->fetch_values<mui::mui_f_wrapper_2D::REAL>(
            std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = static_cast<double>(ret_vals[i]);
    }
    *num_values = static_cast<int>(ret_vals.size());
}

// Temporal sampler: sum
void mui_fetch_values_sum_2f_f(mui_uniface_2f *uniface, const char *attr, float *t,
        mui_chrono_sampler_sum_2f *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_sum_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *t,
        mui_chrono_sampler_sum_2fx *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_sum_2d_f(mui_uniface_2d *uniface, const char *attr, double *t,
        mui_chrono_sampler_sum_2d *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_sum_2dx_f(mui_uniface_2dx *uniface, const char *attr, double *t,
        mui_chrono_sampler_sum_2dx *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_sum_2t_f(mui_uniface_2t *uniface, const char *attr, double *t,
        mui_chrono_sampler_sum_2t *temporal_sampler, double **values, int *num_values) {
    std::vector<mui::mui_f_wrapper_2D::REAL> ret_vals = uniface->fetch_values<mui::mui_f_wrapper_2D::REAL>(
            std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = static_cast<double>(ret_vals[i]);
    }
    *num_values = static_cast<int>(ret_vals.size());
}

/********************************************************************
 * MUI functions for 2d data values only fetch using two time values *
 *********************************************************************/

// Temporal sampler: exact
void mui_fetch_values_exact_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_exact_2f *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_exact_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_exact_2fx *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_exact_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_exact_2d *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_exact_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_exact_2dx *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_exact_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_exact_2t *temporal_sampler, double **values, int *num_values) {
    std::vector<mui::mui_f_wrapper_2D::REAL> ret_vals = uniface->fetch_values<mui::mui_f_wrapper_2D::REAL>(
            std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = static_cast<double>(ret_vals[i]);
    }
    *num_values = static_cast<int>(ret_vals.size());
}

// Temporal sampler: gauss
void mui_fetch_values_gauss_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_gauss_2f *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_gauss_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_gauss_2fx *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_gauss_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_gauss_2d *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_gauss_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_gauss_2dx *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_gauss_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_gauss_2t *temporal_sampler, double **values, int *num_values) {
    std::vector<mui::mui_f_wrapper_2D::REAL> ret_vals = uniface->fetch_values<mui::mui_f_wrapper_2D::REAL>(
            std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = static_cast<double>(ret_vals[i]);
    }
    *num_values = static_cast<int>(ret_vals.size());
}

// Temporal sampler: mean
void mui_fetch_values_mean_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_mean_2f *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_mean_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_mean_2fx *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_mean_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_mean_2d *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_mean_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_mean_2dx *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_mean_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_mean_2t *temporal_sampler, double **values, int *num_values) {
    std::vector<mui::mui_f_wrapper_2D::REAL> ret_vals = uniface->fetch_values<mui::mui_f_wrapper_2D::REAL>(
            std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = static_cast<double>(ret_vals[i]);
    }
    *num_values = static_cast<int>(ret_vals.size());
}

// Temporal sampler: sum
void mui_fetch_values_sum_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_sum_2f *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_sum_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_sum_2fx *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_sum_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_sum_2d *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_sum_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_sum_2dx *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_sum_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_sum_2t *temporal_sampler, double **values, int *num_values) {
    std::vector<mui::mui_f_wrapper_2D::REAL> ret_vals = uniface->fetch_values<mui::mui_f_wrapper_2D::REAL>(
            std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = static_cast<double>(ret_vals[i]);
    }
    *num_values = static_cast<int>(ret_vals.size());
}

/*******************************************
 * MUI functions for single parameter fetch *
 ********************************************/

void mui_fetch_2f_param_f(mui_uniface_2f *uniface, const char *attr, float *return_value) {
    *return_value = uniface->fetch<float>(std::string(attr));
}

void mui_fetch_2fx_param_f(mui_uniface_2fx *uniface, const char *attr, float *return_value) {
    *return_value = uniface->fetch<float>(std::string(attr));
}

void mui_fetch_2d_param_f(mui_uniface_2d *uniface, const char *attr, float *return_value) {
    *return_value = uniface->fetch<double>(std::string(attr));
}

void mui_fetch_2dx_param_f(mui_uniface_2dx *uniface, const char *attr, float *return_value) {
    *return_value = uniface->fetch<double>(std::string(attr));
}

void mui_fetch_2t_param_f(mui_uniface_2t *uniface, const char *attr, float *return_value) {
    mui::mui_f_wrapper_2D::REAL ret_val = uniface->fetch<mui::mui_f_wrapper_2D::REAL>(std::string(attr));
    *return_value = static_cast<double>(ret_val);
}

/******************************************
 * MUI data receive test functions         *
 ******************************************/

// Data ready test using single time value
void mui_is_ready_2f_f(mui_uniface_2f *uniface, const char *attr, float *t, int *return_value) {
    *return_value = uniface->is_ready(std::string(attr), *t);
}

void mui_is_ready_2fx_f(mui_uniface_2fx *uniface, const char *attr, float *t, int *return_value) {
    *return_value = uniface->is_ready(std::string(attr), *t);
}

void mui_is_ready_2d_f(mui_uniface_2d *uniface, const char *attr, double *t, int *return_value) {
    *return_value = uniface->is_ready(std::string(attr), *t);
}

void mui_is_ready_2dx_f(mui_uniface_2dx *uniface, const char *attr, double *t, int *return_value) {
    *return_value = uniface->is_ready(std::string(attr), *t);
}

void mui_is_ready_2t_f(mui_uniface_2t *uniface, const char *attr, double *t, int *return_value) {
    *return_value = uniface->is_ready(std::string(attr), static_cast<mui::mui_f_wrapper_2D::time_type>(*t));
}

// Data ready test using two time values
void mui_is_ready_2f_pair_f(mui_uniface_2f *uniface, const char *attr, float *t_1, float *t_2, int *return_value) {
    *return_value = uniface->is_ready(std::string(attr), *t_1, *t_2);
}

void mui_is_ready_2fx_pair_f(mui_uniface_2fx *uniface, const char *attr, float *t_1, float *t_2, int *return_value) {
    *return_value = uniface->is_ready(std::string(attr), *t_1, *t_2);
}

void mui_is_ready_2d_pair_f(mui_uniface_2d *uniface, const char *attr, double *t_1, double *t_2, int *return_value) {
    *return_value = uniface->is_ready(std::string(attr), *t_1, *t_2);
}

void mui_is_ready_2dx_pair_f(mui_uniface_2dx *uniface, const char *attr, double *t_1, double *t_2, int *return_value) {
    *return_value = uniface->is_ready(std::string(attr), *t_1, *t_2);
}

void mui_is_ready_2t_pair_f(mui_uniface_2t *uniface, const char *attr, double *t_1, double *t_2, int *return_value) {
    *return_value = uniface->is_ready(std::string(attr), static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1),
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2));
}

/******************************************
 * MUI Smart Send functions                *
 ******************************************/

// Send span announce using 2d box geometry
void mui_announce_send_span_2f_box_f(mui_uniface_2f *uniface, float *box_1_1, float *box_1_2, float *box_2_1, float *box_2_2, float *t_start,
        float *t_timeout, int *synchronised) {
    mui::point2f point_1(*box_1_1,*box_1_2);
    mui::point2f point_2(*box_2_1,*box_2_2);
    mui::geometry::box2f bound_box(point_1, point_2);
    uniface->announce_send_span(*t_start, *t_timeout, bound_box, static_cast<bool>(*synchronised));
}

void mui_announce_send_span_2fx_box_f(mui_uniface_2fx *uniface, float *box_1_1, float *box_1_2, float *box_2_1, float *box_2_2, float *t_start,
        float *t_timeout, int *synchronised) {
    mui::point2fx point_1(*box_1_1,*box_1_2);
    mui::point2fx point_2(*box_2_1,*box_2_2);
    mui::geometry::box2fx bound_box(point_1, point_2);
    uniface->announce_send_span(*t_start, *t_timeout, bound_box, static_cast<bool>(*synchronised));
}

void mui_announce_send_span_2d_box_f(mui_uniface_2d *uniface, double *box_1_1, double *box_1_2, double *box_2_1, double *box_2_2, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::point2d point_1(*box_1_1,*box_1_2);
    mui::point2d point_2(*box_2_1,*box_2_2);
    mui::geometry::box2d bound_box(point_1, point_2);
    uniface->announce_send_span(*t_start, *t_timeout, bound_box);
}

void mui_announce_send_span_2dx_box_f(mui_uniface_2dx *uniface, double *box_1_1, double *box_1_2, double *box_2_1, double *box_2_2, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::point2dx point_1(*box_1_1,*box_1_2);
    mui::point2dx point_2(*box_2_1,*box_2_2);
    mui::geometry::box2dx bound_box(point_1, point_2);
    uniface->announce_send_span(*t_start, *t_timeout, bound_box, static_cast<bool>(*synchronised));
}

void mui_announce_send_span_2t_box_f(mui_uniface_2t *uniface, double *box_1_1, double *box_1_2, double *box_2_1, double *box_2_2, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::mui_f_wrapper_2D::point_type point_1(static_cast<mui::mui_f_wrapper_2D::REAL>(*box_1_1,*box_1_2));
    mui::mui_f_wrapper_2D::point_type point_2(static_cast<mui::mui_f_wrapper_2D::REAL>(*box_2_1,*box_2_2));
    mui::geometry::box<mui::mui_f_wrapper_2D> bound_box(point_1, point_2);
    uniface->announce_send_span(static_cast<mui::mui_f_wrapper_2D::time_type>(*t_start),
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_timeout), bound_box, static_cast<bool>(*synchronised));
}

// Send span announce using 2d sphere geometry
void mui_announce_send_span_2f_sphere_f(mui_uniface_2f *uniface, float *centre_1, float *centre_2, float *radius, float *t_start,
        float *t_timeout, int *synchronised) {
    mui::geometry::sphere2f bound_sphere(mui::point2f(*centre_1,*centre_2), *radius);
    uniface->announce_send_span(*t_start, *t_timeout, bound_sphere, static_cast<bool>(*synchronised));
}

void mui_announce_send_span_2fx_sphere_f(mui_uniface_2fx *uniface, float *centre_1, float *centre_2, float *radius, float *t_start,
        float *t_timeout, int *synchronised) {
    mui::geometry::sphere2fx bound_sphere(mui::point2fx(*centre_1,*centre_2), *radius);
    uniface->announce_send_span(*t_start, *t_timeout, bound_sphere, static_cast<bool>(*synchronised));
}

void mui_announce_send_span_2d_sphere_f(mui_uniface_2d *uniface, double *centre_1, double *centre_2, double *radius, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::geometry::sphere2d bound_sphere(mui::point2d(*centre_1,*centre_2), *radius);
    uniface->announce_send_span(*t_start, *t_timeout, bound_sphere, static_cast<bool>(*synchronised));
}

void mui_announce_send_span_2dx_sphere_f(mui_uniface_2dx *uniface, double *centre_1, double *centre_2, double *radius, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::geometry::sphere2dx bound_sphere(mui::point2dx(*centre_1,*centre_2), *radius);
    uniface->announce_send_span(*t_start, *t_timeout, bound_sphere, static_cast<bool>(*synchronised));
}

void mui_announce_send_span_2t_sphere_f(mui_uniface_2t *uniface, double *centre_1, double *centre_2, double *radius, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::geometry::sphere<mui::mui_f_wrapper_2D> bound_sphere(
            mui::mui_f_wrapper_2D::point_type(static_cast<mui::mui_f_wrapper_2D::REAL>(*centre_1,*centre_2)),
            static_cast<mui::mui_f_wrapper_2D::REAL>(*radius));
    uniface->announce_send_span(static_cast<mui::mui_f_wrapper_2D::time_type>(*t_start),
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_timeout), bound_sphere, static_cast<bool>(*synchronised));
}

// Receive span announce using 2d box geometry
void mui_announce_recv_span_2f_box_f(mui_uniface_2f *uniface, float *box_1_1, float *box_1_2, float *box_2_1, float *box_2_2, float *t_start,
        float *t_timeout, int *synchronised) {
    mui::point2f point_1(*box_1_1,*box_1_2);
    mui::point2f point_2(*box_2_1,*box_2_2);
    mui::geometry::box2f bound_box(point_1, point_2);
    uniface->announce_recv_span(*t_start, *t_timeout, bound_box, static_cast<bool>(*synchronised));
}

void mui_announce_recv_span_2fx_box_f(mui_uniface_2fx *uniface, float *box_1_1, float *box_1_2, float *box_2_1, float *box_2_2, float *t_start,
        float *t_timeout, int *synchronised) {
    mui::point2fx point_1(*box_1_1,*box_1_2);
    mui::point2fx point_2(*box_2_1,*box_2_2);
    mui::geometry::box2fx bound_box(point_1, point_2);
    uniface->announce_recv_span(*t_start, *t_timeout, bound_box, static_cast<bool>(*synchronised));
}

void mui_announce_recv_span_2d_box_f(mui_uniface_2d *uniface, double *box_1_1, double *box_1_2, double *box_2_1, double *box_2_2, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::point2d point_1(*box_1_1,*box_1_2);
    mui::point2d point_2(*box_2_1,*box_2_2);
    mui::geometry::box2d bound_box(point_1, point_2);
    uniface->announce_recv_span(*t_start, *t_timeout, bound_box, static_cast<bool>(*synchronised));
}

void mui_announce_recv_span_2dx_box_f(mui_uniface_2dx *uniface, double *box_1_1, double *box_1_2, double *box_2_1, double *box_2_2, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::point2dx point_1(*box_1_1,*box_1_2);
    mui::point2dx point_2(*box_2_1,*box_2_2);
    mui::geometry::box2dx bound_box(point_1, point_2);
    uniface->announce_recv_span(*t_start, *t_timeout, bound_box, static_cast<bool>(*synchronised));
}

void mui_announce_recv_span_2t_box_f(mui_uniface_2t *uniface, double *box_1_1, double *box_1_2, double *box_2_1, double *box_2_2, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::mui_f_wrapper_2D::point_type point_1(static_cast<mui::mui_f_wrapper_2D::REAL>(*box_1_1,*box_1_2));
    mui::mui_f_wrapper_2D::point_type point_2(static_cast<mui::mui_f_wrapper_2D::REAL>(*box_2_1,*box_2_2));
    mui::geometry::box<mui::mui_f_wrapper_2D> bound_box(point_1, point_2);
    uniface->announce_recv_span(static_cast<mui::mui_f_wrapper_2D::time_type>(*t_start),
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_timeout), bound_box, static_cast<bool>(*synchronised));
}

// Receive span announce using 2d sphere geometry
void mui_announce_recv_span_2f_sphere_f(mui_uniface_2f *uniface, float *centre_1, float *centre_2, float *radius, float *t_start,
        float *t_timeout, int *synchronised) {
    mui::geometry::sphere2f bound_sphere(mui::point2f(*centre_1,*centre_2), *radius);
    uniface->announce_recv_span(*t_start, *t_timeout, bound_sphere, static_cast<bool>(*synchronised));
}

void mui_announce_recv_span_2fx_sphere_f(mui_uniface_2fx *uniface, float *centre_1, float *centre_2, float *radius, float *t_start,
        float *t_timeout, int *synchronised) {
    mui::geometry::sphere2fx bound_sphere(mui::point2fx(*centre_1,*centre_2), *radius);
    uniface->announce_recv_span(*t_start, *t_timeout, bound_sphere, static_cast<bool>(*synchronised));
}

void mui_announce_recv_span_2d_sphere_f(mui_uniface_2d *uniface, double *centre_1, double *centre_2, double *radius, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::geometry::sphere2d bound_sphere(mui::point2d(*centre_1,*centre_2), *radius);
    uniface->announce_recv_span(*t_start, *t_timeout, bound_sphere, static_cast<bool>(*synchronised));
}

void mui_announce_recv_span_2dx_sphere_f(mui_uniface_2dx *uniface, double *centre_1, double *centre_2, double *radius, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::geometry::sphere2dx bound_sphere(mui::point2dx(*centre_1,*centre_2), *radius);
    uniface->announce_recv_span(*t_start, *t_timeout, bound_sphere, static_cast<bool>(*synchronised));
}

void mui_announce_recv_span_2t_sphere_f(mui_uniface_2t *uniface, double *centre_1, double *centre_2, double *radius, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::geometry::sphere<mui::mui_f_wrapper_2D> bound_sphere(
            mui::mui_f_wrapper_2D::point_type(static_cast<mui::mui_f_wrapper_2D::REAL>(*centre_1,*centre_2)),
            static_cast<mui::mui_f_wrapper_2D::REAL>(*radius));
    uniface->announce_recv_span(static_cast<mui::mui_f_wrapper_2D::time_type>(*t_start),
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_timeout), bound_sphere, static_cast<bool>(*synchronised));
}

// Send disable announce (local call per MPI rank)
void mui_announce_send_disable_2f_f(mui_uniface_2f *uniface, int *synchronised) {
    uniface->announce_send_disable(static_cast<bool>(*synchronised));
}

void mui_announce_send_disable_2fx_f(mui_uniface_2fx *uniface, int *synchronised) {
    uniface->announce_send_disable(static_cast<bool>(*synchronised));
}

void mui_announce_send_disable_2d_f(mui_uniface_2d *uniface, int *synchronised) {
    uniface->announce_send_disable(static_cast<bool>(*synchronised));
}

void mui_announce_send_disable_2dx_f(mui_uniface_2dx *uniface, int *synchronised) {
    uniface->announce_send_disable(static_cast<bool>(*synchronised));
}

void mui_announce_send_disable_2t_f(mui_uniface_2t *uniface, int *synchronised) {
    uniface->announce_send_disable(static_cast<bool>(*synchronised));
}

// Receive disable announce (local call per MPI rank)
void mui_announce_recv_disable_2f_f(mui_uniface_2f *uniface, int *synchronised) {
    uniface->announce_recv_disable(static_cast<bool>(*synchronised));
}

void mui_announce_recv_disable_2fx_f(mui_uniface_2fx *uniface, int *synchronised) {
    uniface->announce_recv_disable(static_cast<bool>(*synchronised));
}

void mui_announce_recv_disable_2d_f(mui_uniface_2d *uniface, int *synchronised) {
    uniface->announce_recv_disable(static_cast<bool>(*synchronised));
}

void mui_announce_recv_disable_2dx_f(mui_uniface_2dx *uniface, int *synchronised) {
    uniface->announce_recv_disable(static_cast<bool>(*synchronised));
}

void mui_announce_recv_disable_2t_f(mui_uniface_2t *uniface, int *synchronised) {
    uniface->announce_recv_disable(static_cast<bool>(*synchronised));
}

/******************************************
 * MUI barrier functions                   *
 ******************************************/

// Barrier at single time value
void mui_barrier_2f_f(mui_uniface_2f *uniface, float *t) {
    uniface->barrier(*t);
}

void mui_barrier_2fx_f(mui_uniface_2fx *uniface, float *t) {
    uniface->barrier(*t);
}

void mui_barrier_2d_f(mui_uniface_2d *uniface, double *t) {
    uniface->barrier(*t);
}

void mui_barrier_2dx_f(mui_uniface_2dx *uniface, double *t) {
    uniface->barrier(*t);
}

void mui_barrier_2t_f(mui_uniface_2t *uniface, double *t) {
    uniface->barrier(static_cast<mui::mui_f_wrapper_2D::time_type>(*t));
}

// Barrier at two time values
void mui_barrier_2f_pair_f(mui_uniface_2f *uniface, float *t_1, float *t_2) {
    uniface->barrier(*t_1, *t_2);
}

void mui_barrier_2fx_pair_f(mui_uniface_2fx *uniface, float *t_1, float *t_2) {
    uniface->barrier(*t_1, *t_2);
}

void mui_barrier_2d_pair_f(mui_uniface_2d *uniface, double *t_1, double *t_2) {
    uniface->barrier(*t_1, *t_2);
}

void mui_barrier_2dx_pair_f(mui_uniface_2dx *uniface, double *t_1, double *t_2) {
    uniface->barrier(*t_1, *t_2);
}

void mui_barrier_2t_pair_f(mui_uniface_2t *uniface, double *t_1, double *t_2) {
    uniface->barrier(static_cast<mui::mui_f_wrapper_2D::time_type>(*t_1),
            static_cast<mui::mui_f_wrapper_2D::time_type>(*t_2));
}

/******************************************
 * MUI forget functions                    *
 ******************************************/

// Forget log between [-inf, upper]
void mui_forget_upper_2f_f(mui_uniface_2f *uniface, float *upper, int *reset_log) {
    uniface->forget(*upper, static_cast<bool>(*reset_log));
}

void mui_forget_upper_2fx_f(mui_uniface_2fx *uniface, float *upper, int *reset_log) {
    uniface->forget(*upper, static_cast<bool>(*reset_log));
}

void mui_forget_upper_2d_f(mui_uniface_2d *uniface, double *upper, int *reset_log) {
    uniface->forget(*upper, static_cast<bool>(*reset_log));
}

void mui_forget_upper_2dx_f(mui_uniface_2dx *uniface, double *upper, int *reset_log) {
    uniface->forget(*upper, static_cast<bool>(*reset_log));
}

void mui_forget_upper_2t_f(mui_uniface_2t *uniface, double *upper, int *reset_log) {
    uniface->forget(static_cast<mui::mui_f_wrapper_2D::time_type>(*upper), static_cast<bool>(*reset_log));
}

// Forget log between [-inf, -inf], [upper_1, upper_2]
void mui_forget_upper_2f_pair_f(mui_uniface_2f *uniface, float *upper_1, float *upper_2, int *reset_log) {
    std::pair<float, float> forget_time(*upper_1, *upper_2);
    uniface->forget(forget_time, static_cast<bool>(*reset_log));
}

void mui_forget_upper_2fx_pair_f(mui_uniface_2fx *uniface, float *upper_1, float *upper_2, int *reset_log) {
    std::pair<float, float> forget_time(*upper_1, *upper_2);
    uniface->forget(forget_time, static_cast<bool>(*reset_log));
}

void mui_forget_upper_2d_pair_f(mui_uniface_2d *uniface, double *upper_1, double *upper_2, int *reset_log) {
    std::pair<double, double> forget_time(*upper_1, *upper_2);
    uniface->forget(forget_time, static_cast<bool>(*reset_log));
}

void mui_forget_upper_2dx_pair_f(mui_uniface_2dx *uniface, double *upper_1, double *upper_2, int *reset_log) {
    std::pair<double, double> forget_time(*upper_1, *upper_2);
    uniface->forget(forget_time, static_cast<bool>(*reset_log));
}

void mui_forget_upper_2t_pair_f(mui_uniface_2t *uniface, double *upper_1, double *upper_2, int *reset_log) {
    mui::mui_f_wrapper_2D::time_type time_1 = static_cast<mui::mui_f_wrapper_2D::time_type>(*upper_1);
    mui::mui_f_wrapper_2D::time_type time_2 = static_cast<mui::mui_f_wrapper_2D::time_type>(*upper_2);
    std::pair<mui::mui_f_wrapper_2D::time_type, mui::mui_f_wrapper_2D::time_type> forget_time(time_1, time_2);
    uniface->forget(forget_time, static_cast<bool>(*reset_log));
}

// Forget log between [lower, upper]
void mui_forget_lower_upper_2f_f(mui_uniface_2f *uniface, float *lower, float *upper, int *reset_log) {
    uniface->forget(*lower, *upper, static_cast<bool>(*reset_log));
}

void mui_forget_lower_upper_2fx_f(mui_uniface_2fx *uniface, float *lower, float *upper, int *reset_log) {
    uniface->forget(*lower, *upper, static_cast<bool>(*reset_log));
}

void mui_forget_lower_upper_2d_f(mui_uniface_2d *uniface, double *lower, double *upper, int *reset_log) {
    uniface->forget(*lower, *upper, static_cast<bool>(*reset_log));
}

void mui_forget_lower_upper_2dx_f(mui_uniface_2dx *uniface, double *lower, double *upper, int *reset_log) {
    uniface->forget(*lower, *upper, static_cast<bool>(*reset_log));
}

void mui_forget_lower_upper_2t_f(mui_uniface_2t *uniface, double *lower, double *upper, int *reset_log) {
    mui::mui_f_wrapper_2D::time_type forget_time_1 = static_cast<mui::mui_f_wrapper_2D::time_type>(*lower);
    mui::mui_f_wrapper_2D::time_type forget_time_2 = static_cast<mui::mui_f_wrapper_2D::time_type>(*upper);
    uniface->forget(forget_time_1, forget_time_2, static_cast<bool>(*reset_log));
}

// Forget log between [lower_1, lower_2], [upper_1, upper_2]
void mui_forget_lower_upper_2f_pair_f(mui_uniface_2f *uniface, float *lower_1, float *lower_2, float *upper_1, float *upper_2,
        int *reset_log) {
    std::pair<float, float> forget_time_1(*lower_1, *lower_2);
    std::pair<float, float> forget_time_2(*upper_1, *upper_2);
    uniface->forget(forget_time_1, forget_time_2, static_cast<bool>(*reset_log));
}

void mui_forget_lower_upper_2fx_pair_f(mui_uniface_2fx *uniface, float *lower_1, float *lower_2, float *upper_1,
        float *upper_2, int *reset_log) {
    std::pair<float, float> forget_time_1(*lower_1, *lower_2);
    std::pair<float, float> forget_time_2(*upper_1, *upper_2);
    uniface->forget(forget_time_1, forget_time_2, static_cast<bool>(*reset_log));
}

void mui_forget_lower_upper_2d_pair_f(mui_uniface_2d *uniface, double *lower_1, double *lower_2, double *upper_1,
        double *upper_2, int *reset_log) {
    std::pair<double, double> forget_time_1(*lower_1, *lower_2);
    std::pair<double, double> forget_time_2(*upper_1, *upper_2);
    uniface->forget(forget_time_1, forget_time_2, static_cast<bool>(*reset_log));
}

void mui_forget_lower_upper_2dx_pair_f(mui_uniface_2dx *uniface, double *lower_1, double *lower_2, double *upper_1,
        double *upper_2, int *reset_log) {
    std::pair<double, double> forget_time_1(*lower_1, *lower_2);
    std::pair<double, double> forget_time_2(*upper_1, *upper_2);
    uniface->forget(forget_time_1, forget_time_2, static_cast<bool>(*reset_log));
}

void mui_forget_lower_upper_2t_pair_f(mui_uniface_2t *uniface, double *lower_1, double *lower_2, double *upper_1,
        double *upper_2, int *reset_log) {
    mui::mui_f_wrapper_2D::time_type time_1 = static_cast<mui::mui_f_wrapper_2D::time_type>(*lower_1);
    mui::mui_f_wrapper_2D::time_type time_2 = static_cast<mui::mui_f_wrapper_2D::time_type>(*lower_2);
    mui::mui_f_wrapper_2D::time_type time_3 = static_cast<mui::mui_f_wrapper_2D::time_type>(*upper_1);
    mui::mui_f_wrapper_2D::time_type time_4 = static_cast<mui::mui_f_wrapper_2D::time_type>(*upper_2);
    std::pair<mui::mui_f_wrapper_2D::time_type, mui::mui_f_wrapper_2D::time_type> forget_time_1(time_1, time_2);
    std::pair<mui::mui_f_wrapper_2D::time_type, mui::mui_f_wrapper_2D::time_type> forget_time_2(time_3, time_4);
    uniface->forget(forget_time_1, forget_time_2, static_cast<bool>(*reset_log));
}

// Set to forget log between [-inf, current-length] automatically
void mui_set_forget_length_2f_f(mui_uniface_2f *uniface, float *length) {
    uniface->set_memory(*length);
}

void mui_set_forget_length_2fx_f(mui_uniface_2fx *uniface, float *length) {
    uniface->set_memory(*length);
}

void mui_set_forget_length_2d_f(mui_uniface_2d *uniface, double *length) {
    uniface->set_memory(*length);
}

void mui_set_forget_length_2dx_f(mui_uniface_2dx *uniface, double *length) {
    uniface->set_memory(*length);
}

void mui_set_forget_length_2t_f(mui_uniface_2t *uniface, double *length) {
    uniface->set_memory(*length);
}

/******************************************
 * MUI URI functions                      *
 ******************************************/

// Obtain original URI host value from existing interface
void mui_uri_host_2f_f(mui_uniface_2f *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_host().c_str());
}

void mui_uri_host_2fx_f(mui_uniface_2fx *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_host().c_str());
}

void mui_uri_host_2d_f(mui_uniface_2d *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_host().c_str());
}

void mui_uri_host_2dx_f(mui_uniface_2dx *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_host().c_str());
}

void mui_uri_host_2t_f(mui_uniface_2t *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_host().c_str());
}

// Obtain original URI path value from existing interface
void mui_uri_path_2f_f(mui_uniface_2f *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_path().c_str());
}

void mui_uri_path_2fx_f(mui_uniface_2fx *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_path().c_str());
}

void mui_uri_path_2d_f(mui_uniface_2d *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_path().c_str());
}

void mui_uri_path_2dx_f(mui_uniface_2dx *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_path().c_str());
}

void mui_uri_path_2t_f(mui_uniface_2t *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_path().c_str());
}

// Obtain original URI protocol value from existing interface
void mui_uri_protocol_2f_f(mui_uniface_2f *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_protocol().c_str());
}

void mui_uri_protocol_2fx_f(mui_uniface_2fx *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_protocol().c_str());
}

void mui_uri_protocol_2d_f(mui_uniface_2d *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_protocol().c_str());
}

void mui_uri_protocol_2dx_f(mui_uniface_2dx *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_protocol().c_str());
}

void mui_uri_protocol_2t_f(mui_uniface_2t *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_protocol().c_str());
}

}
