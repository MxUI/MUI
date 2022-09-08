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
 * @file mui_f_wrapper_1d.cpp
 * @author S. M. Longshaw (derived from original 3D wrapper by S. Kudo)
 * @date Aug 11, 2021
 * @brief C interface for Fortran wrapper to create and manage 1D MUI interfaces
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
typedef mui::uniface1f mui_uniface_1f;
typedef mui::uniface1fx mui_uniface_1fx;
typedef mui::uniface1d mui_uniface_1d;
typedef mui::uniface1dx mui_uniface_1dx;

// MUI Interface typedefs for template creation (recommended)
typedef mui::uniface<mui::mui_f_wrapper_1D> mui_uniface_1t;

// Exact spatial sampler typedefs for specialism creation
typedef mui::sampler_exact1f<float> mui_sampler_exact_1f;
typedef mui::sampler_exact1fx<float> mui_sampler_exact_1fx;
typedef mui::sampler_exact1d<double> mui_sampler_exact_1d;
typedef mui::sampler_exact1dx<double> mui_sampler_exact_1dx;

// Exact spatial sampler typedef for template creation (recommended)
typedef mui::sampler_exact<mui::mui_f_wrapper_1D> mui_sampler_exact_1t;

// Gaussian spatial sampler typedefs for specialism creation
typedef mui::sampler_gauss1f<float> mui_sampler_gauss_1f;
typedef mui::sampler_gauss1fx<float> mui_sampler_gauss_1fx;
typedef mui::sampler_gauss1d<double> mui_sampler_gauss_1d;
typedef mui::sampler_gauss1dx<double> mui_sampler_gauss_1dx;

// Gaussian spatial sampler typedef for template creation (recommended)
typedef mui::sampler_gauss<mui::mui_f_wrapper_1D> mui_sampler_gauss_1t;

// Moving average spatial sampler typedefs for specialism creation
typedef mui::sampler_moving_average1f<float> mui_sampler_moving_average_1f;
typedef mui::sampler_moving_average1fx<float> mui_sampler_moving_average_1fx;
typedef mui::sampler_moving_average1d<double> mui_sampler_moving_average_1d;
typedef mui::sampler_moving_average1dx<double> mui_sampler_moving_average_1dx;

// Moving average spatial sampler typedef for template creation (recommended)
typedef mui::sampler_moving_average<mui::mui_f_wrapper_1D> mui_sampler_moving_average_1t;

// Nearest neighbor spatial sampler typedefs for specialism creation
typedef mui::sampler_nearest_neighbor1f<float> mui_sampler_nearest_neighbor_1f;
typedef mui::sampler_nearest_neighbor1fx<float> mui_sampler_nearest_neighbor_1fx;
typedef mui::sampler_nearest_neighbor1d<double> mui_sampler_nearest_neighbor_1d;
typedef mui::sampler_nearest_neighbor1dx<double> mui_sampler_nearest_neighbor_1dx;

// Moving average spatial sampler typedef for template creation (recommended)
typedef mui::sampler_nearest_neighbor<mui::mui_f_wrapper_1D> mui_sampler_nearest_neighbor_1t;

// Pseudo-linear n^2 spatial sampler typedefs for specialism creation
typedef mui::sampler_pseudo_n2_linear1f<float> mui_sampler_pseudo_n2_linear_1f;
typedef mui::sampler_pseudo_n2_linear1fx<float> mui_sampler_pseudo_n2_linear_1fx;
typedef mui::sampler_pseudo_n2_linear1d<double> mui_sampler_pseudo_n2_linear_1d;
typedef mui::sampler_pseudo_n2_linear1dx<double> mui_sampler_pseudo_n2_linear_1dx;

// Pseudo-linear n^2 spatial sampler typedef for template creation (recommended)
typedef mui::sampler_pseudo_n2_linear<mui::mui_f_wrapper_1D> mui_sampler_pseudo_n2_linear_1t;

// Pseudo-nearest neighbor spatial sampler typedefs for specialism creation
typedef mui::sampler_pseudo_nearest_neighbor1f<float> mui_sampler_pseudo_nearest_neighbor_1f;
typedef mui::sampler_pseudo_nearest_neighbor1fx<float> mui_sampler_pseudo_nearest_neighbor_1fx;
typedef mui::sampler_pseudo_nearest_neighbor1d<double> mui_sampler_pseudo_nearest_neighbor_1d;
typedef mui::sampler_pseudo_nearest_neighbor1dx<double> mui_sampler_pseudo_nearest_neighbor_1dx;

// Pseudo-nearest neighbor spatial sampler typedef for template creation (recommended)
typedef mui::sampler_pseudo_nearest_neighbor<mui::mui_f_wrapper_1D> mui_sampler_pseudo_nearest_neighbor_1t;

// Shepard interpolation with quintic kernel spatial sampler typedefs for specialism creation
typedef mui::sampler_shepard_quintic1f<float> mui_sampler_shepard_quintic_1f;
typedef mui::sampler_shepard_quintic1fx<float> mui_sampler_shepard_quintic_1fx;
typedef mui::sampler_shepard_quintic1d<double> mui_sampler_shepard_quintic_1d;
typedef mui::sampler_shepard_quintic1dx<double> mui_sampler_shepard_quintic_1dx;

// Shepard interpolation with quintic kernel spatial sampler typedef for template creation (recommended)
typedef mui::sampler_shepard_quintic<mui::mui_f_wrapper_1D> mui_sampler_shepard_quintic_1t;

// Spatial sampler that provides a value at a point using a Smoothed Particle Hydrodynamics (SPH) derived
// interpolation method with a quintic spline kernel spatial sampler typedefs for specialism creation
typedef mui::sampler_sph_quintic1f<float> mui_sampler_sph_quintic_1f;
typedef mui::sampler_sph_quintic1fx<float> mui_sampler_sph_quintic_1fx;
typedef mui::sampler_sph_quintic1d<double> mui_sampler_sph_quintic_1d;
typedef mui::sampler_sph_quintic1dx<double> mui_sampler_sph_quintic_1dx;

// Spatial sampler that provides a value at a point using a Smoothed Particle Hydrodynamics (SPH) derived
// interpolation method with a quintic spline kernel spatial sampler typedefs for template creation (recommended)
typedef mui::sampler_sph_quintic<mui::mui_f_wrapper_1D> mui_sampler_sph_quintic_1t;

// Summation with quintic kernel spatial sampler typedefs for specialism creation
typedef mui::sampler_sum_quintic1f<float> mui_sampler_sum_quintic_1f;
typedef mui::sampler_sum_quintic1fx<float> mui_sampler_sum_quintic_1fx;
typedef mui::sampler_sum_quintic1d<double> mui_sampler_sum_quintic_1d;
typedef mui::sampler_sum_quintic1dx<double> mui_sampler_sum_quintic_1dx;

// Summation with quintic kernel spatial sampler typedef for template creation (recommended)
typedef mui::sampler_sum_quintic<mui::mui_f_wrapper_1D> mui_sampler_sum_quintic_1t;

#ifdef USE_RBF
// Radial Basis Function (RBF) spatial sampler typedefs for specialism creation
typedef mui::sampler_rbf1f<float> mui_sampler_rbf_1f;
typedef mui::sampler_rbf1fx<float> mui_sampler_rbf_1fx;
typedef mui::sampler_rbf1d<double> mui_sampler_rbf_1d;
typedef mui::sampler_rbf1dx<double> mui_sampler_rbf_1dx;

// Radial Basis Function (RBF) spatial sampler typedef for template creation (recommended)
typedef mui::sampler_rbf<mui::mui_f_wrapper_1D> mui_sampler_rbf_1t;
#endif

// Exact temporal sampler typedefs for specialism creation
typedef mui::chrono_sampler_exact1f mui_chrono_sampler_exact_1f;
typedef mui::chrono_sampler_exact1fx mui_chrono_sampler_exact_1fx;
typedef mui::chrono_sampler_exact1d mui_chrono_sampler_exact_1d;
typedef mui::chrono_sampler_exact1dx mui_chrono_sampler_exact_1dx;

// Exact temporal sampler typedefs for template creation (recommended)
typedef mui::chrono_sampler_exact<mui::mui_f_wrapper_1D> mui_chrono_sampler_exact_1t;

// Gaussian temporal sampler typedefs for specialism creation
typedef mui::chrono_sampler_gauss1f mui_chrono_sampler_gauss_1f;
typedef mui::chrono_sampler_gauss1fx mui_chrono_sampler_gauss_1fx;
typedef mui::chrono_sampler_gauss1d mui_chrono_sampler_gauss_1d;
typedef mui::chrono_sampler_gauss1dx mui_chrono_sampler_gauss_1dx;

// Gaussian temporal sampler typedefs for template creation (recommended)
typedef mui::chrono_sampler_gauss<mui::mui_f_wrapper_1D> mui_chrono_sampler_gauss_1t;

// Mean average temporal sampler typedefs for specialism creation
typedef mui::chrono_sampler_mean1f mui_chrono_sampler_mean_1f;
typedef mui::chrono_sampler_mean1fx mui_chrono_sampler_mean_1fx;
typedef mui::chrono_sampler_mean1d mui_chrono_sampler_mean_1d;
typedef mui::chrono_sampler_mean1dx mui_chrono_sampler_mean_1dx;

// Mean average temporal sampler typedefs for template creation (recommended)
typedef mui::chrono_sampler_mean<mui::mui_f_wrapper_1D> mui_chrono_sampler_mean_1t;

// Summation temporal sampler typedefs for specialism creation
typedef mui::chrono_sampler_sum1f mui_chrono_sampler_sum_1f;
typedef mui::chrono_sampler_sum1fx mui_chrono_sampler_sum_1fx;
typedef mui::chrono_sampler_sum1d mui_chrono_sampler_sum_1d;
typedef mui::chrono_sampler_sum1dx mui_chrono_sampler_sum_1dx;

// Summation temporal sampler typedefs for template creation (recommended)
typedef mui::chrono_sampler_sum<mui::mui_f_wrapper_1D> mui_chrono_sampler_sum_1t;

// MUI set of specialism interface on multi-domain coupling
mui_uniface_1f** mui_uniface_multi_1f;
mui_uniface_1fx** mui_uniface_multi_1fx;
mui_uniface_1d** mui_uniface_multi_1d;
mui_uniface_1dx** mui_uniface_multi_1dx;

// MUI set of template interface on multi-domain coupling
mui_uniface_1t** mui_uniface_multi_1t;

/****************************************
 * Create MUI interfaces                 *
 ****************************************/

// 1D interface with float=single and int=int32
void mui_create_uniface_1f_f(mui_uniface_1f **ret, const char *URI) {
    *ret = new mui_uniface_1f(URI);
}

// 1D interface with float=single and int=int64
void mui_create_uniface_1fx_f(mui_uniface_1fx **ret, const char *URI) {
    *ret = new mui_uniface_1fx(URI);
}

// 1D interface with float=double and int=int32
void mui_create_uniface_1d_f(mui_uniface_1d **ret, const char *URI) {
    *ret = new mui_uniface_1d(URI);
}

// 1D interface with float=double and int=int64
void mui_create_uniface_1dx_f(mui_uniface_1dx **ret, const char *URI) {
    *ret = new mui_uniface_1dx(URI);
}

// 1D interface using config from config_f_wrapper.h
void mui_create_uniface_1t_f(mui_uniface_1t **ret, const char *URI) {
    *ret = new mui_uniface_1t(URI);
}

// Set of 1D interfaces with float=single and int=int32
void mui_create_uniface_multi_1f_f(const char *domain, const char **interfaces, int interface_count ) {

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
                std::cerr << "MUI Error [mui_f_wrapper_1d.cpp]: Error MUI Fortran C binding: NULL interfaces at the " << i << "th interface." << std::endl;
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
        std::cerr << "MUI Error [mui_f_wrapper_1d.cpp]: Error MUI Fortran C binding: The size of interfaces array must larger or equals to one." << std::endl;
        std::abort();
    }

    auto created_unifaces = mui::create_uniface<mui::config_1f>(domain, interface_names);

    mui_uniface_multi_1f = new mui_uniface_1f*[created_unifaces.size()];

      for(size_t i=0; i<created_unifaces.size(); i++) {
        mui_uniface_multi_1f[i] = created_unifaces[i].release();
    }
}

// Set of 1D interfaces with float=single and int=int64
void mui_create_uniface_multi_1fx_f(const char *domain, const char **interfaces, int interface_count ) {

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
                std::cerr << "MUI Error [mui_f_wrapper_1d.cpp]: Error MUI Fortran C binding: NULL interfaces at the " << i << "th interface." << std::endl;
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
        std::cerr << "MUI Error [mui_f_wrapper_1d.cpp]: Error MUI Fortran C binding: The size of interfaces array must larger or equals to one." << std::endl;
        std::abort();
    }

    auto created_unifaces = mui::create_uniface<mui::config_1fx>(domain, interface_names);

    mui_uniface_multi_1fx = new mui_uniface_1fx*[created_unifaces.size()];

      for(size_t i=0; i<created_unifaces.size(); i++) {
        mui_uniface_multi_1fx[i] = created_unifaces[i].release();
    }
}

// Set of 1D interfaces with float=double and int=int32
void mui_create_uniface_multi_1d_f(const char *domain, const char **interfaces, int interface_count ) {

    std::vector<std::string> interface_names;
    size_t interface_element_size = sizeof(interfaces)/sizeof(interfaces[0]);
std::cout << "interfaces[0]: " << interfaces[0] << std::endl;
    if (interface_element_size == 1) {
        char interfaces_jointed[strlen(interfaces[0])];
        strcpy(interfaces_jointed, interfaces[0]);
        char delim[] = " ";

        char * pchar = strtok (interfaces_jointed,delim);
        int interface_names_count = 0;
        for(size_t i=0; i<interface_count; i++) {
            if(pchar == NULL) {
                std::cerr << "MUI Error [mui_f_wrapper_1d.cpp]: Error MUI Fortran C binding: NULL interfaces at the " << i << "th interface." << std::endl;
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
        std::cerr << "MUI Error [mui_f_wrapper_1d.cpp]: Error MUI Fortran C binding: The size of interfaces array must larger or equals to one." << std::endl;
        std::abort();
    }

    auto created_unifaces = mui::create_uniface<mui::config_1d>(domain, interface_names);

    mui_uniface_multi_1d = new mui_uniface_1d*[created_unifaces.size()];

      for(size_t i=0; i<created_unifaces.size(); i++) {
        mui_uniface_multi_1d[i] = created_unifaces[i].release();
    }
}

// Set of 1D interfaces with float=double and int=int64
void mui_create_uniface_multi_1dx_f(const char *domain, const char **interfaces, int interface_count ) {

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
                std::cerr << "MUI Error [mui_f_wrapper_1d.cpp]: Error MUI Fortran C binding: NULL interfaces at the " << i << "th interface." << std::endl;
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
        std::cerr << "MUI Error [mui_f_wrapper_1d.cpp]: Error MUI Fortran C binding: The size of interfaces array must larger or equals to one." << std::endl;
        std::abort();
    }

    auto created_unifaces = mui::create_uniface<mui::config_1dx>(domain, interface_names);

    mui_uniface_multi_1dx = new mui_uniface_1dx*[created_unifaces.size()];

      for(size_t i=0; i<created_unifaces.size(); i++) {
        mui_uniface_multi_1dx[i] = created_unifaces[i].release();
    }
}

// Set of 1D interfaces using config from config_c_wrapper.h
void mui_create_uniface_multi_1t_f(const char *domain, const char **interfaces, int interface_count ) {

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
                std::cerr << "MUI Error [mui_f_wrapper_1d.cpp]: Error MUI Fortran C binding: NULL interfaces at the " << i << "th interface." << std::endl;
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
        std::cerr << "MUI Error [mui_f_wrapper_1d.cpp]: Error MUI Fortran C binding: The size of interfaces array must larger or equals to one." << std::endl;
        std::abort();
    }

    auto created_unifaces = mui::create_uniface<mui::mui_f_wrapper_1D>(domain, interface_names);

    mui_uniface_multi_1t = new mui_uniface_1t*[created_unifaces.size()];

      for(size_t i=0; i<created_unifaces.size(); i++) {
        mui_uniface_multi_1t[i] = created_unifaces[i].release();
    }
}

// Access to MUI set of 1D interfaces with float=single and int=int32
mui_uniface_1f* get_mui_uniface_multi_1f_f(int interface_count) {
    return mui_uniface_multi_1f[interface_count-1];
}

// Access to MUI set of 1D interfaces with float=single and int=int64
mui_uniface_1fx* get_mui_uniface_multi_1fx_f(int interface_count) {
    return mui_uniface_multi_1fx[interface_count-1];
}

// Access to MUI set of 1D interfaces with float=double and int=int32
mui_uniface_1d* get_mui_uniface_multi_1d_f(int interface_count) {
    return mui_uniface_multi_1d[interface_count-1];
}

// Access to MUI set of 1D interfaces with float=double and int=int64
mui_uniface_1dx* get_mui_uniface_multi_1dx_f(int interface_count) {
    return mui_uniface_multi_1dx[interface_count-1];
}

// Access to MUI set of 1D interfaces using config from config_f_wrapper.h
mui_uniface_1t* get_mui_uniface_multi_1t_f(int interface_count) {
    return mui_uniface_multi_1t[interface_count-1];
}

/****************************************
 * Destroy MUI interface                 *
 ****************************************/

void mui_destroy_uniface_1f_f(mui_uniface_1f *uniface) {
    delete uniface;
}

void mui_destroy_uniface_1fx_f(mui_uniface_1fx *uniface) {
    delete uniface;
}

void mui_destroy_uniface_1d_f(mui_uniface_1d *uniface) {
    delete uniface;
}

void mui_destroy_uniface_1dx_f(mui_uniface_1dx *uniface) {
    delete uniface;
}

void mui_destroy_uniface_1t_f(mui_uniface_1t *uniface) {
    delete uniface;
}

/******************************************
 * Create 1D spatial samplers              *
 ******************************************/

// Exact sampler
void mui_create_sampler_exact_1f_f(mui_sampler_exact_1f **ret, float *tolerance) {
    *ret = new mui_sampler_exact_1f(*tolerance);
}

void mui_create_sampler_exact_1fx_f(mui_sampler_exact_1fx **ret, float *tolerance) {
    *ret = new mui_sampler_exact_1fx(*tolerance);
}

void mui_create_sampler_exact_1d_f(mui_sampler_exact_1d **ret, double *tolerance) {
    *ret = new mui_sampler_exact_1d(*tolerance);
}

void mui_create_sampler_exact_1dx_f(mui_sampler_exact_1dx** ret, double *tolerance) {
    *ret = new mui_sampler_exact_1dx(*tolerance);
}

void mui_create_sampler_exact_1t_f(mui_sampler_exact_1t** ret, double* tolerance) {
    *ret = new mui_sampler_exact_1t(static_cast<mui::mui_f_wrapper_1D::REAL>(*tolerance));
}

// Gauss sampler
void mui_create_sampler_gauss_1f_f(mui_sampler_gauss_1f **ret, float* r, float* h) {
    *ret = new mui_sampler_gauss_1f(*r, *h);
}

void mui_create_sampler_gauss_1fx_f(mui_sampler_gauss_1fx **ret, float* r, float* h) {
    *ret = new mui_sampler_gauss_1fx(*r, *h);
}

void mui_create_sampler_gauss_1d_f(mui_sampler_gauss_1d **ret, double* r, double* h) {
    *ret = new mui_sampler_gauss_1d(*r, *h);
}

void mui_create_sampler_gauss_1dx_f(mui_sampler_gauss_1dx** ret, double* r, double* h) {
    *ret = new mui_sampler_gauss_1dx(*r, *h);
}

void mui_create_sampler_gauss_1t_f(mui_sampler_gauss_1t** ret, double* r, double* h) {
    *ret = new mui_sampler_gauss_1t(static_cast<mui::mui_f_wrapper_1D::REAL>(*r),
            static_cast<mui::mui_f_wrapper_1D::REAL>(*h));
}

// Moving average sampler
void mui_create_sampler_moving_average_1f_f(mui_sampler_moving_average_1f **ret, float* bbox_1) {
    mui::point1f bbox(*bbox_1);
    *ret = new mui_sampler_moving_average_1f(bbox);
}

void mui_create_sampler_moving_average_1fx_f(mui_sampler_moving_average_1fx **ret, float* bbox_1) {
    mui::point1fx bbox(*bbox_1);
    *ret = new mui_sampler_moving_average_1fx(bbox);
}

void mui_create_sampler_moving_average_1d_f(mui_sampler_moving_average_1d **ret, double* bbox_1) {
    mui::point1d bbox(*bbox_1);
    *ret = new mui_sampler_moving_average_1d(bbox);
}

void mui_create_sampler_moving_average_1dx_f(mui_sampler_moving_average_1dx** ret, double* bbox_1) {
    mui::point1dx bbox(*bbox_1);
    *ret = new mui_sampler_moving_average_1dx(bbox);
}

void mui_create_sampler_moving_average_1t_f(mui_sampler_moving_average_1t** ret, double* bbox_1) {
    mui::mui_f_wrapper_1D::point_type bbox(*bbox_1);
    *ret = new mui_sampler_moving_average_1t(bbox);
}

// Nearest neighbour sampler
void mui_create_sampler_nearest_neighbor_1f_f(mui_sampler_nearest_neighbor_1f **ret) {
    *ret = new mui_sampler_nearest_neighbor_1f();
}

void mui_create_sampler_nearest_neighbor_1fx_f(mui_sampler_nearest_neighbor_1fx **ret) {
    *ret = new mui_sampler_nearest_neighbor_1fx();
}

void mui_create_sampler_nearest_neighbor_1d_f(mui_sampler_nearest_neighbor_1d **ret) {
    *ret = new mui_sampler_nearest_neighbor_1d();
}

void mui_create_sampler_nearest_neighbor_1dx_f(mui_sampler_nearest_neighbor_1dx** ret) {
    *ret = new mui_sampler_nearest_neighbor_1dx();
}

void mui_create_sampler_nearest_neighbor_1t_f(mui_sampler_nearest_neighbor_1t** ret) {
    *ret = new mui_sampler_nearest_neighbor_1t();
}

// Pseudo-linear n^2 interpolation sampler
void mui_create_sampler_pseudo_n2_linear_1f_f(mui_sampler_pseudo_n2_linear_1f **ret, float* r) {
    *ret = new mui_sampler_pseudo_n2_linear_1f(*r);
}

void mui_create_sampler_pseudo_n2_linear_1fx_f(mui_sampler_pseudo_n2_linear_1fx **ret, float* r) {
    *ret = new mui_sampler_pseudo_n2_linear_1fx(*r);
}

void mui_create_sampler_pseudo_n2_linear_1d_f(mui_sampler_pseudo_n2_linear_1d **ret, double* r) {
    *ret = new mui_sampler_pseudo_n2_linear_1d(*r);
}

void mui_create_sampler_pseudo_n2_linear_1dx_f(mui_sampler_pseudo_n2_linear_1dx** ret, double* r) {
    *ret = new mui_sampler_pseudo_n2_linear_1dx(*r);
}

void mui_create_sampler_pseudo_n2_linear_1t_f(mui_sampler_pseudo_n2_linear_1t** ret, double* r) {
    *ret = new mui_sampler_pseudo_n2_linear_1t(static_cast<mui::mui_f_wrapper_1D::REAL>(*r));
}

// Pseudo-linear nearest neighbour interpolation sampler
void mui_create_sampler_pseudo_nearest_neighbor_1f_f(mui_sampler_pseudo_nearest_neighbor_1f **ret, float* h) {
    *ret = new mui_sampler_pseudo_nearest_neighbor_1f(*h);
}

void mui_create_sampler_pseudo_nearest_neighbor_1fx_f(mui_sampler_pseudo_nearest_neighbor_1fx **ret, float* h) {
    *ret = new mui_sampler_pseudo_nearest_neighbor_1fx(*h);
}

void mui_create_sampler_pseudo_nearest_neighbor_1d_f(mui_sampler_pseudo_nearest_neighbor_1d **ret, double* h) {
    *ret = new mui_sampler_pseudo_nearest_neighbor_1d(*h);
}

void mui_create_sampler_pseudo_nearest_neighbor_1dx_f(mui_sampler_pseudo_nearest_neighbor_1dx** ret, double* h) {
    *ret = new mui_sampler_pseudo_nearest_neighbor_1dx(*h);
}

void mui_create_sampler_pseudo_nearest_neighbor_1t_f(mui_sampler_pseudo_nearest_neighbor_1t** ret, double* h) {
    *ret = new mui_sampler_pseudo_nearest_neighbor_1t(static_cast<mui::mui_f_wrapper_1D::REAL>(*h));
}

// Shepard interpolation with a quintic kernel sampler
void mui_create_sampler_shepard_quintic_1f_f(mui_sampler_shepard_quintic_1f **ret, float* r) {
    *ret = new mui_sampler_shepard_quintic_1f(*r);
}

void mui_create_sampler_shepard_quintic_1fx_f(mui_sampler_shepard_quintic_1fx **ret, float* r) {
    *ret = new mui_sampler_shepard_quintic_1fx(*r);
}

void mui_create_sampler_shepard_quintic_1d_f(mui_sampler_shepard_quintic_1d **ret, double* r) {
    *ret = new mui_sampler_shepard_quintic_1d(*r);
}

void mui_create_sampler_shepard_quintic_1dx_f(mui_sampler_shepard_quintic_1dx** ret, double* r) {
    *ret = new mui_sampler_shepard_quintic_1dx(*r);
}

void mui_create_sampler_shepard_quintic_1t_f(mui_sampler_shepard_quintic_1t** ret, double* r) {
    *ret = new mui_sampler_shepard_quintic_1t(static_cast<mui::mui_f_wrapper_1D::REAL>(*r));
}

// SPH derived interpolation method with a quintic spline kernel sampler
void mui_create_sampler_sph_quintic_1f_f(mui_sampler_sph_quintic_1f **ret, float* r) {
    *ret = new mui_sampler_sph_quintic_1f(*r);
}

void mui_create_sampler_sph_quintic_1fx_f(mui_sampler_sph_quintic_1fx **ret, float* r) {
    *ret = new mui_sampler_sph_quintic_1fx(*r);
}

void mui_create_sampler_sph_quintic_1d_f(mui_sampler_sph_quintic_1d **ret, double* r) {
    *ret = new mui_sampler_sph_quintic_1d(*r);
}

void mui_create_sampler_sph_quintic_1dx_f(mui_sampler_sph_quintic_1dx** ret, double* r) {
    *ret = new mui_sampler_sph_quintic_1dx(*r);
}

void mui_create_sampler_sph_quintic_1t_f(mui_sampler_sph_quintic_1t** ret, double* r) {
    *ret = new mui_sampler_sph_quintic_1t(static_cast<mui::mui_f_wrapper_1D::REAL>(*r));
}

// Summation with a quintic kernel sampler
void mui_create_sampler_sum_quintic_1f_f(mui_sampler_sum_quintic_1f **ret, float* r) {
    *ret = new mui_sampler_sum_quintic_1f(*r);
}

void mui_create_sampler_sum_quintic_1fx_f(mui_sampler_sum_quintic_1fx **ret, float* r) {
    *ret = new mui_sampler_sum_quintic_1fx(*r);
}

void mui_create_sampler_sum_quintic_1d_f(mui_sampler_sum_quintic_1d **ret, double* r) {
    *ret = new mui_sampler_sum_quintic_1d(*r);
}

void mui_create_sampler_sum_quintic_1dx_f(mui_sampler_sum_quintic_1dx** ret, double* r) {
    *ret = new mui_sampler_sum_quintic_1dx(*r);
}

void mui_create_sampler_sum_quintic_1t_f(mui_sampler_sum_quintic_1t** ret, double* r) {
    *ret = new mui_sampler_sum_quintic_1t(static_cast<mui::mui_f_wrapper_1D::REAL>(*r));
}

#ifdef USE_RBF
// Radial Basis Function sampler
void mui_create_sampler_rbf_1f_f(mui_sampler_rbf_1f **ret, float* r, float* points_1, int* points_count, int* basis_func,
        int* conservative, int* polynomial, int* smoothFunc, int* readMatrix, int* writeMatrix, const char* file_address,
        float* cutoff, float* cg_solve_tol, int* cg_solve_it, int* pou_size) {
    std::vector<mui::point1f> pts(*points_count);
    for (size_t i = 0; i < *points_count; i++) {
        pts[i][0] = points_1[i];
    }

    *ret = new mui_sampler_rbf_1f(*r, pts, *basis_func, static_cast<bool>(*conservative), static_cast<bool>(*polynomial),
            static_cast<bool>(*smoothFunc), static_cast<bool>(*readMatrix), static_cast<bool>(*writeMatrix), std::string(file_address),
            *cutoff, *cg_solve_tol, *cg_solve_it, *pou_size);
}

void mui_create_sampler_rbf_1fx_f(mui_sampler_rbf_1fx **ret, float* r, float* points_1, int* points_count, int* basis_func,
        int* conservative, int* polynomial, int* smoothFunc, int* readMatrix, int* writeMatrix, const char* file_address,
        float* cutoff, float* cg_solve_tol, int* cg_solve_it, int* pou_size) {
    std::vector<mui::point1fx> pts(*points_count);
    for (size_t i = 0; i < *points_count; i++) {
        pts[i][0] = points_1[i];
    }

    *ret = new mui_sampler_rbf_1fx(*r, pts, *basis_func, static_cast<bool>(*conservative), static_cast<bool>(*polynomial),
            static_cast<bool>(*smoothFunc), static_cast<bool>(*readMatrix), static_cast<bool>(*writeMatrix), std::string(file_address),
            *cutoff, *cg_solve_tol, *cg_solve_it, *pou_size);
}

void mui_create_sampler_rbf_1d_f(mui_sampler_rbf_1d **ret, double* r, double* points_1, int* points_count, int* basis_func,
        int* conservative, int* polynomial, int* smoothFunc, int* readMatrix, int* writeMatrix, const char* file_address,
        double* cutoff, double* cg_solve_tol, int* cg_solve_it, int* pou_size) {
    std::vector<mui::point1d> pts(*points_count);
    for (size_t i = 0; i < *points_count; i++) {
        pts[i][0] = points_1[i];
    }

    *ret = new mui_sampler_rbf_1d(*r, pts, *basis_func, static_cast<bool>(*conservative), static_cast<bool>(*polynomial),
            static_cast<bool>(*smoothFunc), static_cast<bool>(*readMatrix), static_cast<bool>(*writeMatrix), std::string(file_address),
            *cutoff, *cg_solve_tol, *cg_solve_it, *pou_size);
}

void mui_create_sampler_rbf_1dx_f(mui_sampler_rbf_1dx** ret, double* r, double* points_1, int* points_count, int* basis_func,
        int* conservative, int* polynomial, int* smoothFunc, int* readMatrix, int* writeMatrix, const char* file_address,
        double* cutoff, double* cg_solve_tol, int* cg_solve_it, int* pou_size) {
    std::vector<mui::point1dx> pts(*points_count);
    for (size_t i = 0; i < *points_count; i++) {
        pts[i][0] = points_1[i];
    }

    *ret = new mui_sampler_rbf_1dx(*r, pts, *basis_func, static_cast<bool>(*conservative), static_cast<bool>(*polynomial),
            static_cast<bool>(*smoothFunc), static_cast<bool>(*readMatrix), static_cast<bool>(*writeMatrix), std::string(file_address),
            *cutoff, *cg_solve_tol, *cg_solve_it, *pou_size);
}

void mui_create_sampler_rbf_1t_f(mui_sampler_rbf_1t** ret, double* r, double* points_1, int* points_count, int* basis_func,
        int* conservative, int* polynomial, int* smoothFunc, int* readMatrix, int* writeMatrix, const char* file_address,
        double* cutoff, double* cg_solve_tol, int* cg_solve_it, int* pou_size) {
    std::vector<mui::mui_f_wrapper_1D::point_type> pts(*points_count);
    for (size_t i = 0; i < *points_count; i++) {
        pts[i][0] = static_cast<mui::mui_f_wrapper_1D::REAL>(points_1[i]);
    }

    *ret = new mui_sampler_rbf_1t(static_cast<mui::mui_f_wrapper_1D::REAL>(*r), pts, *basis_func,
            static_cast<bool>(*conservative), static_cast<bool>(*polynomial), static_cast<bool>(*smoothFunc),
            static_cast<bool>(*readMatrix), static_cast<bool>(*writeMatrix), std::string(file_address),
            static_cast<mui::mui_f_wrapper_1D::REAL>(*cutoff), static_cast<mui::mui_f_wrapper_1D::REAL>(*cg_solve_tol),
      static_cast<mui::mui_f_wrapper_1D::INT>(*cg_solve_it), static_cast<mui::mui_f_wrapper_1D::INT>(*pou_size));
}
#endif

/*******************************************
 * Destroy 1D spatial samplers              *
 *******************************************/

// Exact sampler
void mui_destroy_sampler_exact_1f_f(mui_sampler_exact_1f* sampler) {
    delete sampler;
}

void mui_destroy_sampler_exact_1fx_f(mui_sampler_exact_1fx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_exact_1d_f(mui_sampler_exact_1d* sampler) {
    delete sampler;
}

void mui_destroy_sampler_exact_1dx_f(mui_sampler_exact_1dx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_exact_1t_f(mui_sampler_exact_1t* sampler) {
    delete sampler;
}

// Gauss sampler
void mui_destroy_sampler_gauss_1f_f(mui_sampler_gauss_1f* sampler) {
    delete sampler;
}

void mui_destroy_sampler_gauss_1fx_f(mui_sampler_gauss_1fx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_gauss_1d_f(mui_sampler_gauss_1d* sampler) {
    delete sampler;
}

void mui_destroy_sampler_gauss_1dx_f(mui_sampler_gauss_1dx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_gauss_1t_f(mui_sampler_gauss_1t* sampler) {
    delete sampler;
}

// Moving average sampler
void mui_destroy_sampler_moving_average_1f_f(mui_sampler_moving_average_1f* sampler) {
    delete sampler;
}

void mui_destroy_sampler_moving_average_1fx_f(mui_sampler_moving_average_1fx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_moving_average_1d_f(mui_sampler_moving_average_1d* sampler) {
    delete sampler;
}

void mui_destroy_sampler_moving_average_1dx_f(mui_sampler_moving_average_1dx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_moving_average_1t_f(mui_sampler_moving_average_1t* sampler) {
    delete sampler;
}

// Nearest neighbour sampler
void mui_destroy_sampler_nearest_neighbor_1f_f(mui_sampler_nearest_neighbor_1f* sampler) {
    delete sampler;
}

void mui_destroy_sampler_nearest_neighbor_1fx_f(mui_sampler_nearest_neighbor_1fx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_nearest_neighbor_1d_f(mui_sampler_nearest_neighbor_1d* sampler) {
    delete sampler;
}

void mui_destroy_sampler_nearest_neighbor_1dx_f(mui_sampler_nearest_neighbor_1dx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_nearest_neighbor_1t_f(mui_sampler_nearest_neighbor_1t* sampler) {
    delete sampler;
}

// Pseudo-linear n^2 interpolation sampler
void mui_destroy_sampler_pseudo_nearest2_linear_1f_f(mui_sampler_pseudo_nearest_neighbor_1f* sampler) {
    delete sampler;
}

void mui_destroy_sampler_pseudo_nearest2_linear_1fx_f(mui_sampler_pseudo_nearest_neighbor_1fx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_pseudo_nearest2_linear_1d_f(mui_sampler_pseudo_nearest_neighbor_1d* sampler) {
    delete sampler;
}

void mui_destroy_sampler_pseudo_nearest2_linear_1dx_f(mui_sampler_pseudo_nearest_neighbor_1dx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_pseudo_nearest2_linear_1t_f(mui_sampler_pseudo_nearest_neighbor_1t* sampler) {
    delete sampler;
}

// Pseudo-linear nearest neighbour interpolation sampler
void mui_destroy_sampler_pseudo_nearest_neighbor_1f_f(mui_sampler_pseudo_nearest_neighbor_1f* sampler) {
    delete sampler;
}

void mui_destroy_sampler_pseudo_nearest_neighbor_1fx_f(mui_sampler_pseudo_nearest_neighbor_1fx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_pseudo_nearest_neighbor_1d_f(mui_sampler_pseudo_nearest_neighbor_1d* sampler) {
    delete sampler;
}

void mui_destroy_sampler_pseudo_nearest_neighbor_1dx_f(mui_sampler_pseudo_nearest_neighbor_1dx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_pseudo_nearest_neighbor_1t_f(mui_sampler_pseudo_nearest_neighbor_1t* sampler) {
    delete sampler;
}

// Shepard interpolation with a quintic kernel sampler
void mui_destroy_sampler_shepard_quintic_1f_f(mui_sampler_shepard_quintic_1f* sampler) {
    delete sampler;
}

void mui_destroy_sampler_shepard_quintic_1fx_f(mui_sampler_shepard_quintic_1fx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_shepard_quintic_1d_f(mui_sampler_shepard_quintic_1d* sampler) {
    delete sampler;
}

void mui_destroy_sampler_shepard_quintic_1dx_f(mui_sampler_shepard_quintic_1dx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_shepard_quintic_1t_f(mui_sampler_shepard_quintic_1t* sampler) {
    delete sampler;
}

// SPH derived interpolation method with a quintic spline kernel sampler
void mui_destroy_sampler_sph_quintic_1f_f(mui_sampler_sph_quintic_1f* sampler) {
    delete sampler;
}

void mui_destroy_sampler_sph_quintic_1fx_f(mui_sampler_sph_quintic_1fx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_sph_quintic_1d_f(mui_sampler_sph_quintic_1d* sampler) {
    delete sampler;
}

void mui_destroy_sampler_sph_quintic_1dx_f(mui_sampler_sph_quintic_1dx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_sph_quintic_1t_f(mui_sampler_sph_quintic_1t* sampler) {
    delete sampler;
}

// Summation with a quintic kernel sampler
void mui_destroy_sampler_sum_quintic_1f_f(mui_sampler_sum_quintic_1f* sampler) {
    delete sampler;
}

void mui_destroy_sampler_sum_quintic_1fx_f(mui_sampler_sum_quintic_1fx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_sum_quintic_1d_f(mui_sampler_sum_quintic_1d* sampler) {
    delete sampler;
}

void mui_destroy_sampler_sum_quintic_1dx_f(mui_sampler_sum_quintic_1dx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_sum_quintic_1t_f(mui_sampler_sum_quintic_1t* sampler) {
    delete sampler;
}

#ifdef USE_RBF
void mui_destroy_sampler_rbf_1f_f(mui_sampler_rbf_1f* sampler) {
    delete sampler;
}

void mui_destroy_sampler_rbf_1fx_f(mui_sampler_rbf_1fx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_rbf_1d_f(mui_sampler_rbf_1d* sampler) {
    delete sampler;
}

void mui_destroy_sampler_rbf_1dx_f(mui_sampler_rbf_1dx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_rbf_1t_f(mui_sampler_rbf_1t* sampler) {
    delete sampler;
}
#endif

/*******************************************
 * Create temporal samplers                 *
 *******************************************/

// Exact temporal sampler
void mui_create_chrono_sampler_exact_1f_f(mui_chrono_sampler_exact_1f **ret, float* tolerance) {
    *ret = new mui_chrono_sampler_exact_1f(*tolerance);
}

void mui_create_chrono_sampler_exact_1fx_f(mui_chrono_sampler_exact_1fx **ret, float* tolerance) {
    *ret = new mui_chrono_sampler_exact_1fx(*tolerance);
}

void mui_create_chrono_sampler_exact_1d_f(mui_chrono_sampler_exact_1d **ret, double* tolerance) {
    *ret = new mui_chrono_sampler_exact_1d(*tolerance);
}

void mui_create_chrono_sampler_exact_1dx_f(mui_chrono_sampler_exact_1dx** ret, double* tolerance) {
    *ret = new mui_chrono_sampler_exact_1dx(*tolerance);
}

void mui_create_chrono_sampler_exact_1t_f(mui_chrono_sampler_exact_1t** ret, double* tolerance) {
    *ret = new mui_chrono_sampler_exact_1t(static_cast<mui::mui_f_wrapper_1D::REAL>(*tolerance));
}

// Gauss temporal sampler
void mui_create_chrono_sampler_gauss_1f_f(mui_chrono_sampler_gauss_1f **ret, float* cutoff, float* sigma) {
    *ret = new mui_chrono_sampler_gauss_1f(*cutoff, *sigma);
}

void mui_create_chrono_sampler_gauss_1fx_f(mui_chrono_sampler_gauss_1fx **ret, float* cutoff, float* sigma) {
    *ret = new mui_chrono_sampler_gauss_1fx(*cutoff, *sigma);
}

void mui_create_chrono_sampler_gauss_1d_f(mui_chrono_sampler_gauss_1d **ret, double* cutoff, double* sigma) {
    *ret = new mui_chrono_sampler_gauss_1d(*cutoff, *sigma);
}

void mui_create_chrono_sampler_gauss_1dx_f(mui_chrono_sampler_gauss_1dx** ret, double* cutoff, double* sigma) {
    *ret = new mui_chrono_sampler_gauss_1dx(*cutoff, *sigma);
}

void mui_create_chrono_sampler_gauss_1t_f(mui_chrono_sampler_gauss_1t** ret, double* cutoff, double* sigma) {
    *ret = new mui_chrono_sampler_gauss_1t(static_cast<mui::mui_f_wrapper_1D::time_type>(*cutoff),
            static_cast<mui::mui_f_wrapper_1D::REAL>(*sigma));
}

// Mean temporal sampler
void mui_create_chrono_sampler_mean_1f_f(mui_chrono_sampler_mean_1f **ret, float* lower, float* upper) {
    *ret = new mui_chrono_sampler_mean_1f(*lower, *upper);
}

void mui_create_chrono_sampler_mean_1fx_f(mui_chrono_sampler_mean_1fx **ret, float* lower, float* upper) {
    *ret = new mui_chrono_sampler_mean_1fx(*lower, *upper);
}

void mui_create_chrono_sampler_mean_1d_f(mui_chrono_sampler_mean_1d **ret, double* lower, double* upper) {
    *ret = new mui_chrono_sampler_mean_1d(*lower, *upper);
}

void mui_create_chrono_sampler_mean_1dx_f(mui_chrono_sampler_mean_1dx** ret, double* lower, double* upper) {
    *ret = new mui_chrono_sampler_mean_1dx(*lower, *upper);
}

void mui_create_chrono_sampler_mean_1t_f(mui_chrono_sampler_mean_1t** ret, double* lower, double* upper) {
    *ret = new mui_chrono_sampler_mean_1t(static_cast<mui::mui_f_wrapper_1D::time_type>(*lower),
            static_cast<mui::mui_f_wrapper_1D::time_type>(*upper));
}

// Sum temporal sampler
void mui_create_chrono_sampler_sum_1f_f(mui_chrono_sampler_sum_1f **ret, float* lower, float* upper) {
    *ret = new mui_chrono_sampler_sum_1f(*lower, *upper);
}

void mui_create_chrono_sampler_sum_1fx_f(mui_chrono_sampler_sum_1fx **ret, float* lower, float* upper) {
    *ret = new mui_chrono_sampler_sum_1fx(*lower, *upper);
}

void mui_create_chrono_sampler_sum_1d_f(mui_chrono_sampler_sum_1d **ret, double* lower, double* upper) {
    *ret = new mui_chrono_sampler_sum_1d(*lower, *upper);
}

void mui_create_chrono_sampler_sum_1dx_f(mui_chrono_sampler_sum_1dx** ret, double* lower, double* upper) {
    *ret = new mui_chrono_sampler_sum_1dx(*lower, *upper);
}

void mui_create_chrono_sampler_sum_1t_f(mui_chrono_sampler_sum_1t** ret, double* lower, double* upper) {
    *ret = new mui_chrono_sampler_sum_1t(static_cast<mui::mui_f_wrapper_1D::time_type>(*lower),
            static_cast<mui::mui_f_wrapper_1D::time_type>(*upper));
}

/*******************************************
 * Destroy temporal samplers                *
 *******************************************/

// Exact temporal sampler
void mui_destroy_chrono_sampler_exact_1f_f(mui_chrono_sampler_exact_1f* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_exact_1fx_f(mui_chrono_sampler_exact_1fx* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_exact_1d_f(mui_chrono_sampler_exact_1d* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_exact_1dx_f(mui_chrono_sampler_exact_1dx* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_exact_1t_f(mui_chrono_sampler_exact_1t* sampler) {
    delete sampler;
}

// Gauss temporal sampler
void mui_destroy_chrono_sampler_gauss_1f_f(mui_chrono_sampler_gauss_1f* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_gauss_1fx_f(mui_chrono_sampler_gauss_1fx* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_gauss_1d_f(mui_chrono_sampler_gauss_1d* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_gauss_1dx_f(mui_chrono_sampler_gauss_1dx* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_gauss_1t_f(mui_chrono_sampler_gauss_1t* sampler) {
    delete sampler;
}

// Mean temporal sampler
void mui_destroy_chrono_sampler_mean_1f_f(mui_chrono_sampler_mean_1f* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_mean_1fx_f(mui_chrono_sampler_mean_1fx* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_mean_1d_f(mui_chrono_sampler_mean_1d* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_mean_1dx_f(mui_chrono_sampler_mean_1dx* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_mean_1t_f(mui_chrono_sampler_mean_1t* sampler) {
    delete sampler;
}

// Sum temporal sampler
void mui_destroy_chrono_sampler_sum_1f_f(mui_chrono_sampler_sum_1f* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_sum_1fx_f(mui_chrono_sampler_sum_1fx* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_sum_1d_f(mui_chrono_sampler_sum_1d* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_sum_1dx_f(mui_chrono_sampler_sum_1dx* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_sum_1t_f(mui_chrono_sampler_sum_1t* sampler) {
    delete sampler;
}

/******************************************
 * MUI functions for data push             *
 ******************************************/

// Standard push functions
void mui_push_1f_f(mui_uniface_1f *uniface, const char* attr, float* point_1, float* value) {
    uniface->push(std::string(attr), mui::point1f(*point_1), *value);
}

void mui_push_1fx_f(mui_uniface_1fx *uniface, const char* attr, float* point_1, float* value) {
    uniface->push(std::string(attr), mui::point1fx(*point_1), *value);
}

void mui_push_1d_f(mui_uniface_1d *uniface, const char* attr, double* point_1, double* value) {
    uniface->push(std::string(attr), mui::point1d(*point_1), *value);
}

void mui_push_1dx_f(mui_uniface_1dx *uniface, const char* attr, double* point_1, double* value) {
    uniface->push(std::string(attr), mui::point1dx(*point_1), *value);
}

void mui_push_1t_f(mui_uniface_1t *uniface, const char* attr, double* point_1, double* value) {
    mui::mui_f_wrapper_1D::point_type push_point(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    uniface->push(std::string(attr), push_point, static_cast<mui::mui_f_wrapper_1D::REAL>(*value));
}

// Single parameter push functions
void mui_push_1f_param_f(mui_uniface_1f *uniface, const char* attr, float* value) {
    uniface->push(std::string(attr), *value);
}

void mui_push_1fx_param_f(mui_uniface_1fx *uniface, const char* attr, float* value) {
    uniface->push(std::string(attr), *value);
}

void mui_push_1d_param_f(mui_uniface_1d *uniface, const char* attr, double* value) {
    uniface->push(std::string(attr), *value);
}

void mui_push_1dx_param_f(mui_uniface_1dx *uniface, const char* attr, double* value) {
    uniface->push(std::string(attr), *value);
}

void mui_push_1t_param_f(mui_uniface_1t *uniface, const char* attr, double* value) {
    uniface->push(std::string(attr), static_cast<mui::mui_f_wrapper_1D::REAL>(*value));
}

/******************************************
 * MUI functions for data commit           *
 ******************************************/

// Commit using one time value
void mui_commit_1f_f(mui_uniface_1f *uniface, float* t) {
    uniface->commit(*t);
}

void mui_commit_1fx_f(mui_uniface_1fx *uniface, float* t) {
    uniface->commit(*t);
}

void mui_commit_1d_f(mui_uniface_1d *uniface, double* t) {
    uniface->commit(*t);
}

void mui_commit_1dx_f(mui_uniface_1dx *uniface, double* t) {
    uniface->commit(*t);
}

void mui_commit_1t_f(mui_uniface_1t *uniface, double* t) {
    uniface->commit(static_cast<mui::mui_f_wrapper_1D::time_type>(*t));
}

// Commit using two time values
void mui_commit_1f_pair_f(mui_uniface_1f *uniface, float* t_1, float* t_2) {
    uniface->commit(*t_1, *t_2);
}

void mui_commit_1fx_pair_f(mui_uniface_1fx *uniface, float* t_1, float* t_2) {
    uniface->commit(*t_1, *t_2);
}

void mui_commit_1d_pair_f(mui_uniface_1d *uniface, double* t_1, double* t_2) {
    uniface->commit(*t_1, *t_2);
}

void mui_commit_1dx_pair_f(mui_uniface_1dx *uniface, double* t_1, double* t_2) {
    uniface->commit(*t_1, *t_2);
}

void mui_commit_1t_pair_f(mui_uniface_1t *uniface, double* t_1, double* t_2) {
    uniface->commit(static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1),
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2));
}

/******************************************
 * MUI functions for data forecast         *
 ******************************************/

// Forecast using one time value
void mui_forecast_1f_f(mui_uniface_1f *uniface, float* t) {
    uniface->forecast(*t);
}

void mui_forecast_1fx_f(mui_uniface_1fx *uniface, float* t) {
    uniface->forecast(*t);
}

void mui_forecast_1d_f(mui_uniface_1d *uniface, double* t) {
    uniface->forecast(*t);
}

void mui_forecast_1dx_f(mui_uniface_1dx *uniface, double* t) {
    uniface->forecast(*t);
}

void mui_forecast_1t_f(mui_uniface_1t *uniface, double* t) {
    uniface->forecast(static_cast<mui::mui_f_wrapper_1D::time_type>(*t));
}

// Forecast using two time values
void mui_forecast_1f_pair_f(mui_uniface_1f *uniface, float* t_1, float* t_2) {
    uniface->forecast(*t_1, *t_2);
}

void mui_forecast_1fx_pair_f(mui_uniface_1fx *uniface, float* t_1, float* t_2) {
    uniface->forecast(*t_1, *t_2);
}

void mui_forecast_1d_pair_f(mui_uniface_1d *uniface, double* t_1, double* t_2) {
    uniface->forecast(*t_1, *t_2);
}

void mui_forecast_1dx_pair_f(mui_uniface_1dx *uniface, double* t_1, double* t_2) {
    uniface->forecast(*t_1, *t_2);
}

void mui_forecast_1t_pair_f(mui_uniface_1t *uniface, double* t_1, double* t_2) {
    uniface->forecast(static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1),
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2));
}

/*******************************************************
 * MUI functions for 1D data fetch using one time value *
 ********************************************************/

// Spatial sampler: exact; temporal sampler: exact
void mui_fetch_exact_exact_1f_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t,
        mui_sampler_exact_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_exact_1fx_f(mui_uniface_1fx *uniface, const char *attr, float* point_1, float* t,
        mui_sampler_exact_1fx *spatial_sampler, mui_chrono_sampler_exact_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_exact_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_exact_1d *spatial_sampler, mui_chrono_sampler_exact_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_exact_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_exact_1dx *spatial_sampler, mui_chrono_sampler_exact_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_exact_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_exact_1t *spatial_sampler, mui_chrono_sampler_exact_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: exact; temporal sampler: gauss
void mui_fetch_exact_gauss_1f_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t,
        mui_sampler_exact_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_gauss_1fx_f(mui_uniface_1fx *uniface, const char *attr, float* point_1, float* t,
        mui_sampler_exact_1fx *spatial_sampler, mui_chrono_sampler_gauss_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_gauss_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_exact_1d *spatial_sampler, mui_chrono_sampler_gauss_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_gauss_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_exact_1dx *spatial_sampler, mui_chrono_sampler_gauss_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_gauss_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_exact_1t *spatial_sampler, mui_chrono_sampler_gauss_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: exact; temporal sampler: mean
void mui_fetch_exact_mean_1f_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t,
        mui_sampler_exact_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_mean_1fx_f(mui_uniface_1fx *uniface, const char *attr, float* point_1, float* t,
        mui_sampler_exact_1fx *spatial_sampler, mui_chrono_sampler_mean_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_mean_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_exact_1d *spatial_sampler, mui_chrono_sampler_mean_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_mean_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_exact_1dx *spatial_sampler, mui_chrono_sampler_mean_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_mean_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_exact_1t *spatial_sampler, mui_chrono_sampler_mean_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: exact; temporal sampler: sum
void mui_fetch_exact_sum_1f_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t,
        mui_sampler_exact_1f *spatial_sampler, mui_chrono_sampler_sum_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_sum_1fx_f(mui_uniface_1fx *uniface, const char *attr, float* point_1, float* t,
        mui_sampler_exact_1fx *spatial_sampler, mui_chrono_sampler_sum_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_sum_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_exact_1d *spatial_sampler, mui_chrono_sampler_sum_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_sum_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_exact_1dx *spatial_sampler, mui_chrono_sampler_sum_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_sum_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_exact_1t *spatial_sampler, mui_chrono_sampler_sum_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: gauss; temporal sampler: exact
void mui_fetch_gauss_exact_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_gauss_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_exact_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_gauss_1fx *spatial_sampler, mui_chrono_sampler_exact_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_exact_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_gauss_1d *spatial_sampler, mui_chrono_sampler_exact_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_exact_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_gauss_1dx *spatial_sampler, mui_chrono_sampler_exact_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_exact_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_gauss_1t *spatial_sampler, mui_chrono_sampler_exact_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: gauss; temporal sampler: gauss
void mui_fetch_gauss_gauss_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_gauss_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_gauss_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_gauss_1fx *spatial_sampler, mui_chrono_sampler_gauss_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_gauss_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_gauss_1d *spatial_sampler, mui_chrono_sampler_gauss_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_gauss_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_gauss_1dx *spatial_sampler, mui_chrono_sampler_gauss_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_gauss_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_gauss_1t *spatial_sampler, mui_chrono_sampler_gauss_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: gauss; temporal sampler: mean
void mui_fetch_gauss_mean_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_gauss_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_mean_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_gauss_1fx *spatial_sampler, mui_chrono_sampler_mean_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_mean_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_gauss_1d *spatial_sampler, mui_chrono_sampler_mean_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_mean_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_gauss_1dx *spatial_sampler, mui_chrono_sampler_mean_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_mean_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_gauss_1t *spatial_sampler, mui_chrono_sampler_mean_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: gauss; temporal sampler: sum
void mui_fetch_gauss_sum_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_gauss_1f *spatial_sampler, mui_chrono_sampler_sum_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_sum_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_gauss_1fx *spatial_sampler, mui_chrono_sampler_sum_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_sum_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_gauss_1d *spatial_sampler, mui_chrono_sampler_sum_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_sum_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_gauss_1dx *spatial_sampler, mui_chrono_sampler_sum_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_sum_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_gauss_1t *spatial_sampler, mui_chrono_sampler_sum_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: moving average; temporal sampler: exact
void mui_fetch_moving_average_exact_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_moving_average_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_exact_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_moving_average_1fx *spatial_sampler, mui_chrono_sampler_exact_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_exact_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_moving_average_1d *spatial_sampler, mui_chrono_sampler_exact_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_exact_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_moving_average_1dx *spatial_sampler, mui_chrono_sampler_exact_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_exact_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_moving_average_1t *spatial_sampler, mui_chrono_sampler_exact_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: moving average; temporal sampler: gauss
void mui_fetch_moving_average_gauss_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_moving_average_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_gauss_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_moving_average_1fx *spatial_sampler, mui_chrono_sampler_gauss_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_gauss_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_moving_average_1d *spatial_sampler, mui_chrono_sampler_gauss_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_gauss_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_moving_average_1dx *spatial_sampler, mui_chrono_sampler_gauss_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_gauss_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_moving_average_1t *spatial_sampler, mui_chrono_sampler_gauss_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: moving average; temporal sampler: mean
void mui_fetch_moving_average_mean_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_moving_average_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_mean_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_moving_average_1fx *spatial_sampler, mui_chrono_sampler_mean_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_mean_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_moving_average_1d *spatial_sampler, mui_chrono_sampler_mean_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_mean_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_moving_average_1dx *spatial_sampler, mui_chrono_sampler_mean_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_mean_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_moving_average_1t *spatial_sampler, mui_chrono_sampler_mean_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: moving average; temporal sampler: sum
void mui_fetch_moving_average_sum_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_moving_average_1f *spatial_sampler, mui_chrono_sampler_sum_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_sum_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_moving_average_1fx *spatial_sampler, mui_chrono_sampler_sum_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_sum_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_moving_average_1d *spatial_sampler, mui_chrono_sampler_sum_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_sum_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_moving_average_1dx *spatial_sampler, mui_chrono_sampler_sum_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_sum_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_moving_average_1t *spatial_sampler, mui_chrono_sampler_sum_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: nearest neighbor; temporal sampler: exact
void mui_fetch_nearest_neighbor_exact_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_nearest_neighbor_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_exact_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_nearest_neighbor_1fx *spatial_sampler, mui_chrono_sampler_exact_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_exact_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_nearest_neighbor_1d *spatial_sampler, mui_chrono_sampler_exact_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_exact_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_nearest_neighbor_1dx *spatial_sampler, mui_chrono_sampler_exact_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_exact_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_nearest_neighbor_1t *spatial_sampler, mui_chrono_sampler_exact_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: nearest neighbor; temporal sampler: gauss
void mui_fetch_nearest_neighbor_gauss_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_nearest_neighbor_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_gauss_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_nearest_neighbor_1fx *spatial_sampler, mui_chrono_sampler_gauss_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_gauss_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_nearest_neighbor_1d *spatial_sampler, mui_chrono_sampler_gauss_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_gauss_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_nearest_neighbor_1dx *spatial_sampler, mui_chrono_sampler_gauss_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_gauss_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_nearest_neighbor_1t *spatial_sampler, mui_chrono_sampler_gauss_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: nearest neighbor; temporal sampler: mean
void mui_fetch_nearest_neighbor_mean_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_nearest_neighbor_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_mean_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_nearest_neighbor_1fx *spatial_sampler, mui_chrono_sampler_mean_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_mean_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_nearest_neighbor_1d *spatial_sampler, mui_chrono_sampler_mean_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_mean_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_nearest_neighbor_1dx *spatial_sampler, mui_chrono_sampler_mean_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_mean_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_nearest_neighbor_1t *spatial_sampler, mui_chrono_sampler_mean_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: nearest neighbor; temporal sampler: sum
void mui_fetch_nearest_neighbor_sum_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_nearest_neighbor_1f *spatial_sampler, mui_chrono_sampler_sum_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_sum_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_nearest_neighbor_1fx *spatial_sampler, mui_chrono_sampler_sum_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_sum_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_nearest_neighbor_1d *spatial_sampler, mui_chrono_sampler_sum_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_sum_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_nearest_neighbor_1dx *spatial_sampler, mui_chrono_sampler_sum_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_sum_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_nearest_neighbor_1t *spatial_sampler, mui_chrono_sampler_sum_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: pseudo nearest neighbor; temporal sampler: exact
void mui_fetch_pseudo_nearest_neighbor_exact_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_pseudo_nearest_neighbor_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_exact_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *point_1,
        float* t, mui_sampler_pseudo_nearest_neighbor_1fx *spatial_sampler,
        mui_chrono_sampler_exact_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_exact_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1,
        double* t, mui_sampler_pseudo_nearest_neighbor_1d *spatial_sampler,
        mui_chrono_sampler_exact_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_exact_1dx(mui_uniface_1dx *uniface, const char *attr, double* point_1,
        double* t, mui_sampler_pseudo_nearest_neighbor_1dx *spatial_sampler,
        mui_chrono_sampler_exact_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_exact_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1,
        double* t, mui_sampler_pseudo_nearest_neighbor_1t *spatial_sampler,
        mui_chrono_sampler_exact_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: nearest neighbor; temporal sampler: gauss
void mui_fetch_pseudo_nearest_neighbor_gauss_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_pseudo_nearest_neighbor_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_gauss_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *point_1,
        float* t, mui_sampler_pseudo_nearest_neighbor_1fx *spatial_sampler,
        mui_chrono_sampler_gauss_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_gauss_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1,
        double* t, mui_sampler_pseudo_nearest_neighbor_1d *spatial_sampler,
        mui_chrono_sampler_gauss_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_gauss_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1,
        double* t, mui_sampler_pseudo_nearest_neighbor_1dx *spatial_sampler,
        mui_chrono_sampler_gauss_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_gauss_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1,
        double* t, mui_sampler_pseudo_nearest_neighbor_1t *spatial_sampler,
        mui_chrono_sampler_gauss_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: nearest neighbor; temporal sampler: mean
void mui_fetch_pseudo_nearest_neighbor_mean_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_pseudo_nearest_neighbor_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_mean_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *point_1,
        float* t, mui_sampler_pseudo_nearest_neighbor_1fx *spatial_sampler,
        mui_chrono_sampler_mean_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_mean_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1,
        double* t, mui_sampler_pseudo_nearest_neighbor_1d *spatial_sampler,
        mui_chrono_sampler_mean_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_mean_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1,
        double* t, mui_sampler_pseudo_nearest_neighbor_1dx *spatial_sampler,
        mui_chrono_sampler_mean_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_mean_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1,
        double* t, mui_sampler_pseudo_nearest_neighbor_1t *spatial_sampler,
        mui_chrono_sampler_mean_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: nearest neighbor; temporal sampler: sum
void mui_fetch_pseudo_nearest_neighbor_sum_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_pseudo_nearest_neighbor_1f *spatial_sampler, mui_chrono_sampler_sum_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_sum_1fx(mui_uniface_1fx *uniface, const char *attr, float *point_1,
        float* t, mui_sampler_pseudo_nearest_neighbor_1fx *spatial_sampler,
        mui_chrono_sampler_sum_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_sum_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_pseudo_nearest_neighbor_1d *spatial_sampler, mui_chrono_sampler_sum_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_sum_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1,
        double* t, mui_sampler_pseudo_nearest_neighbor_1dx *spatial_sampler,
        mui_chrono_sampler_sum_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_sum_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_pseudo_nearest_neighbor_1t *spatial_sampler, mui_chrono_sampler_sum_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: shepard quintic; temporal sampler: exact
void mui_fetch_shepard_quintic_exact_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_shepard_quintic_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_exact_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_shepard_quintic_1fx *spatial_sampler, mui_chrono_sampler_exact_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_exact_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_shepard_quintic_1d *spatial_sampler, mui_chrono_sampler_exact_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_exact_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_shepard_quintic_1dx *spatial_sampler, mui_chrono_sampler_exact_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_exact_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_shepard_quintic_1t *spatial_sampler, mui_chrono_sampler_exact_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: shepard quintic; temporal sampler: gauss
void mui_fetch_shepard_quintic_gauss_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_shepard_quintic_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_gauss_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_shepard_quintic_1fx *spatial_sampler, mui_chrono_sampler_gauss_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_gauss_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_shepard_quintic_1d *spatial_sampler, mui_chrono_sampler_gauss_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_gauss_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_shepard_quintic_1dx *spatial_sampler, mui_chrono_sampler_gauss_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_gauss_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_shepard_quintic_1t *spatial_sampler, mui_chrono_sampler_gauss_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: shepard quintic; temporal sampler: mean
void mui_fetch_shepard_quintic_mean_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_shepard_quintic_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_mean_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_shepard_quintic_1fx *spatial_sampler, mui_chrono_sampler_mean_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_mean_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_shepard_quintic_1d *spatial_sampler, mui_chrono_sampler_mean_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_mean_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_shepard_quintic_1dx *spatial_sampler, mui_chrono_sampler_mean_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_mean_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_shepard_quintic_1t *spatial_sampler, mui_chrono_sampler_mean_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: shepard quintic; temporal sampler: sum
void mui_fetch_shepard_quintic_sum_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_shepard_quintic_1f *spatial_sampler, mui_chrono_sampler_sum_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_sum_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_shepard_quintic_1fx *spatial_sampler, mui_chrono_sampler_sum_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_sum_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_shepard_quintic_1d *spatial_sampler, mui_chrono_sampler_sum_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_sum_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_shepard_quintic_1dx *spatial_sampler, mui_chrono_sampler_sum_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_sum_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_shepard_quintic_1t *spatial_sampler, mui_chrono_sampler_sum_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: sph-derived quintic; temporal sampler: exact
void mui_fetch_sph_quintic_exact_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_sph_quintic_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_exact_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_sph_quintic_1fx *spatial_sampler, mui_chrono_sampler_exact_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_exact_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_sph_quintic_1d *spatial_sampler, mui_chrono_sampler_exact_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_exact_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_sph_quintic_1dx *spatial_sampler, mui_chrono_sampler_exact_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_exact_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_sph_quintic_1t *spatial_sampler, mui_chrono_sampler_exact_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: sph-derived quintic; temporal sampler: gauss
void mui_fetch_sph_quintic_gauss_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_sph_quintic_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_gauss_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_sph_quintic_1fx *spatial_sampler, mui_chrono_sampler_gauss_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_gauss_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_sph_quintic_1d *spatial_sampler, mui_chrono_sampler_gauss_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_gauss_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_sph_quintic_1dx *spatial_sampler, mui_chrono_sampler_gauss_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_gauss_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_sph_quintic_1t *spatial_sampler, mui_chrono_sampler_gauss_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: sph-derived quintic; temporal sampler: mean
void mui_fetch_sph_quintic_mean_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_sph_quintic_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_mean_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_sph_quintic_1fx *spatial_sampler, mui_chrono_sampler_mean_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_mean_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_sph_quintic_1d *spatial_sampler, mui_chrono_sampler_mean_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_mean_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_sph_quintic_1dx *spatial_sampler, mui_chrono_sampler_mean_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_mean_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_sph_quintic_1t *spatial_sampler, mui_chrono_sampler_mean_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: sph-derived quintic; temporal sampler: sum
void mui_fetch_sph_quintic_sum_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_sph_quintic_1f *spatial_sampler, mui_chrono_sampler_sum_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_sum_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_sph_quintic_1fx *spatial_sampler, mui_chrono_sampler_sum_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_sum_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_sph_quintic_1d *spatial_sampler, mui_chrono_sampler_sum_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_sum_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_sph_quintic_1dx *spatial_sampler, mui_chrono_sampler_sum_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_sum_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_sph_quintic_1t *spatial_sampler, mui_chrono_sampler_sum_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: summation quintic; temporal sampler: exact
void mui_fetch_sum_quintic_exact_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_sum_quintic_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_exact_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_sum_quintic_1fx *spatial_sampler, mui_chrono_sampler_exact_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_exact_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_sum_quintic_1d *spatial_sampler, mui_chrono_sampler_exact_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_exact_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_sum_quintic_1dx *spatial_sampler, mui_chrono_sampler_exact_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_exact_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_sum_quintic_1t *spatial_sampler, mui_chrono_sampler_exact_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: summation quintic; temporal sampler: gauss
void mui_fetch_sum_quintic_gauss_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_sum_quintic_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_gauss_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_sum_quintic_1fx *spatial_sampler, mui_chrono_sampler_gauss_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_gauss_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_sum_quintic_1d *spatial_sampler, mui_chrono_sampler_gauss_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_gauss_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_sum_quintic_1dx *spatial_sampler, mui_chrono_sampler_gauss_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_gauss_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_sum_quintic_1t *spatial_sampler, mui_chrono_sampler_gauss_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: summation quintic; temporal sampler: mean
void mui_fetch_sum_quintic_mean_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_sum_quintic_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_mean_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_sum_quintic_1fx *spatial_sampler, mui_chrono_sampler_mean_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_mean_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_sum_quintic_1d *spatial_sampler, mui_chrono_sampler_mean_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_mean_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_sum_quintic_1dx *spatial_sampler, mui_chrono_sampler_mean_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_mean_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_sum_quintic_1t *spatial_sampler, mui_chrono_sampler_mean_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: summation quintic; temporal sampler: sum
void mui_fetch_sum_quintic_sum_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_sum_quintic_1f *spatial_sampler, mui_chrono_sampler_sum_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_sum_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_sum_quintic_1fx *spatial_sampler, mui_chrono_sampler_sum_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_sum_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_sum_quintic_1d *spatial_sampler, mui_chrono_sampler_sum_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_sum_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_sum_quintic_1dx *spatial_sampler, mui_chrono_sampler_sum_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_sum_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_sum_quintic_1t *spatial_sampler, mui_chrono_sampler_sum_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

#ifdef USE_RBF
// Spatial sampler: radial basis function; temporal sampler: exact
void mui_fetch_rbf_exact_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_rbf_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_exact_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_rbf_1fx *spatial_sampler, mui_chrono_sampler_exact_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_exact_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_rbf_1d *spatial_sampler, mui_chrono_sampler_exact_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_exact_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_rbf_1dx *spatial_sampler, mui_chrono_sampler_exact_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_exact_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_rbf_1t *spatial_sampler, mui_chrono_sampler_exact_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: radial basis function; temporal sampler: gauss
void mui_fetch_rbf_gauss_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_rbf_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_gauss_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_rbf_1fx *spatial_sampler, mui_chrono_sampler_gauss_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_gauss_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_rbf_1d *spatial_sampler, mui_chrono_sampler_gauss_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_gauss_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_rbf_1dx *spatial_sampler, mui_chrono_sampler_gauss_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_gauss_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_rbf_1t *spatial_sampler, mui_chrono_sampler_gauss_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: radial basis function; temporal sampler: mean
void mui_fetch_rbf_mean_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_rbf_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_mean_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_rbf_1fx *spatial_sampler, mui_chrono_sampler_mean_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_mean_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_rbf_1d *spatial_sampler, mui_chrono_sampler_mean_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_mean_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_rbf_1dx *spatial_sampler, mui_chrono_sampler_mean_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_mean_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_rbf_1t *spatial_sampler, mui_chrono_sampler_mean_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: radial basis function; temporal sampler: sum
void mui_fetch_rbf_sum_1f_f(mui_uniface_1f *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_rbf_1f *spatial_sampler, mui_chrono_sampler_sum_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_sum_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *point_1, float* t,
        mui_sampler_rbf_1fx *spatial_sampler, mui_chrono_sampler_sum_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_sum_1d_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_rbf_1d *spatial_sampler, mui_chrono_sampler_sum_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_sum_1dx_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_rbf_1dx *spatial_sampler, mui_chrono_sampler_sum_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_sum_1t_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t,
        mui_sampler_rbf_1t *spatial_sampler, mui_chrono_sampler_sum_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}
#endif

/********************************************************
 * MUI functions for 1D data fetch using two time values *
 *********************************************************/

// Spatial sampler: exact; temporal sampler: exact
void mui_fetch_exact_exact_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t_1, float* t_2,
        mui_sampler_exact_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_exact_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_exact_1fx *spatial_sampler, mui_chrono_sampler_exact_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_exact_exact_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_exact_1d *spatial_sampler, mui_chrono_sampler_exact_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_exact_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_exact_1dx *spatial_sampler, mui_chrono_sampler_exact_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_exact_exact_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_exact_1t *spatial_sampler, mui_chrono_sampler_exact_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: exact; temporal sampler: gauss
void mui_fetch_exact_gauss_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t_1, float* t_2,
        mui_sampler_exact_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_gauss_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_exact_1fx *spatial_sampler, mui_chrono_sampler_gauss_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_exact_gauss_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_exact_1d *spatial_sampler, mui_chrono_sampler_gauss_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_gauss_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_exact_1dx *spatial_sampler, mui_chrono_sampler_gauss_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_exact_gauss_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_exact_1t *spatial_sampler, mui_chrono_sampler_gauss_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: exact; temporal sampler: mean
void mui_fetch_exact_mean_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t_1, float* t_2,
        mui_sampler_exact_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_mean_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_exact_1fx *spatial_sampler, mui_chrono_sampler_mean_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_exact_mean_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_exact_1d *spatial_sampler, mui_chrono_sampler_mean_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_mean_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_exact_1dx *spatial_sampler, mui_chrono_sampler_mean_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_exact_mean_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_exact_1t *spatial_sampler, mui_chrono_sampler_mean_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: exact; temporal sampler: sum
void mui_fetch_exact_sum_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t_1, float* t_2,
        mui_sampler_exact_1f *spatial_sampler, mui_chrono_sampler_sum_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_sum_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_exact_1fx *spatial_sampler, mui_chrono_sampler_sum_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_exact_sum_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_exact_1d *spatial_sampler, mui_chrono_sampler_sum_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_sum_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_exact_1dx *spatial_sampler, mui_chrono_sampler_sum_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_exact_sum_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_exact_1t *spatial_sampler, mui_chrono_sampler_sum_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: gauss; temporal sampler: exact
void mui_fetch_gauss_exact_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t_1, float* t_2,
        mui_sampler_gauss_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_exact_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_gauss_1fx *spatial_sampler, mui_chrono_sampler_exact_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_gauss_exact_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_gauss_1d *spatial_sampler, mui_chrono_sampler_exact_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_exact_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_gauss_1dx *spatial_sampler, mui_chrono_sampler_exact_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_gauss_exact_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_gauss_1t *spatial_sampler, mui_chrono_sampler_exact_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: gauss; temporal sampler: gauss
void mui_fetch_gauss_gauss_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t_1, float* t_2,
        mui_sampler_gauss_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_gauss_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_gauss_1fx *spatial_sampler, mui_chrono_sampler_gauss_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_gauss_gauss_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_gauss_1d *spatial_sampler, mui_chrono_sampler_gauss_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_gauss_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_gauss_1dx *spatial_sampler, mui_chrono_sampler_gauss_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_gauss_gauss_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_gauss_1t *spatial_sampler, mui_chrono_sampler_gauss_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: gauss; temporal sampler: mean
void mui_fetch_gauss_mean_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t_1, float* t_2,
        mui_sampler_gauss_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_mean_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_gauss_1fx *spatial_sampler, mui_chrono_sampler_mean_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_gauss_mean_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_gauss_1d *spatial_sampler, mui_chrono_sampler_mean_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_mean_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_gauss_1dx *spatial_sampler, mui_chrono_sampler_mean_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_gauss_mean_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_gauss_1t *spatial_sampler, mui_chrono_sampler_mean_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: moving average; temporal sampler: exact
void mui_fetch_moving_average_exact_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_moving_average_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_exact_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1,
        float* t_1, float* t_2, mui_sampler_moving_average_1fx *spatial_sampler,
        mui_chrono_sampler_exact_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_moving_average_exact_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_moving_average_1d *spatial_sampler, mui_chrono_sampler_exact_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_exact_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_moving_average_1dx *spatial_sampler,
        mui_chrono_sampler_exact_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_moving_average_exact_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_moving_average_1t *spatial_sampler, mui_chrono_sampler_exact_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: moving average; temporal sampler: gauss
void mui_fetch_moving_average_gauss_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_moving_average_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_gauss_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1,
        float* t_1, float* t_2, mui_sampler_moving_average_1fx *spatial_sampler,
        mui_chrono_sampler_gauss_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_moving_average_gauss_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_moving_average_1d *spatial_sampler, mui_chrono_sampler_gauss_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_gauss_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_moving_average_1dx *spatial_sampler,
        mui_chrono_sampler_gauss_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_moving_average_gauss_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_moving_average_1t *spatial_sampler, mui_chrono_sampler_gauss_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: moving average; temporal sampler: mean
void mui_fetch_moving_average_mean_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_moving_average_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_mean_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_moving_average_1fx *spatial_sampler, mui_chrono_sampler_mean_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_moving_average_mean_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_moving_average_1d *spatial_sampler, mui_chrono_sampler_mean_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_mean_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_moving_average_1dx *spatial_sampler,
        mui_chrono_sampler_mean_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_moving_average_mean_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_moving_average_1t *spatial_sampler, mui_chrono_sampler_mean_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: moving average; temporal sampler: sum
void mui_fetch_moving_average_sum_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_moving_average_1f *spatial_sampler, mui_chrono_sampler_sum_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_sum_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_moving_average_1fx *spatial_sampler, mui_chrono_sampler_sum_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_moving_average_sum_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_moving_average_1d *spatial_sampler, mui_chrono_sampler_sum_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_sum_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_moving_average_1dx *spatial_sampler,
        mui_chrono_sampler_sum_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_moving_average_sum_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_moving_average_1t *spatial_sampler, mui_chrono_sampler_sum_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: nearest neighbor; temporal sampler: exact
void mui_fetch_nearest_neighbor_exact_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_nearest_neighbor_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_exact_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1,
        float* t_1, float* t_2, mui_sampler_nearest_neighbor_1fx *spatial_sampler,
        mui_chrono_sampler_exact_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_nearest_neighbor_exact_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_nearest_neighbor_1d *spatial_sampler,
        mui_chrono_sampler_exact_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_exact_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_nearest_neighbor_1dx *spatial_sampler,
        mui_chrono_sampler_exact_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_nearest_neighbor_exact_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_nearest_neighbor_1t *spatial_sampler,
        mui_chrono_sampler_exact_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: nearest neighbor; temporal sampler: gauss
void mui_fetch_nearest_neighbor_gauss_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_nearest_neighbor_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_gauss_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1,
        float* t_1, float* t_2, mui_sampler_nearest_neighbor_1fx *spatial_sampler,
        mui_chrono_sampler_gauss_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_nearest_neighbor_gauss_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_nearest_neighbor_1d *spatial_sampler,
        mui_chrono_sampler_gauss_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_gauss_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_nearest_neighbor_1dx *spatial_sampler,
        mui_chrono_sampler_gauss_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_nearest_neighbor_gauss_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_nearest_neighbor_1t *spatial_sampler,
        mui_chrono_sampler_gauss_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: nearest neighbor; temporal sampler: mean
void mui_fetch_nearest_neighbor_mean_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_nearest_neighbor_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_mean_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1,
        float* t_1, float* t_2, mui_sampler_nearest_neighbor_1fx *spatial_sampler,
        mui_chrono_sampler_mean_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_nearest_neighbor_mean_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_nearest_neighbor_1d *spatial_sampler,
        mui_chrono_sampler_mean_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_mean_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_nearest_neighbor_1dx *spatial_sampler,
        mui_chrono_sampler_mean_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_nearest_neighbor_mean_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_nearest_neighbor_1t *spatial_sampler,
        mui_chrono_sampler_mean_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: nearest neighbor; temporal sampler: sum
void mui_fetch_nearest_neighbor_sum_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_nearest_neighbor_1f *spatial_sampler, mui_chrono_sampler_sum_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_sum_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1,
        float* t_1, float* t_2, mui_sampler_nearest_neighbor_1fx *spatial_sampler,
        mui_chrono_sampler_sum_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_nearest_neighbor_sum_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_nearest_neighbor_1d *spatial_sampler, mui_chrono_sampler_sum_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_sum_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_nearest_neighbor_1dx *spatial_sampler,
        mui_chrono_sampler_sum_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_nearest_neighbor_sum_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_nearest_neighbor_1t *spatial_sampler, mui_chrono_sampler_sum_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: nearest neighbor; temporal sampler: exact
void mui_fetch_pseudo_nearest_neighbor_exact_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1,
        float* t_1, float* t_2, mui_sampler_pseudo_nearest_neighbor_1f *spatial_sampler,
        mui_chrono_sampler_exact_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_exact_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1,
        float* t_1, float* t_2, mui_sampler_pseudo_nearest_neighbor_1fx *spatial_sampler,
        mui_chrono_sampler_exact_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_exact_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_1d *spatial_sampler,
        mui_chrono_sampler_exact_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_exact_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_1dx *spatial_sampler,
        mui_chrono_sampler_exact_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_exact_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_1t *spatial_sampler,
        mui_chrono_sampler_exact_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: nearest neighbor; temporal sampler: gauss
void mui_fetch_pseudo_nearest_neighbor_gauss_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1,
        float* t_1, float* t_2, mui_sampler_pseudo_nearest_neighbor_1f *spatial_sampler,
        mui_chrono_sampler_gauss_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_gauss_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1,
        float* t_1, float* t_2, mui_sampler_pseudo_nearest_neighbor_1fx *spatial_sampler,
        mui_chrono_sampler_gauss_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_gauss_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_1d *spatial_sampler,
        mui_chrono_sampler_gauss_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_gauss_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_1dx *spatial_sampler,
        mui_chrono_sampler_gauss_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_gauss_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_1t *spatial_sampler,
        mui_chrono_sampler_gauss_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: nearest neighbor; temporal sampler: mean
void mui_fetch_pseudo_nearest_neighbor_mean_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1,
        float* t_1, float* t_2, mui_sampler_pseudo_nearest_neighbor_1f *spatial_sampler,
        mui_chrono_sampler_mean_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_mean_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1,
        float* t_1, float* t_2, mui_sampler_pseudo_nearest_neighbor_1fx *spatial_sampler,
        mui_chrono_sampler_mean_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_mean_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_1d *spatial_sampler,
        mui_chrono_sampler_mean_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_mean_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_1dx *spatial_sampler,
        mui_chrono_sampler_mean_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_mean_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_1t *spatial_sampler,
        mui_chrono_sampler_mean_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: nearest neighbor; temporal sampler: sum
void mui_fetch_pseudo_nearest_neighbor_sum_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1,
        float* t_1, float* t_2, mui_sampler_pseudo_nearest_neighbor_1f *spatial_sampler,
        mui_chrono_sampler_sum_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_sum_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1,
        float* t_1, float* t_2, mui_sampler_pseudo_nearest_neighbor_1fx *spatial_sampler,
        mui_chrono_sampler_sum_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_sum_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_1d *spatial_sampler,
        mui_chrono_sampler_sum_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_sum_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_1dx *spatial_sampler,
        mui_chrono_sampler_sum_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_sum_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_1t *spatial_sampler,
        mui_chrono_sampler_sum_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: shepard quintic; temporal sampler: exact
void mui_fetch_shepard_quintic_exact_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_shepard_quintic_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_exact_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1,
        float* t_1, float* t_2, mui_sampler_shepard_quintic_1fx *spatial_sampler,
        mui_chrono_sampler_exact_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_shepard_quintic_exact_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_shepard_quintic_1d *spatial_sampler,
        mui_chrono_sampler_exact_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_exact_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_shepard_quintic_1dx *spatial_sampler,
        mui_chrono_sampler_exact_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_shepard_quintic_exact_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_shepard_quintic_1t *spatial_sampler,
        mui_chrono_sampler_exact_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: shepard quintic; temporal sampler: gauss
void mui_fetch_shepard_quintic_gauss_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_shepard_quintic_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_gauss_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1,
        float* t_1, float* t_2, mui_sampler_shepard_quintic_1fx *spatial_sampler,
        mui_chrono_sampler_gauss_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_shepard_quintic_gauss_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_shepard_quintic_1d *spatial_sampler,
        mui_chrono_sampler_gauss_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_gauss_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_shepard_quintic_1dx *spatial_sampler,
        mui_chrono_sampler_gauss_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_shepard_quintic_gauss_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_shepard_quintic_1t *spatial_sampler,
        mui_chrono_sampler_gauss_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: shepard quintic; temporal sampler: mean
void mui_fetch_shepard_quintic_mean_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_shepard_quintic_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_mean_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1,
        float* t_1, float* t_2, mui_sampler_shepard_quintic_1fx *spatial_sampler,
        mui_chrono_sampler_mean_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_shepard_quintic_mean_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_shepard_quintic_1d *spatial_sampler, mui_chrono_sampler_mean_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_mean_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_shepard_quintic_1dx *spatial_sampler,
        mui_chrono_sampler_mean_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_shepard_quintic_mean_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_shepard_quintic_1t *spatial_sampler, mui_chrono_sampler_mean_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: shepard quintic; temporal sampler: sum
void mui_fetch_shepard_quintic_sum_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_shepard_quintic_1f *spatial_sampler, mui_chrono_sampler_sum_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_sum_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_shepard_quintic_1fx *spatial_sampler, mui_chrono_sampler_sum_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_shepard_quintic_sum_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_shepard_quintic_1d *spatial_sampler, mui_chrono_sampler_sum_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_sum_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1,
        double* t_1, double* t_2, mui_sampler_shepard_quintic_1dx *spatial_sampler,
        mui_chrono_sampler_sum_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_shepard_quintic_sum_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_shepard_quintic_1t *spatial_sampler, mui_chrono_sampler_sum_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: sph-derived quintic; temporal sampler: exact
void mui_fetch_sph_quintic_exact_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_sph_quintic_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_exact_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_sph_quintic_1fx *spatial_sampler, mui_chrono_sampler_exact_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sph_quintic_exact_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_sph_quintic_1d *spatial_sampler, mui_chrono_sampler_exact_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_exact_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_sph_quintic_1dx *spatial_sampler, mui_chrono_sampler_exact_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sph_quintic_exact_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_sph_quintic_1t *spatial_sampler, mui_chrono_sampler_exact_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: sph-derived quintic; temporal sampler: gauss
void mui_fetch_sph_quintic_gauss_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_sph_quintic_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_gauss_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_sph_quintic_1fx *spatial_sampler, mui_chrono_sampler_gauss_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sph_quintic_gauss_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_sph_quintic_1d *spatial_sampler, mui_chrono_sampler_gauss_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_gauss_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_sph_quintic_1dx *spatial_sampler, mui_chrono_sampler_gauss_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sph_quintic_gauss_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_sph_quintic_1t *spatial_sampler, mui_chrono_sampler_gauss_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: sph-derived quintic; temporal sampler: mean
void mui_fetch_sph_quintic_mean_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_sph_quintic_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_mean_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_sph_quintic_1fx *spatial_sampler, mui_chrono_sampler_mean_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sph_quintic_mean_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_sph_quintic_1d *spatial_sampler, mui_chrono_sampler_mean_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_mean_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_sph_quintic_1dx *spatial_sampler, mui_chrono_sampler_mean_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sph_quintic_mean_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_sph_quintic_1t *spatial_sampler, mui_chrono_sampler_mean_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: sph-derived quintic; temporal sampler: sum
void mui_fetch_sph_quintic_sum_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_sph_quintic_1f *spatial_sampler, mui_chrono_sampler_sum_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_sum_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_sph_quintic_1fx *spatial_sampler, mui_chrono_sampler_sum_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sph_quintic_sum_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_sph_quintic_1d *spatial_sampler, mui_chrono_sampler_sum_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_sum_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_sph_quintic_1dx *spatial_sampler, mui_chrono_sampler_sum_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sph_quintic_sum_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_sph_quintic_1t *spatial_sampler, mui_chrono_sampler_sum_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: summation quintic; temporal sampler: exact
void mui_fetch_sum_quintic_exact_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_sum_quintic_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_exact_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_sum_quintic_1fx *spatial_sampler, mui_chrono_sampler_exact_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sum_quintic_exact_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_sum_quintic_1d *spatial_sampler, mui_chrono_sampler_exact_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_exact_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_sum_quintic_1dx *spatial_sampler, mui_chrono_sampler_exact_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sum_quintic_exact_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_sum_quintic_1t *spatial_sampler, mui_chrono_sampler_exact_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: summation quintic; temporal sampler: gauss
void mui_fetch_sum_quintic_gauss_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_sum_quintic_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_gauss_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_sum_quintic_1fx *spatial_sampler, mui_chrono_sampler_gauss_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sum_quintic_gauss_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_sum_quintic_1d *spatial_sampler, mui_chrono_sampler_gauss_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_gauss_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_sum_quintic_1dx *spatial_sampler, mui_chrono_sampler_gauss_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sum_quintic_gauss_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_sum_quintic_1t *spatial_sampler, mui_chrono_sampler_gauss_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: summation quintic; temporal sampler: mean
void mui_fetch_sum_quintic_mean_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_sum_quintic_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_mean_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_sum_quintic_1fx *spatial_sampler, mui_chrono_sampler_mean_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sum_quintic_mean_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_sum_quintic_1d *spatial_sampler, mui_chrono_sampler_mean_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_mean_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_sum_quintic_1dx *spatial_sampler, mui_chrono_sampler_mean_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sum_quintic_mean_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_sum_quintic_1t *spatial_sampler, mui_chrono_sampler_mean_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: summation quintic; temporal sampler: sum
void mui_fetch_sum_quintic_sum_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_sum_quintic_1f *spatial_sampler, mui_chrono_sampler_sum_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_sum_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_sum_quintic_1fx *spatial_sampler, mui_chrono_sampler_sum_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sum_quintic_sum_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_sum_quintic_1d *spatial_sampler, mui_chrono_sampler_sum_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_sum_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_sum_quintic_1dx *spatial_sampler, mui_chrono_sampler_sum_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sum_quintic_sum_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_sum_quintic_1t *spatial_sampler, mui_chrono_sampler_sum_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

#ifdef USE_RBF
// Spatial sampler: radial basis function; temporal sampler: exact
void mui_fetch_rbf_exact_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t_1, float* t_2,
        mui_sampler_rbf_1f *spatial_sampler, mui_chrono_sampler_exact_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_exact_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_rbf_1fx *spatial_sampler, mui_chrono_sampler_exact_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_rbf_exact_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_rbf_1d *spatial_sampler, mui_chrono_sampler_exact_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_exact_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_rbf_1dx *spatial_sampler, mui_chrono_sampler_exact_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_rbf_exact_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_rbf_1t *spatial_sampler, mui_chrono_sampler_exact_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: radial basis function; temporal sampler: gauss
void mui_fetch_rbf_gauss_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t_1, float* t_2,
        mui_sampler_rbf_1f *spatial_sampler, mui_chrono_sampler_gauss_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_gauss_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1, float* t_1,
        float* t_2, mui_sampler_rbf_1fx *spatial_sampler, mui_chrono_sampler_gauss_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_rbf_gauss_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_rbf_1d *spatial_sampler, mui_chrono_sampler_gauss_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_gauss_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_rbf_1dx *spatial_sampler, mui_chrono_sampler_gauss_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_rbf_gauss_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_rbf_1t *spatial_sampler, mui_chrono_sampler_gauss_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: radial basis function; temporal sampler: mean
void mui_fetch_rbf_mean_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t_1, float* t_2,
        mui_sampler_rbf_1f *spatial_sampler, mui_chrono_sampler_mean_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_mean_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1, float* t_1, float* t_2,
        mui_sampler_rbf_1fx *spatial_sampler, mui_chrono_sampler_mean_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_rbf_mean_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t_1, double* t_2,
        mui_sampler_rbf_1d *spatial_sampler, mui_chrono_sampler_mean_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_mean_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_rbf_1dx *spatial_sampler, mui_chrono_sampler_mean_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_rbf_mean_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t_1, double* t_2,
        mui_sampler_rbf_1t *spatial_sampler, mui_chrono_sampler_mean_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: radial basis function; temporal sampler: sum
void mui_fetch_rbf_sum_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float* point_1, float* t_1, float* t_2,
        mui_sampler_rbf_1f *spatial_sampler, mui_chrono_sampler_sum_1f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1f(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_sum_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float* point_1, float* t_1, float* t_2,
        mui_sampler_rbf_1fx *spatial_sampler, mui_chrono_sampler_sum_1fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1fx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_rbf_sum_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double* point_1, double* t_1, double* t_2,
        mui_sampler_rbf_1d *spatial_sampler, mui_chrono_sampler_sum_1d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1d(*point_1), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_sum_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double* point_1, double* t_1,
        double* t_2, mui_sampler_rbf_1dx *spatial_sampler, mui_chrono_sampler_sum_1dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point1dx(*point_1), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_rbf_sum_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double* point_1, double* t_1, double* t_2,
        mui_sampler_rbf_1t *spatial_sampler, mui_chrono_sampler_sum_1t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_1D::point_type point_fetch(static_cast<mui::mui_f_wrapper_1D::REAL>(*point_1));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}
#endif

/*******************************************************************
 * MUI functions for 1D data point only fetch using one time value  *
 ********************************************************************/

// Temporal sampler: exact
void mui_fetch_points_exact_1f_f(mui_uniface_1f *uniface, const char *attr, float *t,
        mui_chrono_sampler_exact_1f *temporal_sampler, float** ret_points_1, int *num_points) {
    std::vector<mui::point1f> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_exact_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *t,
        mui_chrono_sampler_exact_1fx *temporal_sampler, float** ret_points_1, int *num_points) {
    std::vector<mui::point1fx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_exact_1d_f(mui_uniface_1d *uniface, const char *attr, double *t,
        mui_chrono_sampler_exact_1d *temporal_sampler, double** ret_points_1, int *num_points) {
    std::vector<mui::point1d> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_exact_1dx_f(mui_uniface_1dx *uniface, const char *attr, double *t,
        mui_chrono_sampler_exact_1dx *temporal_sampler, double** ret_points_1, int *num_points) {
    std::vector<mui::point1dx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_exact_1t_f(mui_uniface_1t *uniface, const char *attr, double *t,
        mui_chrono_sampler_exact_1t *temporal_sampler, double** ret_points_1, int *num_points) {
    std::vector<mui::mui_f_wrapper_1D::point_type> ret_pts = uniface->fetch_points<float>(std::string(attr), *t,
            *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = static_cast<double>(ret_pts[i][0]);
    }
    *num_points = static_cast<int>(ret_pts.size());
}

// Temporal sampler: gauss
void mui_fetch_points_gauss_1f_f(mui_uniface_1f *uniface, const char *attr, float *t,
        mui_chrono_sampler_gauss_1f *temporal_sampler, float** ret_points_1, int *num_points) {
    std::vector<mui::point1f> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_gauss_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *t,
        mui_chrono_sampler_gauss_1fx *temporal_sampler, float** ret_points_1, int *num_points) {
    std::vector<mui::point1fx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_gauss_1d_f(mui_uniface_1d *uniface, const char *attr, double *t,
        mui_chrono_sampler_gauss_1d *temporal_sampler, double** ret_points_1, int *num_points) {
    std::vector<mui::point1d> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_gauss_1dx_f(mui_uniface_1dx *uniface, const char *attr, double *t,
        mui_chrono_sampler_gauss_1dx *temporal_sampler, double** ret_points_1, int *num_points) {
    std::vector<mui::point1dx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_gauss_1t_f(mui_uniface_1t *uniface, const char *attr, double *t,
        mui_chrono_sampler_gauss_1t *temporal_sampler, double** ret_points_1, int *num_points) {
    std::vector<mui::mui_f_wrapper_1D::point_type> ret_pts = uniface->fetch_points<float>(std::string(attr), *t,
            *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = static_cast<double>(ret_pts[i][0]);
    }
    *num_points = static_cast<int>(ret_pts.size());
}

// Temporal sampler: mean
void mui_fetch_points_mean_1f_f(mui_uniface_1f *uniface, const char *attr, float *t,
        mui_chrono_sampler_mean_1f *temporal_sampler, float** ret_points_1, int *num_points) {
    std::vector<mui::point1f> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_mean_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *t,
        mui_chrono_sampler_mean_1fx *temporal_sampler, float** ret_points_1, int *num_points) {
    std::vector<mui::point1fx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_mean_1d_f(mui_uniface_1d *uniface, const char *attr, double *t,
        mui_chrono_sampler_mean_1d *temporal_sampler, double** ret_points_1, int *num_points) {
    std::vector<mui::point1d> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_mean_1dx_f(mui_uniface_1dx *uniface, const char *attr, double *t,
        mui_chrono_sampler_mean_1dx *temporal_sampler, double** ret_points_1, int *num_points) {
    std::vector<mui::point1dx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_mean_1t_f(mui_uniface_1t *uniface, const char *attr, double *t,
        mui_chrono_sampler_mean_1t *temporal_sampler, double** ret_points_1, int *num_points) {
    std::vector<mui::mui_f_wrapper_1D::point_type> ret_pts = uniface->fetch_points<float>(std::string(attr), *t,
            *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = static_cast<double>(ret_pts[i][0]);
    }
    *num_points = static_cast<int>(ret_pts.size());
}

// Temporal sampler: sum
void mui_fetch_points_sum_1f_f(mui_uniface_1f *uniface, const char *attr, float *t,
        mui_chrono_sampler_sum_1f *temporal_sampler, float** ret_points_1, int *num_points) {
    std::vector<mui::point1f> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_sum_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *t,
        mui_chrono_sampler_sum_1fx *temporal_sampler, float** ret_points_1, int *num_points) {
    std::vector<mui::point1fx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_sum_1d_f(mui_uniface_1d *uniface, const char *attr, double *t,
        mui_chrono_sampler_sum_1d *temporal_sampler, double** ret_points_1, int *num_points) {
    std::vector<mui::point1d> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_sum_1dx_f(mui_uniface_1dx *uniface, const char *attr, double *t,
        mui_chrono_sampler_sum_1dx *temporal_sampler, double** ret_points_1, int *num_points) {
    std::vector<mui::point1dx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_sum_1t_f(mui_uniface_1t *uniface, const char *attr, double *t,
        mui_chrono_sampler_sum_1t *temporal_sampler, double** ret_points_1, int *num_points) {
    std::vector<mui::mui_f_wrapper_1D::point_type> ret_pts = uniface->fetch_points<float>(std::string(attr), *t,
            *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = static_cast<double>(ret_pts[i][0]);
    }
    *num_points = static_cast<int>(ret_pts.size());
}

/*******************************************************************
 * MUI functions for 1D data point only fetch using two time values *
 ********************************************************************/

// Temporal sampler: exact
void mui_fetch_points_exact_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_exact_1f *temporal_sampler, float** ret_points_1, int *num_points) {
    std::vector<mui::point1f> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_exact_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_exact_1fx *temporal_sampler, float** ret_points_1, int *num_points) {
    std::vector<mui::point1fx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_exact_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_exact_1d *temporal_sampler, double** ret_points_1, int *num_points) {
    std::vector<mui::point1d> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_exact_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_exact_1dx *temporal_sampler, double** ret_points_1, int *num_points) {
    std::vector<mui::point1dx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_exact_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_exact_1t *temporal_sampler, double** ret_points_1, int *num_points) {
    std::vector<mui::mui_f_wrapper_1D::point_type> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2,
            *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = static_cast<double>(ret_pts[i][0]);
    }
    *num_points = static_cast<int>(ret_pts.size());
}

// Temporal sampler: gauss
void mui_fetch_points_gauss_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_gauss_1f *temporal_sampler, float** ret_points_1, int *num_points) {
    std::vector<mui::point1f> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_gauss_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_gauss_1fx *temporal_sampler, float** ret_points_1, int *num_points) {
    std::vector<mui::point1fx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_gauss_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_gauss_1d *temporal_sampler, double** ret_points_1, int *num_points) {
    std::vector<mui::point1d> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_gauss_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_gauss_1dx *temporal_sampler, double** ret_points_1, int *num_points) {
    std::vector<mui::point1dx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_gauss_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_gauss_1t *temporal_sampler, double** ret_points_1, int *num_points) {
    std::vector<mui::mui_f_wrapper_1D::point_type> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2,
            *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = static_cast<double>(ret_pts[i][0]);
    }
    *num_points = static_cast<int>(ret_pts.size());
}

// Temporal sampler: mean
void mui_fetch_points_mean_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_mean_1f *temporal_sampler, float** ret_points_1, int *num_points) {
    std::vector<mui::point1f> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_mean_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_mean_1fx *temporal_sampler, float** ret_points_1, int *num_points) {
    std::vector<mui::point1fx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_mean_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_mean_1d *temporal_sampler, double** ret_points_1, int *num_points) {
    std::vector<mui::point1d> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_mean_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_mean_1dx *temporal_sampler, double** ret_points_1, int *num_points) {
    std::vector<mui::point1dx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_mean_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_mean_1t *temporal_sampler, double** ret_points_1, int *num_points) {
    std::vector<mui::mui_f_wrapper_1D::point_type> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2,
            *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = static_cast<double>(ret_pts[i][0]);
    }
    *num_points = static_cast<int>(ret_pts.size());
}

// Temporal sampler: sum
void mui_fetch_points_sum_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_sum_1f *temporal_sampler, float** ret_points_1, int *num_points) {
    std::vector<mui::point1f> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_sum_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_sum_1fx *temporal_sampler, float** ret_points_1, int *num_points) {
    std::vector<mui::point1fx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_sum_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_sum_1d *temporal_sampler, double** ret_points_1, int *num_points) {
    std::vector<mui::point1d> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_sum_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_sum_1dx *temporal_sampler, double** ret_points_1, int *num_points) {
    std::vector<mui::point1dx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_sum_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_sum_1t *temporal_sampler, double** ret_points_1, int *num_points) {
    std::vector<mui::mui_f_wrapper_1D::point_type> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2,
            *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = static_cast<double>(ret_pts[i][0]);
    }
    *num_points = static_cast<int>(ret_pts.size());
}

/*******************************************************************
 * MUI functions for 1D data values only fetch using one time value  *
 ********************************************************************/

// Temporal sampler: exact
void mui_fetch_values_exact_1f_f(mui_uniface_1f *uniface, const char *attr, float *t,
        mui_chrono_sampler_exact_1f *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_exact_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *t,
        mui_chrono_sampler_exact_1fx *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_exact_1d_f(mui_uniface_1d *uniface, const char *attr, double *t,
        mui_chrono_sampler_exact_1d *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_exact_1dx_f(mui_uniface_1dx *uniface, const char *attr, double *t,
        mui_chrono_sampler_exact_1dx *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_exact_1t_f(mui_uniface_1t *uniface, const char *attr, double *t,
        mui_chrono_sampler_exact_1t *temporal_sampler, double **values, int *num_values) {
    std::vector<mui::mui_f_wrapper_1D::REAL> ret_vals = uniface->fetch_values<mui::mui_f_wrapper_1D::REAL>(
            std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = static_cast<double>(ret_vals[i]);
    }
    *num_values = static_cast<int>(ret_vals.size());
}

// Temporal sampler: gauss
void mui_fetch_values_gauss_1f_f(mui_uniface_1f *uniface, const char *attr, float *t,
        mui_chrono_sampler_gauss_1f *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_gauss_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *t,
        mui_chrono_sampler_gauss_1fx *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_gauss_1d_f(mui_uniface_1d *uniface, const char *attr, double *t,
        mui_chrono_sampler_gauss_1d *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_gauss_1dx_f(mui_uniface_1dx *uniface, const char *attr, double *t,
        mui_chrono_sampler_gauss_1dx *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_gauss_1t_f(mui_uniface_1t *uniface, const char *attr, double *t,
        mui_chrono_sampler_gauss_1t *temporal_sampler, double **values, int *num_values) {
    std::vector<mui::mui_f_wrapper_1D::REAL> ret_vals = uniface->fetch_values<mui::mui_f_wrapper_1D::REAL>(
            std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = static_cast<double>(ret_vals[i]);
    }
    *num_values = static_cast<int>(ret_vals.size());
}

// Temporal sampler: mean
void mui_fetch_values_mean_1f_f(mui_uniface_1f *uniface, const char *attr, float *t,
        mui_chrono_sampler_mean_1f *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_mean_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *t,
        mui_chrono_sampler_mean_1fx *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_mean_1d_f(mui_uniface_1d *uniface, const char *attr, double *t,
        mui_chrono_sampler_mean_1d *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_mean_1dx_f(mui_uniface_1dx *uniface, const char *attr, double *t,
        mui_chrono_sampler_mean_1dx *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_mean_1t_f(mui_uniface_1t *uniface, const char *attr, double *t,
        mui_chrono_sampler_mean_1t *temporal_sampler, double **values, int *num_values) {
    std::vector<mui::mui_f_wrapper_1D::REAL> ret_vals = uniface->fetch_values<mui::mui_f_wrapper_1D::REAL>(
            std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = static_cast<double>(ret_vals[i]);
    }
    *num_values = static_cast<int>(ret_vals.size());
}

// Temporal sampler: sum
void mui_fetch_values_sum_1f_f(mui_uniface_1f *uniface, const char *attr, float *t,
        mui_chrono_sampler_sum_1f *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_sum_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *t,
        mui_chrono_sampler_sum_1fx *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_sum_1d_f(mui_uniface_1d *uniface, const char *attr, double *t,
        mui_chrono_sampler_sum_1d *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_sum_1dx_f(mui_uniface_1dx *uniface, const char *attr, double *t,
        mui_chrono_sampler_sum_1dx *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_sum_1t_f(mui_uniface_1t *uniface, const char *attr, double *t,
        mui_chrono_sampler_sum_1t *temporal_sampler, double **values, int *num_values) {
    std::vector<mui::mui_f_wrapper_1D::REAL> ret_vals = uniface->fetch_values<mui::mui_f_wrapper_1D::REAL>(
            std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = static_cast<double>(ret_vals[i]);
    }
    *num_values = static_cast<int>(ret_vals.size());
}

/********************************************************************
 * MUI functions for 1D data values only fetch using two time values *
 *********************************************************************/

// Temporal sampler: exact
void mui_fetch_values_exact_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_exact_1f *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_exact_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_exact_1fx *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_exact_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_exact_1d *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_exact_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_exact_1dx *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_exact_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_exact_1t *temporal_sampler, double **values, int *num_values) {
    std::vector<mui::mui_f_wrapper_1D::REAL> ret_vals = uniface->fetch_values<mui::mui_f_wrapper_1D::REAL>(
            std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = static_cast<double>(ret_vals[i]);
    }
    *num_values = static_cast<int>(ret_vals.size());
}

// Temporal sampler: gauss
void mui_fetch_values_gauss_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_gauss_1f *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_gauss_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_gauss_1fx *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_gauss_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_gauss_1d *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_gauss_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_gauss_1dx *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_gauss_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_gauss_1t *temporal_sampler, double **values, int *num_values) {
    std::vector<mui::mui_f_wrapper_1D::REAL> ret_vals = uniface->fetch_values<mui::mui_f_wrapper_1D::REAL>(
            std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = static_cast<double>(ret_vals[i]);
    }
    *num_values = static_cast<int>(ret_vals.size());
}

// Temporal sampler: mean
void mui_fetch_values_mean_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_mean_1f *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_mean_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_mean_1fx *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_mean_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_mean_1d *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_mean_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_mean_1dx *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_mean_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_mean_1t *temporal_sampler, double **values, int *num_values) {
    std::vector<mui::mui_f_wrapper_1D::REAL> ret_vals = uniface->fetch_values<mui::mui_f_wrapper_1D::REAL>(
            std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = static_cast<double>(ret_vals[i]);
    }
    *num_values = static_cast<int>(ret_vals.size());
}

// Temporal sampler: sum
void mui_fetch_values_sum_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_sum_1f *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_sum_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_sum_1fx *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_sum_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_sum_1d *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_sum_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_sum_1dx *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_sum_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_sum_1t *temporal_sampler, double **values, int *num_values) {
    std::vector<mui::mui_f_wrapper_1D::REAL> ret_vals = uniface->fetch_values<mui::mui_f_wrapper_1D::REAL>(
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

void mui_fetch_1f_param_f(mui_uniface_1f *uniface, const char *attr, float *return_value) {
    *return_value = uniface->fetch<float>(std::string(attr));
}

void mui_fetch_1fx_param_f(mui_uniface_1fx *uniface, const char *attr, float *return_value) {
    *return_value = uniface->fetch<float>(std::string(attr));
}

void mui_fetch_1d_param_f(mui_uniface_1d *uniface, const char *attr, float *return_value) {
    *return_value = uniface->fetch<double>(std::string(attr));
}

void mui_fetch_1dx_param_f(mui_uniface_1dx *uniface, const char *attr, float *return_value) {
    *return_value = uniface->fetch<double>(std::string(attr));
}

void mui_fetch_1t_param_f(mui_uniface_1t *uniface, const char *attr, float *return_value) {
    mui::mui_f_wrapper_1D::REAL ret_val = uniface->fetch<mui::mui_f_wrapper_1D::REAL>(std::string(attr));
    *return_value = static_cast<double>(ret_val);
}

/******************************************
 * MUI data receive test functions         *
 ******************************************/

// Data ready test using single time value
void mui_is_ready_1f_f(mui_uniface_1f *uniface, const char *attr, float *t, int *return_value) {
    *return_value = uniface->is_ready(std::string(attr), *t);
}

void mui_is_ready_1fx_f(mui_uniface_1fx *uniface, const char *attr, float *t, int *return_value) {
    *return_value = uniface->is_ready(std::string(attr), *t);
}

void mui_is_ready_1d_f(mui_uniface_1d *uniface, const char *attr, double *t, int *return_value) {
    *return_value = uniface->is_ready(std::string(attr), *t);
}

void mui_is_ready_1dx_f(mui_uniface_1dx *uniface, const char *attr, double *t, int *return_value) {
    *return_value = uniface->is_ready(std::string(attr), *t);
}

void mui_is_ready_1t_f(mui_uniface_1t *uniface, const char *attr, double *t, int *return_value) {
    *return_value = uniface->is_ready(std::string(attr), static_cast<mui::mui_f_wrapper_1D::time_type>(*t));
}

// Data ready test using two time values
void mui_is_ready_1f_pair_f(mui_uniface_1f *uniface, const char *attr, float *t_1, float *t_2, int *return_value) {
    *return_value = uniface->is_ready(std::string(attr), *t_1, *t_2);
}

void mui_is_ready_1fx_pair_f(mui_uniface_1fx *uniface, const char *attr, float *t_1, float *t_2, int *return_value) {
    *return_value = uniface->is_ready(std::string(attr), *t_1, *t_2);
}

void mui_is_ready_1d_pair_f(mui_uniface_1d *uniface, const char *attr, double *t_1, double *t_2, int *return_value) {
    *return_value = uniface->is_ready(std::string(attr), *t_1, *t_2);
}

void mui_is_ready_1dx_pair_f(mui_uniface_1dx *uniface, const char *attr, double *t_1, double *t_2, int *return_value) {
    *return_value = uniface->is_ready(std::string(attr), *t_1, *t_2);
}

void mui_is_ready_1t_pair_f(mui_uniface_1t *uniface, const char *attr, double *t_1, double *t_2, int *return_value) {
    *return_value = uniface->is_ready(std::string(attr), static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1),
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2));
}

/******************************************
 * MUI Smart Send functions                *
 ******************************************/

// Send span announce using 1D box geometry
void mui_announce_send_span_1f_box_f(mui_uniface_1f *uniface, float *box_1_1, float *box_2_1, float *t_start,
        float *t_timeout, int *synchronised) {
    mui::point1f point_1(*box_1_1);
    mui::point1f point_2(*box_2_1);
    mui::geometry::box1f bound_box(point_1, point_2);
    uniface->announce_send_span(*t_start, *t_timeout, bound_box, static_cast<bool>(*synchronised));
}

void mui_announce_send_span_1fx_box_f(mui_uniface_1fx *uniface, float *box_1_1, float *box_2_1, float *t_start,
        float *t_timeout, int *synchronised) {
    mui::point1fx point_1(*box_1_1);
    mui::point1fx point_2(*box_2_1);
    mui::geometry::box1fx bound_box(point_1, point_2);
    uniface->announce_send_span(*t_start, *t_timeout, bound_box, static_cast<bool>(*synchronised));
}

void mui_announce_send_span_1d_box_f(mui_uniface_1d *uniface, double *box_1_1, double *box_2_1, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::point1d point_1(*box_1_1);
    mui::point1d point_2(*box_2_1);
    mui::geometry::box1d bound_box(point_1, point_2);
    uniface->announce_send_span(*t_start, *t_timeout, bound_box, static_cast<bool>(*synchronised));
}

void mui_announce_send_span_1dx_box_f(mui_uniface_1dx *uniface, double *box_1_1, double *box_2_1, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::point1dx point_1(*box_1_1);
    mui::point1dx point_2(*box_2_1);
    mui::geometry::box1dx bound_box(point_1, point_2);
    uniface->announce_send_span(*t_start, *t_timeout, bound_box, static_cast<bool>(*synchronised));
}

void mui_announce_send_span_1t_box_f(mui_uniface_1t *uniface, double *box_1_1, double *box_2_1, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::mui_f_wrapper_1D::point_type point_1(static_cast<mui::mui_f_wrapper_1D::REAL>(*box_1_1));
    mui::mui_f_wrapper_1D::point_type point_2(static_cast<mui::mui_f_wrapper_1D::REAL>(*box_2_1));
    mui::geometry::box<mui::mui_f_wrapper_1D> bound_box(point_1, point_2);
    uniface->announce_send_span(static_cast<mui::mui_f_wrapper_1D::time_type>(*t_start),
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_timeout), bound_box, static_cast<bool>(*synchronised));
}

// Send span announce using 1D sphere geometry
void mui_announce_send_span_1f_sphere_f(mui_uniface_1f *uniface, float *centre_1, float *radius, float *t_start,
        float *t_timeout, int *synchronised) {
    mui::geometry::sphere1f bound_sphere(mui::point1f(*centre_1), *radius);
    uniface->announce_send_span(*t_start, *t_timeout, bound_sphere, static_cast<bool>(*synchronised));
}

void mui_announce_send_span_1fx_sphere_f(mui_uniface_1fx *uniface, float *centre_1, float *radius, float *t_start,
        float *t_timeout, int *synchronised) {
    mui::geometry::sphere1fx bound_sphere(mui::point1fx(*centre_1), *radius);
    uniface->announce_send_span(*t_start, *t_timeout, bound_sphere, static_cast<bool>(*synchronised));
}

void mui_announce_send_span_1d_sphere_f(mui_uniface_1d *uniface, double *centre_1, double *radius, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::geometry::sphere1d bound_sphere(mui::point1d(*centre_1), *radius);
    uniface->announce_send_span(*t_start, *t_timeout, bound_sphere, static_cast<bool>(*synchronised));
}

void mui_announce_send_span_1dx_sphere_f(mui_uniface_1dx *uniface, double *centre_1, double *radius, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::geometry::sphere1dx bound_sphere(mui::point1dx(*centre_1), *radius);
    uniface->announce_send_span(*t_start, *t_timeout, bound_sphere, static_cast<bool>(*synchronised));
}

void mui_announce_send_span_1t_sphere_f(mui_uniface_1t *uniface, double *centre_1, double *radius, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::geometry::sphere<mui::mui_f_wrapper_1D> bound_sphere(
            mui::mui_f_wrapper_1D::point_type(static_cast<mui::mui_f_wrapper_1D::REAL>(*centre_1)),
            static_cast<mui::mui_f_wrapper_1D::REAL>(*radius));
    uniface->announce_send_span(static_cast<mui::mui_f_wrapper_1D::time_type>(*t_start),
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_timeout), bound_sphere, static_cast<bool>(*synchronised));
}

// Receive span announce using 1D box geometry
void mui_announce_recv_span_1f_box_f(mui_uniface_1f *uniface, float *box_1_1, float *box_2_1, float *t_start,
        float *t_timeout, int *synchronised) {
    mui::point1f point_1(*box_1_1);
    mui::point1f point_2(*box_2_1);
    mui::geometry::box1f bound_box(point_1, point_2);
    uniface->announce_recv_span(*t_start, *t_timeout, bound_box, static_cast<bool>(*synchronised));
}

void mui_announce_recv_span_1fx_box_f(mui_uniface_1fx *uniface, float *box_1_1, float *box_2_1, float *t_start,
        float *t_timeout, int *synchronised) {
    mui::point1fx point_1(*box_1_1);
    mui::point1fx point_2(*box_2_1);
    mui::geometry::box1fx bound_box(point_1, point_2);
    uniface->announce_recv_span(*t_start, *t_timeout, bound_box, static_cast<bool>(*synchronised));
}

void mui_announce_recv_span_1d_box_f(mui_uniface_1d *uniface, double *box_1_1, double *box_2_1, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::point1d point_1(*box_1_1);
    mui::point1d point_2(*box_2_1);
    mui::geometry::box1d bound_box(point_1, point_2);
    uniface->announce_recv_span(*t_start, *t_timeout, bound_box, static_cast<bool>(*synchronised));
}

void mui_announce_recv_span_1dx_box_f(mui_uniface_1dx *uniface, double *box_1_1, double *box_2_1, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::point1dx point_1(*box_1_1);
    mui::point1dx point_2(*box_2_1);
    mui::geometry::box1dx bound_box(point_1, point_2);
    uniface->announce_recv_span(*t_start, *t_timeout, bound_box, static_cast<bool>(*synchronised));
}

void mui_announce_recv_span_1t_box_f(mui_uniface_1t *uniface, double *box_1_1, double *box_2_1, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::mui_f_wrapper_1D::point_type point_1(static_cast<mui::mui_f_wrapper_1D::REAL>(*box_1_1));
    mui::mui_f_wrapper_1D::point_type point_2(static_cast<mui::mui_f_wrapper_1D::REAL>(*box_2_1));
    mui::geometry::box<mui::mui_f_wrapper_1D> bound_box(point_1, point_2);
    uniface->announce_recv_span(static_cast<mui::mui_f_wrapper_1D::time_type>(*t_start),
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_timeout), bound_box, static_cast<bool>(*synchronised));
}

// Receive span announce using 1D sphere geometry
void mui_announce_recv_span_1f_sphere_f(mui_uniface_1f *uniface, float *centre_1, float *radius, float *t_start,
        float *t_timeout, int *synchronised) {
    mui::geometry::sphere1f bound_sphere(mui::point1f(*centre_1), *radius);
    uniface->announce_recv_span(*t_start, *t_timeout, bound_sphere, static_cast<bool>(*synchronised));
}

void mui_announce_recv_span_1fx_sphere_f(mui_uniface_1fx *uniface, float *centre_1, float *radius, float *t_start,
        float *t_timeout, int *synchronised) {
    mui::geometry::sphere1fx bound_sphere(mui::point1fx(*centre_1), *radius);
    uniface->announce_recv_span(*t_start, *t_timeout, bound_sphere, static_cast<bool>(*synchronised));
}

void mui_announce_recv_span_1d_sphere_f(mui_uniface_1d *uniface, double *centre_1, double *radius, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::geometry::sphere1d bound_sphere(mui::point1d(*centre_1), *radius);
    uniface->announce_recv_span(*t_start, *t_timeout, bound_sphere, static_cast<bool>(*synchronised));
}

void mui_announce_recv_span_1dx_sphere_f(mui_uniface_1dx *uniface, double *centre_1, double *radius, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::geometry::sphere1dx bound_sphere(mui::point1dx(*centre_1), *radius);
    uniface->announce_recv_span(*t_start, *t_timeout, bound_sphere, static_cast<bool>(*synchronised));
}

void mui_announce_recv_span_1t_sphere_f(mui_uniface_1t *uniface, double *centre_1, double *radius, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::geometry::sphere<mui::mui_f_wrapper_1D> bound_sphere(
            mui::mui_f_wrapper_1D::point_type(static_cast<mui::mui_f_wrapper_1D::REAL>(*centre_1)),
            static_cast<mui::mui_f_wrapper_1D::REAL>(*radius));
    uniface->announce_recv_span(static_cast<mui::mui_f_wrapper_1D::time_type>(*t_start),
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_timeout), bound_sphere, static_cast<bool>(*synchronised));
}

// Send disable announce (local call per MPI rank)
void mui_announce_send_disable_1f_f(mui_uniface_1f *uniface, int *synchronised) {
    uniface->announce_send_disable(static_cast<bool>(*synchronised));
}

void mui_announce_send_disable_1fx_f(mui_uniface_1fx *uniface, int *synchronised) {
    uniface->announce_send_disable(static_cast<bool>(*synchronised));
}

void mui_announce_send_disable_1d_f(mui_uniface_1d *uniface, int *synchronised) {
    uniface->announce_send_disable(static_cast<bool>(*synchronised));
}

void mui_announce_send_disable_1dx_f(mui_uniface_1dx *uniface, int *synchronised) {
    uniface->announce_send_disable(static_cast<bool>(*synchronised));
}

void mui_announce_send_disable_1t_f(mui_uniface_1t *uniface, int *synchronised) {
    uniface->announce_send_disable(static_cast<bool>(*synchronised));
}

// Receive disable announce (local call per MPI rank)
void mui_announce_recv_disable_1f_f(mui_uniface_1f *uniface, int *synchronised) {
    uniface->announce_recv_disable(static_cast<bool>(*synchronised));
}

void mui_announce_recv_disable_1fx_f(mui_uniface_1fx *uniface, int *synchronised) {
    uniface->announce_recv_disable(static_cast<bool>(*synchronised));
}

void mui_announce_recv_disable_1d_f(mui_uniface_1d *uniface, int *synchronised) {
    uniface->announce_recv_disable(static_cast<bool>(*synchronised));
}

void mui_announce_recv_disable_1dx_f(mui_uniface_1dx *uniface, int *synchronised) {
    uniface->announce_recv_disable(static_cast<bool>(*synchronised));
}

void mui_announce_recv_disable_1t_f(mui_uniface_1t *uniface, int *synchronised) {
    uniface->announce_recv_disable(static_cast<bool>(*synchronised));
}

/******************************************
 * MUI barrier functions                   *
 ******************************************/

// Barrier at single time value
void mui_barrier_1f_f(mui_uniface_1f *uniface, float *t) {
    uniface->barrier(*t);
}

void mui_barrier_1fx_f(mui_uniface_1fx *uniface, float *t) {
    uniface->barrier(*t);
}

void mui_barrier_1d_f(mui_uniface_1d *uniface, double *t) {
    uniface->barrier(*t);
}

void mui_barrier_1dx_f(mui_uniface_1dx *uniface, double *t) {
    uniface->barrier(*t);
}

void mui_barrier_1t_f(mui_uniface_1t *uniface, double *t) {
    uniface->barrier(static_cast<mui::mui_f_wrapper_1D::time_type>(*t));
}

// Barrier at two time values
void mui_barrier_1f_pair_f(mui_uniface_1f *uniface, float *t_1, float *t_2) {
    uniface->barrier(*t_1, *t_2);
}

void mui_barrier_1fx_pair_f(mui_uniface_1fx *uniface, float *t_1, float *t_2) {
    uniface->barrier(*t_1, *t_2);
}

void mui_barrier_1d_pair_f(mui_uniface_1d *uniface, double *t_1, double *t_2) {
    uniface->barrier(*t_1, *t_2);
}

void mui_barrier_1dx_pair_f(mui_uniface_1dx *uniface, double *t_1, double *t_2) {
    uniface->barrier(*t_1, *t_2);
}

void mui_barrier_1t_pair_f(mui_uniface_1t *uniface, double *t_1, double *t_2) {
    uniface->barrier(static_cast<mui::mui_f_wrapper_1D::time_type>(*t_1),
            static_cast<mui::mui_f_wrapper_1D::time_type>(*t_2));
}

/******************************************
 * MUI forget functions                    *
 ******************************************/

// Forget log between [-inf, upper]
void mui_forget_upper_1f_f(mui_uniface_1f *uniface, float *upper, int *reset_log) {
    uniface->forget(*upper, static_cast<bool>(*reset_log));
}

void mui_forget_upper_1fx_f(mui_uniface_1fx *uniface, float *upper, int *reset_log) {
    uniface->forget(*upper, static_cast<bool>(*reset_log));
}

void mui_forget_upper_1d_f(mui_uniface_1d *uniface, double *upper, int *reset_log) {
    uniface->forget(*upper, static_cast<bool>(*reset_log));
}

void mui_forget_upper_1dx_f(mui_uniface_1dx *uniface, double *upper, int *reset_log) {
    uniface->forget(*upper, static_cast<bool>(*reset_log));
}

void mui_forget_upper_1t_f(mui_uniface_1t *uniface, double *upper, int *reset_log) {
    uniface->forget(static_cast<mui::mui_f_wrapper_1D::time_type>(*upper), static_cast<bool>(*reset_log));
}

// Forget log between [-inf, -inf], [upper_1, upper_2]
void mui_forget_upper_1f_pair_f(mui_uniface_1f *uniface, float *upper_1, float *upper_2, int *reset_log) {
    std::pair<float, float> forget_time(*upper_1, *upper_2);
    uniface->forget(forget_time, static_cast<bool>(*reset_log));
}

void mui_forget_upper_1fx_pair_f(mui_uniface_1fx *uniface, float *upper_1, float *upper_2, int *reset_log) {
    std::pair<float, float> forget_time(*upper_1, *upper_2);
    uniface->forget(forget_time, static_cast<bool>(*reset_log));
}

void mui_forget_upper_1d_pair_f(mui_uniface_1d *uniface, double *upper_1, double *upper_2, int *reset_log) {
    std::pair<double, double> forget_time(*upper_1, *upper_2);
    uniface->forget(forget_time, static_cast<bool>(*reset_log));
}

void mui_forget_upper_1dx_pair_f(mui_uniface_1dx *uniface, double *upper_1, double *upper_2, int *reset_log) {
    std::pair<double, double> forget_time(*upper_1, *upper_2);
    uniface->forget(forget_time, static_cast<bool>(*reset_log));
}

void mui_forget_upper_1t_pair_f(mui_uniface_1t *uniface, double *upper_1, double *upper_2, int *reset_log) {
    mui::mui_f_wrapper_1D::time_type time_1 = static_cast<mui::mui_f_wrapper_1D::time_type>(*upper_1);
    mui::mui_f_wrapper_1D::time_type time_2 = static_cast<mui::mui_f_wrapper_1D::time_type>(*upper_2);
    std::pair<mui::mui_f_wrapper_1D::time_type, mui::mui_f_wrapper_1D::time_type> forget_time(time_1, time_2);
    uniface->forget(forget_time, static_cast<bool>(*reset_log));
}

// Forget log between [lower, upper]
void mui_forget_lower_upper_1f_f(mui_uniface_1f *uniface, float *lower, float *upper, int *reset_log) {
    uniface->forget(*lower, *upper, static_cast<bool>(*reset_log));
}

void mui_forget_lower_upper_1fx_f(mui_uniface_1fx *uniface, float *lower, float *upper, int *reset_log) {
    uniface->forget(*lower, *upper, static_cast<bool>(*reset_log));
}

void mui_forget_lower_upper_1d_f(mui_uniface_1d *uniface, double *lower, double *upper, int *reset_log) {
    uniface->forget(*lower, *upper, static_cast<bool>(*reset_log));
}

void mui_forget_lower_upper_1dx_f(mui_uniface_1dx *uniface, double *lower, double *upper, int *reset_log) {
    uniface->forget(*lower, *upper, static_cast<bool>(*reset_log));
}

void mui_forget_lower_upper_1t_f(mui_uniface_1t *uniface, double *lower, double *upper, int *reset_log) {
    mui::mui_f_wrapper_1D::time_type forget_time_1 = static_cast<mui::mui_f_wrapper_1D::time_type>(*lower);
    mui::mui_f_wrapper_1D::time_type forget_time_2 = static_cast<mui::mui_f_wrapper_1D::time_type>(*upper);
    uniface->forget(forget_time_1, forget_time_2, static_cast<bool>(*reset_log));
}

// Forget log between [lower_1, lower_2], [upper_1, upper_2]
void mui_forget_lower_upper_1f_pair_f(mui_uniface_1f *uniface, float *lower_1, float *lower_2, float *upper_1, float *upper_2,
        int *reset_log) {
    std::pair<float, float> forget_time_1(*lower_1, *lower_2);
    std::pair<float, float> forget_time_2(*upper_1, *upper_2);
    uniface->forget(forget_time_1, forget_time_2, static_cast<bool>(*reset_log));
}

void mui_forget_lower_upper_1fx_pair_f(mui_uniface_1fx *uniface, float *lower_1, float *lower_2, float *upper_1,
        float *upper_2, int *reset_log) {
    std::pair<float, float> forget_time_1(*lower_1, *lower_2);
    std::pair<float, float> forget_time_2(*upper_1, *upper_2);
    uniface->forget(forget_time_1, forget_time_2, static_cast<bool>(*reset_log));
}

void mui_forget_lower_upper_1d_pair_f(mui_uniface_1d *uniface, double *lower_1, double *lower_2, double *upper_1,
        double *upper_2, int *reset_log) {
    std::pair<double, double> forget_time_1(*lower_1, *lower_2);
    std::pair<double, double> forget_time_2(*upper_1, *upper_2);
    uniface->forget(forget_time_1, forget_time_2, static_cast<bool>(*reset_log));
}

void mui_forget_lower_upper_1dx_pair_f(mui_uniface_1dx *uniface, double *lower_1, double *lower_2, double *upper_1,
        double *upper_2, int *reset_log) {
    std::pair<double, double> forget_time_1(*lower_1, *lower_2);
    std::pair<double, double> forget_time_2(*upper_1, *upper_2);
    uniface->forget(forget_time_1, forget_time_2, static_cast<bool>(*reset_log));
}

void mui_forget_lower_upper_1t_pair_f(mui_uniface_1t *uniface, double *lower_1, double *lower_2, double *upper_1,
        double *upper_2, int *reset_log) {
    mui::mui_f_wrapper_1D::time_type time_1 = static_cast<mui::mui_f_wrapper_1D::time_type>(*lower_1);
    mui::mui_f_wrapper_1D::time_type time_2 = static_cast<mui::mui_f_wrapper_1D::time_type>(*lower_2);
    mui::mui_f_wrapper_1D::time_type time_3 = static_cast<mui::mui_f_wrapper_1D::time_type>(*upper_1);
    mui::mui_f_wrapper_1D::time_type time_4 = static_cast<mui::mui_f_wrapper_1D::time_type>(*upper_2);
    std::pair<mui::mui_f_wrapper_1D::time_type, mui::mui_f_wrapper_1D::time_type> forget_time_1(time_1, time_2);
    std::pair<mui::mui_f_wrapper_1D::time_type, mui::mui_f_wrapper_1D::time_type> forget_time_2(time_3, time_4);
    uniface->forget(forget_time_1, forget_time_2, static_cast<bool>(*reset_log));
}

// Set to forget log between [-inf, current-length] automatically
void mui_set_forget_length_1f_f(mui_uniface_1f *uniface, float *length) {
    uniface->set_memory(*length);
}

void mui_set_forget_length_1fx_f(mui_uniface_1fx *uniface, float *length) {
    uniface->set_memory(*length);
}

void mui_set_forget_length_1d_f(mui_uniface_1d *uniface, double *length) {
    uniface->set_memory(*length);
}

void mui_set_forget_length_1dx_f(mui_uniface_1dx *uniface, double *length) {
    uniface->set_memory(*length);
}

void mui_set_forget_length_1t_f(mui_uniface_1t *uniface, double *length) {
    uniface->set_memory(*length);
}

/******************************************
 * MUI URI functions                      *
 ******************************************/

// Obtain original URI host value from existing interface
void mui_uri_host_1f_f(mui_uniface_1f *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_host().c_str());
}

void mui_uri_host_1fx_f(mui_uniface_1fx *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_host().c_str());
}

void mui_uri_host_1d_f(mui_uniface_1d *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_host().c_str());
}

void mui_uri_host_1dx_f(mui_uniface_1dx *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_host().c_str());
}

void mui_uri_host_1t_f(mui_uniface_1t *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_host().c_str());
}

// Obtain original URI path value from existing interface
void mui_uri_path_1f_f(mui_uniface_1f *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_path().c_str());
}

void mui_uri_path_1fx_f(mui_uniface_1fx *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_path().c_str());
}

void mui_uri_path_1d_f(mui_uniface_1d *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_path().c_str());
}

void mui_uri_path_1dx_f(mui_uniface_1dx *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_path().c_str());
}

void mui_uri_path_1t_f(mui_uniface_1t *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_path().c_str());
}

// Obtain original URI protocol value from existing interface
void mui_uri_protocol_1f_f(mui_uniface_1f *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_protocol().c_str());
}

void mui_uri_protocol_1fx_f(mui_uniface_1fx *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_protocol().c_str());
}

void mui_uri_protocol_1d_f(mui_uniface_1d *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_protocol().c_str());
}

void mui_uri_protocol_1dx_f(mui_uniface_1dx *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_protocol().c_str());
}

void mui_uri_protocol_1t_f(mui_uniface_1t *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_protocol().c_str());
}

}
