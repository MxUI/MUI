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
 * @file mui_f_wrapper_3d.cpp
 * @author S. M. Longshaw (derived from original 3D wrapper by S. Kudo)
 * @date Dec 08, 2021
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

// MUI interface typedefs for specialism creation
typedef mui::uniface3f mui_uniface_3f;
typedef mui::uniface3fx mui_uniface_3fx;
typedef mui::uniface3d mui_uniface_3d;
typedef mui::uniface3dx mui_uniface_3dx;

// MUI interface typedefs for template creation (recommended)
typedef mui::uniface<mui::mui_f_wrapper_3D> mui_uniface_3t;

// Exact spatial sampler typedefs for specialism creation
typedef mui::sampler_exact3f<float> mui_sampler_exact_3f;
typedef mui::sampler_exact3fx<float> mui_sampler_exact_3fx;
typedef mui::sampler_exact3d<double> mui_sampler_exact_3d;
typedef mui::sampler_exact3dx<double> mui_sampler_exact_3dx;

// Exact spatial sampler typedef for template creation (recommended)
typedef mui::sampler_exact<mui::mui_f_wrapper_3D> mui_sampler_exact_3t;

// Gaussian spatial sampler typedefs for specialism creation
typedef mui::sampler_gauss3f<float> mui_sampler_gauss_3f;
typedef mui::sampler_gauss3fx<float> mui_sampler_gauss_3fx;
typedef mui::sampler_gauss3d<double> mui_sampler_gauss_3d;
typedef mui::sampler_gauss3dx<double> mui_sampler_gauss_3dx;

// Gaussian spatial sampler typedef for template creation (recommended)
typedef mui::sampler_gauss<mui::mui_f_wrapper_3D> mui_sampler_gauss_3t;

// Moving average spatial sampler typedefs for specialism creation
typedef mui::sampler_moving_average3f<float> mui_sampler_moving_average_3f;
typedef mui::sampler_moving_average3fx<float> mui_sampler_moving_average_3fx;
typedef mui::sampler_moving_average3d<double> mui_sampler_moving_average_3d;
typedef mui::sampler_moving_average3dx<double> mui_sampler_moving_average_3dx;

// Moving average spatial sampler typedef for template creation (recommended)
typedef mui::sampler_moving_average<mui::mui_f_wrapper_3D> mui_sampler_moving_average_3t;

// Nearest neighbor spatial sampler typedefs for specialism creation
typedef mui::sampler_nearest_neighbor3f<float> mui_sampler_nearest_neighbor_3f;
typedef mui::sampler_nearest_neighbor3fx<float> mui_sampler_nearest_neighbor_3fx;
typedef mui::sampler_nearest_neighbor3d<double> mui_sampler_nearest_neighbor_3d;
typedef mui::sampler_nearest_neighbor3dx<double> mui_sampler_nearest_neighbor_3dx;

// Moving average spatial sampler typedef for template creation (recommended)
typedef mui::sampler_nearest_neighbor<mui::mui_f_wrapper_3D> mui_sampler_nearest_neighbor_3t;

// Pseudo-linear n^2 spatial sampler typedefs for specialism creation
typedef mui::sampler_pseudo_n2_linear3f<float> mui_sampler_pseudo_n2_linear_3f;
typedef mui::sampler_pseudo_n2_linear3fx<float> mui_sampler_pseudo_n2_linear_3fx;
typedef mui::sampler_pseudo_n2_linear3d<double> mui_sampler_pseudo_n2_linear_3d;
typedef mui::sampler_pseudo_n2_linear3dx<double> mui_sampler_pseudo_n2_linear_3dx;

// Pseudo-linear n^2 spatial sampler typedef for template creation (recommended)
typedef mui::sampler_pseudo_n2_linear<mui::mui_f_wrapper_3D> mui_sampler_pseudo_n2_linear_3t;

// Pseudo-nearest neighbor spatial sampler typedefs for specialism creation
typedef mui::sampler_pseudo_nearest_neighbor3f<float> mui_sampler_pseudo_nearest_neighbor_3f;
typedef mui::sampler_pseudo_nearest_neighbor3fx<float> mui_sampler_pseudo_nearest_neighbor_3fx;
typedef mui::sampler_pseudo_nearest_neighbor3d<double> mui_sampler_pseudo_nearest_neighbor_3d;
typedef mui::sampler_pseudo_nearest_neighbor3dx<double> mui_sampler_pseudo_nearest_neighbor_3dx;

// Pseudo-nearest neighbor spatial sampler typedef for template creation (recommended)
typedef mui::sampler_pseudo_nearest_neighbor<mui::mui_f_wrapper_3D> mui_sampler_pseudo_nearest_neighbor_3t;

// Shepard interpolation with quintic kernel spatial sampler typedefs for specialism creation
typedef mui::sampler_shepard_quintic3f<float> mui_sampler_shepard_quintic_3f;
typedef mui::sampler_shepard_quintic3fx<float> mui_sampler_shepard_quintic_3fx;
typedef mui::sampler_shepard_quintic3d<double> mui_sampler_shepard_quintic_3d;
typedef mui::sampler_shepard_quintic3dx<double> mui_sampler_shepard_quintic_3dx;

// Shepard interpolation with quintic kernel spatial sampler typedef for template creation (recommended)
typedef mui::sampler_shepard_quintic<mui::mui_f_wrapper_3D> mui_sampler_shepard_quintic_3t;

// Spatial sampler that provides a value at a point using a Smoothed Particle Hydrodynamics (SPH) derived
// interpolation method with a quintic spline kernel spatial sampler typedefs for specialism creation
typedef mui::sampler_sph_quintic3f<float> mui_sampler_sph_quintic_3f;
typedef mui::sampler_sph_quintic3fx<float> mui_sampler_sph_quintic_3fx;
typedef mui::sampler_sph_quintic3d<double> mui_sampler_sph_quintic_3d;
typedef mui::sampler_sph_quintic3dx<double> mui_sampler_sph_quintic_3dx;

// Spatial sampler that provides a value at a point using a Smoothed Particle Hydrodynamics (SPH) derived
// interpolation method with a quintic spline kernel spatial sampler typedefs for template creation (recommended)
typedef mui::sampler_sph_quintic<mui::mui_f_wrapper_3D> mui_sampler_sph_quintic_3t;

// Summation with quintic kernel spatial sampler typedefs for specialism creation
typedef mui::sampler_sum_quintic3f<float> mui_sampler_sum_quintic_3f;
typedef mui::sampler_sum_quintic3fx<float> mui_sampler_sum_quintic_3fx;
typedef mui::sampler_sum_quintic3d<double> mui_sampler_sum_quintic_3d;
typedef mui::sampler_sum_quintic3dx<double> mui_sampler_sum_quintic_3dx;

// Summation with quintic kernel spatial sampler typedef for template creation (recommended)
typedef mui::sampler_sum_quintic<mui::mui_f_wrapper_3D> mui_sampler_sum_quintic_3t;

#ifdef USE_RBF
// Radial Basis Function (RBF) spatial sampler typedefs for specialism creation
typedef mui::sampler_rbf3f<float> mui_sampler_rbf_3f;
typedef mui::sampler_rbf3fx<float> mui_sampler_rbf_3fx;
typedef mui::sampler_rbf3d<double> mui_sampler_rbf_3d;
typedef mui::sampler_rbf3dx<double> mui_sampler_rbf_3dx;

// Radial Basis Function (RBF) spatial sampler typedef for template creation (recommended)
typedef mui::sampler_rbf<mui::mui_f_wrapper_3D> mui_sampler_rbf_3t;
#endif

// Exact temporal sampler typedefs for specialism creation
typedef mui::chrono_sampler_exact3f mui_chrono_sampler_exact_3f;
typedef mui::chrono_sampler_exact3fx mui_chrono_sampler_exact_3fx;
typedef mui::chrono_sampler_exact3d mui_chrono_sampler_exact_3d;
typedef mui::chrono_sampler_exact3dx mui_chrono_sampler_exact_3dx;

// Exact temporal sampler typedefs for template creation (recommended)
typedef mui::chrono_sampler_exact<mui::mui_f_wrapper_3D> mui_chrono_sampler_exact_3t;

// Gaussian temporal sampler typedefs for specialism creation
typedef mui::chrono_sampler_gauss3f mui_chrono_sampler_gauss_3f;
typedef mui::chrono_sampler_gauss3fx mui_chrono_sampler_gauss_3fx;
typedef mui::chrono_sampler_gauss3d mui_chrono_sampler_gauss_3d;
typedef mui::chrono_sampler_gauss3dx mui_chrono_sampler_gauss_3dx;

// Gaussian temporal sampler typedefs for template creation (recommended)
typedef mui::chrono_sampler_gauss<mui::mui_f_wrapper_3D> mui_chrono_sampler_gauss_3t;

// Mean average temporal sampler typedefs for specialism creation
typedef mui::chrono_sampler_mean3f mui_chrono_sampler_mean_3f;
typedef mui::chrono_sampler_mean3fx mui_chrono_sampler_mean_3fx;
typedef mui::chrono_sampler_mean3d mui_chrono_sampler_mean_3d;
typedef mui::chrono_sampler_mean3dx mui_chrono_sampler_mean_3dx;

// Mean average temporal sampler typedefs for template creation (recommended)
typedef mui::chrono_sampler_mean<mui::mui_f_wrapper_3D> mui_chrono_sampler_mean_3t;

// Summation temporal sampler typedefs for specialism creation
typedef mui::chrono_sampler_sum3f mui_chrono_sampler_sum_3f;
typedef mui::chrono_sampler_sum3fx mui_chrono_sampler_sum_3fx;
typedef mui::chrono_sampler_sum3d mui_chrono_sampler_sum_3d;
typedef mui::chrono_sampler_sum3dx mui_chrono_sampler_sum_3dx;

// Summation temporal sampler typedefs for template creation (recommended)
typedef mui::chrono_sampler_sum<mui::mui_f_wrapper_3D> mui_chrono_sampler_sum_3t;

// MUI set of specialism interface on multi-domain coupling
mui_uniface_3f** mui_uniface_multi_3f;
mui_uniface_3fx** mui_uniface_multi_3fx;
mui_uniface_3d** mui_uniface_multi_3d;
mui_uniface_3dx** mui_uniface_multi_3dx;

// MUI set of template interface on multi-domain coupling
mui_uniface_3t** mui_uniface_multi_3t;

/****************************************
 * Create MUI interfaces                 *
 ****************************************/

// 3d interface with float=single and int=int32
void mui_create_uniface_3f_f(mui_uniface_3f **ret, const char *URI) {
    *ret = new mui_uniface_3f(URI);
}

// 3d interface with float=single and int=int64
void mui_create_uniface_3fx_f(mui_uniface_3fx **ret, const char *URI) {
    *ret = new mui_uniface_3fx(URI);
}

// 3d interface with float=double and int=int32
void mui_create_uniface_3d_f(mui_uniface_3d **ret, const char *URI) {
    *ret = new mui_uniface_3d(URI);
}

// 3d interface with float=double and int=int64
void mui_create_uniface_3dx_f(mui_uniface_3dx **ret, const char *URI) {
    *ret = new mui_uniface_3dx(URI);
}

// 3d interface using config from config_f_wrapper.h
void mui_create_uniface_3t_f(mui_uniface_3t **ret, const char *URI) {
    *ret = new mui_uniface_3t(URI);
}

// Set of 3D interfaces with float=single and int=int32
void mui_create_uniface_multi_3f_f(const char *domain, const char **interfaces, int interface_count ) {

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
                std::cerr << "MUI Error [mui_f_wrapper_3d.cpp]: Error MUI Fortran C binding: NULL interfaces at the " << i << "th interface." << std::endl;
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
        std::cerr << "MUI Error [mui_f_wrapper_3d.cpp]: Error MUI Fortran C binding: The size of interfaces array must larger or equals to one." << std::endl;
        std::abort();
    }

    auto created_unifaces = mui::create_uniface<mui::config_3f>(domain, interface_names);

    mui_uniface_multi_3f = new mui_uniface_3f*[created_unifaces.size()];

      for(size_t i=0; i<created_unifaces.size(); i++) {
        mui_uniface_multi_3f[i] = created_unifaces[i].release();
    }
}

// Set of 3D interfaces with float=single and int=int64
void mui_create_uniface_multi_3fx_f(const char *domain, const char **interfaces, int interface_count ) {

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
                std::cerr << "MUI Error [mui_f_wrapper_3d.cpp]: Error MUI Fortran C binding: NULL interfaces at the " << i << "th interface." << std::endl;
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
        std::cerr << "MUI Error [mui_f_wrapper_3d.cpp]: Error MUI Fortran C binding: The size of interfaces array must larger or equals to one." << std::endl;
        std::abort();
    }

    auto created_unifaces = mui::create_uniface<mui::config_3fx>(domain, interface_names);

    mui_uniface_multi_3fx = new mui_uniface_3fx*[created_unifaces.size()];

      for(size_t i=0; i<created_unifaces.size(); i++) {
        mui_uniface_multi_3fx[i] = created_unifaces[i].release();
    }
}

// Set of 3D interfaces with float=double and int=int32
void mui_create_uniface_multi_3d_f(const char *domain, const char **interfaces, int interface_count ) {

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
                std::cerr << "MUI Error [mui_f_wrapper_3d.cpp]: Error MUI Fortran C binding: NULL interfaces at the " << i << "th interface." << std::endl;
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
        std::cerr << "MUI Error [mui_f_wrapper_3d.cpp]: Error MUI Fortran C binding: The size of interfaces array must larger or equals to one." << std::endl;
        std::abort();
    }

    auto created_unifaces = mui::create_uniface<mui::config_3d>(domain, interface_names);

    mui_uniface_multi_3d = new mui_uniface_3d*[created_unifaces.size()];

      for(size_t i=0; i<created_unifaces.size(); i++) {
        mui_uniface_multi_3d[i] = created_unifaces[i].release();
    }
}

// Set of 3D interfaces with float=double and int=int64
void mui_create_uniface_multi_3dx_f(const char *domain, const char **interfaces, int interface_count ) {

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
                std::cerr << "MUI Error [mui_f_wrapper_3d.cpp]: Error MUI Fortran C binding: NULL interfaces at the " << i << "th interface." << std::endl;
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
        std::cerr << "MUI Error [mui_f_wrapper_3d.cpp]: Error MUI Fortran C binding: The size of interfaces array must larger or equals to one." << std::endl;
        std::abort();
    }

    auto created_unifaces = mui::create_uniface<mui::config_3dx>(domain, interface_names);

    mui_uniface_multi_3dx = new mui_uniface_3dx*[created_unifaces.size()];

      for(size_t i=0; i<created_unifaces.size(); i++) {
        mui_uniface_multi_3dx[i] = created_unifaces[i].release();
    }
}

// Set of 3D interfaces using config from config_c_wrapper.h
void mui_create_uniface_multi_3t_f(const char *domain, const char **interfaces, int interface_count ) {

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
                std::cerr << "MUI Error [mui_f_wrapper_3d.cpp]: Error MUI Fortran C binding: NULL interfaces at the " << i << "th interface." << std::endl;
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
        std::cerr << "MUI Error [mui_f_wrapper_3d.cpp]: Error MUI Fortran C binding: The size of interfaces array must larger or equals to one." << std::endl;
        std::abort();
    }

    auto created_unifaces = mui::create_uniface<mui::mui_f_wrapper_3D>(domain, interface_names);

    mui_uniface_multi_3t = new mui_uniface_3t*[created_unifaces.size()];

      for(size_t i=0; i<created_unifaces.size(); i++) {
        mui_uniface_multi_3t[i] = created_unifaces[i].release();
    }
}

// Access to MUI set of 3D interfaces with float=single and int=int32
mui_uniface_3f* get_mui_uniface_multi_3f_f(int interface_count) {
    return mui_uniface_multi_3f[interface_count-1];
}

// Access to MUI set of 3D interfaces with float=single and int=int64
mui_uniface_3fx* get_mui_uniface_multi_3fx_f(int interface_count) {
    return mui_uniface_multi_3fx[interface_count-1];
}

// Access to MUI set of 3D interfaces with float=double and int=int32
mui_uniface_3d* get_mui_uniface_multi_3d_f(int interface_count) {
    return mui_uniface_multi_3d[interface_count-1];
}

// Access to MUI set of 3D interfaces with float=double and int=int64
mui_uniface_3dx* get_mui_uniface_multi_3dx_f(int interface_count) {
    return mui_uniface_multi_3dx[interface_count-1];
}

// Access to MUI set of 3D interfaces using config from config_f_wrapper.h
mui_uniface_3t* get_mui_uniface_multi_3t_f(int interface_count) {
    return mui_uniface_multi_3t[interface_count-1];
}

/****************************************
 * Destroy MUI interface                 *
 ****************************************/

void mui_destroy_uniface_3f_f(mui_uniface_3f *uniface) {
    delete uniface;
}

void mui_destroy_uniface_3fx_f(mui_uniface_3fx *uniface) {
    delete uniface;
}

void mui_destroy_uniface_3d_f(mui_uniface_3d *uniface) {
    delete uniface;
}

void mui_destroy_uniface_3dx_f(mui_uniface_3dx *uniface) {
    delete uniface;
}

void mui_destroy_uniface_3t_f(mui_uniface_3t *uniface) {
    delete uniface;
}

/******************************************
 * Create 3d spatial samplers              *
 ******************************************/

// Exact sampler
void mui_create_sampler_exact_3f_f(mui_sampler_exact_3f **ret, float *tolerance) {
    *ret = new mui_sampler_exact_3f(*tolerance);
}

void mui_create_sampler_exact_3fx_f(mui_sampler_exact_3fx **ret, float *tolerance) {
    *ret = new mui_sampler_exact_3fx(*tolerance);
}

void mui_create_sampler_exact_3d_f(mui_sampler_exact_3d **ret, double *tolerance) {
    *ret = new mui_sampler_exact_3d(*tolerance);
}

void mui_create_sampler_exact_3dx_f(mui_sampler_exact_3dx** ret, double *tolerance) {
    *ret = new mui_sampler_exact_3dx(*tolerance);
}

void mui_create_sampler_exact_3t_f(mui_sampler_exact_3t** ret, double* tolerance) {
    *ret = new mui_sampler_exact_3t(static_cast<mui::mui_f_wrapper_3D::REAL>(*tolerance));
}

// Gauss sampler
void mui_create_sampler_gauss_3f_f(mui_sampler_gauss_3f **ret, float* r, float* h) {
    *ret = new mui_sampler_gauss_3f(*r, *h);
}

void mui_create_sampler_gauss_3fx_f(mui_sampler_gauss_3fx **ret, float* r, float* h) {
    *ret = new mui_sampler_gauss_3fx(*r, *h);
}

void mui_create_sampler_gauss_3d_f(mui_sampler_gauss_3d **ret, double* r, double* h) {
    *ret = new mui_sampler_gauss_3d(*r, *h);
}

void mui_create_sampler_gauss_3dx_f(mui_sampler_gauss_3dx** ret, double* r, double* h) {
    *ret = new mui_sampler_gauss_3dx(*r, *h);
}

void mui_create_sampler_gauss_3t_f(mui_sampler_gauss_3t** ret, double* r, double* h) {
    *ret = new mui_sampler_gauss_3t(static_cast<mui::mui_f_wrapper_3D::REAL>(*r),
            static_cast<mui::mui_f_wrapper_3D::REAL>(*h));
}

// Moving average sampler
void mui_create_sampler_moving_average_3f_f(mui_sampler_moving_average_3f **ret, float* bbox_1, float* bbox_2, float* bbox_3) {
    mui::point3f bbox(*bbox_1, *bbox_2, *bbox_3);
    *ret = new mui_sampler_moving_average_3f(bbox);
}

void mui_create_sampler_moving_average_3fx_f(mui_sampler_moving_average_3fx **ret, float* bbox_1, float* bbox_2, float* bbox_3) {
    mui::point3fx bbox(*bbox_1, *bbox_2, *bbox_3);
    *ret = new mui_sampler_moving_average_3fx(bbox);
}

void mui_create_sampler_moving_average_3d_f(mui_sampler_moving_average_3d **ret, double* bbox_1, double* bbox_2, double* bbox_3) {
    mui::point3d bbox(*bbox_1, *bbox_2, *bbox_3);
    *ret = new mui_sampler_moving_average_3d(bbox);
}

void mui_create_sampler_moving_average_3dx_f(mui_sampler_moving_average_3dx** ret, double* bbox_1, double* bbox_2, double* bbox_3) {
    mui::point3dx bbox(*bbox_1, *bbox_2, *bbox_3);
    *ret = new mui_sampler_moving_average_3dx(bbox);
}

void mui_create_sampler_moving_average_3t_f(mui_sampler_moving_average_3t** ret, double* bbox_1, double* bbox_2, double* bbox_3) {
    mui::mui_f_wrapper_3D::point_type bbox(*bbox_1, *bbox_2, *bbox_3);
    *ret = new mui_sampler_moving_average_3t(bbox);
}

// Nearest neighbour sampler
void mui_create_sampler_nearest_neighbor_3f_f(mui_sampler_nearest_neighbor_3f **ret) {
    *ret = new mui_sampler_nearest_neighbor_3f();
}

void mui_create_sampler_nearest_neighbor_3fx_f(mui_sampler_nearest_neighbor_3fx **ret) {
    *ret = new mui_sampler_nearest_neighbor_3fx();
}

void mui_create_sampler_nearest_neighbor_3d_f(mui_sampler_nearest_neighbor_3d **ret) {
    *ret = new mui_sampler_nearest_neighbor_3d();
}

void mui_create_sampler_nearest_neighbor_3dx_f(mui_sampler_nearest_neighbor_3dx** ret) {
    *ret = new mui_sampler_nearest_neighbor_3dx();
}

void mui_create_sampler_nearest_neighbor_3t_f(mui_sampler_nearest_neighbor_3t** ret) {
    *ret = new mui_sampler_nearest_neighbor_3t();
}

// Pseudo-linear n^2 interpolation sampler
void mui_create_sampler_pseudo_n2_linear_3f_f(mui_sampler_pseudo_n2_linear_3f **ret, float* r) {
    *ret = new mui_sampler_pseudo_n2_linear_3f(*r);
}

void mui_create_sampler_pseudo_n2_linear_3fx_f(mui_sampler_pseudo_n2_linear_3fx **ret, float* r) {
    *ret = new mui_sampler_pseudo_n2_linear_3fx(*r);
}

void mui_create_sampler_pseudo_n2_linear_3d_f(mui_sampler_pseudo_n2_linear_3d **ret, double* r) {
    *ret = new mui_sampler_pseudo_n2_linear_3d(*r);
}

void mui_create_sampler_pseudo_n2_linear_3dx_f(mui_sampler_pseudo_n2_linear_3dx** ret, double* r) {
    *ret = new mui_sampler_pseudo_n2_linear_3dx(*r);
}

void mui_create_sampler_pseudo_n2_linear_3t_f(mui_sampler_pseudo_n2_linear_3t** ret, double* r) {
    *ret = new mui_sampler_pseudo_n2_linear_3t(static_cast<mui::mui_f_wrapper_3D::REAL>(*r));
}

// Pseudo-linear nearest neighbour interpolation sampler
void mui_create_sampler_pseudo_nearest_neighbor_3f_f(mui_sampler_pseudo_nearest_neighbor_3f **ret, float* h) {
    *ret = new mui_sampler_pseudo_nearest_neighbor_3f(*h);
}

void mui_create_sampler_pseudo_nearest_neighbor_3fx_f(mui_sampler_pseudo_nearest_neighbor_3fx **ret, float* h) {
    *ret = new mui_sampler_pseudo_nearest_neighbor_3fx(*h);
}

void mui_create_sampler_pseudo_nearest_neighbor_3d_f(mui_sampler_pseudo_nearest_neighbor_3d **ret, double* h) {
    *ret = new mui_sampler_pseudo_nearest_neighbor_3d(*h);
}

void mui_create_sampler_pseudo_nearest_neighbor_3dx_f(mui_sampler_pseudo_nearest_neighbor_3dx** ret, double* h) {
    *ret = new mui_sampler_pseudo_nearest_neighbor_3dx(*h);
}

void mui_create_sampler_pseudo_nearest_neighbor_3t_f(mui_sampler_pseudo_nearest_neighbor_3t** ret, double* h) {
    *ret = new mui_sampler_pseudo_nearest_neighbor_3t(static_cast<mui::mui_f_wrapper_3D::REAL>(*h));
}

// Shepard interpolation with a quintic kernel sampler
void mui_create_sampler_shepard_quintic_3f_f(mui_sampler_shepard_quintic_3f **ret, float* r) {
    *ret = new mui_sampler_shepard_quintic_3f(*r);
}

void mui_create_sampler_shepard_quintic_3fx_f(mui_sampler_shepard_quintic_3fx **ret, float* r) {
    *ret = new mui_sampler_shepard_quintic_3fx(*r);
}

void mui_create_sampler_shepard_quintic_3d_f(mui_sampler_shepard_quintic_3d **ret, double* r) {
    *ret = new mui_sampler_shepard_quintic_3d(*r);
}

void mui_create_sampler_shepard_quintic_3dx_f(mui_sampler_shepard_quintic_3dx** ret, double* r) {
    *ret = new mui_sampler_shepard_quintic_3dx(*r);
}

void mui_create_sampler_shepard_quintic_3t_f(mui_sampler_shepard_quintic_3t** ret, double* r) {
    *ret = new mui_sampler_shepard_quintic_3t(static_cast<mui::mui_f_wrapper_3D::REAL>(*r));
}

// SPH derived interpolation method with a quintic spline kernel sampler
void mui_create_sampler_sph_quintic_3f_f(mui_sampler_sph_quintic_3f **ret, float* r) {
    *ret = new mui_sampler_sph_quintic_3f(*r);
}

void mui_create_sampler_sph_quintic_3fx_f(mui_sampler_sph_quintic_3fx **ret, float* r) {
    *ret = new mui_sampler_sph_quintic_3fx(*r);
}

void mui_create_sampler_sph_quintic_3d_f(mui_sampler_sph_quintic_3d **ret, double* r) {
    *ret = new mui_sampler_sph_quintic_3d(*r);
}

void mui_create_sampler_sph_quintic_3dx_f(mui_sampler_sph_quintic_3dx** ret, double* r) {
    *ret = new mui_sampler_sph_quintic_3dx(*r);
}

void mui_create_sampler_sph_quintic_3t_f(mui_sampler_sph_quintic_3t** ret, double* r) {
    *ret = new mui_sampler_sph_quintic_3t(static_cast<mui::mui_f_wrapper_3D::REAL>(*r));
}

// Summation with a quintic kernel sampler
void mui_create_sampler_sum_quintic_3f_f(mui_sampler_sum_quintic_3f **ret, float* r) {
    *ret = new mui_sampler_sum_quintic_3f(*r);
}

void mui_create_sampler_sum_quintic_3fx_f(mui_sampler_sum_quintic_3fx **ret, float* r) {
    *ret = new mui_sampler_sum_quintic_3fx(*r);
}

void mui_create_sampler_sum_quintic_3d_f(mui_sampler_sum_quintic_3d **ret, double* r) {
    *ret = new mui_sampler_sum_quintic_3d(*r);
}

void mui_create_sampler_sum_quintic_3dx_f(mui_sampler_sum_quintic_3dx** ret, double* r) {
    *ret = new mui_sampler_sum_quintic_3dx(*r);
}

void mui_create_sampler_sum_quintic_3t_f(mui_sampler_sum_quintic_3t** ret, double* r) {
    *ret = new mui_sampler_sum_quintic_3t(static_cast<mui::mui_f_wrapper_3D::REAL>(*r));
}

#ifdef USE_RBF
// Radial Basis Function sampler
void mui_create_sampler_rbf_3f_f(mui_sampler_rbf_3f **ret, float* r, float* points_1, float* points_2, float* points_3, int* points_count, int* basis_func,
        int* conservative, int* polynomial, int* smoothFunc, int* readMatrix, int* writeMatrix, const char* file_address, float* cutoff, float* cg_solve_tol,
        int* cg_solve_it, int* pou_size) {
    std::vector<mui::point3f> pts(*points_count);
    for (size_t i = 0; i < *points_count; i++) {
        pts[i][0] = points_1[i];
        pts[i][1] = points_2[i];
        pts[i][2] = points_3[i];
    }

    *ret = new mui_sampler_rbf_3f(*r, pts, *basis_func, static_cast<bool>(*conservative), static_cast<bool>(*polynomial),
      static_cast<bool>(*smoothFunc), static_cast<bool>(*readMatrix), static_cast<bool>(*writeMatrix), std::string(file_address),
      *cutoff, *cg_solve_tol, *cg_solve_it, *pou_size);
}

void mui_create_sampler_rbf_3fx_f(mui_sampler_rbf_3fx **ret, float* r, float* points_1, float* points_2, float* points_3, int* points_count, int* basis_func,
        int* conservative, int* polynomial, int* smoothFunc, int* readMatrix, int* writeMatrix, const char* file_address, float* cutoff, float* cg_solve_tol,
    int* cg_solve_it, int* pou_size) {
    std::vector<mui::point3fx> pts(*points_count);
    for (size_t i = 0; i < *points_count; i++) {
        pts[i][0] = points_1[i];
        pts[i][1] = points_2[i];
        pts[i][2] = points_3[i];
    }

    *ret = new mui_sampler_rbf_3fx(*r, pts, *basis_func, static_cast<bool>(*conservative), static_cast<bool>(*polynomial),
      static_cast<bool>(*smoothFunc), static_cast<bool>(*readMatrix), static_cast<bool>(*writeMatrix), std::string(file_address),
      *cutoff, *cg_solve_tol, *cg_solve_it, *pou_size);
}

void mui_create_sampler_rbf_3d_f(mui_sampler_rbf_3d **ret, double* r, double* points_1, double* points_2, double* points_3, int* points_count, int* basis_func,
        int* conservative, int* polynomial, int* smoothFunc, int* readMatrix, int* writeMatrix, const char* file_address, double* cutoff, double* cg_solve_tol,
    int* cg_solve_it, int* pou_size) {
    std::vector<mui::point3d> pts(*points_count);
    for (size_t i = 0; i < *points_count; i++) {
        pts[i][0] = points_1[i];
        pts[i][1] = points_2[i];
        pts[i][2] = points_3[i];
    }

    *ret = new mui_sampler_rbf_3d(*r, pts, *basis_func, static_cast<bool>(*conservative), static_cast<bool>(*polynomial),
      static_cast<bool>(*smoothFunc), static_cast<bool>(*readMatrix), static_cast<bool>(*writeMatrix), std::string(file_address),
      *cutoff, *cg_solve_tol, *cg_solve_it, *pou_size);
}

void mui_create_sampler_rbf_3dx_f(mui_sampler_rbf_3dx** ret, double* r, double* points_1, double* points_2, double* points_3, int* points_count, int* basis_func,
        int* conservative, int* polynomial, int* smoothFunc, int* readMatrix, int* writeMatrix, const char* file_address, double* cutoff, double* cg_solve_tol,
    int* cg_solve_it, int* pou_size) {
    std::vector<mui::point3dx> pts(*points_count);
    for (size_t i = 0; i < *points_count; i++) {
        pts[i][0] = points_1[i];
        pts[i][1] = points_2[i];
        pts[i][2] = points_3[i];
    }

    *ret = new mui_sampler_rbf_3dx(*r, pts, *basis_func, static_cast<bool>(*conservative), static_cast<bool>(*polynomial),
      static_cast<bool>(*smoothFunc), static_cast<bool>(*readMatrix), static_cast<bool>(*writeMatrix), std::string(file_address),
      *cutoff, *cg_solve_tol, *cg_solve_it, *pou_size);
}

void mui_create_sampler_rbf_3t_f(mui_sampler_rbf_3t** ret, double* r, double* points_1, double* points_2, double* points_3, int* points_count, int* basis_func,
        int* conservative, int* polynomial, int* smoothFunc, int* readMatrix, int* writeMatrix, const char* file_address, double* cutoff, double* cg_solve_tol,
    int* cg_solve_it, int* pou_size) {
    std::vector<mui::mui_f_wrapper_3D::point_type> pts(*points_count);
    for (size_t i = 0; i < *points_count; i++) {
        pts[i][0] = static_cast<mui::mui_f_wrapper_3D::REAL>(points_1[i]);
        pts[i][1] = static_cast<mui::mui_f_wrapper_3D::REAL>(points_2[i]);
        pts[i][2] = static_cast<mui::mui_f_wrapper_3D::REAL>(points_3[i]);
    }

    *ret = new mui_sampler_rbf_3t(static_cast<mui::mui_f_wrapper_3D::REAL>(*r), pts, *basis_func,
      static_cast<bool>(*conservative), static_cast<bool>(*polynomial), static_cast<bool>(*smoothFunc),
      static_cast<bool>(*readMatrix), static_cast<bool>(*writeMatrix), std::string(file_address),
      static_cast<mui::mui_f_wrapper_3D::REAL>(*cutoff), static_cast<mui::mui_f_wrapper_3D::REAL>(*cg_solve_tol),
      static_cast<mui::mui_f_wrapper_3D::INT>(*cg_solve_it), static_cast<mui::mui_f_wrapper_3D::INT>(*pou_size));
}
#endif

/*******************************************
 * Destroy 3d spatial samplers              *
 *******************************************/

// Exact sampler
void mui_destroy_sampler_exact_3f_f(mui_sampler_exact_3f* sampler) {
    delete sampler;
}

void mui_destroy_sampler_exact_3fx_f(mui_sampler_exact_3fx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_exact_3d_f(mui_sampler_exact_3d* sampler) {
    delete sampler;
}

void mui_destroy_sampler_exact_3dx_f(mui_sampler_exact_3dx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_exact_3t_f(mui_sampler_exact_3t* sampler) {
    delete sampler;
}

// Gauss sampler
void mui_destroy_sampler_gauss_3f_f(mui_sampler_gauss_3f* sampler) {
    delete sampler;
}

void mui_destroy_sampler_gauss_3fx_f(mui_sampler_gauss_3fx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_gauss_3d_f(mui_sampler_gauss_3d* sampler) {
    delete sampler;
}

void mui_destroy_sampler_gauss_3dx_f(mui_sampler_gauss_3dx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_gauss_3t_f(mui_sampler_gauss_3t* sampler) {
    delete sampler;
}

// Moving average sampler
void mui_destroy_sampler_moving_average_3f_f(mui_sampler_moving_average_3f* sampler) {
    delete sampler;
}

void mui_destroy_sampler_moving_average_3fx_f(mui_sampler_moving_average_3fx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_moving_average_3d_f(mui_sampler_moving_average_3d* sampler) {
    delete sampler;
}

void mui_destroy_sampler_moving_average_3dx_f(mui_sampler_moving_average_3dx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_moving_average_3t_f(mui_sampler_moving_average_3t* sampler) {
    delete sampler;
}

// Nearest neighbour sampler
void mui_destroy_sampler_nearest_neighbor_3f_f(mui_sampler_nearest_neighbor_3f* sampler) {
    delete sampler;
}

void mui_destroy_sampler_nearest_neighbor_3fx_f(mui_sampler_nearest_neighbor_3fx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_nearest_neighbor_3d_f(mui_sampler_nearest_neighbor_3d* sampler) {
    delete sampler;
}

void mui_destroy_sampler_nearest_neighbor_3dx_f(mui_sampler_nearest_neighbor_3dx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_nearest_neighbor_3t_f(mui_sampler_nearest_neighbor_3t* sampler) {
    delete sampler;
}

// Pseudo-linear n^2 interpolation sampler
void mui_destroy_sampler_pseudo_nearest2_linear_3f_f(mui_sampler_pseudo_nearest_neighbor_3f* sampler) {
    delete sampler;
}

void mui_destroy_sampler_pseudo_nearest2_linear_3fx_f(mui_sampler_pseudo_nearest_neighbor_3fx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_pseudo_nearest2_linear_3d_f(mui_sampler_pseudo_nearest_neighbor_3d* sampler) {
    delete sampler;
}

void mui_destroy_sampler_pseudo_nearest2_linear_3dx_f(mui_sampler_pseudo_nearest_neighbor_3dx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_pseudo_nearest2_linear_3t_f(mui_sampler_pseudo_nearest_neighbor_3t* sampler) {
    delete sampler;
}

// Pseudo-linear nearest neighbour interpolation sampler
void mui_destroy_sampler_pseudo_nearest_neighbor_3f_f(mui_sampler_pseudo_nearest_neighbor_3f* sampler) {
    delete sampler;
}

void mui_destroy_sampler_pseudo_nearest_neighbor_3fx_f(mui_sampler_pseudo_nearest_neighbor_3fx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_pseudo_nearest_neighbor_3d_f(mui_sampler_pseudo_nearest_neighbor_3d* sampler) {
    delete sampler;
}

void mui_destroy_sampler_pseudo_nearest_neighbor_3dx_f(mui_sampler_pseudo_nearest_neighbor_3dx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_pseudo_nearest_neighbor_3t_f(mui_sampler_pseudo_nearest_neighbor_3t* sampler) {
    delete sampler;
}

// Shepard interpolation with a quintic kernel sampler
void mui_destroy_sampler_shepard_quintic_3f_f(mui_sampler_shepard_quintic_3f* sampler) {
    delete sampler;
}

void mui_destroy_sampler_shepard_quintic_3fx_f(mui_sampler_shepard_quintic_3fx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_shepard_quintic_3d_f(mui_sampler_shepard_quintic_3d* sampler) {
    delete sampler;
}

void mui_destroy_sampler_shepard_quintic_3dx_f(mui_sampler_shepard_quintic_3dx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_shepard_quintic_3t_f(mui_sampler_shepard_quintic_3t* sampler) {
    delete sampler;
}

// SPH derived interpolation method with a quintic spline kernel sampler
void mui_destroy_sampler_sph_quintic_3f_f(mui_sampler_sph_quintic_3f* sampler) {
    delete sampler;
}

void mui_destroy_sampler_sph_quintic_3fx_f(mui_sampler_sph_quintic_3fx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_sph_quintic_3d_f(mui_sampler_sph_quintic_3d* sampler) {
    delete sampler;
}

void mui_destroy_sampler_sph_quintic_3dx_f(mui_sampler_sph_quintic_3dx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_sph_quintic_3t_f(mui_sampler_sph_quintic_3t* sampler) {
    delete sampler;
}

// Summation with a quintic kernel sampler
void mui_destroy_sampler_sum_quintic_3f_f(mui_sampler_sum_quintic_3f* sampler) {
    delete sampler;
}

void mui_destroy_sampler_sum_quintic_3fx_f(mui_sampler_sum_quintic_3fx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_sum_quintic_3d_f(mui_sampler_sum_quintic_3d* sampler) {
    delete sampler;
}

void mui_destroy_sampler_sum_quintic_3dx_f(mui_sampler_sum_quintic_3dx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_sum_quintic_3t_f(mui_sampler_sum_quintic_3t* sampler) {
    delete sampler;
}

#ifdef USE_RBF
void mui_destroy_sampler_rbf_3f_f(mui_sampler_rbf_3f* sampler) {
    delete sampler;
}

void mui_destroy_sampler_rbf_3fx_f(mui_sampler_rbf_3fx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_rbf_3d_f(mui_sampler_rbf_3d* sampler) {
    delete sampler;
}

void mui_destroy_sampler_rbf_3dx_f(mui_sampler_rbf_3dx* sampler) {
    delete sampler;
}

void mui_destroy_sampler_rbf_3t_f(mui_sampler_rbf_3t* sampler) {
    delete sampler;
}
#endif

/*******************************************
 * Create temporal samplers                 *
 *******************************************/

// Exact temporal sampler
void mui_create_chrono_sampler_exact_3f_f(mui_chrono_sampler_exact_3f **ret, float* tolerance) {
    *ret = new mui_chrono_sampler_exact_3f(*tolerance);
}

void mui_create_chrono_sampler_exact_3fx_f(mui_chrono_sampler_exact_3fx **ret, float* tolerance) {
    *ret = new mui_chrono_sampler_exact_3fx(*tolerance);
}

void mui_create_chrono_sampler_exact_3d_f(mui_chrono_sampler_exact_3d **ret, double* tolerance) {
    *ret = new mui_chrono_sampler_exact_3d(*tolerance);
}

void mui_create_chrono_sampler_exact_3dx_f(mui_chrono_sampler_exact_3dx** ret, double* tolerance) {
    *ret = new mui_chrono_sampler_exact_3dx(*tolerance);
}

void mui_create_chrono_sampler_exact_3t_f(mui_chrono_sampler_exact_3t** ret, double* tolerance) {
    *ret = new mui_chrono_sampler_exact_3t(static_cast<mui::mui_f_wrapper_3D::REAL>(*tolerance));
}

// Gauss temporal sampler
void mui_create_chrono_sampler_gauss_3f_f(mui_chrono_sampler_gauss_3f **ret, float* cutoff, float* sigma) {
    *ret = new mui_chrono_sampler_gauss_3f(*cutoff, *sigma);
}

void mui_create_chrono_sampler_gauss_3fx_f(mui_chrono_sampler_gauss_3fx **ret, float* cutoff, float* sigma) {
    *ret = new mui_chrono_sampler_gauss_3fx(*cutoff, *sigma);
}

void mui_create_chrono_sampler_gauss_3d_f(mui_chrono_sampler_gauss_3d **ret, double* cutoff, double* sigma) {
    *ret = new mui_chrono_sampler_gauss_3d(*cutoff, *sigma);
}

void mui_create_chrono_sampler_gauss_3dx_f(mui_chrono_sampler_gauss_3dx** ret, double* cutoff, double* sigma) {
    *ret = new mui_chrono_sampler_gauss_3dx(*cutoff, *sigma);
}

void mui_create_chrono_sampler_gauss_3t_f(mui_chrono_sampler_gauss_3t** ret, double* cutoff, double* sigma) {
    *ret = new mui_chrono_sampler_gauss_3t(static_cast<mui::mui_f_wrapper_3D::time_type>(*cutoff),
            static_cast<mui::mui_f_wrapper_3D::REAL>(*sigma));
}

// Mean temporal sampler
void mui_create_chrono_sampler_mean_3f_f(mui_chrono_sampler_mean_3f **ret, float* lower, float* upper) {
    *ret = new mui_chrono_sampler_mean_3f(*lower, *upper);
}

void mui_create_chrono_sampler_mean_3fx_f(mui_chrono_sampler_mean_3fx **ret, float* lower, float* upper) {
    *ret = new mui_chrono_sampler_mean_3fx(*lower, *upper);
}

void mui_create_chrono_sampler_mean_3d_f(mui_chrono_sampler_mean_3d **ret, double* lower, double* upper) {
    *ret = new mui_chrono_sampler_mean_3d(*lower, *upper);
}

void mui_create_chrono_sampler_mean_3dx_f(mui_chrono_sampler_mean_3dx** ret, double* lower, double* upper) {
    *ret = new mui_chrono_sampler_mean_3dx(*lower, *upper);
}

void mui_create_chrono_sampler_mean_3t_f(mui_chrono_sampler_mean_3t** ret, double* lower, double* upper) {
    *ret = new mui_chrono_sampler_mean_3t(static_cast<mui::mui_f_wrapper_3D::time_type>(*lower),
            static_cast<mui::mui_f_wrapper_3D::time_type>(*upper));
}

// Sum temporal sampler
void mui_create_chrono_sampler_sum_3f_f(mui_chrono_sampler_sum_3f **ret, float* lower, float* upper) {
    *ret = new mui_chrono_sampler_sum_3f(*lower, *upper);
}

void mui_create_chrono_sampler_sum_3fx_f(mui_chrono_sampler_sum_3fx **ret, float* lower, float* upper) {
    *ret = new mui_chrono_sampler_sum_3fx(*lower, *upper);
}

void mui_create_chrono_sampler_sum_3d_f(mui_chrono_sampler_sum_3d **ret, double* lower, double* upper) {
    *ret = new mui_chrono_sampler_sum_3d(*lower, *upper);
}

void mui_create_chrono_sampler_sum_3dx_f(mui_chrono_sampler_sum_3dx** ret, double* lower, double* upper) {
    *ret = new mui_chrono_sampler_sum_3dx(*lower, *upper);
}

void mui_create_chrono_sampler_sum_3t_f(mui_chrono_sampler_sum_3t** ret, double* lower, double* upper) {
    *ret = new mui_chrono_sampler_sum_3t(static_cast<mui::mui_f_wrapper_3D::time_type>(*lower),
            static_cast<mui::mui_f_wrapper_3D::time_type>(*upper));
}

/*******************************************
 * Destroy temporal samplers                *
 *******************************************/

// Exact temporal sampler
void mui_destroy_chrono_sampler_exact_3f_f(mui_chrono_sampler_exact_3f* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_exact_3fx_f(mui_chrono_sampler_exact_3fx* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_exact_3d_f(mui_chrono_sampler_exact_3d* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_exact_3dx_f(mui_chrono_sampler_exact_3dx* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_exact_3t_f(mui_chrono_sampler_exact_3t* sampler) {
    delete sampler;
}

// Gauss temporal sampler
void mui_destroy_chrono_sampler_gauss_3f_f(mui_chrono_sampler_gauss_3f* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_gauss_3fx_f(mui_chrono_sampler_gauss_3fx* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_gauss_3d_f(mui_chrono_sampler_gauss_3d* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_gauss_3dx_f(mui_chrono_sampler_gauss_3dx* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_gauss_3t_f(mui_chrono_sampler_gauss_3t* sampler) {
    delete sampler;
}

// Mean temporal sampler
void mui_destroy_chrono_sampler_mean_3f_f(mui_chrono_sampler_mean_3f* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_mean_3fx_f(mui_chrono_sampler_mean_3fx* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_mean_3d_f(mui_chrono_sampler_mean_3d* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_mean_3dx_f(mui_chrono_sampler_mean_3dx* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_mean_3t_f(mui_chrono_sampler_mean_3t* sampler) {
    delete sampler;
}

// Sum temporal sampler
void mui_destroy_chrono_sampler_sum_3f_f(mui_chrono_sampler_sum_3f* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_sum_3fx_f(mui_chrono_sampler_sum_3fx* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_sum_3d_f(mui_chrono_sampler_sum_3d* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_sum_3dx_f(mui_chrono_sampler_sum_3dx* sampler) {
    delete sampler;
}

void mui_destroy_chrono_sampler_sum_3t_f(mui_chrono_sampler_sum_3t* sampler) {
    delete sampler;
}

/******************************************
 * MUI functions for data push             *
 ******************************************/

// Standard push functions
void mui_push_3f_f(mui_uniface_3f *uniface, const char* attr, float* point_1, float* point_2, float* point_3, float* value) {
    uniface->push(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *value);
}

void mui_push_3fx_f(mui_uniface_3fx *uniface, const char* attr, float* point_1, float* point_2, float* point_3, float* value) {
    uniface->push(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *value);
}

void mui_push_3d_f(mui_uniface_3d *uniface, const char* attr, double* point_1, double* point_2, double* point_3, double* value) {
    uniface->push(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *value);
}

void mui_push_3dx_f(mui_uniface_3dx *uniface, const char* attr, double* point_1, double* point_2, double* point_3, double* value) {
    uniface->push(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *value);
}

void mui_push_3t_f(mui_uniface_3t *uniface, const char* attr, double* point_1, double* point_2, double* point_3, double* value) {
    mui::mui_f_wrapper_3D::point_type push_point(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    uniface->push(std::string(attr), push_point, static_cast<mui::mui_f_wrapper_3D::REAL>(*value));
}

// Single parameter push functions
void mui_push_3f_param_f(mui_uniface_3f *uniface, const char* attr, float* value) {
    uniface->push(std::string(attr), *value);
}

void mui_push_3fx_param_f(mui_uniface_3fx *uniface, const char* attr, float* value) {
    uniface->push(std::string(attr), *value);
}

void mui_push_3d_param_f(mui_uniface_3d *uniface, const char* attr, double* value) {
    uniface->push(std::string(attr), *value);
}

void mui_push_3dx_param_f(mui_uniface_3dx *uniface, const char* attr, double* value) {
    uniface->push(std::string(attr), *value);
}

void mui_push_3t_param_f(mui_uniface_3t *uniface, const char* attr, double* value) {
    uniface->push(std::string(attr), static_cast<mui::mui_f_wrapper_3D::REAL>(*value));
}

/******************************************
 * MUI functions for data commit           *
 ******************************************/

// Commit using one time value
void mui_commit_3f_f(mui_uniface_3f *uniface, float* t) {
    uniface->commit(*t);
}

void mui_commit_3fx_f(mui_uniface_3fx *uniface, float* t) {
    uniface->commit(*t);
}

void mui_commit_3d_f(mui_uniface_3d *uniface, double* t) {
    uniface->commit(*t);
}

void mui_commit_3dx_f(mui_uniface_3dx *uniface, double* t) {
    uniface->commit(*t);
}

void mui_commit_3t_f(mui_uniface_3t *uniface, double* t) {
    uniface->commit(static_cast<mui::mui_f_wrapper_3D::time_type>(*t));
}

// Commit using two time values
void mui_commit_3f_pair_f(mui_uniface_3f *uniface, float* t_1, float* t_2) {
    uniface->commit(*t_1, *t_2);
}

void mui_commit_3fx_pair_f(mui_uniface_3fx *uniface, float* t_1, float* t_2) {
    uniface->commit(*t_1, *t_2);
}

void mui_commit_3d_pair_f(mui_uniface_3d *uniface, double* t_1, double* t_2) {
    uniface->commit(*t_1, *t_2);
}

void mui_commit_3dx_pair_f(mui_uniface_3dx *uniface, double* t_1, double* t_2) {
    uniface->commit(*t_1, *t_2);
}

void mui_commit_3t_pair_f(mui_uniface_3t *uniface, double* t_1, double* t_2) {
    uniface->commit(static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1),
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2));
}

/******************************************
 * MUI functions for data forecast         *
 ******************************************/

// Forecast using one time value
void mui_forecast_3f_f(mui_uniface_3f *uniface, float* t) {
    uniface->forecast(*t);
}

void mui_forecast_3fx_f(mui_uniface_3fx *uniface, float* t) {
    uniface->forecast(*t);
}

void mui_forecast_3d_f(mui_uniface_3d *uniface, double* t) {
    uniface->forecast(*t);
}

void mui_forecast_3dx_f(mui_uniface_3dx *uniface, double* t) {
    uniface->forecast(*t);
}

void mui_forecast_3t_f(mui_uniface_3t *uniface, double* t) {
    uniface->forecast(static_cast<mui::mui_f_wrapper_3D::time_type>(*t));
}

// Forecast using two time values
void mui_forecast_3f_pair_f(mui_uniface_3f *uniface, float* t_1, float* t_2) {
    uniface->forecast(*t_1, *t_2);
}

void mui_forecast_3fx_pair_f(mui_uniface_3fx *uniface, float* t_1, float* t_2) {
    uniface->forecast(*t_1, *t_2);
}

void mui_forecast_3d_pair_f(mui_uniface_3d *uniface, double* t_1, double* t_2) {
    uniface->forecast(*t_1, *t_2);
}

void mui_forecast_3dx_pair_f(mui_uniface_3dx *uniface, double* t_1, double* t_2) {
    uniface->forecast(*t_1, *t_2);
}

void mui_forecast_3t_pair_f(mui_uniface_3t *uniface, double* t_1, double* t_2) {
    uniface->forecast(static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1),
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2));
}

/*******************************************************
 * MUI functions for 3d data fetch using one time value *
 ********************************************************/

// Spatial sampler: exact; temporal sampler: exact
void mui_fetch_exact_exact_3f_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t,
        mui_sampler_exact_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_exact_3fx_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t,
        mui_sampler_exact_3fx *spatial_sampler, mui_chrono_sampler_exact_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_exact_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_exact_3d *spatial_sampler, mui_chrono_sampler_exact_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_exact_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_exact_3dx *spatial_sampler, mui_chrono_sampler_exact_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_exact_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_exact_3t *spatial_sampler, mui_chrono_sampler_exact_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: exact; temporal sampler: gauss
void mui_fetch_exact_gauss_3f_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t,
        mui_sampler_exact_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_gauss_3fx_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t,
        mui_sampler_exact_3fx *spatial_sampler, mui_chrono_sampler_gauss_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_gauss_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_exact_3d *spatial_sampler, mui_chrono_sampler_gauss_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_gauss_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_exact_3dx *spatial_sampler, mui_chrono_sampler_gauss_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_gauss_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_exact_3t *spatial_sampler, mui_chrono_sampler_gauss_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: exact; temporal sampler: mean
void mui_fetch_exact_mean_3f_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t,
        mui_sampler_exact_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_mean_3fx_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t,
        mui_sampler_exact_3fx *spatial_sampler, mui_chrono_sampler_mean_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_mean_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_exact_3d *spatial_sampler, mui_chrono_sampler_mean_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_mean_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_exact_3dx *spatial_sampler, mui_chrono_sampler_mean_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_mean_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_exact_3t *spatial_sampler, mui_chrono_sampler_mean_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: exact; temporal sampler: sum
void mui_fetch_exact_sum_3f_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t,
        mui_sampler_exact_3f *spatial_sampler, mui_chrono_sampler_sum_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_sum_3fx_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t,
        mui_sampler_exact_3fx *spatial_sampler, mui_chrono_sampler_sum_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_sum_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_exact_3d *spatial_sampler, mui_chrono_sampler_sum_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_sum_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_exact_3dx *spatial_sampler, mui_chrono_sampler_sum_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_sum_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_exact_3t *spatial_sampler, mui_chrono_sampler_sum_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: gauss; temporal sampler: exact
void mui_fetch_gauss_exact_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_gauss_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_exact_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_gauss_3fx *spatial_sampler, mui_chrono_sampler_exact_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_exact_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_gauss_3d *spatial_sampler, mui_chrono_sampler_exact_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_exact_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_gauss_3dx *spatial_sampler, mui_chrono_sampler_exact_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_exact_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_gauss_3t *spatial_sampler, mui_chrono_sampler_exact_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: gauss; temporal sampler: gauss
void mui_fetch_gauss_gauss_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_gauss_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_gauss_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_gauss_3fx *spatial_sampler, mui_chrono_sampler_gauss_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_gauss_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_gauss_3d *spatial_sampler, mui_chrono_sampler_gauss_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_gauss_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_gauss_3dx *spatial_sampler, mui_chrono_sampler_gauss_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_gauss_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_gauss_3t *spatial_sampler, mui_chrono_sampler_gauss_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: gauss; temporal sampler: mean
void mui_fetch_gauss_mean_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_gauss_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_mean_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_gauss_3fx *spatial_sampler, mui_chrono_sampler_mean_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_mean_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_gauss_3d *spatial_sampler, mui_chrono_sampler_mean_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_mean_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_gauss_3dx *spatial_sampler, mui_chrono_sampler_mean_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_mean_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_gauss_3t *spatial_sampler, mui_chrono_sampler_mean_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: gauss; temporal sampler: sum
void mui_fetch_gauss_sum_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_gauss_3f *spatial_sampler, mui_chrono_sampler_sum_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_sum_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_gauss_3fx *spatial_sampler, mui_chrono_sampler_sum_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_sum_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_gauss_3d *spatial_sampler, mui_chrono_sampler_sum_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_sum_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_gauss_3dx *spatial_sampler, mui_chrono_sampler_sum_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_sum_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_gauss_3t *spatial_sampler, mui_chrono_sampler_sum_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: moving average; temporal sampler: exact
void mui_fetch_moving_average_exact_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_moving_average_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_exact_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_moving_average_3fx *spatial_sampler, mui_chrono_sampler_exact_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_exact_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_moving_average_3d *spatial_sampler, mui_chrono_sampler_exact_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_exact_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_moving_average_3dx *spatial_sampler, mui_chrono_sampler_exact_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_exact_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_moving_average_3t *spatial_sampler, mui_chrono_sampler_exact_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: moving average; temporal sampler: gauss
void mui_fetch_moving_average_gauss_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_moving_average_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_gauss_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_moving_average_3fx *spatial_sampler, mui_chrono_sampler_gauss_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_gauss_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_moving_average_3d *spatial_sampler, mui_chrono_sampler_gauss_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_gauss_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_moving_average_3dx *spatial_sampler, mui_chrono_sampler_gauss_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_gauss_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_moving_average_3t *spatial_sampler, mui_chrono_sampler_gauss_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: moving average; temporal sampler: mean
void mui_fetch_moving_average_mean_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_moving_average_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_mean_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_moving_average_3fx *spatial_sampler, mui_chrono_sampler_mean_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_mean_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_moving_average_3d *spatial_sampler, mui_chrono_sampler_mean_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_mean_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_moving_average_3dx *spatial_sampler, mui_chrono_sampler_mean_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_mean_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_moving_average_3t *spatial_sampler, mui_chrono_sampler_mean_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: moving average; temporal sampler: sum
void mui_fetch_moving_average_sum_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_moving_average_3f *spatial_sampler, mui_chrono_sampler_sum_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_sum_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_moving_average_3fx *spatial_sampler, mui_chrono_sampler_sum_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_sum_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_moving_average_3d *spatial_sampler, mui_chrono_sampler_sum_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_sum_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_moving_average_3dx *spatial_sampler, mui_chrono_sampler_sum_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_sum_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_moving_average_3t *spatial_sampler, mui_chrono_sampler_sum_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: nearest neighbor; temporal sampler: exact
void mui_fetch_nearest_neighbor_exact_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_nearest_neighbor_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_exact_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_nearest_neighbor_3fx *spatial_sampler, mui_chrono_sampler_exact_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_exact_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_nearest_neighbor_3d *spatial_sampler, mui_chrono_sampler_exact_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_exact_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_nearest_neighbor_3dx *spatial_sampler, mui_chrono_sampler_exact_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_exact_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_nearest_neighbor_3t *spatial_sampler, mui_chrono_sampler_exact_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: nearest neighbor; temporal sampler: gauss
void mui_fetch_nearest_neighbor_gauss_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_nearest_neighbor_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_gauss_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_nearest_neighbor_3fx *spatial_sampler, mui_chrono_sampler_gauss_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_gauss_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_nearest_neighbor_3d *spatial_sampler, mui_chrono_sampler_gauss_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_gauss_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_nearest_neighbor_3dx *spatial_sampler, mui_chrono_sampler_gauss_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_gauss_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_nearest_neighbor_3t *spatial_sampler, mui_chrono_sampler_gauss_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: nearest neighbor; temporal sampler: mean
void mui_fetch_nearest_neighbor_mean_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_nearest_neighbor_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_mean_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_nearest_neighbor_3fx *spatial_sampler, mui_chrono_sampler_mean_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_mean_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_nearest_neighbor_3d *spatial_sampler, mui_chrono_sampler_mean_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_mean_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_nearest_neighbor_3dx *spatial_sampler, mui_chrono_sampler_mean_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_mean_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_nearest_neighbor_3t *spatial_sampler, mui_chrono_sampler_mean_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: nearest neighbor; temporal sampler: sum
void mui_fetch_nearest_neighbor_sum_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_nearest_neighbor_3f *spatial_sampler, mui_chrono_sampler_sum_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_sum_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_nearest_neighbor_3fx *spatial_sampler, mui_chrono_sampler_sum_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_sum_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_nearest_neighbor_3d *spatial_sampler, mui_chrono_sampler_sum_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_sum_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_nearest_neighbor_3dx *spatial_sampler, mui_chrono_sampler_sum_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_sum_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_nearest_neighbor_3t *spatial_sampler, mui_chrono_sampler_sum_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: pseudo nearest neighbor; temporal sampler: exact
void mui_fetch_pseudo_nearest_neighbor_exact_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_pseudo_nearest_neighbor_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_exact_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3,
        float* t, mui_sampler_pseudo_nearest_neighbor_3fx *spatial_sampler,
        mui_chrono_sampler_exact_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_exact_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t, mui_sampler_pseudo_nearest_neighbor_3d *spatial_sampler,
        mui_chrono_sampler_exact_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_exact_3dx(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t, mui_sampler_pseudo_nearest_neighbor_3dx *spatial_sampler,
        mui_chrono_sampler_exact_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_exact_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t, mui_sampler_pseudo_nearest_neighbor_3t *spatial_sampler,
        mui_chrono_sampler_exact_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: nearest neighbor; temporal sampler: gauss
void mui_fetch_pseudo_nearest_neighbor_gauss_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_pseudo_nearest_neighbor_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_gauss_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3,
        float* t, mui_sampler_pseudo_nearest_neighbor_3fx *spatial_sampler,
        mui_chrono_sampler_gauss_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_gauss_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t, mui_sampler_pseudo_nearest_neighbor_3d *spatial_sampler,
        mui_chrono_sampler_gauss_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_gauss_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t, mui_sampler_pseudo_nearest_neighbor_3dx *spatial_sampler,
        mui_chrono_sampler_gauss_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_gauss_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t, mui_sampler_pseudo_nearest_neighbor_3t *spatial_sampler,
        mui_chrono_sampler_gauss_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: nearest neighbor; temporal sampler: mean
void mui_fetch_pseudo_nearest_neighbor_mean_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_pseudo_nearest_neighbor_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_mean_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3,
        float* t, mui_sampler_pseudo_nearest_neighbor_3fx *spatial_sampler,
        mui_chrono_sampler_mean_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_mean_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t, mui_sampler_pseudo_nearest_neighbor_3d *spatial_sampler,
        mui_chrono_sampler_mean_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_mean_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t, mui_sampler_pseudo_nearest_neighbor_3dx *spatial_sampler,
        mui_chrono_sampler_mean_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_mean_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t, mui_sampler_pseudo_nearest_neighbor_3t *spatial_sampler,
        mui_chrono_sampler_mean_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: nearest neighbor; temporal sampler: sum
void mui_fetch_pseudo_nearest_neighbor_sum_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_pseudo_nearest_neighbor_3f *spatial_sampler, mui_chrono_sampler_sum_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_sum_3fx(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3,
        float* t, mui_sampler_pseudo_nearest_neighbor_3fx *spatial_sampler,
        mui_chrono_sampler_sum_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_sum_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_pseudo_nearest_neighbor_3d *spatial_sampler, mui_chrono_sampler_sum_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_sum_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t, mui_sampler_pseudo_nearest_neighbor_3dx *spatial_sampler,
        mui_chrono_sampler_sum_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_sum_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_pseudo_nearest_neighbor_3t *spatial_sampler, mui_chrono_sampler_sum_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: shepard quintic; temporal sampler: exact
void mui_fetch_shepard_quintic_exact_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_shepard_quintic_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_exact_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_shepard_quintic_3fx *spatial_sampler, mui_chrono_sampler_exact_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_exact_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_shepard_quintic_3d *spatial_sampler, mui_chrono_sampler_exact_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_exact_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_shepard_quintic_3dx *spatial_sampler, mui_chrono_sampler_exact_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_exact_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_shepard_quintic_3t *spatial_sampler, mui_chrono_sampler_exact_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: shepard quintic; temporal sampler: gauss
void mui_fetch_shepard_quintic_gauss_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_shepard_quintic_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_gauss_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_shepard_quintic_3fx *spatial_sampler, mui_chrono_sampler_gauss_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_gauss_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_shepard_quintic_3d *spatial_sampler, mui_chrono_sampler_gauss_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_gauss_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_shepard_quintic_3dx *spatial_sampler, mui_chrono_sampler_gauss_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_gauss_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_shepard_quintic_3t *spatial_sampler, mui_chrono_sampler_gauss_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: shepard quintic; temporal sampler: mean
void mui_fetch_shepard_quintic_mean_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_shepard_quintic_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_mean_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_shepard_quintic_3fx *spatial_sampler, mui_chrono_sampler_mean_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_mean_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_shepard_quintic_3d *spatial_sampler, mui_chrono_sampler_mean_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_mean_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_shepard_quintic_3dx *spatial_sampler, mui_chrono_sampler_mean_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_mean_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_shepard_quintic_3t *spatial_sampler, mui_chrono_sampler_mean_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: shepard quintic; temporal sampler: sum
void mui_fetch_shepard_quintic_sum_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_shepard_quintic_3f *spatial_sampler, mui_chrono_sampler_sum_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_sum_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_shepard_quintic_3fx *spatial_sampler, mui_chrono_sampler_sum_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_sum_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_shepard_quintic_3d *spatial_sampler, mui_chrono_sampler_sum_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_sum_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_shepard_quintic_3dx *spatial_sampler, mui_chrono_sampler_sum_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_sum_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_shepard_quintic_3t *spatial_sampler, mui_chrono_sampler_sum_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: sph-derived quintic; temporal sampler: exact
void mui_fetch_sph_quintic_exact_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_sph_quintic_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_exact_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_sph_quintic_3fx *spatial_sampler, mui_chrono_sampler_exact_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_exact_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_sph_quintic_3d *spatial_sampler, mui_chrono_sampler_exact_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_exact_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_sph_quintic_3dx *spatial_sampler, mui_chrono_sampler_exact_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_exact_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_sph_quintic_3t *spatial_sampler, mui_chrono_sampler_exact_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: sph-derived quintic; temporal sampler: gauss
void mui_fetch_sph_quintic_gauss_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_sph_quintic_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_gauss_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_sph_quintic_3fx *spatial_sampler, mui_chrono_sampler_gauss_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_gauss_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_sph_quintic_3d *spatial_sampler, mui_chrono_sampler_gauss_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_gauss_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_sph_quintic_3dx *spatial_sampler, mui_chrono_sampler_gauss_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_gauss_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_sph_quintic_3t *spatial_sampler, mui_chrono_sampler_gauss_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: sph-derived quintic; temporal sampler: mean
void mui_fetch_sph_quintic_mean_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_sph_quintic_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_mean_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_sph_quintic_3fx *spatial_sampler, mui_chrono_sampler_mean_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_mean_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_sph_quintic_3d *spatial_sampler, mui_chrono_sampler_mean_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_mean_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_sph_quintic_3dx *spatial_sampler, mui_chrono_sampler_mean_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_mean_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_sph_quintic_3t *spatial_sampler, mui_chrono_sampler_mean_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: sph-derived quintic; temporal sampler: sum
void mui_fetch_sph_quintic_sum_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_sph_quintic_3f *spatial_sampler, mui_chrono_sampler_sum_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_sum_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_sph_quintic_3fx *spatial_sampler, mui_chrono_sampler_sum_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_sum_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_sph_quintic_3d *spatial_sampler, mui_chrono_sampler_sum_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_sum_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_sph_quintic_3dx *spatial_sampler, mui_chrono_sampler_sum_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_sum_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_sph_quintic_3t *spatial_sampler, mui_chrono_sampler_sum_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: summation quintic; temporal sampler: exact
void mui_fetch_sum_quintic_exact_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_sum_quintic_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_exact_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_sum_quintic_3fx *spatial_sampler, mui_chrono_sampler_exact_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_exact_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_sum_quintic_3d *spatial_sampler, mui_chrono_sampler_exact_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_exact_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_sum_quintic_3dx *spatial_sampler, mui_chrono_sampler_exact_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_exact_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_sum_quintic_3t *spatial_sampler, mui_chrono_sampler_exact_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: summation quintic; temporal sampler: gauss
void mui_fetch_sum_quintic_gauss_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_sum_quintic_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_gauss_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_sum_quintic_3fx *spatial_sampler, mui_chrono_sampler_gauss_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_gauss_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_sum_quintic_3d *spatial_sampler, mui_chrono_sampler_gauss_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_gauss_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_sum_quintic_3dx *spatial_sampler, mui_chrono_sampler_gauss_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_gauss_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_sum_quintic_3t *spatial_sampler, mui_chrono_sampler_gauss_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: summation quintic; temporal sampler: mean
void mui_fetch_sum_quintic_mean_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_sum_quintic_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_mean_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_sum_quintic_3fx *spatial_sampler, mui_chrono_sampler_mean_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_mean_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_sum_quintic_3d *spatial_sampler, mui_chrono_sampler_mean_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_mean_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_sum_quintic_3dx *spatial_sampler, mui_chrono_sampler_mean_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_mean_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_sum_quintic_3t *spatial_sampler, mui_chrono_sampler_mean_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: summation quintic; temporal sampler: sum
void mui_fetch_sum_quintic_sum_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_sum_quintic_3f *spatial_sampler, mui_chrono_sampler_sum_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_sum_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_sum_quintic_3fx *spatial_sampler, mui_chrono_sampler_sum_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_sum_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_sum_quintic_3d *spatial_sampler, mui_chrono_sampler_sum_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_sum_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_sum_quintic_3dx *spatial_sampler, mui_chrono_sampler_sum_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_sum_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_sum_quintic_3t *spatial_sampler, mui_chrono_sampler_sum_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

#ifdef USE_RBF
// Spatial sampler: radial basis function; temporal sampler: exact
void mui_fetch_rbf_exact_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_rbf_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_exact_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_rbf_3fx *spatial_sampler, mui_chrono_sampler_exact_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_exact_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_rbf_3d *spatial_sampler, mui_chrono_sampler_exact_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_exact_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_rbf_3dx *spatial_sampler, mui_chrono_sampler_exact_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_exact_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_rbf_3t *spatial_sampler, mui_chrono_sampler_exact_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: radial basis function; temporal sampler: gauss
void mui_fetch_rbf_gauss_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_rbf_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_gauss_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_rbf_3fx *spatial_sampler, mui_chrono_sampler_gauss_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_gauss_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_rbf_3d *spatial_sampler, mui_chrono_sampler_gauss_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_gauss_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_rbf_3dx *spatial_sampler, mui_chrono_sampler_gauss_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_gauss_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_rbf_3t *spatial_sampler, mui_chrono_sampler_gauss_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: radial basis function; temporal sampler: mean
void mui_fetch_rbf_mean_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_rbf_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_mean_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_rbf_3fx *spatial_sampler, mui_chrono_sampler_mean_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_mean_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_rbf_3d *spatial_sampler, mui_chrono_sampler_mean_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_mean_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_rbf_3dx *spatial_sampler, mui_chrono_sampler_mean_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_mean_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_rbf_3t *spatial_sampler, mui_chrono_sampler_mean_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}

// Spatial sampler: radial basis function; temporal sampler: sum
void mui_fetch_rbf_sum_3f_f(mui_uniface_3f *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_rbf_3f *spatial_sampler, mui_chrono_sampler_sum_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_sum_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *point_1, float *point_2, float *point_3, float* t,
        mui_sampler_rbf_3fx *spatial_sampler, mui_chrono_sampler_sum_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_sum_3d_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_rbf_3d *spatial_sampler, mui_chrono_sampler_sum_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_sum_3dx_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_rbf_3dx *spatial_sampler, mui_chrono_sampler_sum_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_sum_3t_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t,
        mui_sampler_rbf_3t *spatial_sampler, mui_chrono_sampler_sum_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t), *spatial_sampler, *temporal_sampler));

}
#endif

/********************************************************
 * MUI functions for 3d data fetch using two time values *
 *********************************************************/

// Spatial sampler: exact; temporal sampler: exact
void mui_fetch_exact_exact_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1, float* t_2,
        mui_sampler_exact_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_exact_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_exact_3fx *spatial_sampler, mui_chrono_sampler_exact_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_exact_exact_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_exact_3d *spatial_sampler, mui_chrono_sampler_exact_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_exact_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_exact_3dx *spatial_sampler, mui_chrono_sampler_exact_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_exact_exact_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_exact_3t *spatial_sampler, mui_chrono_sampler_exact_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: exact; temporal sampler: gauss
void mui_fetch_exact_gauss_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1, float* t_2,
        mui_sampler_exact_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_gauss_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_exact_3fx *spatial_sampler, mui_chrono_sampler_gauss_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_exact_gauss_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_exact_3d *spatial_sampler, mui_chrono_sampler_gauss_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_gauss_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_exact_3dx *spatial_sampler, mui_chrono_sampler_gauss_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_exact_gauss_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_exact_3t *spatial_sampler, mui_chrono_sampler_gauss_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: exact; temporal sampler: mean
void mui_fetch_exact_mean_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1, float* t_2,
        mui_sampler_exact_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_mean_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_exact_3fx *spatial_sampler, mui_chrono_sampler_mean_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_exact_mean_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_exact_3d *spatial_sampler, mui_chrono_sampler_mean_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_mean_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_exact_3dx *spatial_sampler, mui_chrono_sampler_mean_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_exact_mean_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_exact_3t *spatial_sampler, mui_chrono_sampler_mean_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: exact; temporal sampler: sum
void mui_fetch_exact_sum_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1, float* t_2,
        mui_sampler_exact_3f *spatial_sampler, mui_chrono_sampler_sum_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_sum_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_exact_3fx *spatial_sampler, mui_chrono_sampler_sum_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_exact_sum_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_exact_3d *spatial_sampler, mui_chrono_sampler_sum_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_exact_sum_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_exact_3dx *spatial_sampler, mui_chrono_sampler_sum_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_exact_sum_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_exact_3t *spatial_sampler, mui_chrono_sampler_sum_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: gauss; temporal sampler: exact
void mui_fetch_gauss_exact_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1, float* t_2,
        mui_sampler_gauss_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_exact_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_gauss_3fx *spatial_sampler, mui_chrono_sampler_exact_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_gauss_exact_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_gauss_3d *spatial_sampler, mui_chrono_sampler_exact_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_exact_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_gauss_3dx *spatial_sampler, mui_chrono_sampler_exact_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_gauss_exact_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_gauss_3t *spatial_sampler, mui_chrono_sampler_exact_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: gauss; temporal sampler: gauss
void mui_fetch_gauss_gauss_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1, float* t_2,
        mui_sampler_gauss_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_gauss_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_gauss_3fx *spatial_sampler, mui_chrono_sampler_gauss_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_gauss_gauss_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_gauss_3d *spatial_sampler, mui_chrono_sampler_gauss_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_gauss_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_gauss_3dx *spatial_sampler, mui_chrono_sampler_gauss_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_gauss_gauss_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_gauss_3t *spatial_sampler, mui_chrono_sampler_gauss_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: gauss; temporal sampler: mean
void mui_fetch_gauss_mean_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1, float* t_2,
        mui_sampler_gauss_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_mean_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_gauss_3fx *spatial_sampler, mui_chrono_sampler_mean_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_gauss_mean_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_gauss_3d *spatial_sampler, mui_chrono_sampler_mean_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_gauss_mean_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_gauss_3dx *spatial_sampler, mui_chrono_sampler_mean_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_gauss_mean_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_gauss_3t *spatial_sampler, mui_chrono_sampler_mean_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: moving average; temporal sampler: exact
void mui_fetch_moving_average_exact_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_moving_average_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_exact_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3,
        float* t_1, float* t_2, mui_sampler_moving_average_3fx *spatial_sampler,
        mui_chrono_sampler_exact_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_moving_average_exact_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_moving_average_3d *spatial_sampler, mui_chrono_sampler_exact_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_exact_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_moving_average_3dx *spatial_sampler,
        mui_chrono_sampler_exact_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_moving_average_exact_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_moving_average_3t *spatial_sampler, mui_chrono_sampler_exact_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: moving average; temporal sampler: gauss
void mui_fetch_moving_average_gauss_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_moving_average_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_gauss_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3,
        float* t_1, float* t_2, mui_sampler_moving_average_3fx *spatial_sampler,
        mui_chrono_sampler_gauss_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_moving_average_gauss_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_moving_average_3d *spatial_sampler, mui_chrono_sampler_gauss_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_gauss_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_moving_average_3dx *spatial_sampler,
        mui_chrono_sampler_gauss_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_moving_average_gauss_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_moving_average_3t *spatial_sampler, mui_chrono_sampler_gauss_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: moving average; temporal sampler: mean
void mui_fetch_moving_average_mean_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_moving_average_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_mean_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_moving_average_3fx *spatial_sampler, mui_chrono_sampler_mean_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_moving_average_mean_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_moving_average_3d *spatial_sampler, mui_chrono_sampler_mean_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_mean_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_moving_average_3dx *spatial_sampler,
        mui_chrono_sampler_mean_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_moving_average_mean_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_moving_average_3t *spatial_sampler, mui_chrono_sampler_mean_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: moving average; temporal sampler: sum
void mui_fetch_moving_average_sum_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_moving_average_3f *spatial_sampler, mui_chrono_sampler_sum_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_sum_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_moving_average_3fx *spatial_sampler, mui_chrono_sampler_sum_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_moving_average_sum_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_moving_average_3d *spatial_sampler, mui_chrono_sampler_sum_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_moving_average_sum_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_moving_average_3dx *spatial_sampler,
        mui_chrono_sampler_sum_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_moving_average_sum_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_moving_average_3t *spatial_sampler, mui_chrono_sampler_sum_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: nearest neighbor; temporal sampler: exact
void mui_fetch_nearest_neighbor_exact_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_nearest_neighbor_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_exact_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3,
        float* t_1, float* t_2, mui_sampler_nearest_neighbor_3fx *spatial_sampler,
        mui_chrono_sampler_exact_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_nearest_neighbor_exact_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_nearest_neighbor_3d *spatial_sampler,
        mui_chrono_sampler_exact_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_exact_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_nearest_neighbor_3dx *spatial_sampler,
        mui_chrono_sampler_exact_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_nearest_neighbor_exact_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_nearest_neighbor_3t *spatial_sampler,
        mui_chrono_sampler_exact_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: nearest neighbor; temporal sampler: gauss
void mui_fetch_nearest_neighbor_gauss_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_nearest_neighbor_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_gauss_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3,
        float* t_1, float* t_2, mui_sampler_nearest_neighbor_3fx *spatial_sampler,
        mui_chrono_sampler_gauss_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_nearest_neighbor_gauss_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_nearest_neighbor_3d *spatial_sampler,
        mui_chrono_sampler_gauss_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_gauss_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_nearest_neighbor_3dx *spatial_sampler,
        mui_chrono_sampler_gauss_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_nearest_neighbor_gauss_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_nearest_neighbor_3t *spatial_sampler,
        mui_chrono_sampler_gauss_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: nearest neighbor; temporal sampler: mean
void mui_fetch_nearest_neighbor_mean_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_nearest_neighbor_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_mean_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3,
        float* t_1, float* t_2, mui_sampler_nearest_neighbor_3fx *spatial_sampler,
        mui_chrono_sampler_mean_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_nearest_neighbor_mean_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_nearest_neighbor_3d *spatial_sampler,
        mui_chrono_sampler_mean_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_mean_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_nearest_neighbor_3dx *spatial_sampler,
        mui_chrono_sampler_mean_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_nearest_neighbor_mean_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_nearest_neighbor_3t *spatial_sampler,
        mui_chrono_sampler_mean_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: nearest neighbor; temporal sampler: sum
void mui_fetch_nearest_neighbor_sum_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_nearest_neighbor_3f *spatial_sampler, mui_chrono_sampler_sum_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_sum_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3,
        float* t_1, float* t_2, mui_sampler_nearest_neighbor_3fx *spatial_sampler,
        mui_chrono_sampler_sum_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_nearest_neighbor_sum_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_nearest_neighbor_3d *spatial_sampler, mui_chrono_sampler_sum_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_nearest_neighbor_sum_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_nearest_neighbor_3dx *spatial_sampler,
        mui_chrono_sampler_sum_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_nearest_neighbor_sum_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_nearest_neighbor_3t *spatial_sampler, mui_chrono_sampler_sum_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: nearest neighbor; temporal sampler: exact
void mui_fetch_pseudo_nearest_neighbor_exact_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3,
        float* t_1, float* t_2, mui_sampler_pseudo_nearest_neighbor_3f *spatial_sampler,
        mui_chrono_sampler_exact_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_exact_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3,
        float* t_1, float* t_2, mui_sampler_pseudo_nearest_neighbor_3fx *spatial_sampler,
        mui_chrono_sampler_exact_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_exact_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_3d *spatial_sampler,
        mui_chrono_sampler_exact_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_exact_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_3dx *spatial_sampler,
        mui_chrono_sampler_exact_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_exact_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_3t *spatial_sampler,
        mui_chrono_sampler_exact_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: nearest neighbor; temporal sampler: gauss
void mui_fetch_pseudo_nearest_neighbor_gauss_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3,
        float* t_1, float* t_2, mui_sampler_pseudo_nearest_neighbor_3f *spatial_sampler,
        mui_chrono_sampler_gauss_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_gauss_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3,
        float* t_1, float* t_2, mui_sampler_pseudo_nearest_neighbor_3fx *spatial_sampler,
        mui_chrono_sampler_gauss_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_gauss_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_3d *spatial_sampler,
        mui_chrono_sampler_gauss_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_gauss_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_3dx *spatial_sampler,
        mui_chrono_sampler_gauss_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_gauss_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_3t *spatial_sampler,
        mui_chrono_sampler_gauss_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: nearest neighbor; temporal sampler: mean
void mui_fetch_pseudo_nearest_neighbor_mean_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3,
        float* t_1, float* t_2, mui_sampler_pseudo_nearest_neighbor_3f *spatial_sampler,
        mui_chrono_sampler_mean_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_mean_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3,
        float* t_1, float* t_2, mui_sampler_pseudo_nearest_neighbor_3fx *spatial_sampler,
        mui_chrono_sampler_mean_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_mean_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_3d *spatial_sampler,
        mui_chrono_sampler_mean_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_mean_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_3dx *spatial_sampler,
        mui_chrono_sampler_mean_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_mean_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_3t *spatial_sampler,
        mui_chrono_sampler_mean_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: nearest neighbor; temporal sampler: sum
void mui_fetch_pseudo_nearest_neighbor_sum_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3,
        float* t_1, float* t_2, mui_sampler_pseudo_nearest_neighbor_3f *spatial_sampler,
        mui_chrono_sampler_sum_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_sum_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3,
        float* t_1, float* t_2, mui_sampler_pseudo_nearest_neighbor_3fx *spatial_sampler,
        mui_chrono_sampler_sum_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_sum_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_3d *spatial_sampler,
        mui_chrono_sampler_sum_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_sum_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_3dx *spatial_sampler,
        mui_chrono_sampler_sum_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_pseudo_nearest_neighbor_sum_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_pseudo_nearest_neighbor_3t *spatial_sampler,
        mui_chrono_sampler_sum_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: shepard quintic; temporal sampler: exact
void mui_fetch_shepard_quintic_exact_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_shepard_quintic_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_exact_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3,
        float* t_1, float* t_2, mui_sampler_shepard_quintic_3fx *spatial_sampler,
        mui_chrono_sampler_exact_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_shepard_quintic_exact_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_shepard_quintic_3d *spatial_sampler,
        mui_chrono_sampler_exact_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_exact_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_shepard_quintic_3dx *spatial_sampler,
        mui_chrono_sampler_exact_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_shepard_quintic_exact_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_shepard_quintic_3t *spatial_sampler,
        mui_chrono_sampler_exact_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: shepard quintic; temporal sampler: gauss
void mui_fetch_shepard_quintic_gauss_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_shepard_quintic_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_gauss_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3,
        float* t_1, float* t_2, mui_sampler_shepard_quintic_3fx *spatial_sampler,
        mui_chrono_sampler_gauss_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_shepard_quintic_gauss_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_shepard_quintic_3d *spatial_sampler,
        mui_chrono_sampler_gauss_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_gauss_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_shepard_quintic_3dx *spatial_sampler,
        mui_chrono_sampler_gauss_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_shepard_quintic_gauss_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_shepard_quintic_3t *spatial_sampler,
        mui_chrono_sampler_gauss_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: shepard quintic; temporal sampler: mean
void mui_fetch_shepard_quintic_mean_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_shepard_quintic_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_mean_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3,
        float* t_1, float* t_2, mui_sampler_shepard_quintic_3fx *spatial_sampler,
        mui_chrono_sampler_mean_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_shepard_quintic_mean_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_shepard_quintic_3d *spatial_sampler, mui_chrono_sampler_mean_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_mean_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_shepard_quintic_3dx *spatial_sampler,
        mui_chrono_sampler_mean_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_shepard_quintic_mean_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_shepard_quintic_3t *spatial_sampler, mui_chrono_sampler_mean_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: shepard quintic; temporal sampler: sum
void mui_fetch_shepard_quintic_sum_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_shepard_quintic_3f *spatial_sampler, mui_chrono_sampler_sum_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_sum_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_shepard_quintic_3fx *spatial_sampler, mui_chrono_sampler_sum_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_shepard_quintic_sum_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_shepard_quintic_3d *spatial_sampler, mui_chrono_sampler_sum_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_shepard_quintic_sum_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3,
        double* t_1, double* t_2, mui_sampler_shepard_quintic_3dx *spatial_sampler,
        mui_chrono_sampler_sum_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_shepard_quintic_sum_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_shepard_quintic_3t *spatial_sampler, mui_chrono_sampler_sum_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: sph-derived quintic; temporal sampler: exact
void mui_fetch_sph_quintic_exact_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_sph_quintic_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_exact_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_sph_quintic_3fx *spatial_sampler, mui_chrono_sampler_exact_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sph_quintic_exact_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_sph_quintic_3d *spatial_sampler, mui_chrono_sampler_exact_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_exact_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_sph_quintic_3dx *spatial_sampler, mui_chrono_sampler_exact_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sph_quintic_exact_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_sph_quintic_3t *spatial_sampler, mui_chrono_sampler_exact_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: sph-derived quintic; temporal sampler: gauss
void mui_fetch_sph_quintic_gauss_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_sph_quintic_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_gauss_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_sph_quintic_3fx *spatial_sampler, mui_chrono_sampler_gauss_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sph_quintic_gauss_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_sph_quintic_3d *spatial_sampler, mui_chrono_sampler_gauss_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_gauss_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_sph_quintic_3dx *spatial_sampler, mui_chrono_sampler_gauss_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sph_quintic_gauss_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_sph_quintic_3t *spatial_sampler, mui_chrono_sampler_gauss_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: sph-derived quintic; temporal sampler: mean
void mui_fetch_sph_quintic_mean_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_sph_quintic_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_mean_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_sph_quintic_3fx *spatial_sampler, mui_chrono_sampler_mean_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sph_quintic_mean_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_sph_quintic_3d *spatial_sampler, mui_chrono_sampler_mean_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_mean_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_sph_quintic_3dx *spatial_sampler, mui_chrono_sampler_mean_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sph_quintic_mean_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_sph_quintic_3t *spatial_sampler, mui_chrono_sampler_mean_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: sph-derived quintic; temporal sampler: sum
void mui_fetch_sph_quintic_sum_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_sph_quintic_3f *spatial_sampler, mui_chrono_sampler_sum_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_sum_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_sph_quintic_3fx *spatial_sampler, mui_chrono_sampler_sum_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sph_quintic_sum_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_sph_quintic_3d *spatial_sampler, mui_chrono_sampler_sum_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sph_quintic_sum_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_sph_quintic_3dx *spatial_sampler, mui_chrono_sampler_sum_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sph_quintic_sum_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_sph_quintic_3t *spatial_sampler, mui_chrono_sampler_sum_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: summation quintic; temporal sampler: exact
void mui_fetch_sum_quintic_exact_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_sum_quintic_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_exact_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_sum_quintic_3fx *spatial_sampler, mui_chrono_sampler_exact_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sum_quintic_exact_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_sum_quintic_3d *spatial_sampler, mui_chrono_sampler_exact_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_exact_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_sum_quintic_3dx *spatial_sampler, mui_chrono_sampler_exact_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sum_quintic_exact_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_sum_quintic_3t *spatial_sampler, mui_chrono_sampler_exact_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: summation quintic; temporal sampler: gauss
void mui_fetch_sum_quintic_gauss_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_sum_quintic_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_gauss_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_sum_quintic_3fx *spatial_sampler, mui_chrono_sampler_gauss_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sum_quintic_gauss_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_sum_quintic_3d *spatial_sampler, mui_chrono_sampler_gauss_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_gauss_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_sum_quintic_3dx *spatial_sampler, mui_chrono_sampler_gauss_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sum_quintic_gauss_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_sum_quintic_3t *spatial_sampler, mui_chrono_sampler_gauss_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: summation quintic; temporal sampler: mean
void mui_fetch_sum_quintic_mean_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_sum_quintic_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_mean_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_sum_quintic_3fx *spatial_sampler, mui_chrono_sampler_mean_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sum_quintic_mean_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_sum_quintic_3d *spatial_sampler, mui_chrono_sampler_mean_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_mean_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_sum_quintic_3dx *spatial_sampler, mui_chrono_sampler_mean_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sum_quintic_mean_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_sum_quintic_3t *spatial_sampler, mui_chrono_sampler_mean_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: summation quintic; temporal sampler: sum
void mui_fetch_sum_quintic_sum_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_sum_quintic_3f *spatial_sampler, mui_chrono_sampler_sum_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_sum_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_sum_quintic_3fx *spatial_sampler, mui_chrono_sampler_sum_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sum_quintic_sum_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_sum_quintic_3d *spatial_sampler, mui_chrono_sampler_sum_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_sum_quintic_sum_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_sum_quintic_3dx *spatial_sampler, mui_chrono_sampler_sum_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_sum_quintic_sum_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_sum_quintic_3t *spatial_sampler, mui_chrono_sampler_sum_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

#ifdef USE_RBF
// Spatial sampler: radial basis function; temporal sampler: exact
void mui_fetch_rbf_exact_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1, float* t_2,
        mui_sampler_rbf_3f *spatial_sampler, mui_chrono_sampler_exact_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_exact_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_rbf_3fx *spatial_sampler, mui_chrono_sampler_exact_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_rbf_exact_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_rbf_3d *spatial_sampler, mui_chrono_sampler_exact_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_exact_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_rbf_3dx *spatial_sampler, mui_chrono_sampler_exact_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_rbf_exact_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_rbf_3t *spatial_sampler, mui_chrono_sampler_exact_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: radial basis function; temporal sampler: gauss
void mui_fetch_rbf_gauss_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1, float* t_2,
        mui_sampler_rbf_3f *spatial_sampler, mui_chrono_sampler_gauss_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_gauss_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1,
        float* t_2, mui_sampler_rbf_3fx *spatial_sampler, mui_chrono_sampler_gauss_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_rbf_gauss_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_rbf_3d *spatial_sampler, mui_chrono_sampler_gauss_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_gauss_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_rbf_3dx *spatial_sampler, mui_chrono_sampler_gauss_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_rbf_gauss_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_rbf_3t *spatial_sampler, mui_chrono_sampler_gauss_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: radial basis function; temporal sampler: mean
void mui_fetch_rbf_mean_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1, float* t_2,
        mui_sampler_rbf_3f *spatial_sampler, mui_chrono_sampler_mean_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_mean_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1, float* t_2,
        mui_sampler_rbf_3fx *spatial_sampler, mui_chrono_sampler_mean_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_rbf_mean_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1, double* t_2,
        mui_sampler_rbf_3d *spatial_sampler, mui_chrono_sampler_mean_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_mean_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_rbf_3dx *spatial_sampler, mui_chrono_sampler_mean_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_rbf_mean_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1, double* t_2,
        mui_sampler_rbf_3t *spatial_sampler, mui_chrono_sampler_mean_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}

// Spatial sampler: radial basis function; temporal sampler: sum
void mui_fetch_rbf_sum_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1, float* t_2,
        mui_sampler_rbf_3f *spatial_sampler, mui_chrono_sampler_sum_3f *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3f(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_sum_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float* point_1, float* point_2, float* point_3, float* t_1, float* t_2,
        mui_sampler_rbf_3fx *spatial_sampler, mui_chrono_sampler_sum_3fx *temporal_sampler, float *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3fx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_rbf_sum_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1, double* t_2,
        mui_sampler_rbf_3d *spatial_sampler, mui_chrono_sampler_sum_3d *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3d(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler, *temporal_sampler);
}

void mui_fetch_rbf_sum_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1,
        double* t_2, mui_sampler_rbf_3dx *spatial_sampler, mui_chrono_sampler_sum_3dx *temporal_sampler, double *return_value) {
    *return_value = uniface->fetch(std::string(attr), mui::point3dx(*point_1,*point_2,*point_3), *t_1, *t_2, *spatial_sampler,
            *temporal_sampler);
}

void mui_fetch_rbf_sum_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double* point_1, double* point_2, double* point_3, double* t_1, double* t_2,
        mui_sampler_rbf_3t *spatial_sampler, mui_chrono_sampler_sum_3t *temporal_sampler, double *return_value) {
    mui::mui_f_wrapper_3D::point_type point_fetch(static_cast<mui::mui_f_wrapper_3D::REAL>(*point_1,*point_2,*point_3));
    *return_value = static_cast<double>(uniface->fetch(std::string(attr), point_fetch,
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2),
            *spatial_sampler, *temporal_sampler));
}
#endif

/*******************************************************************
 * MUI functions for 3d data point only fetch using one time value  *
 ********************************************************************/

// Temporal sampler: exact
void mui_fetch_points_exact_3f_f(mui_uniface_3f *uniface, const char *attr, float *t,
        mui_chrono_sampler_exact_3f *temporal_sampler, float** ret_points_1, float** ret_points_2, float** ret_points_3, int *num_points) {
    std::vector<mui::point3f> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_3 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_exact_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *t,
        mui_chrono_sampler_exact_3fx *temporal_sampler, float** ret_points_1, float** ret_points_2, float** ret_points_3, int *num_points) {
    std::vector<mui::point3fx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_3 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_exact_3d_f(mui_uniface_3d *uniface, const char *attr, double *t,
        mui_chrono_sampler_exact_3d *temporal_sampler, double** ret_points_1, double** ret_points_2, double** ret_points_3, int *num_points) {
    std::vector<mui::point3d> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_3 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_exact_3dx_f(mui_uniface_3dx *uniface, const char *attr, double *t,
        mui_chrono_sampler_exact_3dx *temporal_sampler, double** ret_points_1, double** ret_points_2, double** ret_points_3, int *num_points) {
    std::vector<mui::point3dx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_3 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_exact_3t_f(mui_uniface_3t *uniface, const char *attr, double *t,
        mui_chrono_sampler_exact_3t *temporal_sampler, double** ret_points_1, double** ret_points_2, double** ret_points_3, int *num_points) {
    std::vector<mui::mui_f_wrapper_3D::point_type> ret_pts = uniface->fetch_points<float>(std::string(attr), *t,
            *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_3 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = static_cast<double>(ret_pts[i][0]);
        *ret_points_2[i] = static_cast<double>(ret_pts[i][1]);
        *ret_points_3[i] = static_cast<double>(ret_pts[i][2]);
    }
    *num_points = static_cast<int>(ret_pts.size());
}

// Temporal sampler: gauss
void mui_fetch_points_gauss_3f_f(mui_uniface_3f *uniface, const char *attr, float *t,
        mui_chrono_sampler_gauss_3f *temporal_sampler, float** ret_points_1, float** ret_points_2, float** ret_points_3, int *num_points) {
    std::vector<mui::point3f> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_3 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_gauss_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *t,
        mui_chrono_sampler_gauss_3fx *temporal_sampler, float** ret_points_1, float** ret_points_2, float** ret_points_3, int *num_points) {
    std::vector<mui::point3fx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_3 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_gauss_3d_f(mui_uniface_3d *uniface, const char *attr, double *t,
        mui_chrono_sampler_gauss_3d *temporal_sampler, double** ret_points_1, double** ret_points_2, double** ret_points_3, int *num_points) {
    std::vector<mui::point3d> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_3 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_gauss_3dx_f(mui_uniface_3dx *uniface, const char *attr, double *t,
        mui_chrono_sampler_gauss_3dx *temporal_sampler, double** ret_points_1, double** ret_points_2, double** ret_points_3, int *num_points) {
    std::vector<mui::point3dx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_3 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_gauss_3t_f(mui_uniface_3t *uniface, const char *attr, double *t,
        mui_chrono_sampler_gauss_3t *temporal_sampler, double** ret_points_1, double** ret_points_2, double** ret_points_3, int *num_points) {
    std::vector<mui::mui_f_wrapper_3D::point_type> ret_pts = uniface->fetch_points<float>(std::string(attr), *t,
            *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_3 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = static_cast<double>(ret_pts[i][0]);
        *ret_points_2[i] = static_cast<double>(ret_pts[i][1]);
        *ret_points_3[i] = static_cast<double>(ret_pts[i][2]);
    }
    *num_points = static_cast<int>(ret_pts.size());
}

// Temporal sampler: mean
void mui_fetch_points_mean_3f_f(mui_uniface_3f *uniface, const char *attr, float *t,
        mui_chrono_sampler_mean_3f *temporal_sampler, float** ret_points_1, float** ret_points_2, float** ret_points_3, int *num_points) {
    std::vector<mui::point3f> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_3 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_mean_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *t,
        mui_chrono_sampler_mean_3fx *temporal_sampler, float** ret_points_1, float** ret_points_2, float** ret_points_3, int *num_points) {
    std::vector<mui::point3fx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_3 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_mean_3d_f(mui_uniface_3d *uniface, const char *attr, double *t,
        mui_chrono_sampler_mean_3d *temporal_sampler, double** ret_points_1, double** ret_points_2, double** ret_points_3, int *num_points) {
    std::vector<mui::point3d> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_3 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_mean_3dx_f(mui_uniface_3dx *uniface, const char *attr, double *t,
        mui_chrono_sampler_mean_3dx *temporal_sampler, double** ret_points_1, double** ret_points_2, double** ret_points_3, int *num_points) {
    std::vector<mui::point3dx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_3 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_mean_3t_f(mui_uniface_3t *uniface, const char *attr, double *t,
        mui_chrono_sampler_mean_3t *temporal_sampler, double** ret_points_1, double** ret_points_2, double** ret_points_3, int *num_points) {
    std::vector<mui::mui_f_wrapper_3D::point_type> ret_pts = uniface->fetch_points<float>(std::string(attr), *t,
            *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_3 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = static_cast<double>(ret_pts[i][0]);
        *ret_points_2[i] = static_cast<double>(ret_pts[i][1]);
        *ret_points_3[i] = static_cast<double>(ret_pts[i][2]);
    }
    *num_points = static_cast<int>(ret_pts.size());
}

// Temporal sampler: sum
void mui_fetch_points_sum_3f_f(mui_uniface_3f *uniface, const char *attr, float *t,
        mui_chrono_sampler_sum_3f *temporal_sampler, float** ret_points_1, float** ret_points_2, float** ret_points_3, int *num_points) {
    std::vector<mui::point3f> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_3 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_sum_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *t,
        mui_chrono_sampler_sum_3fx *temporal_sampler, float** ret_points_1, float** ret_points_2, float** ret_points_3, int *num_points) {
    std::vector<mui::point3fx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_3 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_sum_3d_f(mui_uniface_3d *uniface, const char *attr, double *t,
        mui_chrono_sampler_sum_3d *temporal_sampler, double** ret_points_1, double** ret_points_2, double** ret_points_3, int *num_points) {
    std::vector<mui::point3d> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_3 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_sum_3dx_f(mui_uniface_3dx *uniface, const char *attr, double *t,
        mui_chrono_sampler_sum_3dx *temporal_sampler, double** ret_points_1, double** ret_points_2, double** ret_points_3, int *num_points) {
    std::vector<mui::point3dx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_3 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_sum_3t_f(mui_uniface_3t *uniface, const char *attr, double *t,
        mui_chrono_sampler_sum_3t *temporal_sampler, double** ret_points_1, double** ret_points_2, double** ret_points_3, int *num_points) {
    std::vector<mui::mui_f_wrapper_3D::point_type> ret_pts = uniface->fetch_points<float>(std::string(attr), *t,
            *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_3 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = static_cast<double>(ret_pts[i][0]);
        *ret_points_2[i] = static_cast<double>(ret_pts[i][1]);
        *ret_points_3[i] = static_cast<double>(ret_pts[i][2]);
    }
    *num_points = static_cast<int>(ret_pts.size());
}

/*******************************************************************
 * MUI functions for 3d data point only fetch using two time values *
 ********************************************************************/

// Temporal sampler: exact
void mui_fetch_points_exact_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_exact_3f *temporal_sampler, float** ret_points_1, float** ret_points_2, float** ret_points_3, int *num_points) {
    std::vector<mui::point3f> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_3 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_exact_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_exact_3fx *temporal_sampler, float** ret_points_1, float** ret_points_2, float** ret_points_3, int *num_points) {
    std::vector<mui::point3fx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_3 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_exact_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_exact_3d *temporal_sampler, double** ret_points_1, double** ret_points_2, double** ret_points_3, int *num_points) {
    std::vector<mui::point3d> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_3 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_exact_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_exact_3dx *temporal_sampler, double** ret_points_1, double** ret_points_2, double** ret_points_3, int *num_points) {
    std::vector<mui::point3dx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_3 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_exact_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_exact_3t *temporal_sampler, double** ret_points_1, double** ret_points_2, double** ret_points_3, int *num_points) {
    std::vector<mui::mui_f_wrapper_3D::point_type> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2,
            *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_3 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = static_cast<double>(ret_pts[i][0]);
        *ret_points_2[i] = static_cast<double>(ret_pts[i][1]);
        *ret_points_3[i] = static_cast<double>(ret_pts[i][2]);
    }
    *num_points = static_cast<int>(ret_pts.size());
}

// Temporal sampler: gauss
void mui_fetch_points_gauss_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_gauss_3f *temporal_sampler, float** ret_points_1, float** ret_points_2, float** ret_points_3, int *num_points) {
    std::vector<mui::point3f> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_3 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_gauss_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_gauss_3fx *temporal_sampler, float** ret_points_1, float** ret_points_2, float** ret_points_3, int *num_points) {
    std::vector<mui::point3fx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_3 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_gauss_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_gauss_3d *temporal_sampler, double** ret_points_1, double** ret_points_2, double** ret_points_3, int *num_points) {
    std::vector<mui::point3d> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_3 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_gauss_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_gauss_3dx *temporal_sampler, double** ret_points_1, double** ret_points_2, double** ret_points_3, int *num_points) {
    std::vector<mui::point3dx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_3 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_gauss_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_gauss_3t *temporal_sampler, double** ret_points_1, double** ret_points_2, double** ret_points_3, int *num_points) {
    std::vector<mui::mui_f_wrapper_3D::point_type> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2,
            *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_3 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = static_cast<double>(ret_pts[i][0]);
        *ret_points_2[i] = static_cast<double>(ret_pts[i][1]);
        *ret_points_3[i] = static_cast<double>(ret_pts[i][2]);
    }
    *num_points = static_cast<int>(ret_pts.size());
}

// Temporal sampler: mean
void mui_fetch_points_mean_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_mean_3f *temporal_sampler, float** ret_points_1, float** ret_points_2, float** ret_points_3, int *num_points) {
    std::vector<mui::point3f> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_3 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_mean_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_mean_3fx *temporal_sampler, float** ret_points_1, float** ret_points_2, float** ret_points_3, int *num_points) {
    std::vector<mui::point3fx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_3 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_mean_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_mean_3d *temporal_sampler, double** ret_points_1, double** ret_points_2, double** ret_points_3, int *num_points) {
    std::vector<mui::point3d> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_3 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_mean_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_mean_3dx *temporal_sampler, double** ret_points_1, double** ret_points_2, double** ret_points_3, int *num_points) {
    std::vector<mui::point3dx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_3 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_mean_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_mean_3t *temporal_sampler, double** ret_points_1, double** ret_points_2, double** ret_points_3, int *num_points) {
    std::vector<mui::mui_f_wrapper_3D::point_type> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2,
            *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_3 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = static_cast<double>(ret_pts[i][0]);
        *ret_points_2[i] = static_cast<double>(ret_pts[i][1]);
        *ret_points_3[i] = static_cast<double>(ret_pts[i][2]);
    }
    *num_points = static_cast<int>(ret_pts.size());
}

// Temporal sampler: sum
void mui_fetch_points_sum_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_sum_3f *temporal_sampler, float** ret_points_1, float** ret_points_2, float** ret_points_3, int *num_points) {
    std::vector<mui::point3f> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_3 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_sum_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_sum_3fx *temporal_sampler, float** ret_points_1, float** ret_points_2, float** ret_points_3, int *num_points) {
    std::vector<mui::point3fx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_2 = (float*) malloc(ret_pts.size() * sizeof(float));
    *ret_points_3 = (float*) malloc(ret_pts.size() * sizeof(float));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_sum_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_sum_3d *temporal_sampler, double** ret_points_1, double** ret_points_2, double** ret_points_3, int *num_points) {
    std::vector<mui::point3d> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_3 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_sum_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_sum_3dx *temporal_sampler, double** ret_points_1, double** ret_points_2, double** ret_points_3, int *num_points) {
    std::vector<mui::point3dx> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_3 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = ret_pts[i][0];
        *ret_points_2[i] = ret_pts[i][1];
        *ret_points_3[i] = ret_pts[i][2];
    }
    *num_points = static_cast<int>(ret_pts.size());
}

void mui_fetch_points_sum_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_sum_3t *temporal_sampler, double** ret_points_1, double** ret_points_2, double** ret_points_3, int *num_points) {
    std::vector<mui::mui_f_wrapper_3D::point_type> ret_pts = uniface->fetch_points<float>(std::string(attr), *t_1, *t_2,
            *temporal_sampler);
    *ret_points_1 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_2 = (double*) malloc(ret_pts.size() * sizeof(double));
    *ret_points_3 = (double*) malloc(ret_pts.size() * sizeof(double));
    for (size_t i = 0; i < ret_pts.size(); i++) {
        *ret_points_1[i] = static_cast<double>(ret_pts[i][0]);
        *ret_points_2[i] = static_cast<double>(ret_pts[i][1]);
        *ret_points_3[i] = static_cast<double>(ret_pts[i][2]);
    }
    *num_points = static_cast<int>(ret_pts.size());
}

/*******************************************************************
 * MUI functions for 3d data values only fetch using one time value  *
 ********************************************************************/

// Temporal sampler: exact
void mui_fetch_values_exact_3f_f(mui_uniface_3f *uniface, const char *attr, float *t,
        mui_chrono_sampler_exact_3f *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_exact_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *t,
        mui_chrono_sampler_exact_3fx *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_exact_3d_f(mui_uniface_3d *uniface, const char *attr, double *t,
        mui_chrono_sampler_exact_3d *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_exact_3dx_f(mui_uniface_3dx *uniface, const char *attr, double *t,
        mui_chrono_sampler_exact_3dx *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_exact_3t_f(mui_uniface_3t *uniface, const char *attr, double *t,
        mui_chrono_sampler_exact_3t *temporal_sampler, double **values, int *num_values) {
    std::vector<mui::mui_f_wrapper_3D::REAL> ret_vals = uniface->fetch_values<mui::mui_f_wrapper_3D::REAL>(
            std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = static_cast<double>(ret_vals[i]);
    }
    *num_values = static_cast<int>(ret_vals.size());
}

// Temporal sampler: gauss
void mui_fetch_values_gauss_3f_f(mui_uniface_3f *uniface, const char *attr, float *t,
        mui_chrono_sampler_gauss_3f *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_gauss_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *t,
        mui_chrono_sampler_gauss_3fx *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_gauss_3d_f(mui_uniface_3d *uniface, const char *attr, double *t,
        mui_chrono_sampler_gauss_3d *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_gauss_3dx_f(mui_uniface_3dx *uniface, const char *attr, double *t,
        mui_chrono_sampler_gauss_3dx *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_gauss_3t_f(mui_uniface_3t *uniface, const char *attr, double *t,
        mui_chrono_sampler_gauss_3t *temporal_sampler, double **values, int *num_values) {
    std::vector<mui::mui_f_wrapper_3D::REAL> ret_vals = uniface->fetch_values<mui::mui_f_wrapper_3D::REAL>(
            std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = static_cast<double>(ret_vals[i]);
    }
    *num_values = static_cast<int>(ret_vals.size());
}

// Temporal sampler: mean
void mui_fetch_values_mean_3f_f(mui_uniface_3f *uniface, const char *attr, float *t,
        mui_chrono_sampler_mean_3f *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_mean_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *t,
        mui_chrono_sampler_mean_3fx *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_mean_3d_f(mui_uniface_3d *uniface, const char *attr, double *t,
        mui_chrono_sampler_mean_3d *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_mean_3dx_f(mui_uniface_3dx *uniface, const char *attr, double *t,
        mui_chrono_sampler_mean_3dx *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_mean_3t_f(mui_uniface_3t *uniface, const char *attr, double *t,
        mui_chrono_sampler_mean_3t *temporal_sampler, double **values, int *num_values) {
    std::vector<mui::mui_f_wrapper_3D::REAL> ret_vals = uniface->fetch_values<mui::mui_f_wrapper_3D::REAL>(
            std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = static_cast<double>(ret_vals[i]);
    }
    *num_values = static_cast<int>(ret_vals.size());
}

// Temporal sampler: sum
void mui_fetch_values_sum_3f_f(mui_uniface_3f *uniface, const char *attr, float *t,
        mui_chrono_sampler_sum_3f *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_sum_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *t,
        mui_chrono_sampler_sum_3fx *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_sum_3d_f(mui_uniface_3d *uniface, const char *attr, double *t,
        mui_chrono_sampler_sum_3d *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_sum_3dx_f(mui_uniface_3dx *uniface, const char *attr, double *t,
        mui_chrono_sampler_sum_3dx *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_sum_3t_f(mui_uniface_3t *uniface, const char *attr, double *t,
        mui_chrono_sampler_sum_3t *temporal_sampler, double **values, int *num_values) {
    std::vector<mui::mui_f_wrapper_3D::REAL> ret_vals = uniface->fetch_values<mui::mui_f_wrapper_3D::REAL>(
            std::string(attr), *t, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = static_cast<double>(ret_vals[i]);
    }
    *num_values = static_cast<int>(ret_vals.size());
}

/********************************************************************
 * MUI functions for 3d data values only fetch using two time values *
 *********************************************************************/

// Temporal sampler: exact
void mui_fetch_values_exact_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_exact_3f *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_exact_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_exact_3fx *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_exact_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_exact_3d *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_exact_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_exact_3dx *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_exact_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_exact_3t *temporal_sampler, double **values, int *num_values) {
    std::vector<mui::mui_f_wrapper_3D::REAL> ret_vals = uniface->fetch_values<mui::mui_f_wrapper_3D::REAL>(
            std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = static_cast<double>(ret_vals[i]);
    }
    *num_values = static_cast<int>(ret_vals.size());
}

// Temporal sampler: gauss
void mui_fetch_values_gauss_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_gauss_3f *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_gauss_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_gauss_3fx *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_gauss_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_gauss_3d *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_gauss_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_gauss_3dx *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_gauss_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_gauss_3t *temporal_sampler, double **values, int *num_values) {
    std::vector<mui::mui_f_wrapper_3D::REAL> ret_vals = uniface->fetch_values<mui::mui_f_wrapper_3D::REAL>(
            std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = static_cast<double>(ret_vals[i]);
    }
    *num_values = static_cast<int>(ret_vals.size());
}

// Temporal sampler: mean
void mui_fetch_values_mean_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_mean_3f *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_mean_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_mean_3fx *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_mean_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_mean_3d *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_mean_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_mean_3dx *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_mean_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_mean_3t *temporal_sampler, double **values, int *num_values) {
    std::vector<mui::mui_f_wrapper_3D::REAL> ret_vals = uniface->fetch_values<mui::mui_f_wrapper_3D::REAL>(
            std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = static_cast<double>(ret_vals[i]);
    }
    *num_values = static_cast<int>(ret_vals.size());
}

// Temporal sampler: sum
void mui_fetch_values_sum_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_sum_3f *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_sum_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float *t_1, float *t_2,
        mui_chrono_sampler_sum_3fx *temporal_sampler, float **values, int *num_values) {
    std::vector<float> ret_vals = uniface->fetch_values<float>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (float*) malloc(ret_vals.size() * sizeof(float));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_sum_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_sum_3d *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_sum_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_sum_3dx *temporal_sampler, double **values, int *num_values) {
    std::vector<double> ret_vals = uniface->fetch_values<double>(std::string(attr), *t_1, *t_2, *temporal_sampler);
    *values = (double*) malloc(ret_vals.size() * sizeof(double));
    for (size_t i = 0; i < ret_vals.size(); i++) {
        *values[i] = ret_vals[i];
    }
    *num_values = static_cast<int>(ret_vals.size());
}

void mui_fetch_values_sum_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double *t_1, double *t_2,
        mui_chrono_sampler_sum_3t *temporal_sampler, double **values, int *num_values) {
    std::vector<mui::mui_f_wrapper_3D::REAL> ret_vals = uniface->fetch_values<mui::mui_f_wrapper_3D::REAL>(
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

void mui_fetch_3f_param_f(mui_uniface_3f *uniface, const char *attr, float *return_value) {
    *return_value = uniface->fetch<float>(std::string(attr));
}

void mui_fetch_3fx_param_f(mui_uniface_3fx *uniface, const char *attr, float *return_value) {
    *return_value = uniface->fetch<float>(std::string(attr));
}

void mui_fetch_3d_param_f(mui_uniface_3d *uniface, const char *attr, float *return_value) {
    *return_value = uniface->fetch<double>(std::string(attr));
}

void mui_fetch_3dx_param_f(mui_uniface_3dx *uniface, const char *attr, float *return_value) {
    *return_value = uniface->fetch<double>(std::string(attr));
}

void mui_fetch_3t_param_f(mui_uniface_3t *uniface, const char *attr, float *return_value) {
    mui::mui_f_wrapper_3D::REAL ret_val = uniface->fetch<mui::mui_f_wrapper_3D::REAL>(std::string(attr));
    *return_value = static_cast<double>(ret_val);
}

/******************************************
 * MUI data receive test functions         *
 ******************************************/

// Data ready test using single time value
void mui_is_ready_3f_f(mui_uniface_3f *uniface, const char *attr, float *t, int *return_value) {
    *return_value = uniface->is_ready(std::string(attr), *t);
}

void mui_is_ready_3fx_f(mui_uniface_3fx *uniface, const char *attr, float *t, int *return_value) {
    *return_value = uniface->is_ready(std::string(attr), *t);
}

void mui_is_ready_3d_f(mui_uniface_3d *uniface, const char *attr, double *t, int *return_value) {
    *return_value = uniface->is_ready(std::string(attr), *t);
}

void mui_is_ready_3dx_f(mui_uniface_3dx *uniface, const char *attr, double *t, int *return_value) {
    *return_value = uniface->is_ready(std::string(attr), *t);
}

void mui_is_ready_3t_f(mui_uniface_3t *uniface, const char *attr, double *t, int *return_value) {
    *return_value = uniface->is_ready(std::string(attr), static_cast<mui::mui_f_wrapper_3D::time_type>(*t));
}

// Data ready test using two time values
void mui_is_ready_3f_pair_f(mui_uniface_3f *uniface, const char *attr, float *t_1, float *t_2, int *return_value) {
    *return_value = uniface->is_ready(std::string(attr), *t_1, *t_2);
}

void mui_is_ready_3fx_pair_f(mui_uniface_3fx *uniface, const char *attr, float *t_1, float *t_2, int *return_value) {
    *return_value = uniface->is_ready(std::string(attr), *t_1, *t_2);
}

void mui_is_ready_3d_pair_f(mui_uniface_3d *uniface, const char *attr, double *t_1, double *t_2, int *return_value) {
    *return_value = uniface->is_ready(std::string(attr), *t_1, *t_2);
}

void mui_is_ready_3dx_pair_f(mui_uniface_3dx *uniface, const char *attr, double *t_1, double *t_2, int *return_value) {
    *return_value = uniface->is_ready(std::string(attr), *t_1, *t_2);
}

void mui_is_ready_3t_pair_f(mui_uniface_3t *uniface, const char *attr, double *t_1, double *t_2, int *return_value) {
    *return_value = uniface->is_ready(std::string(attr), static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1),
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2));
}

/******************************************
 * MUI Smart Send functions                *
 ******************************************/

// Send span announce using 3d box geometry
void mui_announce_send_span_3f_box_f(mui_uniface_3f *uniface, float *box_1_1, float *box_1_2, float *box_1_3, float *box_2_1, float *box_2_2, float *box_2_3, float *t_start,
        float *t_timeout, int *synchronised) {
    mui::point3f point_1(*box_1_1,*box_1_2,*box_1_3);
    mui::point3f point_2(*box_2_1,*box_2_2,*box_2_3);
    mui::geometry::box3f bound_box(point_1, point_2);
    uniface->announce_send_span(*t_start, *t_timeout, bound_box, static_cast<bool>(*synchronised));
}

void mui_announce_send_span_3fx_box_f(mui_uniface_3fx *uniface, float *box_1_1, float *box_1_2, float *box_1_3, float *box_2_1, float *box_2_2, float *box_2_3, float *t_start,
        float *t_timeout, int *synchronised) {
    mui::point3fx point_1(*box_1_1,*box_1_2,*box_1_3);
    mui::point3fx point_2(*box_2_1,*box_2_2,*box_2_3);
    mui::geometry::box3fx bound_box(point_1, point_2);
    uniface->announce_send_span(*t_start, *t_timeout, bound_box, static_cast<bool>(*synchronised));
}

void mui_announce_send_span_3d_box_f(mui_uniface_3d *uniface, double *box_1_1, double *box_1_2, double *box_1_3, double *box_2_1, double *box_2_2, double *box_2_3, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::point3d point_1(*box_1_1,*box_1_2,*box_1_3);
    mui::point3d point_2(*box_2_1,*box_2_2,*box_2_3);
    mui::geometry::box3d bound_box(point_1, point_2);
    uniface->announce_send_span(*t_start, *t_timeout, bound_box, static_cast<bool>(*synchronised));
}

void mui_announce_send_span_3dx_box_f(mui_uniface_3dx *uniface, double *box_1_1, double *box_1_2, double *box_1_3, double *box_2_1, double *box_2_2, double *box_2_3, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::point3dx point_1(*box_1_1,*box_1_2,*box_1_3);
    mui::point3dx point_2(*box_2_1,*box_2_2,*box_2_3);
    mui::geometry::box3dx bound_box(point_1, point_2);
    uniface->announce_send_span(*t_start, *t_timeout, bound_box, static_cast<bool>(*synchronised));
}

void mui_announce_send_span_3t_box_f(mui_uniface_3t *uniface, double *box_1_1, double *box_1_2, double *box_1_3, double *box_2_1, double *box_2_2, double *box_2_3, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::mui_f_wrapper_3D::point_type point_1(static_cast<mui::mui_f_wrapper_3D::REAL>(*box_1_1,*box_1_2,*box_1_3));
    mui::mui_f_wrapper_3D::point_type point_2(static_cast<mui::mui_f_wrapper_3D::REAL>(*box_2_1,*box_2_2,*box_2_3));
    mui::geometry::box<mui::mui_f_wrapper_3D> bound_box(point_1, point_2);
    uniface->announce_send_span(static_cast<mui::mui_f_wrapper_3D::time_type>(*t_start),
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_timeout), bound_box, static_cast<bool>(*synchronised));
}

// Send span announce using 3d sphere geometry
void mui_announce_send_span_3f_sphere_f(mui_uniface_3f *uniface, float *centre_1, float *centre_2, float *centre_3, float *radius, float *t_start,
        float *t_timeout, int *synchronised) {
    mui::geometry::sphere3f bound_sphere(mui::point3f(*centre_1,*centre_2,*centre_3), *radius);
    uniface->announce_send_span(*t_start, *t_timeout, bound_sphere, static_cast<bool>(*synchronised));
}

void mui_announce_send_span_3fx_sphere_f(mui_uniface_3fx *uniface, float *centre_1, float *centre_2, float *centre_3, float *radius, float *t_start,
        float *t_timeout, int *synchronised) {
    mui::geometry::sphere3fx bound_sphere(mui::point3fx(*centre_1,*centre_2,*centre_3), *radius);
    uniface->announce_send_span(*t_start, *t_timeout, bound_sphere, static_cast<bool>(*synchronised));
}

void mui_announce_send_span_3d_sphere_f(mui_uniface_3d *uniface, double *centre_1, double *centre_2, double *centre_3, double *radius, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::geometry::sphere3d bound_sphere(mui::point3d(*centre_1,*centre_2,*centre_3), *radius);
    uniface->announce_send_span(*t_start, *t_timeout, bound_sphere, static_cast<bool>(*synchronised));
}

void mui_announce_send_span_3dx_sphere_f(mui_uniface_3dx *uniface, double *centre_1, double *centre_2, double *centre_3, double *radius, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::geometry::sphere3dx bound_sphere(mui::point3dx(*centre_1,*centre_2,*centre_3), *radius);
    uniface->announce_send_span(*t_start, *t_timeout, bound_sphere, static_cast<bool>(*synchronised));
}

void mui_announce_send_span_3t_sphere_f(mui_uniface_3t *uniface, double *centre_1, double *centre_2, double *centre_3, double *radius, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::geometry::sphere<mui::mui_f_wrapper_3D> bound_sphere(
            mui::mui_f_wrapper_3D::point_type(static_cast<mui::mui_f_wrapper_3D::REAL>(*centre_1,*centre_2,*centre_3)),
            static_cast<mui::mui_f_wrapper_3D::REAL>(*radius));
    uniface->announce_send_span(static_cast<mui::mui_f_wrapper_3D::time_type>(*t_start),
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_timeout), bound_sphere, static_cast<bool>(*synchronised));
}

// Receive span announce using 3d box geometry
void mui_announce_recv_span_3f_box_f(mui_uniface_3f *uniface, float *box_1_1, float *box_1_2, float *box_1_3, float *box_2_1, float *box_2_2, float *box_2_3, float *t_start,
        float *t_timeout, int *synchronised) {
    mui::point3f point_1(*box_1_1,*box_1_2,*box_1_3);
    mui::point3f point_2(*box_2_1,*box_2_2,*box_2_3);
    mui::geometry::box3f bound_box(point_1, point_2);
    uniface->announce_recv_span(*t_start, *t_timeout, bound_box, static_cast<bool>(*synchronised));
}

void mui_announce_recv_span_3fx_box_f(mui_uniface_3fx *uniface, float *box_1_1, float *box_1_2, float *box_1_3, float *box_2_1, float *box_2_2, float *box_2_3, float *t_start,
        float *t_timeout, int *synchronised) {
    mui::point3fx point_1(*box_1_1,*box_1_2,*box_1_3);
    mui::point3fx point_2(*box_2_1,*box_2_2,*box_2_3);
    mui::geometry::box3fx bound_box(point_1, point_2);
    uniface->announce_recv_span(*t_start, *t_timeout, bound_box, static_cast<bool>(*synchronised));
}

void mui_announce_recv_span_3d_box_f(mui_uniface_3d *uniface, double *box_1_1, double *box_1_2, double *box_1_3, double *box_2_1, double *box_2_2, double *box_2_3, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::point3d point_1(*box_1_1,*box_1_2,*box_1_3);
    mui::point3d point_2(*box_2_1,*box_2_2,*box_2_3);
    mui::geometry::box3d bound_box(point_1, point_2);
    uniface->announce_recv_span(*t_start, *t_timeout, bound_box, static_cast<bool>(*synchronised));
}

void mui_announce_recv_span_3dx_box_f(mui_uniface_3dx *uniface, double *box_1_1, double *box_1_2, double *box_1_3, double *box_2_1, double *box_2_2, double *box_2_3, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::point3dx point_1(*box_1_1,*box_1_2,*box_1_3);
    mui::point3dx point_2(*box_2_1,*box_2_2,*box_2_3);
    mui::geometry::box3dx bound_box(point_1, point_2);
    uniface->announce_recv_span(*t_start, *t_timeout, bound_box, static_cast<bool>(*synchronised));
}

void mui_announce_recv_span_3t_box_f(mui_uniface_3t *uniface, double *box_1_1, double *box_1_2, double *box_1_3, double *box_2_1, double *box_2_2, double *box_2_3, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::mui_f_wrapper_3D::point_type point_1(static_cast<mui::mui_f_wrapper_3D::REAL>(*box_1_1,*box_1_2,*box_1_3));
    mui::mui_f_wrapper_3D::point_type point_2(static_cast<mui::mui_f_wrapper_3D::REAL>(*box_2_1,*box_2_2,*box_2_3));
    mui::geometry::box<mui::mui_f_wrapper_3D> bound_box(point_1, point_2);
    uniface->announce_recv_span(static_cast<mui::mui_f_wrapper_3D::time_type>(*t_start),
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_timeout), bound_box, static_cast<bool>(*synchronised));
}

// Receive span announce using 3d sphere geometry
void mui_announce_recv_span_3f_sphere_f(mui_uniface_3f *uniface, float *centre_1, float *centre_2, float *centre_3, float *radius, float *t_start,
        float *t_timeout, int *synchronised) {
    mui::geometry::sphere3f bound_sphere(mui::point3f(*centre_1,*centre_2,*centre_3), *radius);
    uniface->announce_recv_span(*t_start, *t_timeout, bound_sphere, static_cast<bool>(*synchronised));
}

void mui_announce_recv_span_3fx_sphere_f(mui_uniface_3fx *uniface, float *centre_1, float *centre_2, float *centre_3, float *radius, float *t_start,
        float *t_timeout, int *synchronised) {
    mui::geometry::sphere3fx bound_sphere(mui::point3fx(*centre_1,*centre_2,*centre_3), *radius);
    uniface->announce_recv_span(*t_start, *t_timeout, bound_sphere, static_cast<bool>(*synchronised));
}

void mui_announce_recv_span_3d_sphere_f(mui_uniface_3d *uniface, double *centre_1, double *centre_2, double *centre_3, double *radius, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::geometry::sphere3d bound_sphere(mui::point3d(*centre_1,*centre_2,*centre_3), *radius);
    uniface->announce_recv_span(*t_start, *t_timeout, bound_sphere);
}

void mui_announce_recv_span_3dx_sphere_f(mui_uniface_3dx *uniface, double *centre_1, double *centre_2, double *centre_3, double *radius, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::geometry::sphere3dx bound_sphere(mui::point3dx(*centre_1,*centre_2,*centre_3), *radius);
    uniface->announce_recv_span(*t_start, *t_timeout, bound_sphere, static_cast<bool>(*synchronised));
}

void mui_announce_recv_span_3t_sphere_f(mui_uniface_3t *uniface, double *centre_1, double *centre_2, double *centre_3, double *radius, double *t_start,
        double *t_timeout, int *synchronised) {
    mui::geometry::sphere<mui::mui_f_wrapper_3D> bound_sphere(
            mui::mui_f_wrapper_3D::point_type(static_cast<mui::mui_f_wrapper_3D::REAL>(*centre_1,*centre_2,*centre_3)),
            static_cast<mui::mui_f_wrapper_3D::REAL>(*radius));
    uniface->announce_recv_span(static_cast<mui::mui_f_wrapper_3D::time_type>(*t_start),
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_timeout), bound_sphere, static_cast<bool>(*synchronised));
}

// Send disable announce (local call per MPI rank)
void mui_announce_send_disable_3f_f(mui_uniface_3f *uniface, int *synchronised) {
    uniface->announce_send_disable(static_cast<bool>(*synchronised));
}

void mui_announce_send_disable_3fx_f(mui_uniface_3fx *uniface, int *synchronised) {
    uniface->announce_send_disable(static_cast<bool>(*synchronised));
}

void mui_announce_send_disable_3d_f(mui_uniface_3d *uniface, int *synchronised) {
    uniface->announce_send_disable(static_cast<bool>(*synchronised));
}

void mui_announce_send_disable_3dx_f(mui_uniface_3dx *uniface, int *synchronised) {
    uniface->announce_send_disable(static_cast<bool>(*synchronised));
}

void mui_announce_send_disable_3t_f(mui_uniface_3t *uniface, int *synchronised) {
    uniface->announce_send_disable(static_cast<bool>(*synchronised));
}

// Receive disable announce (local call per MPI rank)
void mui_announce_recv_disable_3f_f(mui_uniface_3f *uniface, int *synchronised) {
    uniface->announce_recv_disable(static_cast<bool>(*synchronised));
}

void mui_announce_recv_disable_3fx_f(mui_uniface_3fx *uniface, int *synchronised) {
    uniface->announce_recv_disable(static_cast<bool>(*synchronised));
}

void mui_announce_recv_disable_3d_f(mui_uniface_3d *uniface, int *synchronised) {
    uniface->announce_recv_disable(static_cast<bool>(*synchronised));
}

void mui_announce_recv_disable_3dx_f(mui_uniface_3dx *uniface, int *synchronised) {
    uniface->announce_recv_disable(static_cast<bool>(*synchronised));
}

void mui_announce_recv_disable_3t_f(mui_uniface_3t *uniface, int *synchronised) {
    uniface->announce_recv_disable(static_cast<bool>(*synchronised));
}

/******************************************
 * MUI barrier functions                   *
 ******************************************/

// Barrier at single time value
void mui_barrier_3f_f(mui_uniface_3f *uniface, float *t) {
    uniface->barrier(*t);
}

void mui_barrier_3fx_f(mui_uniface_3fx *uniface, float *t) {
    uniface->barrier(*t);
}

void mui_barrier_3d_f(mui_uniface_3d *uniface, double *t) {
    uniface->barrier(*t);
}

void mui_barrier_3dx_f(mui_uniface_3dx *uniface, double *t) {
    uniface->barrier(*t);
}

void mui_barrier_3t_f(mui_uniface_3t *uniface, double *t) {
    uniface->barrier(static_cast<mui::mui_f_wrapper_3D::time_type>(*t));
}

// Barrier at two time values
void mui_barrier_3f_pair_f(mui_uniface_3f *uniface, float *t_1, float *t_2) {
    uniface->barrier(*t_1, *t_2);
}

void mui_barrier_3fx_pair_f(mui_uniface_3fx *uniface, float *t_1, float *t_2) {
    uniface->barrier(*t_1, *t_2);
}

void mui_barrier_3d_pair_f(mui_uniface_3d *uniface, double *t_1, double *t_2) {
    uniface->barrier(*t_1, *t_2);
}

void mui_barrier_3dx_pair_f(mui_uniface_3dx *uniface, double *t_1, double *t_2) {
    uniface->barrier(*t_1, *t_2);
}

void mui_barrier_3t_pair_f(mui_uniface_3t *uniface, double *t_1, double *t_2) {
    uniface->barrier(static_cast<mui::mui_f_wrapper_3D::time_type>(*t_1),
            static_cast<mui::mui_f_wrapper_3D::time_type>(*t_2));
}

/******************************************
 * MUI forget functions                    *
 ******************************************/

// Forget log between [-inf, upper]
void mui_forget_upper_3f_f(mui_uniface_3f *uniface, float *upper, int *reset_log) {
    uniface->forget(*upper, static_cast<bool>(*reset_log));
}

void mui_forget_upper_3fx_f(mui_uniface_3fx *uniface, float *upper, int *reset_log) {
    uniface->forget(*upper, static_cast<bool>(*reset_log));
}

void mui_forget_upper_3d_f(mui_uniface_3d *uniface, double *upper, int *reset_log) {
    uniface->forget(*upper, static_cast<bool>(*reset_log));
}

void mui_forget_upper_3dx_f(mui_uniface_3dx *uniface, double *upper, int *reset_log) {
    uniface->forget(*upper, static_cast<bool>(*reset_log));
}

void mui_forget_upper_3t_f(mui_uniface_3t *uniface, double *upper, int *reset_log) {
    uniface->forget(static_cast<mui::mui_f_wrapper_3D::time_type>(*upper), static_cast<bool>(*reset_log));
}

// Forget log between [-inf, -inf], [upper_1, upper_2]
void mui_forget_upper_3f_pair_f(mui_uniface_3f *uniface, float *upper_1, float *upper_2, int *reset_log) {
    std::pair<float, float> forget_time(*upper_1, *upper_2);
    uniface->forget(forget_time, static_cast<bool>(*reset_log));
}

void mui_forget_upper_3fx_pair_f(mui_uniface_3fx *uniface, float *upper_1, float *upper_2, int *reset_log) {
    std::pair<float, float> forget_time(*upper_1, *upper_2);
    uniface->forget(forget_time, static_cast<bool>(*reset_log));
}

void mui_forget_upper_3d_pair_f(mui_uniface_3d *uniface, double *upper_1, double *upper_2, int *reset_log) {
    std::pair<double, double> forget_time(*upper_1, *upper_2);
    uniface->forget(forget_time, static_cast<bool>(*reset_log));
}

void mui_forget_upper_3dx_pair_f(mui_uniface_3dx *uniface, double *upper_1, double *upper_2, int *reset_log) {
    std::pair<double, double> forget_time(*upper_1, *upper_2);
    uniface->forget(forget_time, static_cast<bool>(*reset_log));
}

void mui_forget_upper_3t_pair_f(mui_uniface_3t *uniface, double *upper_1, double *upper_2, int *reset_log) {
    mui::mui_f_wrapper_3D::time_type time_1 = static_cast<mui::mui_f_wrapper_3D::time_type>(*upper_1);
    mui::mui_f_wrapper_3D::time_type time_2 = static_cast<mui::mui_f_wrapper_3D::time_type>(*upper_2);
    std::pair<mui::mui_f_wrapper_3D::time_type, mui::mui_f_wrapper_3D::time_type> forget_time(time_1, time_2);
    uniface->forget(forget_time, static_cast<bool>(*reset_log));
}

// Forget log between [lower, upper]
void mui_forget_lower_upper_3f_f(mui_uniface_3f *uniface, float *lower, float *upper, int *reset_log) {
    uniface->forget(*lower, *upper, static_cast<bool>(*reset_log));
}

void mui_forget_lower_upper_3fx_f(mui_uniface_3fx *uniface, float *lower, float *upper, int *reset_log) {
    uniface->forget(*lower, *upper, static_cast<bool>(*reset_log));
}

void mui_forget_lower_upper_3d_f(mui_uniface_3d *uniface, double *lower, double *upper, int *reset_log) {
    uniface->forget(*lower, *upper, static_cast<bool>(*reset_log));
}

void mui_forget_lower_upper_3dx_f(mui_uniface_3dx *uniface, double *lower, double *upper, int *reset_log) {
    uniface->forget(*lower, *upper, static_cast<bool>(*reset_log));
}

void mui_forget_lower_upper_3t_f(mui_uniface_3t *uniface, double *lower, double *upper, int *reset_log) {
    mui::mui_f_wrapper_3D::time_type forget_time_1 = static_cast<mui::mui_f_wrapper_3D::time_type>(*lower);
    mui::mui_f_wrapper_3D::time_type forget_time_2 = static_cast<mui::mui_f_wrapper_3D::time_type>(*upper);
    uniface->forget(forget_time_1, forget_time_2, static_cast<bool>(*reset_log));
}

// Forget log between [lower_1, lower_2], [upper_1, upper_2]
void mui_forget_lower_upper_3f_pair_f(mui_uniface_3f *uniface, float *lower_1, float *lower_2, float *upper_1, float *upper_2,
        int *reset_log) {
    std::pair<float, float> forget_time_1(*lower_1, *lower_2);
    std::pair<float, float> forget_time_2(*upper_1, *upper_2);
    uniface->forget(forget_time_1, forget_time_2, static_cast<bool>(*reset_log));
}

void mui_forget_lower_upper_3fx_pair_f(mui_uniface_3fx *uniface, float *lower_1, float *lower_2, float *upper_1,
        float *upper_2, int *reset_log) {
    std::pair<float, float> forget_time_1(*lower_1, *lower_2);
    std::pair<float, float> forget_time_2(*upper_1, *upper_2);
    uniface->forget(forget_time_1, forget_time_2, static_cast<bool>(*reset_log));
}

void mui_forget_lower_upper_3d_pair_f(mui_uniface_3d *uniface, double *lower_1, double *lower_2, double *upper_1,
        double *upper_2, int *reset_log) {
    std::pair<double, double> forget_time_1(*lower_1, *lower_2);
    std::pair<double, double> forget_time_2(*upper_1, *upper_2);
    uniface->forget(forget_time_1, forget_time_2, static_cast<bool>(*reset_log));
}

void mui_forget_lower_upper_3dx_pair_f(mui_uniface_3dx *uniface, double *lower_1, double *lower_2, double *upper_1,
        double *upper_2, int *reset_log) {
    std::pair<double, double> forget_time_1(*lower_1, *lower_2);
    std::pair<double, double> forget_time_2(*upper_1, *upper_2);
    uniface->forget(forget_time_1, forget_time_2, static_cast<bool>(*reset_log));
}

void mui_forget_lower_upper_3t_pair_f(mui_uniface_3t *uniface, double *lower_1, double *lower_2, double *upper_1,
        double *upper_2, int *reset_log) {
    mui::mui_f_wrapper_3D::time_type time_1 = static_cast<mui::mui_f_wrapper_3D::time_type>(*lower_1);
    mui::mui_f_wrapper_3D::time_type time_2 = static_cast<mui::mui_f_wrapper_3D::time_type>(*lower_2);
    mui::mui_f_wrapper_3D::time_type time_3 = static_cast<mui::mui_f_wrapper_3D::time_type>(*upper_1);
    mui::mui_f_wrapper_3D::time_type time_4 = static_cast<mui::mui_f_wrapper_3D::time_type>(*upper_2);
    std::pair<mui::mui_f_wrapper_3D::time_type, mui::mui_f_wrapper_3D::time_type> forget_time_1(time_1, time_2);
    std::pair<mui::mui_f_wrapper_3D::time_type, mui::mui_f_wrapper_3D::time_type> forget_time_2(time_3, time_4);
    uniface->forget(forget_time_1, forget_time_2, static_cast<bool>(*reset_log));
}

// Set to forget log between [-inf, current-length] automatically
void mui_set_forget_length_3f_f(mui_uniface_3f *uniface, float *length) {
    uniface->set_memory(*length);
}

void mui_set_forget_length_3fx_f(mui_uniface_3fx *uniface, float *length) {
    uniface->set_memory(*length);
}

void mui_set_forget_length_3d_f(mui_uniface_3d *uniface, double *length) {
    uniface->set_memory(*length);
}

void mui_set_forget_length_3dx_f(mui_uniface_3dx *uniface, double *length) {
    uniface->set_memory(*length);
}

void mui_set_forget_length_3t_f(mui_uniface_3t *uniface, double *length) {
    uniface->set_memory(*length);
}

/******************************************
 * MUI URI functions                      *
 ******************************************/

// Obtain original URI host value from existing interface
void mui_uri_host_3f_f(mui_uniface_3f *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_host().c_str());
}

void mui_uri_host_3fx_f(mui_uniface_3fx *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_host().c_str());
}

void mui_uri_host_3d_f(mui_uniface_3d *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_host().c_str());
}

void mui_uri_host_3dx_f(mui_uniface_3dx *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_host().c_str());
}

void mui_uri_host_3t_f(mui_uniface_3t *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_host().c_str());
}

// Obtain original URI path value from existing interface
void mui_uri_path_3f_f(mui_uniface_3f *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_path().c_str());
}

void mui_uri_path_3fx_f(mui_uniface_3fx *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_path().c_str());
}

void mui_uri_path_3d_f(mui_uniface_3d *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_path().c_str());
}

void mui_uri_path_3dx_f(mui_uniface_3dx *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_path().c_str());
}

void mui_uri_path_3t_f(mui_uniface_3t *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_path().c_str());
}

// Obtain original URI protocol value from existing interface
void mui_uri_protocol_3f_f(mui_uniface_3f *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_protocol().c_str());
}

void mui_uri_protocol_3fx_f(mui_uniface_3fx *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_protocol().c_str());
}

void mui_uri_protocol_3d_f(mui_uniface_3d *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_protocol().c_str());
}

void mui_uri_protocol_3dx_f(mui_uniface_3dx *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_protocol().c_str());
}

void mui_uri_protocol_3t_f(mui_uniface_3t *uniface, char **return_val) {
    *return_val = const_cast<char*>(uniface->uri_protocol().c_str());
}

}
