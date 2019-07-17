/*
Multiscale Universal Interface Code Coupling Library

Copyright (C) 2017 Y. H. Tang, S. Kudo, X. Bian, Z. Li, G. E. Karniadakis

This software is jointly licensed under the Apache License, Version 2.0
and the GNU General Public License version 3, you may use it according
to either.

** Apache License, version 2.0 **

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

** GNU General Public License, version 3 **

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

** File Details **

Filename: mui_3df.cpp
Created: Jan 21, 2015
Author: S. Kudo
Description: Fortran wrapper to create 3D MUI uniface.
*/

#include "../../mui.h"

extern "C" {

using namespace mui;

typedef uniface3d                        mui_uniface3d;
typedef sampler_gauss3d<double>          mui_sampler_gauss3d;
typedef sampler_moving_average3d<double> mui_sampler_moving_average3d;
typedef chrono_sampler_exact3d           mui_chrono_sampler_exact3d;
typedef chrono_sampler_mean3d            mui_chrono_sampler_mean3d;

// allocator
void mui_create_uniface3d_f( mui_uniface3d** ret, const char *URI ) {
	*ret = new mui_uniface3d( URI );
}

void mui_create_sampler_3d_f_( mui_sampler_gauss3d** ret, double* r, double* h ) {
	*ret = new mui_sampler_gauss3d( *r, *h );
}

void mui_create_sampler_moving_average3d_f_( mui_sampler_moving_average3d** ret, double* dx, double* dy, double* dz ) {
	*ret = new mui_sampler_moving_average3d( point3d(*dx,*dy,*dz) );
}

void mui_create_chrono_sampler_exact3d_f_(mui_chrono_sampler_exact3d** ret) {
	*ret = new mui_chrono_sampler_exact3d;
}

void mui_create_chrono_sampler_mean3d_f_( mui_chrono_sampler_mean3d** ret, double* past, double* future ) {
	*ret = new mui_chrono_sampler_mean3d( *past, *future );
}

// deallocator
void mui_destroy_uniface3d_f( mui_uniface3d *uniface ) {
	delete uniface;
}

void mui_destroy_sampler_3d_f_( mui_sampler_gauss3d* sampler ) {
	delete sampler;
}

void mui_destroy_sampler_moving_average3d_f_( mui_sampler_moving_average3d* sampler ) {
	delete sampler;
}

void mui_destroy_chrono_sampler_exact3d_f_( mui_chrono_sampler_exact3d* sampler ) {
	delete sampler;
}

void mui_destroy_chrono_sampler_mean3d_f_( mui_chrono_sampler_mean3d* sampler ) {
	delete sampler;
}

// push
void mui_push_f( mui_uniface3d* uniface, const char *attr, double* x, double* y, double* z, double* value ) {
	uniface->push( std::string(attr), point3d(*x,*y,*z), *value );
}

// spatial sampler: gaussian
// temporal sampler: exact point
double mui_fetch_gaussian_exact_f_( mui_uniface3d* uniface, const char *attr, double* x, double* y, double* z, double* t, mui_sampler_gauss3d *spatial, mui_chrono_sampler_exact3d *temporal ) {
	return uniface->fetch( std::string(attr), point3d(*x,*y,*z), *t, *spatial, *temporal );
}

// spatial sampler: gaussian
// temporal sampler: mean
double mui_fetch_gaussian_mean_f_( mui_uniface3d* uniface, const char *attr, double* x, double* y, double* z, double* t, mui_sampler_gauss3d *spatial, mui_chrono_sampler_mean3d *temporal ) {
	return uniface->fetch( std::string(attr), point3d(*x,*y,*z), *t, *spatial, *temporal );
}

// spatial sampler: moving average
// temporal sampler: exact point
double mui_fetch_moving_average_exact_f_( mui_uniface3d* uniface, const char *attr, double* x, double* y, double* z, double* t, mui_sampler_moving_average3d *spatial, mui_chrono_sampler_exact3d *temporal ) {
	return uniface->fetch( std::string(attr), point3d(*x,*y,*z), *t, *spatial, *temporal );
}

// spatial sampler: moving_average
// temporal sampler: mean
double mui_fetch_moving_average_mean_f_( mui_uniface3d* uniface, const char *attr, double* x, double* y, double* z, double* t, mui_sampler_moving_average3d *spatial, mui_chrono_sampler_mean3d *temporal ) {
	return uniface->fetch( std::string(attr), point3d(*x,*y,*z), *t, *spatial, *temporal );
}

// commit all data in buffer
void mui_commit_f_( mui_uniface3d* uniface, double* t ) {
	uniface->commit( *t );
}

// wait for peers
void mui_barrier_f_( mui_uniface3d* uniface, double* t ) {
	uniface->barrier( *t );
}

// remove obsolete data
void mui_forget_f_( mui_uniface3d* uniface, double* first, double* last ) {
	uniface->forget( *first, *last );
}

// set automatic deletion
void mui_set_memory_f_( mui_uniface3d* uniface, double* length ) {
	return uniface->set_memory( *length );
}

}
