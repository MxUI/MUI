/*
 * mui_3d.cpp
 *
 *  Created on: Jan 20, 2015
 *      Author: ytang
 */
#include "../mui.h"

extern "C" {

using namespace mui;

typedef uniface3d                        mui_uniface3d;
typedef sampler_gauss3d<double>          mui_sampler_gauss3d;
typedef sampler_moving_average3d<double> mui_sampler_moving_average3d;
typedef chrono_sampler_exact3d           mui_chrono_sampler_exact3d;
typedef chrono_sampler_mean3d            mui_chrono_sampler_mean3d;

// allocator
mui_uniface3d* mui_create_uniface3d( const char *URI ) {
	return new mui_uniface3d( URI );
}

mui_sampler_gauss3d* mui_create_sampler_3d( double r, double h ) {
	return new mui_sampler_gauss3d( r, h );
}

mui_sampler_moving_average3d* mui_create_sampler_moving_average3d( double dx, double dy, double dz ) {
	return new mui_sampler_moving_average3d( point3d(dx,dy,dz) );
}

mui_chrono_sampler_exact3d* mui_create_chrono_sampler_exact3d() {
	return new mui_chrono_sampler_exact3d;
}

mui_chrono_sampler_mean3d* mui_create_chrono_sampler_mean3d( double past, double future ) {
	return new mui_chrono_sampler_mean3d( past, future );
}

// deallocator
void mui_destroy_uniface3d( mui_uniface3d *uniface ) {
	delete uniface;
}

void mui_destroy_sampler_3d( mui_sampler_gauss3d* sampler ) {
	delete sampler;
}

void mui_destroy_sampler_moving_average3d( mui_sampler_moving_average3d* sampler ) {
	delete sampler;
}

void mui_destroy_chrono_sampler_exact3d( mui_chrono_sampler_exact3d* sampler ) {
	delete sampler;
}

void mui_destroy_chrono_sampler_mean3d( mui_chrono_sampler_mean3d* sampler ) {
	delete sampler;
}

// push
void mui_push( mui_uniface3d* uniface, const char *attr, double x, double y, double z, double value ) {
	uniface->push( std::string(attr), point3d(x,y,z), value );
}

// spatial sampler: gaussian
// temporal sampler: exact point
double mui_fetch_gaussian_exact( mui_uniface3d* uniface, const char *attr, double x, double y, double z, double t, mui_sampler_gauss3d *spatial, mui_chrono_sampler_exact3d *temporal ) {
	return uniface->fetch( std::string(attr), point3d(x,y,z), t, *spatial, *temporal );
}

// spatial sampler: gaussian
// temporal sampler: mean
double mui_fetch_gaussian_mean( mui_uniface3d* uniface, const char *attr, double x, double y, double z, double t, mui_sampler_gauss3d *spatial, mui_chrono_sampler_mean3d *temporal ) {
	return uniface->fetch( std::string(attr), point3d(x,y,z), t, *spatial, *temporal );
}

// spatial sampler: moving average
// temporal sampler: exact point
double mui_fetch_moving_average_exact( mui_uniface3d* uniface, const char *attr, double x, double y, double z, double t, mui_sampler_moving_average3d *spatial, mui_chrono_sampler_exact3d *temporal ) {
	return uniface->fetch( std::string(attr), point3d(x,y,z), t, *spatial, *temporal );
}

// spatial sampler: moving_average
// temporal sampler: mean
double mui_fetch_moving_average_mean( mui_uniface3d* uniface, const char *attr, double x, double y, double z, double t, mui_sampler_moving_average3d *spatial, mui_chrono_sampler_mean3d *temporal ) {
	return uniface->fetch( std::string(attr), point3d(x,y,z), t, *spatial, *temporal );
}

// commit all data in buffer
void mui_commit( mui_uniface3d* uniface, double t ) {
	uniface->commit( t );
}

// wait for peers
void mui_barrier( mui_uniface3d* uniface, double t ) {
	uniface->barrier( t );
}

// remove obsolete data
void mui_forget( mui_uniface3d* uniface, double first, double last ) {
	uniface->forget( first, last );
}

// set automatic deletion
void mui_set_memory( mui_uniface3d* uniface, double length ) {
	return uniface->set_memory( length );
}

}
