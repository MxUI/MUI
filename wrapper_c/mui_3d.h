/*
 * mui_3d.h
 *
 *  Created on: Jan 19, 2015
 *      Author: ytang
 */

/*
This is a C wrapper of MUI using mui::config3d
Only minimal functionality is covered here
*/
#ifndef MUI_3D_H_
#define MUI_3D_H_

typedef struct mui_uniface3d                mui_uniface3d;
typedef struct mui_sampler_gauss3d          mui_sampler_gauss3d;
typedef struct mui_sampler_moving_average3d mui_sampler_moving_average3d;
typedef struct mui_chrono_sampler_exact3d   mui_chrono_sampler_exact3d;
typedef struct mui_chrono_sampler_mean3d    mui_chrono_sampler_mean3d;

/* allocator */
mui_uniface3d* mui_create_uniface3d( const char *URI );
mui_sampler_gauss3d* mui_create_sampler_3d( double r, double h );
mui_sampler_moving_average3d* mui_create_sampler_moving_average3d( double dx, double dy, double dz );
mui_chrono_sampler_exact3d* mui_create_chrono_sampler_exact3d();
mui_chrono_sampler_mean3d* mui_create_chrono_sampler_mean3d( double past, double future );

/* deallocator */
void mui_destroy_uniface3d( mui_uniface3d *uniface );
void mui_destroy_sampler_3d( mui_sampler_gauss3d* sampler );
void mui_destroy_sampler_moving_average3d( mui_sampler_moving_average3d* sampler );
void mui_destroy_chrono_sampler_exact3d( mui_chrono_sampler_exact3d* sampler );
void mui_destroy_chrono_sampler_mean3d( mui_chrono_sampler_mean3d* sampler );

/* push */
void mui_push( mui_uniface3d* uniface, const char *attr, double x, double y, double z, double t );

/*  spatial sampler: gaussian */
/*  temporal sampler: exact point */
double mui_fetch_gaussian_exact( mui_uniface3d* uniface, const char *attr, double x, double y, double z, double t, mui_sampler_gauss3d *spatial, mui_chrono_sampler_exact3d *temporal );

/*  spatial sampler: gaussian */
/*  temporal sampler: mean */
double mui_fetch_gaussian_mean( mui_uniface3d* uniface, const char *attr, double x, double y, double z, double t, mui_sampler_gauss3d *spatial, mui_chrono_sampler_mean3d *temporal );

/*  spatial sampler: moving average */
/*  temporal sampler: exact point */
double mui_fetch_moving_average_exact( mui_uniface3d* uniface, const char *attr, double x, double y, double z, double t, mui_sampler_moving_average3d *spatial, mui_chrono_sampler_exact3d *temporal );

/*  spatial sampler: moving_average */
/*  temporal sampler: mean */
double mui_fetch_moving_average_mean( mui_uniface3d* uniface, const char *attr, double x, double y, double z, double t, mui_sampler_moving_average3d *spatial, mui_chrono_sampler_mean3d *temporal );

/*  commit all data in buffer */
void mui_commit( mui_uniface3d*, double t );

/*  wait for peers */
void mui_barrier( mui_uniface3d*, double t );

/*  remove obsolete data */
void mui_forget( mui_uniface3d*, double first, double last );

/*  set automatic deletion */
void mui_set_memory( mui_uniface3d*, double length );

#endif /* MUI_3D_H_ */
