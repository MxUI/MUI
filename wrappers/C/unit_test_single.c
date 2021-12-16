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
 * @file unit_test_single.c
 * @author S. M. Longshaw (derived from original example by Y. H. Tang)
 * @date 20 October 2021
 * @brief Unit test for MUI C wrapper - creates 1D/2D/3D interface examples
 */

// Standard C includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// Include general non-dimensional MUI functions
#include "mui_c_wrapper_general.h"
// Include 1D MUI functions
#include "mui_c_wrapper_1d.h"
// Include 2D MUI functions
#include "mui_c_wrapper_2d.h"
// Include 3D MUI functions
#include "mui_c_wrapper_3d.h"

int main(int argc, char **argv) {

	// No need to call mui_mpi_split_by_app() function first as with the _multi example
	// as direct uniface creation also called MPI_Init()

	/***********************
	 * Single 1D interface *
	 ***********************/

	// Create URI in the format mpi://[domain]_1d/[interface_name]
	char *uri1d = (char*) malloc(strlen(argv[1]) + strlen(argv[2]) + 11);
	strcpy(uri1d, "mpi://");
	strcat(uri1d, argv[1]);
	strcat(uri1d, "_1d/");
	strcat(uri1d, argv[2]);
	strcat(uri1d, "\0");

	// Create MUI interface
	mui_uniface_1d *uniface1d = mui_create_uniface_1d((const char*) uri1d);

	// Define point for push
	mui_point_1d push_point1d = { 0 };

	// Push value=1 through the MUI interface with the tag "position" at location=push_point
	mui_push_1d(uniface1d, "position", push_point1d, 1);

	// Commit (transmit) the pushed value at time=0
	mui_commit_1d(uniface1d, 0);

	// Create spatial and temporal samplers for fetch operation
	mui_sampler_exact_1d *spatial_sampler1d = mui_create_sampler_exact_1d(1e-37);
	mui_chrono_sampler_exact_1d *temporal_sampler1d = mui_create_chrono_sampler_exact_1d(1e-37);

	// Define point for fetch
	mui_point_1d fetch_point1d = { 0 };

	// Fetch the value for tag "position" at fetch_point at time=0
	double result1d = mui_fetch_exact_exact_1d(uniface1d, "position", fetch_point1d, 0, spatial_sampler1d,
			temporal_sampler1d);

	printf("Fetched 1D interface value = %f \n", result1d);

	// Destroy created 1D MUI objects
	mui_destroy_sampler_exact_1d(spatial_sampler1d);
	mui_destroy_chrono_sampler_exact_1d(temporal_sampler1d);
	free(uri1d);

	/**********************
	 * 2D interface       *
	 *********************/

	// Create URI in the format mpi://[domain]_2d/interface_name
	char *uri2d = (char*) malloc(strlen(argv[1]) + strlen(argv[2]) + 11);
	strcpy(uri2d, "mpi://");
	strcat(uri2d, argv[1]);
	strcat(uri2d, "_2d/");
	strcat(uri2d, argv[2]);
	strcat(uri2d, "\0");

	// Create MUI interface
	mui_uniface_2d *uniface2d = mui_create_uniface_2d((const char*) uri2d);

	// Define point for push
	mui_point_2d push_point2d = { 0, 0 };

	// Push value=1 through the MUI interface with the tag "position" at location=push_point
	mui_push_2d(uniface2d, "position", push_point2d, 1);

	// Commit (transmit) the pushed value at time=0
	mui_commit_2d(uniface2d, 0);

	// Create spatial and temporal samplers for fetch operation
	mui_sampler_exact_2d *spatial_sampler2d = mui_create_sampler_exact_2d(1e-37);
	mui_chrono_sampler_exact_2d *temporal_sampler2d = mui_create_chrono_sampler_exact_2d(1e-37);

	// Define point for fetch
	mui_point_2d fetch_point2d = { 0, 0 };

	// Fetch the value for tag "position" at fetch_point at time=0
	double result2d = mui_fetch_exact_exact_2d(uniface2d, "position", fetch_point2d, 0, spatial_sampler2d,
			temporal_sampler2d);

	printf("Fetched 2D interface value = %f \n", result2d);

	// Destroy created 2D MUI objects
	mui_destroy_sampler_exact_2d(spatial_sampler2d);
	mui_destroy_chrono_sampler_exact_2d(temporal_sampler2d);
	free(uri2d);

	/*********************
	 * 3D interface       *
	 *********************/

	// Create URI in the format mpi://[domain]_3d/interface_name
	char *uri3d = (char*) malloc(strlen(argv[1]) + strlen(argv[2]) + 11);
	strcpy(uri3d, "mpi://");
	strcat(uri3d, argv[1]);
	strcat(uri3d, "_3d/");
	strcat(uri3d, argv[2]);
	strcat(uri3d, "\0");

	// Create MUI interface
	mui_uniface_3d *uniface3d = mui_create_uniface_3d((const char*) uri3d);

	// Define point for push
	mui_point_3d push_point3d = { 0, 0, 0 };

	// Push value=1 through the MUI interface with the tag "position" at location=push_point
	mui_push_3d(uniface3d, "position", push_point3d, 1);

	// Commit (transmit) the pushed value at time=0
	mui_commit_3d(uniface3d, 0);

	// Create spatial and temporal samplers for fetch operation
	mui_sampler_exact_3d *spatial_sampler3d = mui_create_sampler_exact_3d(1e-37);
	mui_chrono_sampler_exact_3d *temporal_sampler3d = mui_create_chrono_sampler_exact_3d(1e-37);

	// Define point for fetch
	mui_point_3d fetch_point3d = { 0, 0, 0 };

	// Fetch the value for tag "position" at fetch_point at time=0
	double result3d = mui_fetch_exact_exact_3d(uniface3d, "position", fetch_point3d, 0, spatial_sampler3d,
			temporal_sampler3d);

	printf("Fetched 3D interface value = %f \n", result3d);

	// Destroy created 3D MUI objects
	mui_destroy_sampler_exact_3d(spatial_sampler3d);
	mui_destroy_chrono_sampler_exact_3d(temporal_sampler3d);

	// Destroy created MUI interfaces note: calls MPI_Finalize(), so need to do last
	mui_destroy_uniface_1d(uniface1d);
	mui_destroy_uniface_2d(uniface2d);
	mui_destroy_uniface_3d(uniface3d);

	return 0;
}
