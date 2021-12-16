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
 * @file unit_test_multi.c
 * @author S. M. Longshaw (derived from original example by Y. H. Tang)
 * @date 23 November 2021
 * @brief Unit test for MUI C wrapper - creates 1D/2D/3D interface examples
 * 		  using the helper function for multiple interfaces for a single domain
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
	//Initialise MPI (returns an MPI comm world we might use in a local MPI_Init() call)
	mui_mpi_split_by_app();

	// Get the number of interfaces to create
	int num_interfaces = atoi(argv[3]);

	/***********************
	 * Multi 1D interfaces *
	 ***********************/

	char *domain1d = (char*) malloc(strlen(argv[1]) + 4);
	strcpy(domain1d, argv[1]);
	strcat(domain1d, "_1d");
	strcat(domain1d, "\0");

	// Define interface names in the format [interface_name]_[1,2,3...]
	char** interfaces1d = (char**) malloc(sizeof(char*) * num_interfaces);
	for(int i=0; i<num_interfaces; i++) {
		interfaces1d[i] = (char*) malloc(strlen(argv[2] + 3));
		strcpy(interfaces1d[i], argv[2]);
		strcat(interfaces1d[i], "_");
		char tmp[11];
		sprintf(tmp, "%d", i);
		strcat(interfaces1d[i], tmp);
		strcat(interfaces1d[i], "\0");
	}

	// Define return structure to hold created unifaces
	mui_uniface_1d** unifaces1d = (mui_uniface_1d**) malloc(sizeof(mui_uniface_1d*) * num_interfaces);

	// Create MUI interfaces using helper function
	unifaces1d = mui_create_uniface_multi_1d( (const char*)domain1d, (const char**)interfaces1d, num_interfaces );

	// Free char memory after uniface creation
	for(int i=0; i<num_interfaces; i++) {
		free(interfaces1d[i]);
	}

	free(domain1d);

	// Define point for push
	mui_point_1d push_point1d = { 0 };

	// Push value=1 through the MUI interfaces with the tag "position" at location=push_point
	for(int i=0; i<num_interfaces; i++) {
		mui_push_1d(unifaces1d[i], "position", push_point1d, 1);
	}

	// Commit (transmit) the pushed value at time=0
	for(int i=0; i<num_interfaces; i++) {
		mui_commit_1d(unifaces1d[i], 0);
	}

	// Create spatial and temporal samplers for fetch operation
	mui_sampler_exact_1d *spatial_sampler1d = mui_create_sampler_exact_1d(1e-37);
	mui_chrono_sampler_exact_1d *temporal_sampler1d = mui_create_chrono_sampler_exact_1d(1e-37);

	// Define point for fetch
	mui_point_1d fetch_point1d = { 0 };

	double* results1d = (double*) malloc(sizeof(double) * num_interfaces);

	// Fetch the value for tag "position" at fetch_point at time=0
	for(int i=0; i<num_interfaces; i++) {
		results1d[i] = mui_fetch_exact_exact_1d(unifaces1d[i], "position", fetch_point1d, 0, spatial_sampler1d,
				temporal_sampler1d);

		printf("Fetched 1D interface_%d value = %f \n", i, results1d[i]);
	}

	// Destroy created 1D MUI objects
	mui_destroy_sampler_exact_1d(spatial_sampler1d);
	mui_destroy_chrono_sampler_exact_1d(temporal_sampler1d);

	free(results1d);

	/***********************
	 * Multi 2D interfaces *
	 ***********************/

	char *domain2d = (char*) malloc(strlen(argv[1]) + 4);
	strcpy(domain2d, argv[1]);
	strcat(domain2d, "_2d");
	strcat(domain2d, "\0");

	// Define interface names in the format [interface_name]_[1,2,3...]
	char** interfaces2d = (char**) malloc(sizeof(char*) * num_interfaces);
	for(int i=0; i<num_interfaces; i++) {
		interfaces2d[i] = (char*) malloc(strlen(argv[2] + 3));
		strcpy(interfaces2d[i], argv[2]);
		strcat(interfaces2d[i], "_");
		char tmp[11];
		sprintf(tmp, "%d", i);
		strcat(interfaces2d[i], tmp);
		strcat(interfaces2d[i], "\0");
	}

	// Define return structure to hold created unifaces
	mui_uniface_2d** unifaces2d = (mui_uniface_2d**) malloc(sizeof(mui_uniface_2d*) * num_interfaces);

	// Create MUI interfaces using helper function
	unifaces2d = mui_create_uniface_multi_2d( (const char*)domain2d, (const char**)interfaces2d, num_interfaces );

	// Free char memory after uniface creation
	for(int i=0; i<num_interfaces; i++) {
		free(interfaces2d[i]);
	}

	free(domain2d);

	// Define point for push
	mui_point_2d push_point2d = { 0, 0 };

	// Push value=1 through the MUI interfaces with the tag "position" at location=push_point
	for(int i=0; i<num_interfaces; i++) {
		mui_push_2d(unifaces2d[i], "position", push_point2d, 1);
	}

	// Commit (transmit) the pushed value at time=0
	for(int i=0; i<num_interfaces; i++) {
		mui_commit_2d(unifaces2d[i], 0);
	}

	// Create spatial and temporal samplers for fetch operation
	mui_sampler_exact_2d *spatial_sampler2d = mui_create_sampler_exact_2d(1e-37);
	mui_chrono_sampler_exact_2d *temporal_sampler2d = mui_create_chrono_sampler_exact_2d(1e-37);

	// Define point for fetch
	mui_point_2d fetch_point2d = { 0, 0 };

	double* results2d = (double*) malloc(sizeof(double) * num_interfaces);

	// Fetch the value for tag "position" at fetch_point at time=0
	for(int i=0; i<num_interfaces; i++) {
		results2d[i] = mui_fetch_exact_exact_2d(unifaces2d[i], "position", fetch_point2d, 0, spatial_sampler2d,
				temporal_sampler2d);

		printf("Fetched 2D interface_%d value = %f \n", i, results2d[i]);
	}

	// Destroy created 2D MUI objects
	mui_destroy_sampler_exact_2d(spatial_sampler2d);
	mui_destroy_chrono_sampler_exact_2d(temporal_sampler2d);

	free(results2d);

	/***********************
	 * Multi 3D interfaces *
	 ***********************/

	char *domain3d = (char*) malloc(strlen(argv[1]) + 4);
	strcpy(domain3d, argv[1]);
	strcat(domain3d, "_2d");
	strcat(domain3d, "\0");

	// Define interface names in the format [interface_name]_[1,2,3...]
	char** interfaces3d = (char**) malloc(sizeof(char*) * num_interfaces);
	for(int i=0; i<num_interfaces; i++) {
		interfaces3d[i] = (char*) malloc(strlen(argv[2] + 3));
		strcpy(interfaces3d[i], argv[2]);
		strcat(interfaces3d[i], "_");
		char tmp[11];
		sprintf(tmp, "%d", i);
		strcat(interfaces3d[i], tmp);
		strcat(interfaces3d[i], "\0");
	}

	// Define return structure to hold created unifaces
	mui_uniface_3d** unifaces3d = (mui_uniface_3d**) malloc(sizeof(mui_uniface_3d*) * num_interfaces);

	// Create MUI interfaces using helper function
	unifaces3d = mui_create_uniface_multi_3d( (const char*)domain3d, (const char**)interfaces3d, num_interfaces );

	// Free char memory after uniface creation
	for(int i=0; i<num_interfaces; i++) {
		free(interfaces3d[i]);
	}

	free(domain3d);

	// Define point for push
	mui_point_3d push_point3d = { 0, 0, 0 };

	// Push value=1 through the MUI interfaces with the tag "position" at location=push_point
	for(int i=0; i<num_interfaces; i++) {
		mui_push_3d(unifaces3d[i], "position", push_point3d, 1);
	}

	// Commit (transmit) the pushed value at time=0
	for(int i=0; i<num_interfaces; i++) {
		mui_commit_3d(unifaces3d[i], 0);
	}

	// Create spatial and temporal samplers for fetch operation
	mui_sampler_exact_3d *spatial_sampler3d = mui_create_sampler_exact_3d(1e-37);
	mui_chrono_sampler_exact_3d *temporal_sampler3d = mui_create_chrono_sampler_exact_3d(1e-37);

	// Define point for fetch
	mui_point_3d fetch_point3d = { 0, 0, 0 };

	double* results3d = (double*) malloc(sizeof(double) * num_interfaces);

	// Fetch the value for tag "position" at fetch_point at time=0
	for(int i=0; i<num_interfaces; i++) {
		results3d[i] = mui_fetch_exact_exact_3d(unifaces3d[i], "position", fetch_point3d, 0, spatial_sampler3d,
				temporal_sampler3d);

		printf("Fetched 3D interface_%d value = %f \n", i, results3d[i]);
	}

	// Destroy created 2D MUI objects
	mui_destroy_sampler_exact_3d(spatial_sampler3d);
	mui_destroy_chrono_sampler_exact_3d(temporal_sampler3d);

	free(results3d);

	// Destroy created MUI interfaces note: calls MPI_Finalize(), so need to do last
	for(int i=0; i<num_interfaces; i++) {
		mui_destroy_uniface_1d(unifaces1d[i]);
		mui_destroy_uniface_2d(unifaces2d[i]);
		mui_destroy_uniface_3d(unifaces3d[i]);
	}

	return 0;
}
