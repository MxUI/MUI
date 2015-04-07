/*
 * unit_test.c
 *
 *  Created on: Jan 20, 2015
 *      Author: ytang
 */

/* test whehter the wrapper works */

#include "mui_3d.h"

int main( int argc, char **argv )
{
	mui_uniface3d *uniface = mui_create_uniface3d( argv[1] );

	mui_push( uniface, "position", 0, 0, 0, 0 );

	mui_destroy_uniface3d( uniface );

	return 0;
}



