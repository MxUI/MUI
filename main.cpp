/*
 * main.cpp
 *
 *  Created on: Jan 28, 2014
 *      Author: ytang
 */

//#include "mui.h"
#include "sampler_gauss.h"

int main()
{
//	mui::uniface<> interface("mpi://atomistic/a2c");

	mui::sampler_gauss<double> s1;
	mui::chrono_sampler_gauss<> s2;
	return 0;
}

