/*
 * comm_factory.h
 *
 *  Created on: Mar 14, 2014
 *      Author: ytang
 */

#ifndef COMM_FACTORY_H_
#define COMM_FACTORY_H_

#include "util.h"
#include "exception.h"
#include "lib_uri.h"
#include "lib_dispatcher.h"
#include "lib_singleton.h"
#include "comm.h"

namespace mui {

struct comm_factory: public singleton<dispatcher<std::string, std::function<communicator *(const char [])> > >
{
	static communicator *create_comm( const char URI[] ) {
		if ( !instance().exist(uri(URI).protocol()) ) {
			exception_segv( "Unknown communicator type ", uri(URI).protocol() );
		}
		return instance()[uri(URI).protocol()]( URI );
	}
};

}

#endif /* COMM_FACTORY_H_ */
