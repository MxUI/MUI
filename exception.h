/*
 * exception.h
 *
 *  Created on: Mar 10, 2014
 *      Author: ytang
 */

#ifndef EXCEPTION_H_
#define EXCEPTION_H_

#include "util.h"

namespace mui {

struct exception_segv {
	exception_segv( const std::exception& err ) {
		std::cerr << err.what() << std::endl;
		raise( SIGSEGV );
	}
	exception_segv( std::string err=std::string() ) {
		std::cerr << err << std::endl;
		raise( SIGSEGV );
	}
};

struct exception_abort {
	exception_abort( const std::exception& except ) {
		std::cerr << except.what() << std::endl;
		raise( SIGABRT );
	}
};

struct exception_throw {
	exception_throw( const std::exception& ) {
		throw;
	}
};

}
#endif /* EXCEPTION_H_ */
