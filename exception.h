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

inline std::ostream& operator << ( std::ostream &out, std::exception const& err ) {
	return ( out << err.what() );
}

struct exception_segv {
	explicit inline exception_segv() {
		handle();
	}
	template<typename HEAD, typename ...TAIL>
	inline exception_segv( HEAD const &head, TAIL const &... tail ) {
		handle( head, tail... );
	}
	template<typename HEAD, typename ...TAIL>
	inline void handle( HEAD const &head, TAIL const &... tail ) {
		std::cerr << head;
		handle( tail... );
	}
	inline void handle( void ) {
		raise( SIGSEGV );
	}
};

struct exception_abort {
	inline exception_abort( const std::exception& except ) {
		std::cerr << except.what() << std::endl;
		raise( SIGABRT );
	}
};

struct exception_throw {
	inline exception_throw( const std::exception& ) {
		throw;
	}
};

}
#endif /* EXCEPTION_H_ */
