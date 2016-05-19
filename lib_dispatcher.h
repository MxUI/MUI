/*
 * lib_dispatcher.h
 *
 *  Created on: Mar 10, 2014
 *	  Author: ytang
 */

#ifndef LIB_DISPATCHER_H_
#define LIB_DISPATCHER_H_

#include "util.h"
#include "exception.h"
#include <unordered_map>

namespace mui
{

template<
	typename UUID,
	class    FPTR,
	class    EXCEPTION=exception_segv>
struct dispatcher
{
	FPTR dispatch( const UUID &id ) {
		auto i = dtable_.find(id);
		if ( i == dtable_.end() ) EXCEPTION();
		return i->second;
	}
	bool exist( const UUID &id ) {
		return dtable_.find(id) != dtable_.end();
	}
	FPTR operator [] ( const UUID &id ) {
		return dispatch(id);
	}
	bool link( const UUID &id, FPTR parser ) {
		return dtable_.insert( std::make_pair(id,parser) ).second;
	}
	bool unlink( const UUID &id ) {
		return dtable_.erase(id) == 1;
	}
protected:
	using assoc_table = std::unordered_map<UUID,FPTR>;
	assoc_table dtable_;
};

}

#endif /* LIB_DISPATCHER_H_ */
