/*
 * lib_factory.h
 *
 *  Created on: Mar 10, 2014
 *	  Author: ytang
 */

#ifndef LIB_FACTORY_H_
#define LIB_FACTORY_H_

#include "util.h"
#include "exception.h"

namespace mui
{

template<
	typename TAG,
	class    CREATOR,
	class    EXCEPTION=exception_segv>
struct factory
{
	CREATOR& create( const TAG id ) {
		assoc_iter i = dtable_.find(id);
		if ( i == dtable_.end() ) EXCEPTION();
		return i->second;
	}
	bool link( const TAG id, CREATOR creator ) {
		return dtable_.insert( make_pair(id,creator) ).second;
	}
	bool unlink( const TAG id ) {
		return dtable_.erase(id) == 1;
	}
protected:
	using assoc_table = std::unordered_map<TAG,CREATOR>;
	using assoc_iter  = typename assoc_table::iterator;
	assoc_table dtable_;
};

}

#endif /* LIB_FACTORY_H_ */
