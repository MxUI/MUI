/*
 * lib_singleton.h
 *
 *  Created on: Mar 14, 2014
 *      Author: ytang
 */

#ifndef LIB_SINGLETON_H_
#define LIB_SINGLETON_H_

#include "util.h"

namespace mui {

template<class T>
class singleton
{
public:
	static T& instance() {
		static T instance_;
		return instance_;
	}
private:
	// make all constructors etc. private to guarantee uniqueness
	singleton();
	singleton(const singleton&);
	singleton& operator=(const singleton&);
	~singleton();
};

}

#endif /* LIB_SINGLETON_H_ */
