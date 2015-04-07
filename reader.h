/*
 * reader.h
 *
 *  Created on: Mar 16, 2014
 *      Author: ytang
 */

#ifndef READER_H_
#define READER_H_

#include "util.h"
#include "message.h"

namespace mui {

struct reader {
	reader( const std::function<void(message)>& functor ): functor_(functor) {}
	void scan( message msg ) { functor_(msg); };
private:
	std::function<void(message)> functor_;
};

template<class CLASS, class FPTR>
reader make_reader( CLASS &obj, const FPTR &f ) {
	return reader( std::bind( std::mem_fn(f), &obj, std::placeholders::_1 ) );
}

}

#endif /* READER_H_ */
