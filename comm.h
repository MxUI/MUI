/* -*- c++ -*-
 * comm.h
 *
 *  Created on: Feb 10, 2014
 *      Author: skudo
 */

#ifndef COMM_H_
#define COMM_H_

#include <vector>
#include "message.h"

namespace mui {
/*
 * Communication Interface class
 * Comm is the base class of communication libraries.
 * This class is designed for internal-use.
 */
class communicator {
private:
	// Comm is not copyable.
	communicator( const communicator& ) = delete;
	communicator& operator=( const communicator& ) = delete;
public:
	communicator() {}
	virtual ~communicator() {}

	virtual int local_rank() const { return 0; }
	virtual int local_size() const { return 1; }
	virtual int remote_size() const { return 1; }

	// send message
	void send( message msg, const std::vector<bool> &is_sending ) {
		if( is_sending.size() == remote_size() )
			return send_impl_(std::move(msg), is_sending);
		else {
			std::vector<bool> dest = is_sending;
			dest.resize(remote_size(), true);
			return send_impl_(std::move(msg), dest);
		}
	}
	void send( message msg ) {
		std::vector<bool> is_sending(remote_size(),true);
		return send_impl_(std::move(msg),is_sending);
	}

	// recv messages
	message recv() {
		return recv_impl_();
	}


protected:
	virtual void send_impl_( message msg, const std::vector<bool> &is_sending ) = 0;
	virtual message recv_impl_() = 0;
};
}

#endif
