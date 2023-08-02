/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library                       *
*                                                                            *
* Copyright (C) 2019 Y. H. Tang, S. Kudo, X. Bian, Z. Li, G. E. Karniadakis  *
*                                                                            *
* This software is jointly licensed under the Apache License, Version 2.0    *
* and the GNU General Public License version 3, you may use it according     *
* to either.                                                                 *
*                                                                            *
* ** Apache License, version 2.0 **                                          *
*                                                                            *
* Licensed under the Apache License, Version 2.0 (the "License");            *
* you may not use this file except in compliance with the License.           *
* You may obtain a copy of the License at                                    *
*                                                                            *
* http://www.apache.org/licenses/LICENSE-2.0                                 *
*                                                                            *
* Unless required by applicable law or agreed to in writing, software        *
* distributed under the License is distributed on an "AS IS" BASIS,          *
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
* See the License for the specific language governing permissions and        *
* limitations under the License.                                             *
*                                                                            *
* ** GNU General Public License, version 3 **                                *
*                                                                            *
* This program is free software: you can redistribute it and/or modify       *
* it under the terms of the GNU General Public License as published by       *
* the Free Software Foundation, either version 3 of the License, or          *
* (at your option) any later version.                                        *
*                                                                            *
* This program is distributed in the hope that it will be useful,            *
* but WITHOUT ANY WARRANTY; without even the implied warranty of             *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
* GNU General Public License for more details.                               *
*                                                                            *
* You should have received a copy of the GNU General Public License          *
* along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
*****************************************************************************/

/**
 * @file comm.h
 * @author S. Kudo
 * @date 10 February 2014
 * @brief File containing class definition of communication interface.
 * This is the base class for all other communication related classes.
 */

#ifndef COMM_H_
#define COMM_H_

#include <vector>
#include "message/message.h"

namespace mui {
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
	virtual std::string uri_host() const { return std::string(); }
	virtual std::string uri_path() const { return std::string(); }
	virtual std::string uri_protocol() const { return std::string(); }

	// send message
	void send( message msg, const std::vector<bool> &is_sending ) {
		if( is_sending.size() == static_cast<size_t>(remote_size()) )
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
