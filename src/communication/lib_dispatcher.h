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
 * @file lib_dispatcher.h
 * @author Y. H. Tang
 * @date 10 March 2014
 * @brief Structure for communicator used in comm_factory.h
 */

#ifndef LIB_DISPATCHER_H_
#define LIB_DISPATCHER_H_

#include "../general/util.h"
#include "../general/exception.h"
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
