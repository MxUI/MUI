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
 * @file lib_uri.h
 * @author Y. H. Tang
 * @date 14 March 2014
 * @brief Base class to contain and manipulate a unique URI (Uniform Resource
 * Identifier).
 */

#ifndef LIB_URI_H_
#define LIB_URI_H_

#include "../general/util.h"

namespace mui {

class uri {
public:
	uri(const std::string& url_s) {
		parse(url_s);
	}
	uri(const char url_c[]) {
		parse(url_c);
	}
	const std::string& protocol() const { return protocol_; }
	const std::string& host()     const { return host_; }
	const std::string& path()     const { return path_; }

	uri( const uri &another ) = delete;
	uri& operator = ( const uri &another ) = delete;
private:
	void parse(const std::string& url_s) {
		// "__protocol__://__host__/__path__"
		std::size_t prot_end = url_s.find("://");
		protocol_ = url_s.substr(0,prot_end);
		std::size_t host_end = url_s.find("/",prot_end+3);
		host_ = url_s.substr(prot_end+3,host_end-prot_end-3);
		path_ = url_s.substr(host_end+1);

		std::transform(protocol_.begin(), protocol_.end(), protocol_.begin(), ::tolower);
		std::transform(host_.begin(), host_.end(), host_.begin(), ::tolower);
	}

	std::string protocol_, host_, path_;
};

}

#endif /* LIB_URI_H_ */
