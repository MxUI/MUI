/*
Multiscale Universal Interface Code Coupling Library

Copyright (C) 2017 Y. H. Tang, S. Kudo, X. Bian, Z. Li, G. E. Karniadakis

This software is jointly licensed under the Apache License, Version 2.0
and the GNU General Public License version 3, you may use it according
to either.

** Apache License, version 2.0 **

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

** GNU General Public License, version 3 **

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

** File Details **

Filename: reader.h
Created: Mar 16, 2014
Author: Y. H. Tang
Description:
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
