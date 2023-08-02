/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library                       *
*                                                                            *
* Copyright (C) 2023 W. Liu                                                  *
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
 * @file linalg_util.h
 * @author W. Liu
 * @date 17 mAY 2023
 * @brief Utility functions for mui::linalg.
 */

#ifndef MUI_LINALG_UTIL_H_
#define MUI_LINALG_UTIL_H_

#include <sstream>
#include <cassert>
#include <algorithm>
#include <cctype>

namespace mui {
namespace linalg {

// Function to left trim a string - helper function on matrix file I/O
inline std::string ltrim(const std::string &s) {
    const std::string WHITESPACE = " \n\r\t\f\v";
    std::string::size_type start = s.find_first_not_of(WHITESPACE);
    return (start == std::string::npos) ? "" : s.substr(start);
}

// Function to right trim a string - helper function on matrix file I/O
inline std::string rtrim(const std::string &s) {
    const std::string WHITESPACE = " \n\r\t\f\v";
    std::string::size_type end = s.find_last_not_of(WHITESPACE);
    return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

// Function to trim a string on both sides - helper function on matrix file I/O
inline std::string trim(const std::string &s) {
    return rtrim(ltrim(s));
}

// Function to convert the input string to all lowercase characters - helper function on matrix file I/O
inline std::string string_to_lower(const std::string &s) {
    std::string lower;
    std::transform(s.begin(), s.end(), std::back_inserter(lower), [](unsigned char c){ return std::tolower(c); });
    return lower;
}

// Function to convert the input string to all uppercase characters - helper function on matrix file I/O
inline std::string string_to_upper(const std::string &s) {
    std::string upper;
    std::transform(s.begin(), s.end(), std::back_inserter(upper), [](unsigned char c){ return std::toupper(c); });
    return upper;
}

} // linalg
} // mui

#endif /* MUI_LINALG_UTIL_H_ */
