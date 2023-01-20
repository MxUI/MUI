/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library                       *
*                                                                            *
* Copyright (C) 2019 Y. H. Tang, S. Kudo, X. Bian, Z. Li, G. E. Karniadakis, *
*                    R. W. Nash                                              *
*                                                                            *
* (* The University of Edinburgh)                                            *
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
 * @file endian_traits.h
 * @author R. W. Nash <r.nash@epcc.ed.ac.uk>
 * @date 14 November 2018
 * @brief Support for dealing with big and little endian platforms.
 *
 * Currently this only supports big and little endian platforms, if
 * you're on something odder, then this wont work.
 *
 * The primary export of this header is a template struct
 * endian_traits - see below for a details
 *
 * You can configure endianness and conversion by defining one of the following
 *
 *  1. no endian options, in which case this header will try to figure
 *     it out but that probably only works on GCC
 *
 *  2. MUI_IGNORE_ENDIAN - this mode will never reorder the bytes of a
 *     multibyte number
 *
 *  3. (MUI_INT_BIG_ENDIAN ^ MUI_INT_LITTLE_ENDIAN) && (MUI_FLOAT_BIG_ENDIAN ^ MUI_FLOAT_LITTLE_ENDIAN)
 *
 *     where:
 *      - MUI_INT_BIG_ENDIAN - this will **assume** that the host has a big endian representation of integers
 *
 *      - MUI_INT_LITTLE_ENDIAN - this will **assume** that the host has a little endian representation of integers
 *
 *      - MUI_FLOAT_BIG_ENDIAN - this will **assume** that the host has a big endian representation of floating values
 *
 *      - MUI_FLOAT_LITTLE_ENDIAN - this will **assume** that the host has a little endian representation of floating values
 */

#ifndef MUI_ENDIAN_TRAITS_H_
#define MUI_ENDIAN_TRAITS_H_

#include <cstdint>

#ifdef __APPLE__
// Provide mappings to Apple specific function/macros
#  include <machine/endian.h>
#  include <libkern/OSByteOrder.h>

#  define htobe16(x) OSSwapHostToBigInt16(x)
#  define htole16(x) OSSwapHostToLittleInt16(x)
#  define be16toh(x) OSSwapBigToHostInt16(x)
#  define le16toh(x) OSSwapLittleToHostInt16(x)

#  define htobe32(x) OSSwapHostToBigInt32(x)
#  define htole32(x) OSSwapHostToLittleInt32(x)
#  define be32toh(x) OSSwapBigToHostInt32(x)
#  define le32toh(x) OSSwapLittleToHostInt32(x)

#  define htobe64(x) OSSwapHostToBigInt64(x)
#  define htole64(x) OSSwapHostToLittleInt64(x)
#  define be64toh(x) OSSwapBigToHostInt64(x)
#  define le64toh(x) OSSwapLittleToHostInt64(x)
#else
#  include <endian.h>
#endif

#  define MUI_POSITIVE true
#  define MUI_NEGATIVE false

// Convert input definitions (see above) to two macros returning bool
// (MUI_CONVERT_FLOAT & MUI_CONVERT_INT) that tell the implementation
// what to do
#ifdef MUI_IGNORE_ENDIAN
// Sanity check
#  if defined(MUI_INT_BIG_ENDIAN) || defined(MUI_INT_LITTLE_ENDIAN) || defined(MUI_FLOAT_BIG_ENDIAN) || defined(MUI_FLOAT_LITTLE_ENDIAN)
#    error "MUI Error [endian_traits.h]: Must set no other MUI endian options with MUI_IGNORE_ENDIAN"
#  endif

// Never convert
#  define MUI_CONVERT_INT false
#  define MUI_CONVERT_FLOAT false

#else // Not ignoring

// Integer
#  ifdef MUI_INT_BIG_ENDIAN
#    ifdef MUI_INT_LITTLE_ENDIAN
#      error "MUI Error [endian_traits.h]: Both MUI_INT_BIG_ENDIAN and MUI_INT_LITTLE_ENDIAN defined"
#    else
#      define MUI_INT_DEFINED MUI_POSITIVE
#      define MUI_CONVERT_INT false
#    endif
#  else
#    ifdef MUI_INT_LITTLE_ENDIAN
#      define MUI_INT_DEFINED MUI_POSITIVE
#      define MUI_CONVERT_INT true
#    else
// We have no endian options - try to figure it out.
#      if defined(__BYTE_ORDER__)
#        define MUI_INT_DEFINED MUI_POSITIVE
#        define MUI_CONVERT_INT (__BYTE_ORDER__ != __ORDER_BIG_ENDIAN__)
#      else
#        define MUI_INT_DEFINED MUI_NEGATIVE
#        error "MUI Error [endian_traits.h]: Cannot auto-detect integer endianness of platform - please set at compilation (-DMUI_INT_LITTLE_ENDIAN or -DMUI_INT_BIG_ENDIAN) or set to ignore (-DMUI_IGNORE_ENDIAN)"
#      endif
#    endif
#  endif

// Floating point
#  ifdef MUI_FLOAT_BIG_ENDIAN
#    ifdef MUI_FLOAT_LITTLE_ENDIAN
#      error "MUI Error [endian_traits.h]: Both MUI_FLOAT_BIG_ENDIAN and MUI_FLOAT_LITTLE_ENDIAN defined"
#    else
#      define MUI_CONVERT_FLOAT false
#    endif
#  else
#    ifdef MUI_FLOAT_LITTLE_ENDIAN
#      define MUI_CONVERT_FLOAT true
#    else
// We have no endian options - try to figure it out.
#      if defined(__FLOAT_WORD_ORDER__)
#        define MUI_CONVERT_FLOAT (__FLOAT_WORD_ORDER__ != __ORDER_BIG_ENDIAN__)
#      else // Likely using something other than GNU
#        if MUI_INT_DEFINED == MUI_POSITIVE
#          define MUI_CONVERT_FLOAT MUI_CONVERT_INT
#          warning "MUI Warning [endian_traits.h]: Cannot auto-detect float endianness of platform - please set at compilation (-DMUI_FLOAT_LITTLE_ENDIAN or -DMUI_FLOAT_BIG_ENDIAN) or set to ignore (-DMUI_IGNORE_ENDIAN), integer endianness will be followed by default"
#        else
#          error "MUI Error [endian_traits.h]: Cannot auto-detect float endianness of platform - please set at compilation (-DMUI_FLOAT_LITTLE_ENDIAN or -DMUI_FLOAT_BIG_ENDIAN) or set to ignore (-DMUI_IGNORE_ENDIAN)"
#        endif
#      endif
#    endif
#  endif

#endif // End definition processing

namespace mui {
  
  namespace detail {
    // Metafunction to get an unsigned integer type of the specified
    // size in bytes as member type named "type".
    template<size_t size_bytes>
    struct uint;
    // Specialisations for 8, 16, 32, 64 bits
    template<>
    struct uint<1> {
      using type = uint8_t;
    };   
    template<>
    struct uint<2> {
      using type = uint16_t;
    };
    template<>
    struct uint<4> {
      using type = uint32_t;
    };
    template<>
    struct uint<8> {
      using type = uint64_t;
    };
    // Helper class to convert betweem host and big-endian
    // (i.e. network) ordering using the htobe/betoh family of
    // functions.
    template<size_t size_bytes>
    struct endian_converter {
      union data_t {
	char buf[size_bytes];
	typename uint<size_bytes>::type val;
      };
      data_t data;
      void htobe();
      void betoh();
    };
    // Specialisations
    template<>
    inline void endian_converter<2>::htobe() {
      data.val = htobe16(data.val);
    }
    template<>
    inline void endian_converter<2>::betoh() {
      data.val = be16toh(data.val);
    }
    template<>
    inline void endian_converter<4>::htobe() {
      data.val = htobe32(data.val);
    }
    template<>
    inline void endian_converter<4>::betoh() {
      data.val = be32toh(data.val);
    }
    template<>
    inline void endian_converter<8>::htobe() {
      data.val = htobe64(data.val);
    }
    template<>
    inline void endian_converter<8>::betoh() {
      data.val = be64toh(data.val);
    }
  }
  
  // Traits class for controlling MUI's behaviour about endianness.
  //
  // Uses SFINAE to select the appropriate specialisation based on the
  // type.
  // 
  // It defines a bool constexpr member "convert" that indicates if
  // the type should undergo byte reordering as it is (de-)serialised.
  //
  // Currently will handle integral (signed or unsigned, widths of 8,
  // 16, 32, 64 b) and floating (32 or 64 b) types.
  template<typename T, typename enable = void>
  struct endian_traits;

  // Integers
  template<typename T>
  struct endian_traits<T, typename std::enable_if<std::is_integral<T>::value>::type>
  {
    static constexpr bool convert = (sizeof(T) > 1) && MUI_CONVERT_INT;
  };

  // Floats
  template<typename T>
  struct endian_traits<T, typename std::enable_if<std::is_floating_point<T>::value>::type>
  {
    static constexpr bool convert = MUI_CONVERT_FLOAT;
  };

}

#endif // Include guard
