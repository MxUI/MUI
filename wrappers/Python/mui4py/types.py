"""
##############################################################################
# Multiscale Universal Interface Code Coupling Library                       #
#                                                                            #
# Copyright (C) 2023 E. R. Fernandez, W. Liu                                 #
#                                                                            #
# This software is jointly licensed under the Apache License, Version 2.0    #
# and the GNU General Public License version 3, you may use it according     #
# to either.                                                                 #
#                                                                            #
# ** Apache License, version 2.0 **                                          #
#                                                                            #
# Licensed under the Apache License, Version 2.0 (the "License");            #
# you may not use this file except in compliance with the License.           #
# You may obtain a copy of the License at                                    #
#                                                                            #
# http://www.apache.org/licenses/LICENSE-2.0                                 #
#                                                                            #
# Unless required by applicable law or agreed to in writing, software        #
# distributed under the License is distributed on an "AS IS" BASIS,          #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   #
# See the License for the specific language governing permissions and        #
# limitations under the License.                                             #
#                                                                            #
# ** GNU General Public License, version 3 **                                #
#                                                                            #
# This program is free software: you can redistribute it and/or modify       #
# it under the terms of the GNU General Public License as published by       #
# the Free Software Foundation, either version 3 of the License, or          #
# (at your option) any later version.                                        #
#                                                                            #
# This program is distributed in the hope that it will be useful,            #
# but WITHOUT ANY WARRANTY; without even the implied warranty of             #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
# GNU General Public License for more details.                               #
#                                                                            #
# You should have received a copy of the GNU General Public License          #
# along with this program.  If not, see <http://www.gnu.org/licenses/>.      #
##############################################################################
#
# @file types.py
# @author E. R. Fernandez, W. Liu
# @date 25 January 2019
# @brief Data type functions for the MUI Python wrapper.
#
"""

import numpy as np

# Type enumeration
STRING = str
BOOL = bool
INT32 = np.int32
INT64 = np.int64
INT = int
UINT32 = np.uint32
UINT64 = np.uint64
UINT = UINT64
FLOAT32 = np.float32
FLOAT64 = np.float64
FLOAT = float

map_type = {STRING: STRING, BOOL: BOOL,
            INT32: INT32, INT64: INT64,
            INT: INT64,
            UINT32: UINT32, UINT64: UINT64,
            UINT: UINT64,
            FLOAT32: FLOAT32, FLOAT64: FLOAT64,
            FLOAT: FLOAT64,
            None: None}

__int_size = int(str(np.iinfo(int).dtype)[-2:])
if __int_size == 32:
    map_type[INT] = INT32
__uint_size = int(str(np.iinfo(int).dtype)[-2:])
if __uint_size == 32:
    _map_type[UINT] = UINT32
__float_size = int(str(np.finfo(float).dtype)[-2:])
if __float_size == 32:
    map_type[FLOAT] = FLOAT32
assert __int_size in [32, 64]
assert __float_size in [32, 64]
__io_float_map = {64: "double", 32: "float"}
__io_int_map = {64: "int64_t", 32: "int32_t"}
__io_uint_map = {64: "uint64_t", 32: "uint32_t"}

# Types allowed for configuring the library
ALLOWED_INT_TYPES = {INT32: "i32", INT64: "i64"}  # , int: "i%d" % __int_size}
ALLOWED_UINT_TYPES = {UINT32: "u32", UINT64: "u64"} #, uint: "i%d" % __int_size}
ALLOWED_FLOAT_TYPES = {FLOAT32: "f32", FLOAT64: "f64"}  # , float: "f%d" % __float_size}

# Types allowed to be pushed/fetched
ALLOWED_IO_TYPES = {FLOAT32: "float", FLOAT64: "double",
                    INT32: "int32_t", INT64: "int64_t",
                    UINT32: "uint32_t", UINT64: "uint64_t",
                    FLOAT: __io_float_map[__float_size],
                    UINT: __io_uint_map[__int_size],
                    INT: __io_int_map[__int_size],
                    str: "string"}


def get_int_type_str(typein):
    if typein in ALLOWED_INT_TYPES.keys():
        return ALLOWED_INT_TYPES[typein]
    else:
        raise Exception("Integer type '{}' not supported. Supported types : [int, np.int32, np.int64]".format(typein))

def get_uint_type_str(typein):
    if typein in ALLOWED_UINT_TYPES.keys():
        return ALLOWED_UINT_TYPES[typein]
    else:
        pass

def get_float_type_str(typein):
    if typein in ALLOWED_FLOAT_TYPES.keys():
        return ALLOWED_FLOAT_TYPES[typein]
    else:
        raise Exception("Float type '{}' not supported. "
                        "Supported types : [float, np.float32, np.float64]".format(typein))

def get_io_type_str(typein):
    if typein in ALLOWED_IO_TYPES.keys():
        return ALLOWED_IO_TYPES[typein]
    else:
        raise Exception("Float type '{}' not supported. "
                        "Supported types : [float, np.float32, np.float64]".format(typein))


def safe_cast(value_type, value):
    if not isinstance(value, str):  # Only for numerics types
        if not np.can_cast(value, value_type):
            raise Exception("Value '{}' cannot be safely casted to type '{}'.".format(value, value_type.__name__))
    return value_type(value)
