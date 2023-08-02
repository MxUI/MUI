"""
##############################################################################
# Multiscale Universal Interface Code Coupling Library                       #
#                                                                            #
# Copyright (C) 2023 E. R. Fernandez                                         #
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
# @file config.py
# @author E. R. Fernandez
# @date 25 January 2019
# @brief Configure functions for MUI Python wrapper.
#
"""

from mui4py.types import map_type, UINT
__default_config = None


def set_default_config(config):
    global __default_config
    __default_config = config


def get_default_config():
    if __default_config is None:
        raise Exception("Default configuration not defined.")
    return __default_config


class Config:
    def __init__(self, dim=None, float_type=float, force_casting=True):
        # NOTE: int_type is fixed to Python int size at compile time
        self.int_type = map_type[int]
        self.uint_type = map_type[UINT]
        self.float_type = map_type[float_type]
        self.dim = dim
        self._check_types()
        self.force_casting = force_casting

    def _check_types(self):
        pass
