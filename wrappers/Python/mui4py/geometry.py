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
# @file geometry.py
# @author E. R. Fernandez
# @date 25 January 2019
# @brief Geometry functions for MUI Python wrapper.
#
"""

from mui4py.common import CppClass, _Point
from mui4py.config import Config


class Geometry(CppClass):
    def __init__(self, args=(), kwargs={}):
        super(Geometry, self).__init__(Config(), args, kwargs)
        self.namespace = "geometry"

    def bbox(self):
        return self.raw.bbox()


def collide(shape1, shape2):
    pass


class Box(Geometry):
    def __init__(self, x1, x2):
        super(Box, self).__init__(args=(_Point(x1), _Point(x2)))


class Sphere(Geometry):
    def __init__(self, x0, r):
        super(Sphere, self).__init__(args=(_Point(x0), r))


class Point(Geometry):
    def __init__(self, x0):
        super(Point, self).__init__(args=(_Point(x0), ))
