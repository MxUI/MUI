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
# @file temporal_sampler.py
# @author E. R. Fernandez, W. Liu
# @date 18 March 2019
# @brief Temporal sampler functions for the MUI Python wrapper.
#
"""

from mui4py.common import CppClass
from mui4py.config import Config
from mui4py.types import UINT, UINT32, UINT64, FLOAT, FLOAT32, FLOAT64, INT, INT32, INT64, STRING, BOOL

# Interface for TemporalSampler
def sampler_fetch_signature(self):
    sig = self._split_class_name(title=False)
    sig = sig.replace("_sampler", "")
    return sig.replace("sampler_", "")

# Temporal samplers
class TemporalSampler(CppClass):
    def __init__(self, args=(), kwargs={}):
        # Empty config to not trigger error in default config.
        super(TemporalSampler, self).__init__(Config(), args, kwargs)


TemporalSampler.fetch_signature = sampler_fetch_signature


class TemporalSamplerExact(TemporalSampler):
    def __init__(self, tol=None):
        super(TemporalSamplerExact, self).__init__(kwargs={"tol": tol})
        self._ALLOWED_IO_TYPES = [INT, INT32, INT64, UINT, UINT32, UINT64, FLOAT, FLOAT32, FLOAT64, STRING, BOOL]


class TemporalSamplerGauss(TemporalSampler):
    def __init__(self, cutoff, sigma):
        super(TemporalSamplerGauss, self).__init__(args=(cutoff, sigma))
        self._ALLOWED_IO_TYPES = [UINT, UINT32, UINT64, INT, INT32, INT64, FLOAT, FLOAT32, FLOAT64]


class TemporalSamplerMean(TemporalSampler):
    def __init__(self, newleft=None, newright=None):
        super(TemporalSamplerMean, self).__init__(kwargs={"newleft": newleft, "newright": newright})
        self._ALLOWED_IO_TYPES = [UINT, UINT32, UINT64, INT, INT32, INT64, FLOAT, FLOAT32, FLOAT64]


class TemporalSamplerSum(TemporalSampler):
    def __init__(self, newleft=None, newright=None):
        super(TemporalSamplerSum, self).__init__(kwargs={"newleft": newleft, "newright": newright})
        self._ALLOWED_IO_TYPES = [UINT, UINT32, UINT64, INT, INT32, INT64, FLOAT, FLOAT32, FLOAT64]
