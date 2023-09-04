"""
##############################################################################
# Multiscale Universal Interface Code Coupling Library                       #
#                                                                            #
# Copyright (C) 2023 W. Liu                                                  #
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
# @file algorithms.py
# @author W. Liu
# @date 18 March 2019
# @brief Coupling algorithms for the MUI Python wrapper.
#
"""

from mui4py.common import CppClass
from mui4py.config import Config
from mui4py.types import UINT, UINT32, UINT64, FLOAT, FLOAT32, FLOAT64, INT, INT32, INT64, STRING, BOOL

# Interface for TemporalSampler
def algorithm_signature(self):
    sig = self._split_class_name(title=False)
    return sig

# Algorithms
class Algorithm(CppClass):
    def __init__(self, config=None, args=(), kwargs={}):
        if config is not None:
            super(Algorithm, self).__init__(config, args, kwargs)
            self.configure(self.config)
        else:
            # Empty config to not trigger error in default config.
            super(Algorithm, self).__init__(Config(), args, kwargs)

Algorithm.fetch_signature = algorithm_signature

class AlgorithmFixedRelaxation(Algorithm):
    def __init__(self, under_relaxation_factor=None, local_comm=None, pts_value_init=None, config=None):
        if config is None:
            self.config=Config()
        else:
            self.config=config
        if pts_value_init is not None:
            super(AlgorithmFixedRelaxation, self).__init__(config, args=(under_relaxation_factor, local_comm, pts_value_init))
        else:
            if local_comm is not None:
                super(AlgorithmFixedRelaxation, self).__init__(config, args=(under_relaxation_factor, local_comm, []))
            else:
                if under_relaxation_factor is not None:
                    super(AlgorithmFixedRelaxation, self).__init__(config, args=(under_relaxation_factor, None, []))
                else:
                    super(AlgorithmFixedRelaxation, self).__init__(1.0, None, [])
        self._ALLOWED_IO_TYPES = [INT, INT32, INT64, UINT, UINT32, UINT64, FLOAT, FLOAT32, FLOAT64]

    def get_under_relaxation_factor(self, t1, t2=None):
        if t2 is not None:
            return self.raw.get_under_relaxation_factor(t1, t2)
        else:
            return self.raw.get_under_relaxation_factor(t1)

    def get_residual_L2_Norm(self, t1, t2=None):
        if t2 is not None:
            return self.raw.get_residual_L2_Norm(t1, t2)
        else:
            return self.raw.get_residual_L2_Norm(t1)

class AlgorithmAitken(Algorithm):
    def __init__(self, under_relaxation_factor=None, under_relaxation_factor_max=None, local_comm=None, pts_vlu_init=None, res_l2_norm_nm1=None, config=None):
        if config is None:
            self.config=Config()
        else:
            self.config=config
        if res_l2_norm_nm1 is not None:
            super(AlgorithmAitken, self).__init__(config, args=(under_relaxation_factor, under_relaxation_factor_max, local_comm, pts_vlu_init, res_l2_norm_nm1))
        else:
            if pts_vlu_init is not None:
                super(AlgorithmAitken, self).__init__(config, args=(under_relaxation_factor, under_relaxation_factor_max, local_comm, pts_vlu_init, 0.0))
            else:
                if local_comm is not None:
                    super(AlgorithmAitken, self).__init__(config, args=(under_relaxation_factor, under_relaxation_factor_max, local_comm, [], 0.0))
                else:
                    if under_relaxation_factor_max is not None:
                        super(AlgorithmAitken, self).__init__(config, args=(under_relaxation_factor, under_relaxation_factor_max, None, [], 0.0))
                    else:
                        if under_relaxation_factor is not None:
                            super(AlgorithmAitken, self).__init__(config, args=(under_relaxation_factor, 1.0, None, [], 0.0))
                        else:
                            super(AlgorithmAitken, self).__init__(config, args=(1.0, 1.0, None, [], 0.0))
        self._ALLOWED_IO_TYPES = [INT, INT32, INT64, UINT, UINT32, UINT64, FLOAT, FLOAT32, FLOAT64]

    def get_under_relaxation_factor(self, t1, t2=None):
        if t2 is not None:
            return self.raw.get_under_relaxation_factor(t1, t2)
        else:
            return self.raw.get_under_relaxation_factor(t1)

    def get_residual_L2_Norm(self, t1, t2=None):
        if t2 is not None:
            return self.raw.get_residual_L2_Norm(t1, t2)
        else:
            return self.raw.get_residual_L2_Norm(t1)
