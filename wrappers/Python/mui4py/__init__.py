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
# @file __init__.py
# @author E. R. Fernandez, W. Liu
# @date 25 January 2019
# @brief Marker file denotes the directory containing the file as the MUI Python wrapper package.
#
"""

# flake8: noqa
from mui4py.mui4py import Uniface, mpi_split_by_app, set_quiet,\
                          set_data_types_unifaces, create_unifaces,\
                          get_mpi_version, get_compiler_version, get_compiler_config
from mui4py.samplers import SamplerExact, SamplerGauss, SamplerMovingAverage,\
                            SamplerNearestNeighbor, SamplerPseudoNearest2Linear,\
                            SamplerPseudoNearestNeighbor, SamplerSherpardQuintic,\
                            SamplerSphQuintic, SamplerSumQuintic, SamplerRbf
from mui4py.temporal_samplers import TemporalSamplerExact, TemporalSamplerGauss,\
                            TemporalSamplerMean, TemporalSamplerSum
from mui4py.algorithms import AlgorithmFixedRelaxation, AlgorithmAitken
from mui4py.types import STRING, INT32, INT64, INT, UINT32, UINT64, UINT, FLOAT32, FLOAT64, FLOAT, BOOL
from mui4py.config import Config, set_default_config, get_default_config
import mui4py.geometry
