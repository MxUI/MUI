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
# @file mui4py.py
# @author E. R. Fernandez, W. Liu
# @date 25 January 2019
# @brief main file for the MUI Python wrapper.
#
"""

import mui4py.mui4py_mod as mui4py_mod
from mui4py.common import CppClass, get_cpp_name, array2Point
from mui4py.config import Config
from mui4py.types import safe_cast, map_type, ALLOWED_IO_TYPES, STRING, BOOL, INT32, INT64, INT, UINT32, UINT64, UINT, FLOAT32, FLOAT64, FLOAT
from mui4py.geometry import Geometry
from mui4py.samplers import Sampler
from mui4py.temporal_samplers import TemporalSampler
from mui4py.algorithms import Algorithm,AlgorithmAitken,AlgorithmFixedRelaxation
import copy


def create_unifaces(domain, ifaces_names, config, world=None):
    assert type(ifaces_names) == list
    assert type(domain) == str
    assert issubclass(config.__class__, Config)
    ifaces_out = {}
    cpp_obj_name = get_cpp_name("create_uniface", config.dim,
                                config.float_type, config.int_type)
    if world is not None:
        ifaceraw = getattr(mui4py_mod, cpp_obj_name)(domain, ifaces_names, world)
    else:
        ifaceraw = getattr(mui4py_mod, cpp_obj_name)(domain, ifaces_names)
    for i, obj in enumerate(ifaceraw):
        ifaces_out[ifaces_names[i]] = Uniface(config=config, cpp_obj=obj)
    return ifaces_out


def get_mpi_version():
    return mui4py_mod.get_mpi_version()


def get_compiler_version():
    return mui4py_mod.get_compiler_version()


def get_compiler_config():
    return mui4py_mod.get_compiler_config()


def mpi_split_by_app():
    return mui4py_mod.mpi_split_by_app()


def set_quiet(q):
    mui4py_mod.set_quiet(q)


def set_data_types_unifaces(ifaces, data):
    for iface_name, iface in ifaces.items():
        iface.set_data_types(data[iface_name])


# MUI Classes
class Uniface(CppClass):
    def __init__(self, uri=None, cpp_obj=None, config=None):
        super(Uniface, self).__init__(config, args=(uri,))
        self.uri = uri
        self.configure(self.config, cpp_obj=cpp_obj)
        self.tags_type = {}
        self._tags_spatial_samplers = {}
        self._tags_temporal_samplers = {}
        self._tags_coupling_algorithms = {}
        self._tags_fetch = {}
        self._tags_fetch_points = {}
        self._ALLOWED_PROTOCOLS = ["mpi"]
        self._aitken_under_relaxation_factor = []
        self._aitken_residual_l2_norm = []
        self._latest_aitken_under_relaxation_factor = 0.0
        self._latest_aitken_residual_l2_norm = 0.0
        self._fixed_relaxation_under_relaxation_factor = []
        self._fixed_relaxation_residual_l2_norm = []
        self._latest_fixed_relaxation_under_relaxation_factor = 0.0
        self._latest_fixed_relaxation_residual_l2_norm = 0.0

    def _get_tag_type(self, tag):
        try:
            return self.tags_type[tag]
        except KeyError:
            raise Exception("A type has not been defined for ata tag '{}'. Use 'Uniface.set_data_type()'.".format(tag))

    def set_data_types(self, data):
        for tag, data_type in data.items():
            try:
                self._set_data_type(tag, data_type)
            except KeyError:
                raise Exception("Uniface does not exist.")

    def _set_data_type(self, tag, data_type):
        if data_type not in ALLOWED_IO_TYPES.keys():
            raise Exception("Data type not supported. Supported types: {..}")
        try:
            data_type_stored = self.tags_type[tag]
        except KeyError:
            self.tags_type[tag] = data_type
            self._tags_spatial_samplers[tag] = {}
            self._tags_temporal_samplers[tag] = {}
            self._tags_coupling_algorithms[tag] = {}
            self._tags_fetch[tag] = {}
            self._tags_fetch_points[tag] = {}
        else:
            raise Exception("Type '{}' has already been defined for tag '{}'.".format(data_type_stored.__name__, tag))

    def _get_pushfname(self, fname_root, tag, val=None, type_in=None):
        assert val is not None or type is not None
        stored_data_type = map_type[self._get_tag_type(tag)]
        if self.config.force_casting:
            data_type = stored_data_type
        else:
            if val is not None:
                data_type = map_type[type(val)]
            elif type is not None:
                data_type = type_in
        if stored_data_type != data_type:
            raise Exception("Data type set for tag '{}' do not match with the "
                            "data type of the value provided.".format(tag))
        return (fname_root + ALLOWED_IO_TYPES[data_type], data_type)

    def push(self, *args, **kwargs):
        if len(args) == 1:
            loc = array2Point(args[0], self.config, self.raw_point)
            push_fname = "push_"
            pargs = (loc, )
        elif len(args) == 2:
            tag = args[0]
            val = args[1]
            push_fname, data_type = self._get_pushfname("push_", tag, val=val)
            pargs = (tag, safe_cast(data_type, val))
        elif len(args) == 3:
            tag = args[0]
            loc = array2Point(args[1], self.config, self.raw_point)
            val = args[2]
            push_fname, data_type = self._get_pushfname("push_", tag, val=val)
            try:
                pargs = (tag, loc, safe_cast(data_type, val))
            except ValueError:
                raise Exception("Forced type casting failed in push.")
        else:
            raise Exception("MUI Error [mui4py.py]: Push function accept 1, 2 or 3 parameters.")

        push = getattr(self.raw, push_fname)
        push(*pargs)

    def push_many(self, tag, points, values):
        # TODO: Try to apply safe_cast
        push_fname, data_type = self._get_pushfname("push_many_", tag, type_in=values.dtype.type)
        getattr(self.raw, push_fname)(tag, points, values)

    def commit(self, t1, t2=None):
        if t2 is not None:
            return self.raw.commit(t1, t2)
        else:
            return self.raw.commit(t1, mui4py_mod.numeric_limits_uint)

    def forecast(self, t1, t2=None):
        if t2 is not None:
            self.raw.forecast(t1, t2)
        else:
            self.raw.forecast(t1, mui4py_mod.numeric_limits_uint)

    def is_ready(self, attr, t1, t2=None):
        if t2 is not None:
            return self.raw.is_ready(attr, t1, t2)
        else:
            return self.raw.is_ready(attr, t1)

    def barrier(self, t1, t2=None):
        if t2 is not None:
            self.raw.barrier(t1, t2)
        else:
            self.raw.barrier(t1)

    # In case to pass (std::pair<Ttime,Titer>, bool) as args for `forget()` function:
    #    ```
    #    t1 = (2.5, 5)
    #    forget(t1, True)
    #    ```
    def forget(self, t1, t2=None, reset_log=True):
        if t2 is not None:
            self.raw.forget(t1, t2, reset_log)
        else:
            self.raw.forget(t1, reset_log)

    def set_memory(self, t):
        self.raw.set_memory(t)

    def uri_host(self):
        return self.raw.uri_host()

    def uri_path(self):
        return self.raw.uri_path()

    def uri_protocol(self):
        return self.raw.uri_protocol()

    def announce_send_span(self, tinit, timeout, geometry, synchronised):
        assert issubclass(geometry.__class__, Geometry)
        geometry.configure(self.config)
        self.raw.announce_send_span(tinit, timeout, geometry.raw, synchronised)

    def announce_recv_span(self, tinit, timeout, geometry, synchronised):
        assert issubclass(geometry.__class__, Geometry)
        geometry.configure(self.config)
        self.raw.announce_recv_span(tinit, timeout, geometry.raw, synchronised)

    def announce_send_disable(self, synchronised=False):
        self.raw.announce_send_disable(synchronised)

    def announce_recv_disable(self, synchronised=False):
        self.raw.announce_recv_disable(synchronised)

    def update_smart_send(self, t1):
        self.raw.update_smart_send(t1)

    def barrier_ss_send(self):
        self.raw.barrier_ss_send()

    def barrier_ss_recv(self):
        self.raw.barrier_ss_recv()

    def assign(self, tag, val):
        data_type = map_type[self._get_tag_type(tag)]
        assign = getattr(self.raw, "assign_" + ALLOWED_IO_TYPES[data_type])
        assign(tag, safe_cast(data_type, val))

    def _get_fetch_args(self, fname_root, tag, data_type, spatial_sampler, temporal_sampler):
        assert issubclass(spatial_sampler.__class__, Sampler)
        assert issubclass(temporal_sampler.__class__, TemporalSampler)
        ss = None
        cs = None
        rehash_fetch = False
        try:
            ss = self._tags_spatial_samplers[tag][spatial_sampler.signature]
        except KeyError:
            ss = copy.copy(spatial_sampler)
            ss.configure(self.config, data_type)
            self._tags_spatial_samplers[tag][ss.signature] = ss
            rehash_fetch = True

        try:
            cs = self._tags_temporal_samplers[tag][temporal_sampler.signature]
        except KeyError:
            cs = copy.copy(temporal_sampler)
            cs.configure(self.config, data_type, onlycheck=True)
            self._tags_temporal_samplers[tag][cs.signature] = cs
            rehash_fetch = True
        if rehash_fetch:
            self._tags_fetch[tag][("fetch", cs.signature, ss.signature)] = \
                "fetch_{}_{}_{}".format(ALLOWED_IO_TYPES[data_type],
                                        ss.fetch_signature(),
                                        cs.fetch_signature())
            self._tags_fetch[tag][("fetch_many", cs.signature, ss.signature)] = \
                "fetch_many_{}_{}_{}".format(ALLOWED_IO_TYPES[data_type],
                                             ss.fetch_signature(),
                                             cs.fetch_signature())
            self._tags_fetch[tag][("fetch_dual", cs.signature, ss.signature)] = \
                    "{}_{}_{}_{}".format("fetch_dual",
                                         ALLOWED_IO_TYPES[data_type],
                                         ss.fetch_signature(),
                                         cs.fetch_signature())
            self._tags_fetch[tag][("fetch_many_dual", cs.signature, ss.signature)] = \
                "{}_{}_{}_{}".format("fetch_many_dual",
                                     ALLOWED_IO_TYPES[data_type],
                                     ss.fetch_signature(),
                                     cs.fetch_signature())
        return self._tags_fetch[tag][(fname_root, cs.signature, ss.signature)], ss, cs

    def _get_fetch_algo_args(self, fname_root, tag, data_type, spatial_sampler, temporal_sampler, coupling_algorithm):
        assert issubclass(spatial_sampler.__class__, Sampler)
        assert issubclass(temporal_sampler.__class__, TemporalSampler)
        assert issubclass(coupling_algorithm.__class__, Algorithm)
        ss = None
        cs = None
        ca = None
        rehash_fetch = False
        try:
            ss = self._tags_spatial_samplers[tag][spatial_sampler.signature]
        except KeyError:
            ss = copy.copy(spatial_sampler)
            ss.configure(self.config, data_type)
            self._tags_spatial_samplers[tag][ss.signature] = ss
            rehash_fetch = True

        try:
            cs = self._tags_temporal_samplers[tag][temporal_sampler.signature]
        except KeyError:
            cs = copy.copy(temporal_sampler)
            cs.configure(self.config, data_type, onlycheck=True)
            self._tags_temporal_samplers[tag][cs.signature] = cs
            rehash_fetch = True

        try:
            ca = self._tags_coupling_algorithms[tag][coupling_algorithm.signature]
        except KeyError:
            ca = copy.copy(coupling_algorithm)
            ca.configure(self.config, data_type, onlycheck=True)
            self._tags_coupling_algorithms[tag][ca.signature] = ca
            rehash_fetch = True
        if rehash_fetch:
            self._tags_fetch[tag][("fetch_algorithm", cs.signature, ss.signature, ca.signature)] = \
                "fetch_{}_{}_{}_{}".format(ALLOWED_IO_TYPES[data_type],
                                        ss.fetch_signature(),
                                        cs.fetch_signature(),
                                        ca.fetch_signature())
            self._tags_fetch[tag][("fetch_algorithm_many", cs.signature, ss.signature, ca.signature)] = \
                "fetch_many_{}_{}_{}_{}".format(ALLOWED_IO_TYPES[data_type],
                                             ss.fetch_signature(),
                                             cs.fetch_signature(),
                                             ca.fetch_signature())
            self._tags_fetch[tag][("fetch_algorithm_dual", cs.signature, ss.signature, ca.signature)] = \
                    "{}_{}_{}_{}_{}".format("fetch_dual",
                                         ALLOWED_IO_TYPES[data_type],
                                         ss.fetch_signature(),
                                         cs.fetch_signature(),
                                         ca.fetch_signature())
            self._tags_fetch[tag][("fetch_algorithm_many_dual", cs.signature, ss.signature, ca.signature)] = \
                "{}_{}_{}_{}_{}".format("fetch_many_dual",
                                     ALLOWED_IO_TYPES[data_type],
                                     ss.fetch_signature(),
                                     cs.fetch_signature(),
                                     ca.fetch_signature())
        return self._tags_fetch[tag][(fname_root, cs.signature, ss.signature, ca.signature)], ss, cs, ca

    def _get_fetch_points_values_args(self, fname_root, tag, data_type, temporal_sampler):
        assert issubclass(temporal_sampler.__class__, TemporalSampler)
        cs = None
        rehash_fetch = False

        try:
            cs = self._tags_temporal_samplers[tag][temporal_sampler.signature]
        except KeyError:
            cs = copy.copy(temporal_sampler)
            cs.configure(self.config, data_type, onlycheck=True)
            self._tags_temporal_samplers[tag][cs.signature] = cs
            rehash_fetch = True
        if rehash_fetch:
            self._tags_fetch_points[tag][("fetch_points", cs.signature)] = \
                "fetch_points_{}_{}".format(ALLOWED_IO_TYPES[data_type],
                                        cs.fetch_signature())
            self._tags_fetch_points[tag][("fetch_points_dual", cs.signature)] = \
                "fetch_points_dual_{}_{}".format(ALLOWED_IO_TYPES[data_type],
                                        cs.fetch_signature())
            self._tags_fetch_points[tag][("fetch_values", cs.signature)] = \
                "fetch_values_{}_{}".format(ALLOWED_IO_TYPES[data_type],
                                        cs.fetch_signature())
            self._tags_fetch_points[tag][("fetch_values_dual", cs.signature)] = \
                "fetch_values_dual_{}_{}".format(ALLOWED_IO_TYPES[data_type],
                                        cs.fetch_signature())
        return self._tags_fetch_points[tag][(fname_root, cs.signature)], cs

    def fetch_points(self, *args, **kwargs):
        tag = args[0]
        time = args[1]
        if len(args) == 5:
            temporal_sampler = args[2]
            barrier_enabled = args[3]
            test_value = args[4]
            fetch_Pname, cs = self._get_fetch_points_values_args("fetch_points", tag, type(test_value), temporal_sampler)
            fetch_points = getattr(self.raw, fetch_Pname)
            return fetch_points(tag, time, cs.raw, barrier_enabled, test_value)
        elif len(args) == 6:
            iterator = args[2]
            temporal_sampler = args[3]
            barrier_enabled = args[4]
            test_value = args[5]
            fetch_Pname, cs = self._get_fetch_points_values_args("fetch_points_dual", tag, type(test_value), temporal_sampler)
            fetch_points = getattr(self.raw, fetch_Pname)
            return fetch_points(tag, time, iterator, cs.raw, barrier_enabled, test_value)
        else:
            raise Exception("MUI Error [mui4py.py]: fetch_points() can only take 5 or 6 args")

    def fetch_values(self, *args, **kwargs):
        tag = args[0]
        time = args[1]
        if len(args) == 5:
            temporal_sampler = args[2]
            barrier_enabled = args[3]
            test_value = args[4]
            fetch_Pname, cs = self._get_fetch_points_values_args("fetch_values", tag, type(test_value), temporal_sampler)
            fetch_values = getattr(self.raw, fetch_Pname)
            return fetch_values(tag, time, cs.raw, barrier_enabled)
        elif len(args) == 6:
            iterator = args[2]
            temporal_sampler = args[3]
            barrier_enabled = args[4]
            test_value = args[5]
            fetch_Pname, cs = self._get_fetch_points_values_args("fetch_values_dual", tag, type(test_value), temporal_sampler)
            fetch_values = getattr(self.raw, fetch_Pname)
            return fetch_values(tag, time, iterator, cs.raw, barrier_enabled)
        else:
            raise Exception("MUI Error [mui4py.py]: fetch_values() can only take 5 or 6 args")

    def fetch_many(self, *args, **kwargs):
        tag = args[0]
        points = args[1]
        time = args[2]
        if len(args) == 5:
            spatial_sampler = args[3]
            temporal_sampler = args[4]
            fetch_fname, ss, cs = self._get_fetch_args("fetch_many", tag, points.dtype.type,
                                                        spatial_sampler, temporal_sampler)
            barrier_enabled = True
            fetch = getattr(self.raw, fetch_fname)
            return fetch(tag, points, time, ss.raw, cs.raw, barrier_enabled)
        elif len(args) == 6:
            if isinstance(args[3], (FLOAT, UINT, UINT32, UINT64, INT, INT32, INT64, FLOAT32, FLOAT64)):
              time2 = args[3]
              spatial_sampler = args[4]
              temporal_sampler = args[5]
              fetch_fname, ss, cs = self._get_fetch_args("fetch_many_dual", tag, points.dtype.type,
                                                    spatial_sampler, temporal_sampler)
              barrier_enabled = True
              fetch = getattr(self.raw, fetch_fname)
              return fetch(tag, points, time, time2, ss.raw, cs.raw, barrier_enabled)
            else:
              spatial_sampler = args[3]
              temporal_sampler = args[4]
              coupling_algorithm = args[5]
              fetch_fname, ss, cs, ca = self._get_fetch_algo_args("fetch_algorithm_many", tag, points.dtype.type, spatial_sampler, temporal_sampler, coupling_algorithm)
              if(isinstance(ca, AlgorithmAitken)):
                  self._aitken_under_relaxation_factor.append(((time, 0), ca.get_under_relaxation_factor(time)))
                  self._aitken_residual_l2_norm.append(((time, 0), ca.get_residual_L2_Norm(time)))
                  self._latest_aitken_under_relaxation_factor = ca.get_under_relaxation_factor(time)
                  self._latest_aitken_residual_l2_norm = ca.get_residual_L2_Norm(time)
              elif(isinstance(ca, AlgorithmFixedRelaxation)):
                  self._fixed_relaxation_under_relaxation_factor.append(((time, 0), ca.get_under_relaxation_factor(time)))
                  self._fixed_relaxation_residual_l2_norm.append(((time, 0), ca.get_residual_L2_Norm(time)))
                  self._latest_fixed_relaxation_under_relaxation_factor = ca.get_under_relaxation_factor(time)
                  self._latest_fixed_relaxation_residual_l2_norm = ca.get_residual_L2_Norm(time)
              barrier_enabled = True
              fetch = getattr(self.raw, fetch_fname)
              return fetch(tag, points, time, ss.raw, cs.raw, ca.raw, barrier_enabled)
        elif len(args) == 7:
          time2 = args[3]
          spatial_sampler = args[4]
          temporal_sampler = args[5]
          coupling_algorithm = args[6]
          fetch_fname, ss, cs, ca = self._get_fetch_algo_args("fetch_algorithm_many_dual", tag, points.dtype.type, spatial_sampler, temporal_sampler, coupling_algorithm)
          if(isinstance(ca, AlgorithmAitken)):
              self._aitken_under_relaxation_factor.append(((time, time2), ca.get_under_relaxation_factor(time, time2)))
              self._aitken_residual_l2_norm.append(((time, time2), ca.get_residual_L2_Norm(time, time2)))
              self._latest_aitken_under_relaxation_factor = ca.get_under_relaxation_factor(time, time2)
              self._latest_aitken_residual_l2_norm = ca.get_residual_L2_Norm(time, time2)
          elif(isinstance(ca, AlgorithmFixedRelaxation)):
              self._fixed_relaxation_under_relaxation_factor.append(((time, time2), ca.get_under_relaxation_factor(time, time2)))
              self._fixed_relaxation_residual_l2_norm.append(((time, time2), ca.get_residual_L2_Norm(time, time2)))
              self._latest_fixed_relaxation_under_relaxation_factor = ca.get_under_relaxation_factor(time, time2)
              self._latest_fixed_relaxation_residual_l2_norm = ca.get_residual_L2_Norm(time, time2)
          barrier_enabled = True
          fetch = getattr(self.raw, fetch_fname)
          return fetch(tag, points, time, time2, ss.raw, cs.raw, ca.raw, barrier_enabled)
        else:
            raise Exception("MUI Error [mui4py.py]: fetch_many() can only take 5, 6 or 7 args")

    def fetch(self, *args, **kwargs):
        tag = args[0]
        data_type = map_type[self._get_tag_type(tag)]
        if len(args) == 1:
            fetch_fname = "fetch_single_" + ALLOWED_IO_TYPES[data_type]
            fargs = (tag,)
        elif len(args) == 5:
            loc = array2Point(args[1], self.config, self.raw_point)
            time = args[2]
            spatial_sampler = args[3]
            temporal_sampler = args[4]
            fetch_fname, ss, cs = self._get_fetch_args("fetch", tag, data_type, spatial_sampler, temporal_sampler)
            barrier_enabled = True
            if type(time).__name__ == 'float':
                barrier_time = mui4py_mod.numeric_limits_real
            elif type(time).__name__ == 'int':
                barrier_time = mui4py_mod.numeric_limits_int
            else:
                raise Exception("Unrecognized time type '{}'.".format(type(time).__name__))
            fargs = (tag, loc, time, ss.raw, cs.raw, barrier_enabled)
        elif len(args) == 6:
            if isinstance(args[3], (FLOAT, UINT, UINT32, UINT64, INT, INT32, INT64, FLOAT32, FLOAT64)):
                loc = array2Point(args[1], self.config, self.raw_point)
                time1 = args[2]
                time2 = args[3]
                spatial_sampler = args[4]
                temporal_sampler = args[5]
                fetch_fname, ss, cs = self._get_fetch_args("fetch_dual", tag, data_type, spatial_sampler, temporal_sampler)
                barrier_enabled = True
                if type(time1).__name__ == 'float':
                    barrier_time = mui4py_mod.numeric_limits_real
                elif type(time1).__name__ == 'int':
                    barrier_time = mui4py_mod.numeric_limits_int
                else:
                    raise Exception("Unrecognized time1 type '{}'.".format(type(time1).__name__))
                fargs = (tag, loc, time1, time2, ss.raw, cs.raw, barrier_enabled)
            else:
                loc = array2Point(args[1], self.config, self.raw_point)
                time1 = args[2]
                spatial_sampler = args[3]
                temporal_sampler = args[4]
                coupling_algorithm = args[5]
                fetch_fname, ss, cs, ca = self._get_fetch_algo_args("fetch_algorithm", tag, data_type, spatial_sampler, temporal_sampler, coupling_algorithm)
                if(isinstance(ca, AlgorithmAitken)):
                    self._aitken_under_relaxation_factor.append(((time1, 0), ca.get_under_relaxation_factor(time1)))
                    self._aitken_residual_l2_norm.append(((time1, 0), ca.get_residual_L2_Norm(time1)))
                    self._latest_aitken_under_relaxation_factor = ca.get_under_relaxation_factor(time1)
                    self._latest_aitken_residual_l2_norm = ca.get_residual_L2_Norm(time1)
                if(isinstance(ca, AlgorithmFixedRelaxation)):
                    self._fixed_relaxation_under_relaxation_factor.append(((time1, 0), ca.get_under_relaxation_factor(time1)))
                    self._fixed_relaxation_residual_l2_norm.append(((time1, 0), ca.get_residual_L2_Norm(time1)))
                    self._latest_fixed_relaxation_under_relaxation_factor = ca.get_under_relaxation_factor(time1)
                    self._latest_fixed_relaxation_residual_l2_norm = ca.get_residual_L2_Norm(time1)
                barrier_enabled = True
                if type(time1).__name__ == 'float':
                    barrier_time = mui4py_mod.numeric_limits_real
                elif type(time1).__name__ == 'int':
                    barrier_time = mui4py_mod.numeric_limits_int
                else:
                    raise Exception("Unrecognized time1 type '{}'.".format(type(time1).__name__))
                fargs = (tag, loc, time1, ss.raw, cs.raw, ca.raw, barrier_enabled)
        elif len(args) == 7:
            loc = array2Point(args[1], self.config, self.raw_point)
            time1 = args[2]
            time2 = args[3]
            spatial_sampler = args[4]
            temporal_sampler = args[5]
            coupling_algorithm = args[6]
            fetch_fname, ss, cs, ca = self._get_fetch_algo_args("fetch_algorithm_dual", tag, data_type, spatial_sampler, temporal_sampler, coupling_algorithm)
            if(isinstance(ca, AlgorithmAitken)):
                self._aitken_under_relaxation_factor.append(((time1, time2), ca.get_under_relaxation_factor(time1, time2)))
                self._aitken_residual_l2_norm.append(((time1, time2), ca.get_residual_L2_Norm(time1, time2)))
                self._latest_aitken_under_relaxation_factor = ca.get_under_relaxation_factor(time1, time2)
                self._latest_aitken_residual_l2_norm = ca.get_residual_L2_Norm(time1, time2)
            if(isinstance(ca, AlgorithmFixedRelaxation)):
                self._fixed_relaxation_under_relaxation_factor.append(((time1, time2), ca.get_under_relaxation_factor(time1, time2)))
                self._fixed_relaxation_residual_l2_norm.append(((time1, time2), ca.get_residual_L2_Norm(time1, time2)))
                self._latest_fixed_relaxation_under_relaxation_factor = ca.get_under_relaxation_factor(time1, time2)
                self._latest_fixed_relaxation_residual_l2_norm = ca.get_residual_L2_Norm(time1, time2)
            barrier_enabled = True
            if type(time1).__name__ == 'float':
                barrier_time = mui4py_mod.numeric_limits_real
            elif type(time1).__name__ == 'int':
                barrier_time = mui4py_mod.numeric_limits_int
            else:
                raise Exception("Unrecognized time1 type '{}'.".format(type(time1).__name__))
            fargs = (tag, loc, time1, time2, ss.raw, cs.raw, ca.raw, barrier_enabled)
        else:
            raise Exception("MUI Error [mui4py.py]: fetch() can only take 1, 5, 6 or 7 args")
        fetch = getattr(self.raw, fetch_fname)
        return safe_cast(self._get_tag_type(tag), fetch(*fargs))

    def Point(self, points):
        return array2Point(points, self.config, self.raw_point)

    def get_aitken_under_relaxation_factor(self, t, iter=None):
        return_value = None
        for element in self._aitken_under_relaxation_factor:
            if element[0][0] == t:
                if iter is not None:
                    if element[0][1] == iter:
                        return_value = element[1]
                        break
                else:
                    return_value = element[1]
                    break
        if return_value is None:
            raise Exception("MUI Error [mui4py.py]: No match found for given t and iter in get_aitken_under_relaxation_factor")
        else:
            return return_value

    def get_aitken_residual_l2_norm(self, t, iter=None):
        return_value = None
        for element in self._aitken_residual_l2_norm:
            if element[0][0] == t:
                if iter is not None:
                    if element[0][1] == iter:
                        return_value = element[1]
                        break
                else:
                    return_value = element[1]
                    break
        if return_value is None:
            raise Exception("MUI Error [mui4py.py]: No match found for given t and iter in get_aitken_residual_l2_norm")
        else:
            return return_value

    def get_latest_aitken_under_relaxation_factor(self):
        return self._latest_aitken_under_relaxation_factor

    def get_latest_aitken_residual_l2_norm(self):
        return self._latest_aitken_residual_l2_norm

    def get_fixed_relaxation_under_relaxation_factor(self, t, iter=None):
        return_value = None
        for element in self._fixed_relaxation_under_relaxation_factor:
            if element[0][0] == t:
                if iter is not None:
                    if element[0][1] == iter:
                        return_value = element[1]
                        break
                else:
                    return_value = element[1]
                    break
        if return_value is None:
            raise Exception("MUI Error [mui4py.py]: No match found for given t and iter in get_fixed_relaxation_under_relaxation_factor")
        else:
            return return_value

    def get_fixed_relaxation_residual_l2_norm(self, t, iter=None):
        return_value = None
        for element in self._fixed_relaxation_residual_l2_norm:
            if element[0][0] == t:
                if iter is not None:
                    if element[0][1] == iter:
                        return_value = element[1]
                        break
                else:
                    return_value = element[1]
                    break
        if return_value is None:
            raise Exception("MUI Error [mui4py.py]: No match found for given t and iter in get_fixed_relaxation_residual_l2_norm")
        else:
            return return_value

    def get_latest_fixed_relaxation_under_relaxation_factor(self):
        return self._latest_fixed_relaxation_under_relaxation_factor

    def get_latest_fixed_relaxation_residual_l2_norm(self):
        return self._latest_fixed_relaxation_residual_l2_norm