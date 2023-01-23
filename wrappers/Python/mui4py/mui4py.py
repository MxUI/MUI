import mui4py.mui4py_mod as mui4py_mod
from mui4py.common import CppClass, get_cpp_name, array2Point
from mui4py.config import Config
from mui4py.types import safe_cast, map_type, ALLOWED_IO_TYPES
from mui4py.geometry import Geometry
from mui4py.samplers import Sampler, ChronoSampler
import copy


# TODO: Add communicator parameter
def create_unifaces(domain, ifaces_names, config):
    assert type(ifaces_names) == list
    assert type(domain) == str
    assert issubclass(config.__class__, Config)
    ifaces_out = {}
    cpp_obj_name = get_cpp_name("create_uniface", config.dim,
                                config.float_type, config.int_type)
    ifaceraw = getattr(mui4py_mod, cpp_obj_name)(domain, ifaces_names)
    for i, obj in enumerate(ifaceraw):
        ifaces_out[ifaces_names[i]] = Uniface(config=config,
                                              cpp_obj=obj)
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
        self._tags_chrono_samplers = {}
        self._tags_fetch = {}
        # self._tags_push = {}
        self._ALLOWED_PROTOCOLS = ["mpi"]

    def _check_uri(self):
        # protocol://domain/interface
        # if not re.match("[mpi]+//[
        pass

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
            self._tags_chrono_samplers[tag] = {}
            self._tags_fetch[tag] = {}
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
            raise Exception("Push function accept 1, 2 or 3 parameters.")

        push = getattr(self.raw, push_fname)
        push(*pargs)

    def push_many(self, tag, points, values):
        # TODO: Try to apply safe_cast
        push_fname, data_type = self._get_pushfname("push_many_", tag, type_in=values.dtype.type)
        getattr(self.raw, push_fname)(tag, points, values)

    def commit(self, tstamp):
        return self.raw.commit(tstamp)

    def barrier(self, t1, t2=None):
        if t2 is not None:
            self.raw.barrier(t1, t2)
        else:
            self.raw.barrier(t1)

    def forget(self, tend, tbegin=0.0):
        self.raw.forget(tend, True)

    def forecast(self, timestamp):
        self.raw.forecast(timestamp)

    def is_ready(self, attr, t1, t2=None):
        if t2 is not None:
            self.raw.is_ready(attr, t1, t2)
        else:
            self.raw.is_ready(attr, t1)

    def set_memory(self, t):
        self.raw.set_memmory(t)

    def assign(self, tag, val):
        data_type = map_type[self._get_tag_type(tag)]
        assign = getattr(self.raw, "assign_" + ALLOWED_IO_TYPES[data_type])
        assign(tag, safe_cast(data_type, val))

    def announce_recv_span(self, tinit, timeout, geometry, synchronised):
        assert issubclass(geometry.__class__, Geometry)
        geometry.configure(self.config)
        self.raw.announce_recv_span(tinit, timeout, geometry.raw, synchronised)

    def announce_send_span(self, tinit, timeout, geometry, synchronised):
        assert issubclass(geometry.__class__, Geometry)
        geometry.configure(self.config)
        self.raw.announce_send_span(tinit, timeout, geometry.raw, synchronised)

    def announce_recv_disable(self):
        self.raw.announce_recv_disable()

    def announce_send_disable(self):
        self.raw.announce_send_disable()

    def _get_fetch_5args(self, fname_root, tag, data_type, spatial_sampler, chrono_sampler):
        assert issubclass(spatial_sampler.__class__, Sampler)
        assert issubclass(chrono_sampler.__class__, ChronoSampler)
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
            cs = self._tags_chrono_samplers[tag][chrono_sampler.signature]
        except KeyError:
            cs = copy.copy(chrono_sampler)
            cs.configure(self.config, data_type, onlycheck=True)
            self._tags_chrono_samplers[tag][cs.signature] = cs
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
        return self._tags_fetch[tag][(fname_root, cs.signature, ss.signature)], ss, cs

    def _get_fetch_6args(self, fname_root, tag, data_type, spatial_sampler, chrono_sampler):
        assert issubclass(spatial_sampler.__class__, Sampler)
        assert issubclass(chrono_sampler.__class__, ChronoSampler)
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
            cs = self._tags_chrono_samplers[tag][chrono_sampler.signature]
        except KeyError:
            cs = copy.copy(chrono_sampler)
            cs.configure(self.config, data_type, onlycheck=True)
            self._tags_chrono_samplers[tag][cs.signature] = cs
            rehash_fetch = True
        if rehash_fetch:
            self._tags_fetch[tag][("fetch6", cs.signature, ss.signature)] = \
                    "{}_{}_{}_{}".format("fetch",
                                         ALLOWED_IO_TYPES[data_type],
                                         ss.fetch_signature(),
                                         cs.fetch_signature())
            self._tags_fetch[tag][("fetch_many6", cs.signature, ss.signature)] = \
                "{}_{}_{}_{}".format("fetch_many6",
                                     ALLOWED_IO_TYPES[data_type],
                                     ss.fetch_signature(),
                                     cs.fetch_signature())

        return self._tags_fetch[tag][(fname_root, cs.signature, ss.signature)], ss, cs

    def fetch_points(self, tag, time):
        data_type = map_type[self._get_tag_type(tag)]
        fetch_points = getattr(self.raw, "fetch_points_" + ALLOWED_IO_TYPES[data_type])
        return fetch_points(tag, time)

    def fetch_many(self, tag, points, time, spatial_sampler, chrono_sampler):
        fetch_fname, ss, cs = self._get_fetch_5args("fetch_many", tag, points.dtype.type,
                                                    spatial_sampler, chrono_sampler)
        fetch = getattr(self.raw, fetch_fname)
        return fetch(tag, points, time, ss.raw, cs.raw)

    def fetch_many6(self, tag, points, time1, time2, spatial_sampler, chrono_sampler):
        fetch_fname, ss, cs = self._get_fetch_6args("fetch_many6", tag, points.dtype.type,
                                                    spatial_sampler, chrono_sampler)
        fetch = getattr(self.raw, fetch_fname)
        return fetch(tag, points, time1, time2, ss.raw, cs.raw)

    def fetch(self, *args, **kwargs):
        tag = args[0]
        data_type = map_type[self._get_tag_type(tag)]
        if len(args) == 1:
            fetch_fname = "fetch_" + ALLOWED_IO_TYPES[data_type]
            fargs = (tag, )
        if len(args) == 5:
            loc = array2Point(args[1], self.config, self.raw_point)
            time = args[2]
            spatial_sampler = args[3]
            chrono_sampler = args[4]
            fetch_fname, ss, cs = self._get_fetch_5args("fetch", tag, data_type, spatial_sampler, chrono_sampler)
            barrier_enabled = True
            if type(time).__name__ == 'float':
                barrier_time = mui4py_mod.numeric_limits_real
            elif type(time).__name__ == 'int':
                barrier_time = mui4py_mod.numeric_limits_int
            else:
                raise Exception("Unrecognized time type '{}'.".format(type(time).__name__))
            fargs = (tag, loc, time, ss.raw, cs.raw, barrier_enabled)
        if len(args) == 6:
            loc = array2Point(args[1], self.config, self.raw_point)
            time1 = args[2]
            time2 = args[3]
            spatial_sampler = args[4]
            chrono_sampler = args[5]
            fetch_fname, ss, cs = self._get_fetch_6args("fetch", tag, data_type, spatial_sampler, chrono_sampler)
            barrier_enabled = True
            if type(time1).__name__ == 'float':
                barrier_time = mui4py_mod.numeric_limits_real
            elif type(time1).__name__ == 'int':
                barrier_time = mui4py_mod.numeric_limits_int
            else:
                raise Exception("Unrecognized time1 type '{}'.".format(type(time1).__name__))
            if type(time1).__name__ != type(time2).__name__:
                raise Exception("time1 type '{}'. doesn't same as time2 type".format(type(time1).__name__))
            fargs = (tag, loc, time1, time2, ss.raw, cs.raw, barrier_enabled)
        fetch = getattr(self.raw, fetch_fname)
        return safe_cast(self._get_tag_type(tag), fetch(*fargs))

    def Point(self, points):
        return array2Point(points, self.config, self.raw_point)
