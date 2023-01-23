from mui4py.common import CppClass
from mui4py.config import Config
from mui4py.types import INT, FLOAT, FLOAT32, FLOAT64, INT32, INT64, STRING


# Inteface for Sampler and ChronoSampler
def sampler_fetch_signature(self):
    sig = self._split_class_name(title=False)
    sig = sig.replace("_sampler", "")
    return sig.replace("sampler_", "")


# Spatial samplers
class Sampler(CppClass):
    def __init__(self, args=(), kwargs={}):
        # Empty config to not trigger error in default config.
        super(Sampler, self).__init__(Config(), args=args, kwargs=kwargs)


Sampler.fetch_signature = sampler_fetch_signature


class SamplerExact(Sampler):
    def __init__(self, tol=None):
        super(SamplerExact, self).__init__(kwargs={"tol": tol})
        self._ALLOWED_IO_TYPES = [FLOAT, INT, INT32, INT64, FLOAT32, FLOAT64, STRING]


class SamplerGauss(Sampler):
    def __init__(self, r, h):
        super(SamplerGauss, self).__init__(args=(r, h))
        self._ALLOWED_IO_TYPES = [INT32, INT64, FLOAT32, FLOAT64]


class SamplerMovingAverage(Sampler):
    def __init__(self, bbox):
        super(SamplerMovingAverage, self).__init__(args=(_Point(bbox), ))
        self._ALLOWED_IO_TYPES = [INT32, INT64, FLOAT32, FLOAT64]


class SamplerNearestNeighbor(Sampler):
    def __init__(self, r, h):
        super(SamplerNearestNeighbor, self).__init__(args=(r, h))
        self._ALLOWED_IO_TYPES = [INT32, INT64, FLOAT32, FLOAT64]


class SamplerPseudoNearest2Linear(Sampler):
    def __init__(self, h):
        super(SamplerPseudoNearest2Linear, self).__init__(args=(h, ))
        self._ALLOWED_IO_TYPES = [INT32, INT64, FLOAT32, FLOAT64]


class SamplerPseudoNearestNeighbor(Sampler):
    def __init__(self, h):
        super(SamplerPseudoNearestNeighbor, self).__init__(args=(h,))
        self._ALLOWED_IO_TYPES = [INT32, INT64, FLOAT32, FLOAT64]


class SamplerSherpardQuintic(Sampler):
    def __init__(self, r):
        super(SamplerSherpardQuintic, self).__init__(args=(r,))
        self._ALLOWED_IO_TYPES = [INT32, INT64, FLOAT32, FLOAT64]


class SamplerSphQuintic(Sampler):
    def __init__(self, r):
        super(SamplerSphQuintic, self).__init__(args=(r,))
        self._ALLOWED_IO_TYPES = [INT32, INT64, FLOAT32, FLOAT64]


class SamplerSumQuintic(Sampler):
    def __init__(self, r):
        super(SamplerSumQuintic, self).__init__(args=(r,))
        self._ALLOWED_IO_TYPES = [INT32, INT64, FLOAT32, FLOAT64]


class SamplerRbf(Sampler, CppClass):
    def __init__(self, r, pointvect, basisFunc, conservative, polynomial,
                 smoothFunc, readMatrix, writeMatrix, fileAddress, cutoff, cgSolveTol, cgMaxIter, pouSize):
        super(SamplerRbf, self).__init__(args=(r, pointvect, basisFunc,
                                         conservative, polynomial, smoothFunc,
                                         readMatrix, writeMatrix, fileAddress,
                                         cutoff, cgSolveTol, cgMaxIter, pouSize, ))
        self._ALLOWED_IO_TYPES = [INT32, INT64, FLOAT32, FLOAT64]


# Chrono samplers
class ChronoSampler(CppClass):
    def __init__(self, args=(), kwargs={}):
        # Empty config to not trigger error in default config.
        super(ChronoSampler, self).__init__(Config(), args, kwargs)


ChronoSampler.fetch_signature = sampler_fetch_signature


class ChronoSamplerExact(ChronoSampler):
    def __init__(self, tol=None):
        super(ChronoSamplerExact, self).__init__(kwargs={"tol": tol})
        self._ALLOWED_IO_TYPES = [INT, INT32, INT64, FLOAT, FLOAT32, FLOAT64, STRING]


class ChronoSamplerGauss(ChronoSampler):
    def __init__(self, cutoff, sigma):
        super(ChronoSamplerGauss, self).__init__(args=(cutoff, sigma))
        self._ALLOWED_IO_TYPES = [INT, INT32, INT64, FLOAT, FLOAT32, FLOAT64]


class ChronoSamplerMean(ChronoSampler):
    def __init__(self, newleft=None, newright=None):
        super(ChronoSamplerMean, self).__init__(kwargs={"newleft": newleft, "newright": newright})
        self._ALLOWED_IO_TYPES = [INT, INT32, INT64, FLOAT, FLOAT32, FLOAT64]


class ChronoSamplerSum(ChronoSampler):
    def __init__(self, newleft=None, newright=None):
        super(ChronoSamplerSum, self).__init__(kwargs={"newleft": newleft, "newright": newright})
        self._ALLOWED_IO_TYPES = [INT, INT32, INT64, FLOAT, FLOAT32, FLOAT64]
