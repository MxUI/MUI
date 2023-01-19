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
