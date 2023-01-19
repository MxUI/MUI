from mui4py.types import map_type
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
        self.float_type = map_type[float_type]
        self.dim = dim
        self._check_types()
        self.force_casting = force_casting

    def _check_types(self):
        pass
