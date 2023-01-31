import sys
import os
import shlex
from setuptools import setup, Extension
import sysconfig
import subprocess
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            _ = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: "
                               + ", ".join(e.name for e in self.extensions))

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = shlex.split(os.environ.get("CMAKE_ARGS", ""))
        cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir + "/mui4py",
                       '-DPython3_EXECUTABLE=' + sys.executable,
                       f'-DPython3_LIBRARIES={sysconfig.get_config_var("LIBDEST")}',
                       f'-DPython3_INCLUDE_DIRS={sysconfig.get_config_var("INCLUDEPY")}']

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]
        cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]

        env = os.environ.copy()
        # default to 3 build threads
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in env:
            env["CMAKE_BUILD_PARALLEL_LEVEL"] = "3"

        import pybind11
        env['pybind11_DIR'] = pybind11.get_cmake_dir()

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp, env=env)


setup(
      ext_modules=[CMakeExtension('mui4py_mod')],
      cmdclass=dict(build_ext=CMakeBuild),
      name='mui4py',
      version='1.2.4',
      description='Python bindings for MUI coupling library.',
      url='http://mxui.github.io',
      author='Eduardo Ramos Fernandez',
      author_email='eduardo.rf159@gmail.com',
      license='GPLv3 or Apache v2.0',
      packages=['mui4py', 'mui4py.cpp'],
      install_requires=[
            'mpi4py',
            'numpy>=1.21',
      ],
      setup_requires=["pybind11"],
      include_package_data=True,
      zip_safe=False)
