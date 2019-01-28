import sys
import os
import pybind11
import mpi4py
from setuptools import setup, Extension
import setuptools.command.build_py
import sysconfig
import subprocess

#TODO: Maybe produce a meta-package which download the library, mui4py branch and produce the binary package with some
#      extra configuration like MPI compiler and flags.

# os.environ["CC"] = "mpic++"
# os.environ["LD"] = "mpic++"
# extra_compile_args = sysconfig.get_config_var('CFLAGS').split()
# extra_compile_args += ["-Wall", "-std=c++11", "-O3"]
# extra_link_args = sysconfig.get_config_var('LDFLAGS').split()
# extra_link_args += ["-Wl,-undefined dynamic_lookup"]
# includedir_mpi4py = os.path.dirname(sys.modules['mpi4py'].__file__)
# includedir_mpi4py = os.path.join(includedir_mpi4py, "include")
# includedir_pybind = pybind11.get_include()
# mui4py_mod = Extension('mui4py_mod',
# 		    # Do this for MAC/LINUX/compiiler
#                     define_macros = [('MAJOR_VERSION', '1'),
#                                      ('MINOR_VERSION', '0')],
#                     include_dirs = [includedir_mpi4py, includedir_pybind],
#                     extra_compile_args = extra_compile_args,
#                     extra_link_args = extra_link_args,
#                     sources = ['mui4py/mui4py.cpp'])
#
# class Build(setuptools.command.build_py.build_py):
#     """Customized setuptools build command - builds mui_mod on build."""
#     def run(self):
#         protoc_command = ["make", "mui4py_mod"]
#         if subprocess.call(protoc_command) != 0:
#             sys.exit(-1)
#         setuptools.command.build_py.build_py.run(self)
#
setup(
      #  cmdclass={
	#        'build_py': Build,
      # },
      # ext_modules = [mui4py_mod],
      name='mui4py',
      version='0.1',
      description='Python bindings for MUI coupling library.',
      url='',
      author='Eduardo Ramos Fernandez',
      author_email='eduardo.rf159@gmail.com',
      license='Apache v2',
      packages=['mui4py'],
      install_requires=[
            'mpi4py',
            'numpy',
      ],
      include_package_data=True,
      zip_safe=False)
