"""
##############################################################################
# Multiscale Universal Interface Code Coupling Library                       #
#                                                                            #
# Copyright (C) 2023 C. Richardson, E. R. Fernandez                          #
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
# @file setup.py
# @author C. Richardson, E. R. Fernandez
# @date 20 January 2023
# @brief Setup file for MUI Python wrapper.
#
"""

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
        cfg = 'Debug' if self.debug else 'Release'
        cmake_args = shlex.split(os.environ.get("CMAKE_ARGS", ""))
        cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir + "/mui4py",
                       '-DPython3_EXECUTABLE=' + sys.executable,
                       '-DPython3_INCLUDE_DIRS={sysconfig.get_config_var("INCLUDEPY")}',
                       '-DCMAKE_BUILD_TYPE=' + cfg]
        build_args = ['--config', cfg]

        env = os.environ.copy()
        # default to 2 build threads
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in env:
            env["CMAKE_BUILD_PARALLEL_LEVEL"] = "1"

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
      version='2.0',
      description='Python bindings for MUI coupling library.',
      url='http://mxui.github.io',
      author='Eduardo Ramos Fernandez, Chris Richardson, Wendi Liu',
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
