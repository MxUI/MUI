COMPILER ?= GNU
MPI ?= OPENMPI
ifeq ('$(MPI)','OPENMPI')
	CXX = mpic++
else ifeq ('$(MPI)','INTEL')
	CXX = mpicxx
else ifeq ('$(MPI)','MPICH')
	CXX = mpicxx
endif
CXXFLAGS = -Wall -std=c++11 -O3 -Wno-unused-local-typedefs -Wno-comment -Wno-reorder
LDFLAGS = -shared
USE_RBF ?= 0
include $(COMPILER).make
INC_EIGEN =
INC_MPI4PY = -I$(shell python3 -c "import mpi4py;print(mpi4py.get_include())")
INC_PYBIND = $(shell python3 -m pybind11 --includes)
CXX_VERSION_STR = '"$(shell $(CXX) --version | head -n 1)"'
MPI_VERSION_STR = '"$(shell bash mpi_version.sh)"'
CXX_CONFIG_STR = '"$(shell $(CXX) -show) $(CXXFLAGS) $(INC_MPI4PY) $(INC_PYBIND)"'
PYTHON_INT_SIZE = $(shell python3 -c "import numpy; print(str(numpy.iinfo(int).dtype)[-2:])")
ifeq ($(PYTHON_INT_SIZE),64)
	CXXFLAGS_D=-DPYTHON_INT_64
else ifeq ($(PYTHON_INT_SIZE),32)
	CXXFLAGS_D = -DPYTHON_INT_32
endif
CXXFLAGS_D += -DCOMPILER_VERSION_STR=$(CXX_VERSION_STR) -DMPI_VERSION_STR=$(MPI_VERSION_STR) -DCOMPILER_CONFIG_STR=$(CXX_CONFIG_STR) -DPYTHON_BINDINGS
ifeq ($(USE_RBF),1)
	CXXFLAGS += -fpermissive
	CXXFLAGS_D += -DUSE_RBF
endif
TARGET_LIB = mui4py/mui4py_mod.so
SRCS = mui4py/mui4py.cpp
OBJS = $(SRCS:.cpp=.o)

.PHONY: mui4py_mod
none:
	@echo "No default Makefile option. Choose from: 'mui4py_mod', 'package', 'meta-package', 'clean-all', 'clean-package', 'clean-mod' or 'pip-install'."

mui4py_mod: $(TARGET_LIB)

$(TARGET_LIB): $(OBJS)
	$(CXX) $(LDFLAGS) $^ -o $@

%.o: %.cpp
ifeq ($(USE_RBF),1)
	$(CXX) -c $(CXXFLAGS) $(CXXFLAGS_D) $(INC_MPI4PY) $(INC_PYBIND) -I$(INC_EIGEN) $< -o $@
else
	$(CXX) -c $(CXXFLAGS) $(CXXFLAGS_D) $(INC_MPI4PY) $(INC_PYBIND) $< -o $@
endif

.PHONY: package meta-package clean-package clean-mod clean-all install pip-install clean-obj

package: mui4py_mod
	@python3 setup.py sdist
	rm -rf mui4py.egg-info

meta-package:
	@echo "Not implemented yet"

clean-package:
	rm -rf mui4py.egg-info
	rm -rf dist
	rm -rf build

clean-mod:
	rm -f $(TARGET_LIB)
	rm -rf mui4py/__pycache__

clean-obj:
	rm -f mui4py/*.o

clean-all: clean-package clean-mod clean-obj

install:
	@python3 setup.py install

pip-install:
	@pip3 install dist/mui4py*.tar.gz
