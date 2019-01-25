#!/bin/bash                                                                                                                                                                                                                                                                                                             
#NOTE: Only detects MPICH and OpenMPI/Spectrum implementations
#TODO: Expand this to other implementations (CRAY, INTEL, etc.)
if type "mpichversion" &> /dev/null; then
    echo $(mpichversion | head -n 1 |  tr -s ' ' | cut -d " " -f 1,3)
elif type "ompi_info" &> /dev/null; then
    echo $(mpicc -showme:version | head -n 1 |  tr -s ' ' | cut -d " " -f 2,3,4)
fi 
