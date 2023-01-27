#!/bin/bash

CC	= mpic++
CFLAGS	= -std=c++11 -O3
INC_EIGEN =-I../eigen-3.4.0

SCR = $(wildcard *.cpp)
EXE = $(SCR:.cpp=)

default: $(EXE)

% : %.cpp
	$(CC) $(CFLAGS) $(INC_EIGEN) $< -o $@

clean:
	rm -f $(EXE)
