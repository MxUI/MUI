#!/bin/bash

CC	= mpic++
CFLAGS	= -std=c++11 -O3
EIGEN_PATH ?= ./

SCR = $(wildcard *.cpp)
EXE = $(SCR:.cpp=)

default: $(EXE)

% : %.cpp
	$(CC) $(CFLAGS) -I$(EIGEN_PATH) $< -o $@

clean:
	rm -f $(EXE)
