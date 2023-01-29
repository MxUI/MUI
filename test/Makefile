#!/bin/bash

CC	= mpic++
CFLAGS	= -std=c++11 -O3

SCR = $(wildcard *.cpp)
EXE = $(SCR:.cpp=)

default: $(EXE)

% : %.cpp
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(EXE) *.csv
