# Compile settings
CC=g++
CFLAGS=-Wall -Wextra -std=c++98 -g -Ofast
# CFLAGS=-Wall -Wextra -std=c++11 -g -Ofast
LFLAGS=
LIBS=-lm

CFLAGS:=$(CFLAGS) $(shell gsl-config --cflags)
LIBS:=$(LIBS) $(shell gsl-config --libs)

# SRC=$(wildcard *.cpp)
SRC=main.cpp

# Targets
all: release

release: $(SRC)
	$(CC) $(CFLAGS) $^ $(LFLAGS) $(LIBS)


