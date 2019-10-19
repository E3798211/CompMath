# Compile settings
CC=gcc
CFLAGS=-Wall -Wextra -std=c99 -Ofast
LFLAGS=
LIBS=-lm

SRC=$(wildcard *.c)

# Targets
all: release

release: $(SRC)
	$(CC) $(CFLAGS) $^ $(LFLAGS) $(LIBS)


