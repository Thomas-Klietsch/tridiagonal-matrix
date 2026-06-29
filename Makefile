# Simple makefile

CC := g++ -std=c++11 -O3

CCW := -Wall -Werror -Wextra

SOURCES := ./src/main.cpp

.DEFAULT_GOAL := compile

compile:
	mkdir -p bin
	$(CC) $(CCW) -o ./bin/main $(SOURCES)

clean:
	rm -rf ./bin/main

all: clean main
