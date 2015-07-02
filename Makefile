CC = g++
CFLAGS = -Wall -std=c++11 -fdiagnostics-color=auto
CFLAGS+= -Wextra -g #debug flags
CFLAGS+=$(shell pkg-config --cflags eigen3)

LIBS =-lboost_program_options
LIBS+=$(shell pkg-config --libs eigen3)

CPP_FILES:=$(wildcard src/*.cpp)
OBJ_FILES:=$(addprefix obj/, $(notdir $(CPP_FILES:.cpp=.o)))

.PHONY:all
all:emrg

emrg:$(OBJ_FILES)
	$(CC) $(LIBS) -o $@ $^

obj/%.o: src/%.cpp | obj
	$(CC) $(CFLAGS) -c -o $@ $<

obj:
	mkdir -p $@

.PHONY:clean
clean:
	rm -rf emrg obj *.dat
