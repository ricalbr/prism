CC=gcc
CXX=g++
RM=rm -f
# CPPFLAGS=-g $(shell root-config --cflags)
CPPFLAGS=-g -fopenmp -O3 -lyaml-cpp -march=native -mavx2 -mfma -finline-limit=1000000 -ffp-contract=fast
# LDFLAGS=-g $(shell root-config --ldflags)
# LDLIBS=$(shell root-config --libs)

SRCS=click.cpp
OBJS=$(subst .cc,.o,$(SRCS))

all: click

click: $(OBJS)
	$(CXX) -o click.o $(OBJS) $(CPPFLAGS)

clean:
	$(RM) $(OBJS)

distclean: clean
	$(RM) click
