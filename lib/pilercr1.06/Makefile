CFLAGS = -O3 -DNDEBUG=1
LDLIBS = -lm -static
# LDLIBS = -lm

OBJ = .o
EXE =

RM = rm -f
CP = cp

GPP = g++
LD = $(GPP) $(CFLAGS)
CPP = $(GPP) -c $(CFLAGS) 

all: pilercr

CPPSRC = $(sort $(wildcard *.cpp))
CPPOBJ	= $(subst .cpp,.o,$(CPPSRC))

$(CPPOBJ): %.o: %.cpp
	$(CPP) $< -o $@

pilercr: $(CPPOBJ)
	$(LD) -o pilercr $(CPPOBJ) $(LDLIBS)
	strip pilercr
