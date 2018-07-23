OBJDIR = $(GARFIELD_HOME)/Object
SRCDIR = $(GARFIELD_HOME)/Source
INCDIR = $(GARFIELD_HOME)/Include
HEEDDIR = $(GARFIELD_HOME)/Heed
LIBDIR = $(GARFIELD_HOME)/Library

# Compiler flags
CFLAGS = -Wall -Wextra -Wno-long-long \
	`root-config --cflags` \
	-O3 -fno-common -c \
	-I$(INCDIR) -I$(HEEDDIR)

# Debug flags
# CFLAGS += -g

LDFLAGS = -L$(LIBDIR) -lGarfield
LDFLAGS += `root-config --glibs` -lGeom -lgfortran -lm
# LDFLAGS += -g

scriptisto: scriptisto.C 
	$(CXX) $(CFLAGS) scriptisto.C
	$(CXX) `root-config --cflags` -o scriptisto scriptisto.o $(LDFLAGS)
	rm scriptisto.o

gem2: gem2.C 
	$(CXX) $(CFLAGS) gem2.C
	$(CXX) `root-config --cflags` -o gem2 gem2.o $(LDFLAGS)
	rm gem2.o
