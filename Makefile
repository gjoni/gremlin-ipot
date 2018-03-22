
CXX=g++
CXXFLAGS += -Wall -Wno-unused-result -pedantic -Ofast -march=native -std=c++0x -g -ggdb -fopenmp
CXXLIBS = -lgsl -lgslcblas
INCDIRS = -I./src
LIBDIRS =

OBJDIR = obj
SRCDIR = src

DEPS = $(shell find $(SRCDIR) -name '*.h')
SRCS = $(shell find $(SRCDIR) -name '*.cpp')
OBJS = $(patsubst $(SRCDIR)%.cpp, $(OBJDIR)%.o, $(SRCS))

fo_gr3=obj/gremlin3c.o obj/map_align.o
fo_map=obj/gremlin3.o obj/gremlin3c.o

all: $(OBJDIR) gremlin3

$(OBJDIR):
	mkdir -p $(OBJDIR)

gremlin3: $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(INCDIRS) $(LIBDIRS) -o gremlin3 $(filter-out $(fo_gr3), $(OBJS)) $(CXXLIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCDIRS) -c -o $@ $<

clean:
	rm $(OBJS) gremlin3
	rmdir $(OBJDIR)
