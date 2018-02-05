
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

#objgr3=obj/gremlin3.o
#objgr3c=obj/gremlin3c.o
#objmap=obj/map_align.o
fo_gr3=obj/gremlin3c.o obj/map_align.o
fo_map=obj/gremlin3.o obj/gremlin3c.o

all: $(OBJDIR) gremlin3 map_align

$(OBJDIR):
	mkdir -p $(OBJDIR)

gremlin3: $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(INCDIRS) $(LIBDIRS) -o gremlin3 $(filter-out $(fo_gr3), $(OBJS)) $(CXXLIBS)

map_align: $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(INCDIRS) $(LIBDIRS) -o map_align $(filter-out $(fo_map), $(OBJS)) $(CXXLIBS)

#score_gramm: $(OBJS) $(DEPS)
#	$(CXX) $(CXXFLAGS) $(INCDIRS) $(LIBDIRS) -o score_gramm $(filter-out $(objgr), $(OBJS)) $(CXXLIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCDIRS) -c -o $@ $<

clean:
	rm $(OBJS) gremlin3 map_align
	rmdir $(OBJDIR)
