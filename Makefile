
CXX=g++
CXXFLAGS += -Wall -Wno-unused-result -pedantic -Ofast -march=native -std=c++0x -g -ggdb
CXXLIBS = -lgsl -lgslcblas
INCDIRS = -I./src
LIBDIRS =

OBJDIR = obj
SRCDIR = src

DEPS = $(shell find $(SRCDIR) -name '*.h')
SRCS = $(shell find $(SRCDIR) -name '*.cpp')
OBJS = $(patsubst $(SRCDIR)%.cpp, $(OBJDIR)%.o, $(SRCS))

objgr3=obj/gremlin3.o
objgr3c=obj/gremlin3c.o

all: $(OBJDIR) gremlin3

$(OBJDIR):
	mkdir -p $(OBJDIR)

gremlin3: $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(INCDIRS) $(LIBDIRS) -o gremlin3 $(filter-out $(objgr3c), $(OBJS)) $(CXXLIBS)

#score_gramm: $(OBJS) $(DEPS)
#	$(CXX) $(CXXFLAGS) $(INCDIRS) $(LIBDIRS) -o score_gramm $(filter-out $(objgr), $(OBJS)) $(CXXLIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCDIRS) -c -o $@ $<

clean:
	rm $(OBJS) gremlin3
	rmdir $(OBJDIR)
