
CXX=g++
CXXFLAGS += -Wall -pedantic -O3 -mtune=native -std=c++0x -g -ggdb
CXXLIBS = -lgsl -lgslcblas
INCDIRS =
LIBDIRS =

OBJDIR = obj
SRCDIR = src

DEPS = $(shell find $(SRCDIR) -name '*.h')
SRCS = $(shell find $(SRCDIR) -name '*.cpp')
OBJS = $(patsubst $(SRCDIR)%.cpp, $(OBJDIR)%.o, $(SRCS))

objgr=obj/gremlin3.o
objsg=obj/score_gramm.o

all: $(OBJDIR) gremlin3

$(OBJDIR):
	mkdir -p $(OBJDIR)

gremlin3: $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(INCDIRS) $(LIBDIRS) -o gremlin3 $(filter-out $(objsg), $(OBJS)) $(CXXLIBS)

#score_gramm: $(OBJS) $(DEPS)
#	$(CXX) $(CXXFLAGS) $(INCDIRS) $(LIBDIRS) -o score_gramm $(filter-out $(objgr), $(OBJS)) $(CXXLIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCDIRS) -c -o $@ $<

clean:
	rm $(OBJS) gremlin3 score_gramm
	rmdir $(OBJDIR)
