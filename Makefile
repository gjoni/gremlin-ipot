
CXX=g++
CXXFLAGS += -Wall -Wno-unused-result -pedantic -Ofast -march=native -std=c++0x -g -ggdb -fopenmp
CXXLIBS =
INCDIRS = -I./src
LIBDIRS =

OBJDIR = obj
SRCDIR = src
BINDIR = bin

DEPS = $(shell find $(SRCDIR) -name '*.h')
SRCS = $(shell find $(SRCDIR) -name '*.cpp')
OBJS = $(patsubst $(SRCDIR)%.cpp, $(OBJDIR)%.o, $(SRCS))

fo_gr3=obj/gremlin3c.o obj/neff.o obj/rstgen.o
fo_map=obj/gremlin3.o obj/gremlin3c.o

all: $(OBJDIR) $(BINDIR) gremlin3 neff rstgen

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

gremlin3: $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(INCDIRS) $(LIBDIRS) -o $(BINDIR)/gremlin3 $(filter-out $(fo_gr3), $(OBJS)) $(CXXLIBS)

neff: $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(INCDIRS) $(LIBDIRS) -o $(BINDIR)/neff obj/neff.o obj/MSAclass.o $(CXXLIBS)

rstgen: $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(INCDIRS) $(LIBDIRS) -o $(BINDIR)/rstgen obj/rstgen.o obj/MSAclass.o $(CXXLIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCDIRS) -c -o $@ $<

clean:
	rm $(OBJS)
	rm $(BINDIR)/gremlin3
	rm $(BINDIR)/neff
	rm $(BINDIR)/rstgen
	rmdir $(OBJDIR)
	rmdir $(BINDIR)
