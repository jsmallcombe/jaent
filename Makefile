
LIB = $(shell pwd)
INCLUDE = $(shell pwd)/include

ROOT_LIBS = `root-config --glibs`
ROOT_GCC_FLAGS = `root-config --cflags`

EXTERNAL_LIBS = -ljames_phys -lJanalysistools
CC = g++
CFLAGS = -std=c++11 -g -fPIC -Wall $(ROOT_GCC_FLAGS) -I$(INCLUDE) 
LIBRS = -L$(INCLUDE) $(EXTERNAL_LIBS) $(ROOT_LIBS)

INCLUDE = include
SOURCEDIR = src
BINDIR = bin

TARG = bin/jaent_simm_one
MAIN = src/wrapper.cpp

EXPCORE = $(wildcard $(SOURCEDIR)/exp_core_*)
EXPCOROBJ = $(patsubst $(SOURCEDIR)/%.cxx,$(BINDIR)/%.o,$(EXPCORE))

$(TARG):$(MAIN) $(BINDIR)/exp_core_total.o $(BINDIR)/auto_setup.o  $(BINDIR)/detector_class.o
	$(CC) -o $@ $(CFLAGS) $(MAIN) $(LIBRS) $(BINDIR)/detector_class.o $(BINDIR)/exp_core_total.o $(BINDIR)/auto_setup.o
# 	-lstdc++
# 	$(CC) $(CFLAGS) -o $@ $(LIBRS) bin/DictOutput.cxx -I. $(OBJECTS)

$(BINDIR)/exp_core_total.o: $(EXPCOROBJ)
	ld -r -o  $@ $(EXPCOROBJ)

$(BINDIR)/auto_setup.o: $(SOURCEDIR)/auto_setup.cxx $(INCLUDE)/auto_setup.h $(BINDIR)/exp_core_total.o
	$(CC) $(CFLAGS) -o $@ -c $< $(LIBRS)
	
$(BINDIR)/detector_class.o: $(SOURCEDIR)/detector_class.cxx $(INCLUDE)/detector_class.h
	$(CC) $(CFLAGS) -o $@ -c $< $(LIBRS)	

$(BINDIR)/%.o: $(SOURCEDIR)/%.cxx $(INCLUDE)/exp_core.h $(BINDIR)/detector_class.o
	$(CC) $(CFLAGS) -o $@ -c $< $(LIBRS)
	
clean: 
	rm $(BINDIR)/*.o
	rm $(TARG)

examples/%: $(TARG) FORCE
	$(CC) -o bin/$(subst examples/,,$(subst .cpp,,$@)) $(CFLAGS) $@ $(LIBRS) $(BINDIR)/detector_class.o $(BINDIR)/exp_core_total.o $(BINDIR)/auto_setup.o

FORCE:
