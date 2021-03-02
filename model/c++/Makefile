#configuration file
include config.mk

#-------------------------------------------------------------------------------
#local directories

#source code directory
dirs=src
#built object directory
diro=obj
#built executable directory
dirb=bin
#flags to include libode
odesrc=-I$(odepath)/src
odelib=-L$(odepath)/bin -lode

#-------------------------------------------------------------------------------
#stuff to compile

obj=$(diro)/util.o \
	$(diro)/linalg.o \
	$(diro)/newton.o \
	$(diro)/msc_gcv.o

all: libodemake $(obj) $(dirb)/single.exe $(dirb)/sweep.exe

#-------------------------------------------------------------------------------
#compilation rules

libodemake:
	$(MAKE) -C $(odepath)

$(diro)/util.o: $(dirs)/util.cc $(dirs)/util.h
	$(cxx) $(flags) -o $@ -c $< -I$(dirs)

$(diro)/linalg.o: $(dirs)/linalg.cc $(dirs)/linalg.h
	$(cxx) $(flags) -o $@ -c $< -I$(dirs)

$(diro)/newton.o: $(dirs)/newton.cc $(dirs)/newton.h $(diro)/linalg.o
	$(cxx) $(flags) -o $@ -c $< -I$(dirs)

$(diro)/msc_gcv.o: $(dirs)/msc_gcv.cc $(dirs)/msc_gcv.h $(diro)/util.o $(diro)/newton.o
	$(cxx) $(flags) -o $@ -c $< -I$(dirs) $(odesrc) $(odelib)

$(dirb)/single.exe: $(dirs)/main_single.cc $(diro)/msc_gcv.o
	$(cxx) $(flags) $(omp) -o $@ $< $(obj) $(odesrc) $(odelib)

$(dirb)/sweep.exe: $(dirs)/main_sweep.cc $(diro)/msc_gcv.o
	$(cxx) $(flags) $(omp) -o $@ $< $(obj) $(odesrc) $(odelib)
