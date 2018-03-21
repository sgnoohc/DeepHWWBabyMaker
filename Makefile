#
# stuff to make
#
ifndef ROOTSYS
all:
	@echo "ROOTSYS is not set. Please set ROOT environment properly"; echo
else

all: 	build
help:
	@echo "Available Targets:";\
	cat Makefile | perl -ne 'printf("\t%-15s %s\n",$$1,$$2) if(/^(\S+):[^#]+(#.*)$$/)'

ifndef VERBOSE
  QUIET := @
endif

PACKAGES = rooutil coreutil

LIBDIR = libs
OBJDIR = objs

CC = g++
CMSROOT = ./
INCLUDE = $(shell root-config --cflags) -I$(CMSROOT) -I$(CMSROOT)/CORE
CFLAGS = -Wall -Wno-unused-function -g -O2 -fPIC $(INCLUDE) $(EXTRACFLAGS) -fno-var-tracking
ROOTLIBS = $(shell root-config --ldflags --cflags --libs) -lTMVA #-lEG -lGenVector

DICTINCLUDE = $(ROOTSYS)/include/Math/QuantFuncMathCore.h $(ROOTSYS)/include/TLorentzVector.h $(ROOTSYS)/include/Math/Vector4D.h

LINKER = g++
LINKERFLAGS = $(shell root-config --ldflags --libs) -lEG -lGenVector -lTMVA -O2

DIR = ./

SOURCES = $(wildcard *.cc)
OBJECTS = $(SOURCES:%.cc=$(OBJDIR)/%.o)
LIB = $(LIBDIR)/libScanChain.so

CORESOURCES = $(wildcard CORE/*.cc)
COREOBJECTS = $(CORESOURCES:.cc=.o) 
CORELIB = $(LIBDIR)/libBabymakerCORE.so

TOOLSSOURCES = $(wildcard CORE/Tools/*.cc) $(wildcard CORE/Tools/MT2/*.cc) $(wildcard CORE/Tools/btagsf/*.cc) $(wildcard CORE/Tools/datasetinfo/*.cc)
TOOLSOBJECTS = $(TOOLSSOURCES:.cc=.o) 
TOOLSLIB = $(LIBDIR)/libBabymakerTools.so

#DICT = LinkDef_out.o
DICT = 

LIBS = $(LIB) $(CORELIB) $(TOOLSLIB) $(FWLIB) 

EXE = processBaby

.PHONY: all help compile clean cms2env $(PACKAGES)

#
# how to make it
#

libs:	$(LIBS)

$(CORELIB): $(DICT) $(COREOBJECTS)
	$(QUIET) echo "Linking $@"; \
	echo "$(LINKER) -shared -o $@ $(COREOBJECTS) $(DICT) $(LINKERFLAGS)"; \
	$(LINKER) -shared -o $@ $(COREOBJECTS) $(DICT) $(LINKERFLAGS)

$(TOOLSLIB): $(DICT) $(TOOLSOBJECTS)
	$(QUIET) echo "Linking $@"; \
	echo "$(LINKER) -shared -o $@ $(TOOLSOBJECTS) $(DICT) $(LINKERFLAGS)"; \
	$(LINKER) -shared -o $@ $(TOOLSOBJECTS) $(DICT) $(LINKERFLAGS)

$(LIB):	$(DICT) $(OBJECTS)
	$(QUIET) echo "Linking $@"; \
	echo "$(LINKER) -shared -o $@ $(OBJECTS) $(DICT) $(LINKERFLAGS)"; \
	$(LINKER) -shared -o $@ $(OBJECTS) $(DICT) $(LINKERFLAGS)

$(PACKAGES):
	make -j 15 -C $@
#	cp $@/lib$@.so $(LIBDIR)/

# the option "-Wl,-rpath,./" adds ./ to the runtime search path for libraries
$(EXE): $(PACKAGES) makedir $(LIBS)
	$(QUIET) echo "Building $@"; \
	echo "$(CC) -o $@ $(LIBS) $(ROOTLIBS) -Wl,-rpath,./ $(addprefix -L,$(PACKAGES)) $(addprefix -l,$(PACKAGES))"; \
	$(CC) -o $@ $(LIBS) $(ROOTLIBS) -Wl,-rpath,./  $(addprefix -L,$(PACKAGES)) $(addprefix -l,$(PACKAGES))

%.o: 	%.cc %.h
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

$(OBJDIR)/%.o: 	%.C %.h
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

$(OBJDIR)/%.o: 	%.cc
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

$(OBJDIR)/%.o: 	%.C
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

$(OBJDIR)/%.o:    %.cxx 
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

libs:	$(LIBS)

build:  $(EXE)

makedir:
	mkdir -p $(LIBDIR) $(OBJDIR)

b: build

clean:
	rm -f \
	$(OBJECTS) \
	$(LIBS) \
	$(EXE) \
	CORE/*.o \
	CORE/*.d \
	CORE/*.so \
	CORE/Tools/*.o \
	CORE/Tools/*.d \
	CORE/Tools/*.so \
	CORE/Tools/MT2/*.o \
	CORE/Tools/MT2/*.d \
	CORE/Tools/MT2/*.so 

localclean:
	rm -f \
	LinkDef_out* \
	*.o \
	*.so 

endif
