###########################################################
#
# Example GNUmakefile to build the LCFI Vertex Package
#
# @author Frank Gaede, Ben Jeffery
#
###########################################################

# requires path to MARLIN and Marlin v00-09 or higher
# can be called from top level makefile


#----------------------------------------------------------------------------
#                user variables 
#----------------------------------------------------------------------------


# name of the program/library
PROGNAME = LCFIVertex

#fg: changed this to make it an independent marlin package
VERTEX_LCFI = $(MARLINWORKDIR)/packages/LCFIVertex/vertex_lcfi


#----  additional user include/libraries
#-- optionally can be specified in the to level GNUmakefile of $MARLIN

#ifndef USERINCLUDES
# USERINCLUDES = -I $(VERTEX_LCFI) -I $(BOOSTINC)
#endif
VPINCLUDES = -I $(VERTEX_LCFI) -I $(MARLINWORKDIR)/packages/LCFIVertex/boost

ifndef USERLIBS
 USERLIBS =
endif

%.o: %.cc
	$(CXX) $(CPPFLAGS) $(USERINCLUDES) $(VPINCLUDES) -g -ansi -o $@ -c $?

$(VERTEX_LCFI)/src/%.o: $(VERTEX_LCFI)/src/%.cpp
	$(CXX) $(CPPFLAGS) $(USERINCLUDES) $(VPINCLUDES) -g -ansi -o $@ -c $?
COBJS = $(patsubst %.cpp,%.o, $(wildcard $(VERTEX_LCFI)/src/*.cpp) )

$(VERTEX_LCFI)/algo/src/%.o: $(VERTEX_LCFI)/algo/src/%.cpp
	$(CXX) $(CPPFLAGS) $(USERINCLUDES) $(VPINCLUDES) -g -ansi -o $@ -c $?
COBJS += $(patsubst %.cpp,%.o, $(wildcard $(VERTEX_LCFI)/algo/src/*.cpp) )

$(VERTEX_LCFI)/util/src/%.o: $(VERTEX_LCFI)/util/src/%.cpp
	$(CXX) $(CPPFLAGS) $(USERINCLUDES) $(VPINCLUDES) -g -ansi -o $@ -c $?
COBJS += $(patsubst %.cpp,%.o, $(wildcard $(VERTEX_LCFI)/util/src/*.cpp) )

$(VERTEX_LCFI)/zvtop/src/%.o: $(VERTEX_LCFI)/zvtop/src/%.cpp
	$(CXX) $(CPPFLAGS) $(USERINCLUDES) $(VPINCLUDES) -g -ansi -o $@ -c $?
COBJS += $(patsubst %.cpp,%.o, $(wildcard $(VERTEX_LCFI)/zvtop/src/*.cpp) )


#This next line has marlin includes for the tinyxml header. This directory should
#already be included but do it again because it's not always.
$(VERTEX_LCFI)/nnet/src/%.o: $(VERTEX_LCFI)/nnet/src/%.cpp
	$(CXX) $(USERINCLUDES) $(VPINCLUDES) -I $(VERTEX_LCFI)/nnet/inc -I $(MARLIN)/include -g -ansi -pedantic -o $@ -c $?
COBJS += $(patsubst %.cpp,%.o, $(wildcard $(VERTEX_LCFI)/nnet/src/*.cpp) )


#----------------------------------------------------------------------------
#        do not modify anything below in order to use top level makefile
#----------------------------------------------------------------------------

# the base directory of your application
# this is where you  would normaly have the lib, bin and source directories
BASEDIR = ..
INSTALL = $(BASEDIR)/lib
BINDIR = $(BASEDIR)/bin

ifdef MARLINDEBUG
 CPPFLAGS = -g -Wall -ansi
else
 CPPFLAGS = -O3 -D NDEBUG -Wall -ansi
endif


# check for private compiler definition
ifdef MY_CXX
  CXX = $(MY_CXX)
else
  CXX = g++
endif

AR = ar

#--------------------------------------------------------------------------------------
LIBS +=   -Wl,-whole-archive  -L $(INSTALL)  -l$(PROGNAME) -Wl,-no-whole-archive
#--------------------------------------------------------------------------------------

INCLUDE +=  -I $(BASEDIR)/include
INCLUDE += $(USERINCLUDES)
CPPFLAGS += $(INCLUDE)
CPPFLAGS += `. $(MARLIN)/bin/marlin_includes.sh`

LIBS += $(USERLIBS)
LIBS += `. $(MARLIN)/bin/marlin_libs.sh`

bin:  lib $(BINDIR)/$(PROGNAME)

$(BINDIR)/$(PROGNAME): $(INSTALL)/lib$(PROGNAME).a
	mkdir -p $(BINDIR)
	$(CXX) -o $(BINDIR)/$(PROGNAME) $(CPPFLAGS) $(MARLIN)/src/Marlin.cc $(LIBS)
	
objects = $(patsubst %.cc,%.o, $(wildcard *.cc)  )
objects += $(COBJS)

$(INSTALL)/lib$(PROGNAME).a: $(objects)
	mkdir -p $(INSTALL)
	$(AR) cr $(INSTALL)/lib$(PROGNAME).a  $?

lib: $(INSTALL)/lib$(PROGNAME).a

clean:
	rm -f *.o $(INSTALL)/lib$(PROGNAME).a $(BINDIR)/$(PROGNAME) $(VERTEX_LCFI)/src/*.o $(VERTEX_LCFI)/algo/src/*.o $(VERTEX_LCFI)/util/src/*.o $(VERTEX_LCFI)/nnet/src/*.o $(VERTEX_LCFI)/zvtop/src/*.o

doc:
	cd $(BASEDIR)/doc ; doxygen ;


