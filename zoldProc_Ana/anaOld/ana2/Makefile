# Makefile for phaar analysis code
# Heavily based off Thomas Lindner's makefile for analyzer library
# Edited by Austin dsc, #ac means commented out by austin

### executable files	.cpp
### source files			src/.cxx
### headers						include/.h
### objecsts					objs/.o

CXX = g++ -std=c++11			#what compiler and c++ version we want
CXXFLAGS = -O2 -Wall -lm -Wno-unused-variable	#-g is debugger, -02 is optimizer, -Wall gives compiler warnings -lm is basic math library

# ROOT libraries
ifdef ROOTSYS
#if ROOTSYS is defined in the bash rc profile
# the equal sign := means set ROOTCFLAGS right away (expanded variable)
ROOTCFLAGS := $(shell root-config --cflags) 	#this gives root/include
ROOTFLAGS := $(shell root-config --libs)		#root library stuff?
ROOTGLIBS := $(shell root-config --glibs) 		#this gives lthread library

#ac decided to put rootcflags in actual lines below instead of in CXXFLAGS
CXXFLAGS  += -DHAVE_ROOT $(ROOTCFLAGS) $(ROOTFLAGS)
endif


# CXXFLAGS+=-I /usr/local/include/

# ### our headers ###
# CXXFLAGS+= -I./include

LIBS+= -lz -lpthread -lutil -lgsl


# this is a print statement
#$(info Printthis [$(LIBS)])
# $(info xxROOTCFLAGS [$(ROOTGLIBS)])
#$(info )
$(info )
# $(info xxROOTGLibs [$(ROOTGLIBS)])

### pass objects to executables (just pass all to be 'safe') ###
all: test.exe # gslTest.exe

test.exe: compGas2_betaFit.cpp
		$(CXX) -o $@ $(ROOTCFLAGS) $(CXXFLAGS) $^ $(ROOTGLIBS) $(LIBS)

# gslTest.exe: testGSL.cpp
# 		$(CXX) -o $@ $(CXXFLAGS) $^ $(LIBS)

# 'make clean' deletes all these objects
clean::
	-rm -f *.o objs/*.o *.a
	-rm -f *.exe
	-rm -rf *.exe.dSYM

# end
