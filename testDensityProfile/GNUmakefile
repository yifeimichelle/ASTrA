#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.

HOME = ./..

# We've made a common set of makefile definitions and rules in the home directory.
# By including that file, you have all the definitions you need.

DEBUGFLAGS = -DDEBUG
include $(HOME)/GNUmakefile
#CFLAGS += -O3
CFLAGS += -g

testDensityProfile: GNUmakefile $(OBJS) testDensityProfile.cpp
	  $(CXX) $(CPPFLAGS) $(CFLAGS) -o testDensityProfile.exe testDensityProfile.cpp  $(OBJS)

