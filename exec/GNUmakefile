#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.

HOME = ./..

# We've made a common set of makefile definitions and rules in the home directory.
# By including that file, you have all the definitions you need.

CHARGE=1 # 1 if with charges, 0 if without
#DEBUGFLAGS = -DDEBUG
include $(HOME)/GNUmakefile
#CFLAGS += -O3
CFLAGS += -g

# ON NERSC:
# change name of executable to reflect machine:
# .cori for cori, .exe for edison
trajAnalyzer: GNUmakefile $(OBJS) trajMain.cpp
	  $(CXX) $(CPPFLAGS) $(CFLAGS) -o trajAnalyzerChg.exe trajMain.cpp  $(OBJS)

