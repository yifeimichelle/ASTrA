# Place commonly changed variables
# at the top of the makefile

DIM = 3
TRAJ_SRC = $(HOME)/src/trajectory
CFLAGS = -Wall -I$(TRAJ_SRC) -std=c++11

#CXX = g++
CXX = clang++
CPPFLAGS = -D DIM=$(DIM)

TRAJ_SRCFILES:=$(wildcard $(TRAJ_SRC)/*.cpp)
TRAJ_OBJS:=$(patsubst %.cpp, %.o, $(TRAJ_SRCFILES))
OBJS = $(TRAJ_OBJS)

%.o:%.cpp GNUmakefile
	$(CXX) -c $(CPPFLAGS) $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*.d

clean:
	rm -r $(OBJS) $(OBJS:.o=.d) *.exe *.dSYM

#-include $(OBJS:-o=.d)
