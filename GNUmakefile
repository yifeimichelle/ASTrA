# Place commonly changed variables
# at the top of the makefile

DIM = 3
TRAJ_SRC = $(HOME)/src/trajectory
ANALYSIS_SRC = $(HOME)/src/analysis
WRITER_SRC = $(HOME)/src/writer
CFLAGS = -Wall -I$(TRAJ_SRC) -I$(ANALYSIS_SRC) -I$(WRITER_SRC) -std=c++11

#CXX = g++
CXX = clang++
CPPFLAGS = -D DIM=$(DIM)

TRAJ_SRCFILES:=$(wildcard $(TRAJ_SRC)/*.cpp)
TRAJ_OBJS:=$(patsubst %.cpp, %.o, $(TRAJ_SRCFILES))
ANALYSIS_SRCFILES:=$(wildcard $(ANALYSIS_SRC)/*.cpp)
ANALYSIS_OBJS:=$(patsubst %.cpp, %.o, $(ANALYSIS_SRCFILES))
WRITER_SRCFILES:=$(wildcard $(WRITER_SRC)/*.cpp)
WRITER_OBJS:=$(patsubst %.cpp, %.o, $(WRITER_SRCFILES))
OBJS = $(TRAJ_OBJS) $(ANALYSIS_OBJS) $(WRITER_OBJS)

%.o:%.cpp GNUmakefile
	$(CXX) -c $(DEBUGFLAGS) $(CPPFLAGS) $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*.d

clean:
	rm -r $(OBJS) $(OBJS:.o=.d) *.exe *.dSYM

#-include $(OBJS:-o=.d)
