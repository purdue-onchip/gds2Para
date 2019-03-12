# Custom makefile for LimboInterface
# Directories and Names
LIB_PREFIX = gds
LIMBO_ROOT_DIR = $(realpath ../Limbo/)
PARSER_SPEF_ROOT_DIR = $(realpath ../Parser-SPEF/)
EIGEN_ROOT_DIR = $(realpath ../eigen-git-mirror/)
MKL_ROOT_DIR = ${MKL_DIR}
OBJDIR = $(realpath ./)/obj
LIBDIR = $(LIMBO_ROOT_DIR)/lib
MKDIR = if [ ! -d $(@D) ]; then mkdir -p $(@D); fi
HYPRE_DIR = $(realpath ../hypre/src/hypre)

# Compilation Flags
DBG =0 # Off by default
#include $(LIMBO_ROOT_DIR)/Include.mk # Include environ config

ifeq ($(DBG), 1)
	CXXFLAGS =$(CXXFLAGS_DEBUG) -DDEBUG_GDSREADER -DDEBUG_GDSWRITER
else
	CXXFLAGS =$(CXXFLAGS_RELEASE) -std=c++17 -O2 -lstdc++fs
endif

ifdef false #ZLIB_DIR # Compression support
ifdef BOOST_DIR # Boost library support
	CXXFLAGS += -DZLIB=1
endif
endif

MKLFLAGS =-Wl,--start-group $(MKL_ROOT_DIR)/lib/intel64/libmkl_intel_lp64.a $(MKL_ROOT_DIR)/lib/intel64/libmkl_core.a $(MKL_ROOT_DIR)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lm -pthread -ldl



# Special Libraries to Include
INCLUDE =-I $(LIMBO_ROOT_DIR) -I $(PARSER_SPEF_ROOT_DIR) -I $(EIGEN_ROOT_DIR) -I $(MKL_ROOT_DIR)/include -I $(HYPRE_DIR)/include
LIB =-L $(LIBDIR) -l$(LIB_PREFIX)parser -L $(HYPRE_DIR)/lib -lHYPRE -lm

ifdef false #ZLIB_DIR
ifdef BOOST_DIR
	INCLUDE += -I $(BOOST_DIR)/include $(ZLIB_INCLUDE_FLAG)
	LIB += -L $(BOOST_DIR)/lib -lboost_iostreams \
		   $(ZLIB_LINK_FLAG) -lz
endif
endif

# Standard make Settings and Recipes
SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:%.cpp=$(OBJDIR)/%.o)
DEPS = $(OBJS:%.o=%.d) # Dependency file for each source

all: LayoutAnalyzer

$(OBJDIR)/%.d: %.cpp
	@$(MKDIR)
	$(CXX) $(CXXFLAGS) $< -MM -MT $(@:%.d=%.o) >$@ $(INCLUDE)
-include %(DEPS)

%(OBJDIR)/%.o: %.cpp
	@$(MKDIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDE)

LimboInterface: $(OBJS) $(LIBDIR)/lib$(LIB_PREFIX)parser.a
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LIB) -l$(LIB_PREFIX)parser $(INCLUDE)

LayoutAnalyzer: TestMain.o mesh.o matrixCon.o hypreSolve.o
	mpicxx $(CXXFLAGS) -o $@ TestMain.o mesh.o matrixCon.o hypreSolve.o $(LIB) $(INCLUDE) $(MKLFLAGS) $(LFLAGS)

explicit: TestMain.o mesh.o matrixCon.o
	g++ -std=c++17 -g -lstdc++fs -o Test_$@ TestMain.o mesh.o matrixCon.o hypreSolve.o -L $(LIBDIR) -l$(LIB_PREFIX)parser -I $(LIMBO_ROOT_DIR)/limbo/parsers/gdsii/stream/ -I $(PARSER_SPEF_ROOT_DIR) -I $(EIGEN_ROOT_DIR) -I$(MKL_ROOT_DIR)/include -Wl,--start-group $(MKL_ROOT_DIR)/lib/intel64/libmkl_intel_lp64.a $(MKL_ROOT_DIR)/lib/intel64/libmkl_core.a $(MKL_ROOT_DIR)/lib/intel64/libmkl_sequential.a -Wl,--end-group  -lm -pthread -ldl -L $(HYPRE_DIR)/lib

TestMain.o: TestMain.cpp fdtd.h limboint.h solnoutclass.h #$(LIMBO_ROOT_DIR)/limbo/parsers/gdsii/stream/GdsReader.h $(LIMBO_ROOT_DIR)/limbo/parsers/gdsii/stream/GdsWriter.h $(PARSER_SPEF_ROOT_DIR)/parser-spef/parser-spef.hpp $(EIGEN_ROOT_DIR)/Eigen/Sparse
	g++ -std=c++17 -g -lstdc++fs -c TestMain.cpp -L $(LIBDIR) -l$(LIB_PREFIX)parser -I $(LIMBO_ROOT_DIR) -I $(LIMBO_ROOT_DIR)/limbo/parsers/gdsii/stream/ -I $(PARSER_SPEF_ROOT_DIR) -I $(EIGEN_ROOT_DIR) -I $(MKL_ROOT_DIR)/include -L $(HYPRE_DIR)/lib -lHYPRE -lm
	
mesh.o: mesh.cpp fdtd.h
	g++ -c -g mesh.cpp -I $(MKL_ROOT_DIR)/include

matrixCon.o: matrixCon.cpp fdtd.h
	mpicxx -c -g matrixCon.cpp -I $(MKL_ROOT_DIR)/include -I $(HYPRE_DIR)/include -L $(HYPRE_DIR)/lib -lHYPRE -lm $(LFLAGS)

hypreSolve.o: hypreSolve.cpp fdtd.h vis.c
	mpicxx -c -g hypreSolve.cpp -I $(MKL_ROOT_DIR)/include -I $(HYPRE_DIR)/include -L $(HYPRE_DIR)/lib -lHYPRE -lm $(LFLAGS)


.PHONY: clean
clean: cleandep
	rm -f LimboInterface core

.PHONY: cleandep
cleandep:
	rm -f $(DEPS)
