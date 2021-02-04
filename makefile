# Custom makefile for gds2Para LayoutAnalyzer
# Directories and Names
LIB_PREFIX = gds
BOOST_ROOT_DIR = ${BOOST_DIR}
LIMBO_ROOT_DIR = ${LIMBO_DIR}
LIMBO_LIB_DIR = ${LIMBO_DIR}/lib
PARSER_SPEF_ROOT_DIR = ${PARSER_SPEF_DIR}
EIGEN_ROOT_DIR = ${EIGEN_DIR}
MKL_ROOT_DIR = ${MKL_DIR}
INTEL_LIB_DIR = ${INTEL_MATH_LIB}
HYPRE_LIB_DIR = ${HYPRE_DIR}/src/hypre/lib
HYPRE_HEAD_DIR = ${HYPRE_DIR}/src/hypre/include
SRCDIR = $(realpath ./)/src
OBJDIR = $(realpath ./)/obj
MKDIR = if [ ! -d $(@D) ]; then mkdir -p $(@D); fi

# Compilation Flags (MKL is serial xor threaded with Intel BLACS and debug options)
MKL_LINK_FLAGS =-Wl,--start-group $(MKL_ROOT_DIR)/lib/intel64/libmkl_intel_ilp64.a $(MKL_ROOT_DIR)/lib/intel64/libmkl_core.a $(MKL_ROOT_DIR)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lm -pthread -ldl
#MKL_LINK_FLAGS =-Wl,--start-group $(MKL_ROOT_DIR)/lib/intel64/libmkl_intel_ilp64.a $(MKL_ROOT_DIR)/lib/intel64/libmkl_intel_thread.a $(MKL_ROOT_DIR)/lib/intel64/libmkl_core.a $(MKL_ROOT_DIR)/lib/intel64/libmkl_blacs_openmpi_ilp64.a -Wl,--end-group -L $(INTEL_LIB_DIR) -liomp5 -lpthread -lm -ldl
MKL_COMP_FLAGS =-m64 -I $(MKL_ROOT_DIR)/include
DBG =0 # Off by default
#include $(LIMBO_LIB_DIR)/../Include.mk # Include environ config

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

# Flags -I -L for Limbo DEF and LEF parsers
LIMBO_DEF_LEF_FLAGS = -L $(LIMBO_LIB_DIR) -ldefparseradapt -llefparseradapt -I $(LIMBO_ROOT_DIR)

# Special Libraries to Include
INCLUDE =-I $(LIMBO_ROOT_DIR) -I $(PARSER_SPEF_ROOT_DIR) -I $(EIGEN_ROOT_DIR) -I $(HYPRE_HEAD_DIR) $(MKL_COMP_FLAGS)
LIB =-L $(LIMBO_LIB_DIR) -l$(LIB_PREFIX)parser -ldefparseradapt -llefparseradapt -L $(HYPRE_LIB_DIR) -lHYPRE -lm

ifdef false #ZLIB_DIR
ifdef BOOST_DIR
	INCLUDE += -I $(BOOST_DIR) $(ZLIB_INCLUDE_FLAG)
	LIB += -L $(BOOST_DIR)/../libs -lboost_iostreams \
		   $(ZLIB_LINK_FLAG) -lz
endif
endif

# Standard make Settings and Recipes
SRCS = $(wildcard $(SRCDIR)/*.cpp)
OBJS = $(SRCS:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
DEPS = $(OBJS:%.o=%.d) # Dependency file for each source

all: LayoutAnalyzer

LayoutAnalyzer: $(OBJS)
	mpicxx -w $(CXXFLAGS) -o $@ $(OBJS) $(LIB) $(INCLUDE) $(MKL_LINK_FLAGS) $(LFLAGS)

$(OBJDIR)/TestMain.o: $(SRCDIR)/TestMain.cpp $(SRCDIR)/fdtd.hpp $(SRCDIR)/limboint.hpp $(SRCDIR)/solnoutclass.hpp
	@$(MKDIR)
	mpicxx -w -std=c++17 -g -lstdc++fs -O0 -c $(SRCDIR)/TestMain.cpp -o $(OBJDIR)/TestMain.o -L $(LIMBO_LIB_DIR) -l$(LIB_PREFIX)parser -ldefparseradapt -llefparseradapt -I $(LIMBO_ROOT_DIR) -I $(PARSER_SPEF_ROOT_DIR) -I $(EIGEN_ROOT_DIR) $(MKL_COMP_FLAGS)

$(OBJDIR)/mesh.o: $(SRCDIR)/mesh.cpp $(SRCDIR)/fdtd.hpp
	@$(MKDIR)
	mpicxx -w -g -O1 -c $(SRCDIR)/mesh.cpp -o $(OBJDIR)/mesh.o $(MKL_COMP_FLAGS)

$(OBJDIR)/matrixCon.o: $(SRCDIR)/matrixCon.cpp $(SRCDIR)/fdtd.hpp $(SRCDIR)/hypreSolver.h
	@$(MKDIR)
	mpicxx -w -g -O1 -c $(SRCDIR)/matrixCon.cpp -o $(OBJDIR)/matrixCon.o $(MKL_COMP_FLAGS) -I $(HYPRE_HEAD_DIR) -L $(HYPRE_LIB_DIR) -lHYPRE -lm $(LFLAGS)

$(OBJDIR)/hypreSolve.o: $(SRCDIR)/hypreSolve.cpp $(SRCDIR)/fdtd.hpp $(SRCDIR)/hypreSolver.h
	@$(MKDIR)
	mpicxx -w -g -O1 -c $(SRCDIR)/hypreSolve.cpp -o $(OBJDIR)/hypreSolve.o -I $(MKL_COMP_FLAGS) -I $(HYPRE_HEAD_DIR) -L $(HYPRE_LIB_DIR) -lHYPRE -lm $(LFLAGS)

$(OBJDIR)/generateStiff.o: $(SRCDIR)/generateStiff.cpp $(SRCDIR)/fdtd.hpp
	@$(MKDIR)
	mpicxx -w -g -O1 -c $(SRCDIR)/generateStiff.cpp -o $(OBJDIR)/generateStiff.o $(MKL_COMP_FLAGS)

$(OBJDIR)/findVh.o: $(SRCDIR)/findVh.cpp $(SRCDIR)/fdtd.hpp
	@$(MKDIR)
	mpicxx -w -g -O1 -c $(SRCDIR)/findVh.cpp -o $(OBJDIR)/findVh.o $(MKL_COMP_FLAGS)

$(OBJDIR)/autoPortFromDefLef.o: $(SRCDIR)/autoPortFromDefLef.cpp $(SRCDIR)/autoPortFromDefLef.hpp
	@$(MKDIR)
	mpicxx -w -g -O1 -c $(SRCDIR)/autoPortFromDefLef.cpp -o $(OBJDIR)/autoPortFromDefLef.o $(LIMBO_DEF_LEF_FLAGS)


.PHONY: clean
clean: cleandep
	rm -f LayoutAnalyzer

.PHONY: cleandep
cleandep:
	rm -f $(OBJS)
