# Custom makefile for LimboInterface
# Directories and Names
LIB_PREFIX = gds
LIMBO_ROOT_DIR = $(realpath ../Limbo/)
PARSER_SPEF_ROOT_DIR = $(realpath ../Parser-SPEF/)
OBJDIR = $(realpath ./)/obj
LIBDIR = $(LIMBO_ROOT_DIR)/lib
MKDIR = if [ ! -d $(@D) ]; then mkdir -p $(@D); fi

# Compilation Flags
DBG = 0 # Off by default
#include $(LIMBO_ROOT_DIR)/Include.mk # Include environ config

ifeq ($(DBG), 1)
	CXXFLAGS = $(CXXFLAGS_DEBUG) -DDEBUG_GDSREADER -DDEBUG_GDSWRITER
else
	CXXFLAGS = $(CXXFLAGS_RELEASE) -std=C++17 -g -lstdc++fs
endif

ifdef ZLIB_DIR # Compression support
ifdef BOOST_DIR # Boost library support
	CXXFLAGS += -DZLIB=1
endif
endif

# Special Libraries to Include
INCLUDE = -I $(LIMBO_ROOT_DIR) -I $(PARSER_SPEF_ROOT_DIR)

ifdef ZLIB_DIR
ifdef BOOST_DIR
	INCLUDE += -I $(BOOST_DIR)/include $(ZLIB_INCLUDE_FLAG)
	LIB += -L $(BOOST_DIR)/lib -lboost_iostreams \
		   $(ZLIB_LINK_FLAG) -lz
endif
endif
LIB += $(LIBDIR)

# Standard make Settings and Recipes
SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:%.cpp=$(OBJDIR)/%.o)
DEPS = $(OBJS:%.o=%.d) # Dependency file for each source

all: $(SRCS) LimboInterface

$(OBJDIR)/%.d: %.cpp
	@$(MKDIR)
	$(CXX) $(CXXFLAGS) $< -MM -MT $(@:%.d=%.o) >$@ $(INCLUDE)
-include %(DEPS)

%(OBJDIR)/%.o: %.cpp
	@$(MKDIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDE)

LimboInterface: $(OBJS) $(LIBDIR)/lib$(LIB_PREFIX)parser.a
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LIB) -l$(LIB_PREFIX)parser $(INCLUDE)

explicit: TestLimboInterface.cpp
	g++ -std=c++17 -g -lstdc++fs -o Test_$@ TestLimboInterface.cpp -L $(LIBDIR) -l$(LIB_PREFIX)parser -I $(LIMBO_ROOT_DIR) -I $(PARSER_SPEF_ROOT_DIR)
			
.PHONY: clean
clean: cleandep
	rm -f LimboInterface core

.PHONY: cleandep
cleandep:
	rm -f $(DEPS)
