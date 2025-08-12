#  Update 2025-08-12:
#  o Add flags for codecov
#                                                            Moon yying7@jh.edu
#
#  Update 2020-01-29: 
#  o Uses required argument of serial, omp, or mpi. 
#  o Use VPATH for finding cpp file in different directories -- this simplifies rules
#  o Abort if gsl-config isn't available
#  o Fixed (INTEL) compiler search 0=found | 1=notfound ; make conditional simple (ifeq 0|1)
#  o Also use conditional for GCC
#  o For objects, use basename to get base file name
#  o Clean up directory prefix (shorten variable names and group)
#  o Simplified obj and bin rule logic and readability.
#  o put rules in canonical order
#  o Now has PROF for profiling. (This is by default overrided with empty PROF.)
#  o Now uses INCS. CXXFLAGS is used for C++ specific options.
#  o Make executables with suffixes ( continuum_membrane_serial | continuum_membrane_mpi | continuum_membrane_omp).
#  --  a bit cleaner                                            Kent milfeld@tacc.utexas.edu
#
# TODO: use function to create VPATH
# TODO: Make rules for *.hpp's
#
# Set terminal width to 220 to avoid viewing wrapped lines in output. A width of 200 avoids most wrapping.


# LDFLAGS = -larmadillo
#  -- armadillo is not a dependency any more. (changed 2022)

VPATH = \
	src/ \
	src/io \
	src/math \
	src/parser \
	src/energy_force \
	src/mesh \
	src/parameters \
	src/diffusion \
	src/vector_math \
	src/dynamics \
	src/gsl_matrix_methods \
	src/linear_algebra \
	src/model

BDIR   = bin
EDIR   = EXEs
ODIR   = obj
SDIR   = src
EPDIR  = EXE_PAR
ECDIR = EXE_CLUSTER

# Add test to path if enabled
ifeq (test,$(MAKECMDGOALS))
	VPATH += tests
	SDIR += tests
endif

# Set version to C++14
CXXFLAGS = -std=c++14

# Enable coverage for "coverage" and "test" targets
# --- Coverage toggle (OFF by default)
COVERAGE ?= 0

# GCC/g++ coverage (works on Linux easily)
ifeq ($(COVERAGE),1)
  CFLAGS   += -O0 -g --coverage
  CXXFLAGS += -O0 -g --coverage
  LDFLAGS  += --coverage
  # Optional: ensure gcov is linked explicitly
  # LIBS += -lgcov
endif

# (Optional) If you use Clang/LLVM for coverage instead of GCC,
# use COVERAGE=llvm and uncomment this block:
ifeq ($(COVERAGE),llvm)
  CFLAGS   += -O0 -g -fprofile-instr-generate -fcoverage-mapping
  CXXFLAGS += -O0 -g -fprofile-instr-generate -fcoverage-mapping
  LDFLAGS  += -fprofile-instr-generate
endif


PROF   = 

.PHONY: any

#               REQUIREMENTS: gls and directories

hasGSL = $(shell type gsl-config >/dev/null 2>&1; echo $$?)
ifeq ($(hasGSL),1)
$(error " GSL must be installed, and gsl-config must be in path.")
else
$(shell mkdir -p bin)
$(shell mkdir -p obj)
endif

#               EXECUTABLE SETUP for serial, MPI, OpenMP (omp), cluster
#
ifeq (serial,$(MAKECMDGOALS))
	_EXEC = continuum_membrane
endif

ifeq (dyna,$(MAKECMDGOALS))
	_EXEC = membrane_dynamics
endif

ifeq (mpi,$(MAKECMDGOALS))
	_EXEC = continuum_membrane
         DEFS = -DMPI
endif

ifeq (omp,$(MAKECMDGOALS))
	_EXEC  = continuum_membrane
         DEFS  = -DOMP
         PLANG = -fopenmp
endif

ifeq (multi,$(MAKECMDGOALS))
	_EXEC = continuum_membrane_multithreading
endif

ifeq (dyna_omp,$(MAKECMDGOALS))
	_EXEC = membrane_dynamics
		 DEFS  = -DOMP
		 PLANG = -fopenmp
endif

ifeq (dyna_multi,$(MAKECMDGOALS))
	_EXEC = membrane_dynamics_multithreading
endif

ifeq (test,$(MAKECMDGOALS))
	_EXEC = test_main
endif

ifeq (clean,$(MAKECMDGOALS))
	MAKECMDGOALS = dummy
endif

EXEC  = $(patsubst %,$(BDIR)/%,$(_EXEC))


OS    := $(shell uname)
INTEL  = $(shell type icpc  >/dev/null 2>&1; echo $$?)
GCC    = $(shell type g++   >/dev/null 2>&1; echo $$?)

INCS    = $(shell gsl-config --cflags) -Iinclude -Iinclude/*
LIBS     = $(shell gsl-config --libs)

# Add pthread to library if multithreading (embarrassingly parallel) is needed
ifeq (multi,$(MAKECMDGOALS))
	LIBS += -pthread
endif

ifeq (dyna_multi,$(MAKECMDGOALS))
	LIBS += -pthread
endif

# Add gtest to library if running unittest
ifeq (test,$(MAKECMDGOALS))
	LIBS += -lgtest -lgtest_main -pthread
	CXXFLAGS += -I./tests
	INCS += -Itests
endif

#---------------COMPILER SETUP

#               comment out next statement to allow profiling with gprof
override PROF   = 

ifeq ($(GCC),0)          # Will use GCC. (See Intel below.)
	CC      = g++
	MPCC    = mpicxx
	CFLAGS  = -O3
	PROF    = -pg -g
#	MPCFLAG = -I /cm/shared/mpi/openmpi/2.1/intel/17.0/include
endif


ifeq ($(INTEL),0)        # Will use Intel, overrides GCC if both present.
	CC      = icpc
	MPCC    = mpicxx
	CFLAGS  = -O3
	PROF    = -pg -g
endif


#---------------OBJECT FILES

ifeq ($(OS),Linux)
       _OBJS = $(shell find $(SDIR) -name "*.cpp" | xargs -n 1 basename | sed -r 's/(\.cc|.cpp)/.o/')
else
       _OBJS = $(shell find $(SDIR) -name "*.cpp" | xargs -n 1 basename | sed -E 's/(\.cc|.cpp)/.o/')
endif

        OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))


#---------------RULES

syntax:
	@echo "------------------------------------"
	@printf '\033[31m%s\033[0m\n' "   USAGE: make serial|mpi|omp"
	@echo "------------------------------------"
	exit 0

#             Rules: for $(MAKECMDGOALS)  serial,     mpi, or            omp            build 
#                        $(EXEC)          bin/continuum_membrane, bin/continuum_membrane_mpi or /bincontinuum_membrane_omp
# Build the executable with or without tests

$(MAKECMDGOALS):$(EXEC)
	@echo "Finished making (re-)building $(MAKECMDGOALS) version, $(EXEC)."

$(EXEC): $(OBJS) $(TEST_OBJ)
	@echo "Compiling $(EDIR)/$(@F).cpp"
	$(CXX) $(CFLAGS) $(CXXFLAGS) $(INCS) $(PROF) $(LDFLAGS) -o $@ $(EDIR)/$(@F).cpp $(OBJS) $(LIBS) $(PLANG)
	@echo "------------"

obj/%.o: %.cpp
	@echo "Compiling $< at $(<F) $(<D)"
	$(CXX) $(CFLAGS) $(CXXFLAGS) $(INCS) $(PROF) $(LDFLAGS) -c $< -o $@ $(PLANG) $(DEFS)
	@echo "------------"

clean:
	rm -rf $(ODIR) bin
	rm -rf *.gcno *.gcda *.gcov coverage/ coverage.info

# Reference: https://www.gnu.org/software/make/manual/html_node/Quick-Reference.html
#            https://www.gnu.org/software/make/
#            https://www.cmcrossroads.com/article/basics-vpath-and-vpath
#            https://www.gnu.org/software/make/manual/html_node/Implicit-Variables.html
