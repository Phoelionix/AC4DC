
CPP := g++-10

LIB := -fopenmp

INC := -I$(HOME)/Programming/include -Iinclude

SRCDIR := src
BUILDDIR := build
TARGET := bin/ac4dc

SRCEXT := cpp

SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT) )

OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
CXXFLAGS := -std=c++17 -fopenmp -MD -g -Wall

debug: CXXFLAGS += -DDEBUG -Wpedantic
release: CXXFLAGS += -O3 -DNDEBUG

# For tests:
TESTS := $(shell find tests -type f -name *.$(SRCEXT) )
TINC := $(INC) -I./

$(TARGET): $(OBJECTS)
	@echo " Linking $(TARGET)... ";
	@echo " $(CPP) $^ $(LIB) -o $(TARGET) "; $(CPP) $^ $(LIB) -o $(TARGET)

all: $(TARGET)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR) bin $(BUILDDIR)/Wigner $(BUILDDIR)/$(MAINSUBDIR)
	@echo " $(CPP) $(CXXFLAGS) $(INC) -c -o $@ $<"; $(CPP) $(CXXFLAGS) $(INC) -c -o $@ $<

test: $(TESTS)
	@mkdir -p bin/tests tmp testOutput
	# $(CPP) -g -std=c++11 tests/abm_verif.cpp $(TINC) -o bin/tests/abm_verif
	# $(CPP) -g tests/binsearch.cpp -o bin/tests/binsearch
	# $(CPP) -g tests/spline_check.cpp src/SplineBasis.cpp $(TINC) -o bin/tests/splinecheck
	# $(CPP) -g tests/test_inverse.cpp src/Constant.cpp $(TINC) -o bin/tests/rate_inverse
	# $(CPP) -g $(CXXFLAGS) $(TINC) -DDEBUG tests/q_eii.cpp src/Constant.cpp src/FreeDistribution.cpp src/SplineIntegral.cpp src/SplineBasis.cpp src/Dipole.cpp -o bin/tests/q_eii
	# $(CPP) -g $(CXXFLAGS) $(TINC) -DDEBUG tests/q_tbr.cpp src/Constant.cpp src/FreeDistribution.cpp src/SplineIntegral.cpp src/SplineBasis.cpp src/Dipole.cpp -o bin/tests/q_tbr
	# $(CPP) -g $(CXXFLAGS) $(TINC) -DDEBUG tests/q_ee.cpp src/Constant.cpp src/FreeDistribution.cpp src/SplineIntegral.cpp src/SplineBasis.cpp src/Dipole.cpp -o bin/tests/q_ee
	$(CPP) -g $(CXXFLAGS) $(TINC) -DDEBUG tests/ee_dynamics.cpp src/Constant.cpp src/FreeDistribution.cpp src/SplineIntegral.cpp src/SplineBasis.cpp src/Dipole.cpp -o bin/tests/ee_dynamics
	# $(CPP) -g $(CXXFLAGS) $(TINC) -DDEBUG tests/basis_checker.cpp src/Constant.cpp src/FreeDistribution.cpp src/SplineIntegral.cpp src/SplineBasis.cpp src/Dipole.cpp -o bin/tests/basis_checker
	# $(CPP) -g $(TINC) tests/sigma_test.cpp src/Dipole.cpp -o bin/tests/sigma_test
	# $(CPP) -g $(TINC) -DDEBUG_VERBOSE tests/rate_io_test.cpp src/Constant.cpp -o bin/tests/rate_io

debug: all
release: all

scrub:
	$(RM) build/FreeDistribution.o build/ElectronSolver.o build/SplineBasis.o build/SplineIntegral.o

clean:
	@echo " Cleaning...";
	@echo " $(RM) -r $(BUILDDIR)"; $(RM) -r $(BUILDDIR)

.PHONY: clean, test, all
-include $(OBJECTS:.o=.d)
