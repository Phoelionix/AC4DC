
CPP := g++-9

LIB := -fopenmp

INC := -I$(HOME)/Programming/include -Iinclude

SRCDIR := src
BUILDDIR := build
TARGET := bin/ac4dc2

SRCEXT := cpp

SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT) )

OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
CXXFLAGS := -std=c++11 -fopenmp -g

debug: CXXFLAGS += -DDEBUG
release: CXXFLAGS += -O3 -DNDEBUG

# For tests:
TESTS := $(shell find tests -type f -name *.$(SRCEXT) )
TINC := $(INC) -Isrc

$(TARGET): $(OBJECTS)
	@echo " Linking $(TARGET)... ";
	@echo " $(CPP) $^ $(LIB) -o $(TARGET) "; $(CPP) $^ $(LIB) -o $(TARGET)

all: $(TARGET)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR) bin $(BUILDDIR)/Wigner $(BUILDDIR)/$(MAINSUBDIR)
	@echo " $(CPP) $(CXXFLAGS) $(INC) -c -o $@ $<"; $(CPP) $(CXXFLAGS) $(INC) -c -o $@ $<

test: $(TESTS)
	@mkdir -p bin/tests tmp testOutput
	@echo " Compiling tests/abm_verif.cpp "
	$(CPP) -g -std=c++11 tests/abm_verif.cpp $(TINC) -o bin/tests/abm_verif
	$(CPP) -g -std=c++11 tests/integral_verif.cpp src/RateSystem.cpp src/FreeDistribution.cpp src/Dipole.cpp $(TINC) -o bin/tests/integral_verif

debug: all
release: all

scrub:
	$(RM) build/FreeDistribution.o build/ElectronSolver.o build/SplineBasis.o

clean:
	@echo " Cleaning...";
	@echo " $(RM) -r $(BUILDDIR)"; $(RM) -r $(BUILDDIR)

.PHONY: clean, test, all
