
CPP := g++-9

LIB := -fopenmp -L/usr/local/opt/llvm/lib -L/Users/alaric-mba/Programming/lib
INC := -I/usr/local/opt/llvm/include -I/Users/alaric-mba/Programming/include


SRCDIR := src
BUILDDIR := build
TARGET := bin/ac4dc

SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
CPPFLAGS := -std=c++11 -Xpreprocessor -fopenmp -g -O3

$(TARGET): $(OBJECTS)
	@echo " Linking... "
	@echo " $(CPP) $^ -o $(TARGET) $(LIB) "; $(CPP) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR) bin $(BUILDDIR)/Wigner
	@echo " $(CPP) $(CPPFLAGS) $(INC) -c -o $@ $<"; $(CPP) $(CPPFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning...";
	@echo " $(RM) -r $(BUILDDIR)"; $(RM) -r $(BUILDDIR)

.PHONY: clean
