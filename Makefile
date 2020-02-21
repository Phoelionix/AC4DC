
CPP := /usr/local/opt/llvm/bin/clang

SRCDIR := src
BUILDDIR := build
TARGET := bin/plasMD

HOME := $(shell echo $HOME)

SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
CPPFLAGS := -O3 -std=c++11 -fopenmp
LIB := -L/usr/local/opt/llvm/lib -L$(HOME)/Programming/lib
INC := -I/usr/local/opt/llvm/include -I$(HOME)/Programming/include -I/usr/local/opt/llvm/include

$(TARGET): $(OBJECTS)
	@echo " Linking..."
	@echo " $(CPP) $^ -o $(TARGET) $(LIB) "; $(CC) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR) bin
	@echo " $(CPP) $(CPPFLAGS) $(INC) -c -o $@ $<"; $(CPP) $(CPPFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning...";
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

# Tests
tester:
	$(CC) $(CFLAGS) test/tester.cpp $(INC) $(LIB) -o bin/tester

.PHONY: clean
